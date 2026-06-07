"""Normalized signed Laplacian smoothing — direction-aware, topology-controlled.

A drop-in alternative to ``run_bayesian_inference`` that addresses the two
structural failures of the BP+gate framework exposed by perm-null
analysis (April 2026):

  1. **Topology dominance** — BP messages accumulate at hub reactions
     regardless of label content. We use the symmetric normalized
     Laplacian L_sym = I - D^{-1/2} A D^{-1/2}, so each node's
     contribution is divided by sqrt(degree) on both sides → hubs
     can't broadcast.

  2. **Direction loss** — BP couplings can't distinguish substrate UP
     (mass-action evidence reaction is UP) from product UP (also
     evidence reaction is UP, by accumulation). The signed adjacency
     A_signed encodes mass-action edges with role-aware weights so
     direction propagates through the linear smoothing solve instead
     of the BP message passing.

Algorithm
---------

Build signed weighted adjacency from graph edges:

  - SUBSTRATE  (metabolite → reaction): weight = +stoichiometry
        substrate concordance — substrate UP supports reaction UP
  - PRODUCT    (reaction → metabolite): weight = +stoichiometry
        product concordance — product UP supports reaction UP
  - MODIFIER   (gene → reaction): weight = +1
        capacity — enzyme available supports reaction UP
  - skip currency-flagged metabolites and `skip_edges` from
    conditional-currency analysis

Symmetrize: A := (A + A^T) / 2

Degree-normalize: L_sym = I - D^{-1/2} A D^{-1/2}

Evidence vector: y[node] = sum(effect_size * confidence over records)

Solve: (I + α L_sym) x = y for α > 0 (smoothness regularization)

Map x → posterior triples for downstream compatibility:
  p_up   = sigmoid(x / scale) * (1 - normal_floor)
  p_down = sigmoid(-x / scale) * (1 - normal_floor)
  p_normal = 1 - p_up - p_down

Why this should work
--------------------

- The signed evidence vector y carries direction at every observed node.
- The signed edges propagate that direction with the right sign
  (substrate UP and product UP both support reaction UP, so all three
  hang together; if substrate is UP and product is DOWN, reaction
  prediction will be near-zero — "inconsistent observation" rather than
  "average UP").
- D^{-1/2} normalization means a hub with 100 evidence-bearing neighbors
  doesn't get 100x more signal than a node with 5 neighbors of equal
  effect — the per-neighbor contribution is normalized.
- α controls how diffuse the smoothing is; small α keeps signal local
  to evidence, large α propagates farther. Roughly equivalent to
  hops=1..n in BP but with a single closed-form solve.

Complexity
----------

For sparse A with E nonzeros and N nodes, the linear solve via
scipy.sparse.linalg.cg is O(E · iter), typically ≪ BP's iterative
message passing.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy
import scipy.sparse as sp
from scipy.sparse.linalg import cg

from gizmo.schema import EdgeRole

_CG_TOL_KWARG = (
    "rtol" if tuple(int(x) for x in scipy.__version__.split(".")[:2]) >= (1, 12)
    else "tol"
)


@dataclass
class LaplacianConfig:
    """Hyperparameters for normalized signed Laplacian smoothing."""
    alpha: float = 1.0
    """Smoothness regularization. Larger α → more graph diffusion.
    α=0 returns y unchanged on observed nodes, zero elsewhere.
    α=∞ converges to a constant (graph mean)."""

    normal_floor: float = 0.05
    """Minimum p_normal in the output triple (regularizer for the
    softmax-like mapping)."""

    score_scale: float = 0.5
    """Scale factor in the sigmoid mapping x → p_up/p_down. Tunes
    how aggressively the signed score is converted to up/down
    probability mass. Smaller = sharper, larger = softer."""

    cg_tol: float = 1e-6
    """Conjugate gradient tolerance for the linear solve."""

    cg_maxiter: int = 200
    """Conjugate gradient iteration cap."""

    skip_currency: bool = True
    """Skip metabolite nodes flagged as currency (water, ATP, etc.)."""

    skip_edges: Optional[set] = None
    """Additional (u, v) tuples to skip — typically from
    compute_conditional_currency_edges."""

    substrate_weight: float = 1.0
    """Multiplier on substrate→reaction edges (sign-concordant)."""

    product_weight: float = 1.0
    """Multiplier on reaction→product edges (sign-concordant)."""

    modifier_weight: float = 1.0
    """Multiplier on gene→reaction (catalysis) edges. Lower means
    transcript signal contributes less to reaction state."""

    normalize_per_omic: bool = True
    """When True, L2-normalize the metabolite-channel and gene-channel
    evidence vectors separately before summing into the joint vector.

    Why this matters: trans panels typically have ~14k genes with
    nonzero effect_size while metab panels have ~30-1000 metabolites.
    Without normalization the larger-magnitude trans channel swamps
    metab in the linear solve, so ``both_directional`` ends up
    dominantly inheriting trans-only's top reactions and metab-only's
    top reactions get pushed down the ranking. With normalization,
    each channel contributes equal energy to the integrated solve and
    the surfaced novel biology is genuinely the union of both omics."""


@dataclass
class LaplacianResult:
    """Mirrors BP `Result` interface so callers don't need to change.

    posteriors : node_id → (p_down, p_normal, p_up)
    raw_scores : node_id → signed smoothed score x (real-valued)
    """
    posteriors: dict[str, tuple[float, float, float]]
    raw_scores: dict[str, float]


def _build_signed_adjacency(mg, cfg: LaplacianConfig,
                             gene_capacity: Optional[dict[str, float]] = None):
    """Construct the signed weighted adjacency matrix and node index.

    If ``gene_capacity`` is provided (gene_node_id → multiplier), MODIFIER
    edges from those genes are scaled by the multiplier — the capacity
    gate from the asymmetric design implemented natively in the
    adjacency rather than as a post-hoc reweighting.

    Returns (A_csr, idx, nodes_list).
    """
    g = mg.graph
    skip_edges = cfg.skip_edges or set()
    cap = gene_capacity or {}

    nodes = list(g.nodes)
    idx = {nid: i for i, nid in enumerate(nodes)}
    n = len(nodes)

    rows, cols, data = [], [], []
    for u, v, ed in g.edges(data=True):
        # Currency skip
        if cfg.skip_currency:
            if g.nodes[u].get("is_currency") or g.nodes[v].get("is_currency"):
                continue
        if (u, v) in skip_edges or (v, u) in skip_edges:
            continue

        role = ed.get("role")
        if hasattr(role, "value"):
            role_s = role.value
        else:
            role_s = str(role).lower() if role else ""

        stoich = float(ed.get("stoichiometry", 1.0)) or 1.0

        if role_s == EdgeRole.SUBSTRATE.value:
            w = cfg.substrate_weight * stoich
        elif role_s == EdgeRole.PRODUCT.value:
            w = cfg.product_weight * stoich
        elif role_s == EdgeRole.MODIFIER.value:
            w = cfg.modifier_weight
            # Asymmetric capacity gate: scale gene→reaction edges by
            # per-gene multiplier. u or v could be the gene depending on
            # edge orientation; check both.
            if u in cap:
                w *= cap[u]
            elif v in cap:
                w *= cap[v]
        else:
            edge_type = ed.get("edge_type", "").lower()
            if "catalysis" in edge_type or "gene" in edge_type:
                w = cfg.modifier_weight
                if u in cap:
                    w *= cap[u]
                elif v in cap:
                    w *= cap[v]
            else:
                continue

        i, j = idx[u], idx[v]
        rows.append(i); cols.append(j); data.append(w)

    A = sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=np.float64)
    A = (A + A.T) * 0.5
    return A.tocsr(), idx, nodes


def _build_evidence_vector(ctx, idx, normalize_per_omic: bool = True):
    """Aggregate evidence records into a length-N signed vector.

    Multiple records on the same node are summed (effect_size *
    confidence). Unobserved nodes are zero.

    When ``normalize_per_omic`` is True, the metabolite-channel and
    gene-channel evidence vectors are L2-normalized separately before
    being summed. This prevents the larger-cardinality channel (usually
    transcriptomics with thousands of features) from dominating the
    integrated solve over a sparser metabolomics channel.
    """
    n = len(idx)
    if not normalize_per_omic:
        y = np.zeros(n, dtype=np.float64)
        for rec in ctx.records():
            i = idx.get(rec.node_id)
            if i is None:
                continue
            c = rec.confidence if rec.confidence is not None else 1.0
            y[i] += float(rec.effect_size) * float(c)
        return y

    # Per-omic separated channels — L2-normalize each, then sum.
    # Both channels keep their signed direction; only their total
    # energy is rescaled to be comparable.
    METABOLITE_NODE_TYPES = {"metabolite"}
    GENE_NODE_TYPES = {"gene"}
    y_metab = np.zeros(n, dtype=np.float64)
    y_gene  = np.zeros(n, dtype=np.float64)
    y_other = np.zeros(n, dtype=np.float64)
    for rec in ctx.records():
        i = idx.get(rec.node_id)
        if i is None:
            continue
        c = rec.confidence if rec.confidence is not None else 1.0
        v = float(rec.effect_size) * float(c)
        nt = (rec.node_type or "").lower()
        if nt in METABOLITE_NODE_TYPES:
            y_metab[i] += v
        elif nt in GENE_NODE_TYPES:
            y_gene[i] += v
        else:
            y_other[i] += v

    def _norm(vec):
        n2 = float(np.linalg.norm(vec))
        return vec / n2 if n2 > 1e-12 else vec

    return _norm(y_metab) + _norm(y_gene) + _norm(y_other)


def _normalized_laplacian(A: sp.csr_matrix) -> sp.csr_matrix:
    """L_sym = I - D^{-1/2} |A| D^{-1/2}, but applied to signed A.

    We use |A| for the degree (so degree counts edge presence, not
    signed cancellation). The signed entries propagate sign through
    the solve.
    """
    abs_A = abs(A)
    deg = np.asarray(abs_A.sum(axis=1)).ravel()
    deg_safe = np.where(deg > 0, deg, 1.0)
    d_inv_sqrt = sp.diags(1.0 / np.sqrt(deg_safe))
    return sp.eye(A.shape[0]) - d_inv_sqrt @ A @ d_inv_sqrt


def _solve(L: sp.csr_matrix, y: np.ndarray, cfg: LaplacianConfig) -> np.ndarray:
    """Solve (I + α L) x = y via conjugate gradient.

    The system matrix (I + α L) is symmetric positive definite when
    L is SPSD (which the normalized Laplacian is), so CG converges.
    """
    M = sp.eye(L.shape[0]) + cfg.alpha * L
    x, info = cg(M, y, **{_CG_TOL_KWARG: cfg.cg_tol}, maxiter=cfg.cg_maxiter)
    return x


def _to_posterior(x: float, cfg: LaplacianConfig) -> tuple[float, float, float]:
    """Map signed smoothed score → (p_down, p_normal, p_up)."""
    s = x / cfg.score_scale
    p_up = 1.0 / (1.0 + np.exp(-s))
    p_down = 1.0 / (1.0 + np.exp(+s))
    # symmetric tail probabilities sum to 1; keep mass for "normal"
    pert = p_up + p_down
    if pert <= 1.0 - cfg.normal_floor:
        p_normal = 1.0 - pert
    else:
        # Renormalize so p_normal == normal_floor
        scale = (1.0 - cfg.normal_floor) / pert
        p_up *= scale; p_down *= scale
        p_normal = cfg.normal_floor
    return float(p_down), float(p_normal), float(p_up)


def run_laplacian_inference(
    mg, ctx, cfg: Optional[LaplacianConfig] = None,
    *,
    gene_capacity: Optional[dict[str, float]] = None,
) -> LaplacianResult:
    """Smooth signed evidence over the metabolic graph via normalized
    signed Laplacian regularization.

    ``gene_capacity`` (optional) maps gene_node_id → multiplicative
    capacity multiplier ∈ (0, 1]. Used by the asymmetric design to
    dampen gene→reaction edge weights when the gene is transcriptionally
    silent. Genes not in the dict get the default weight (1.0).

    Drop-in API analogue of ``gizmo.inference.run_bayesian_inference``.
    Returns posterior triples + raw signed scores per node.
    """
    cfg = cfg or LaplacianConfig()

    A, idx, nodes = _build_signed_adjacency(mg, cfg, gene_capacity)
    y = _build_evidence_vector(ctx, idx,
                                normalize_per_omic=cfg.normalize_per_omic)
    L = _normalized_laplacian(A)
    x = _solve(L, y, cfg)

    posteriors: dict[str, tuple[float, float, float]] = {}
    raw: dict[str, float] = {}
    for nid, i in idx.items():
        xi = float(x[i])
        raw[nid] = xi
        posteriors[nid] = _to_posterior(xi, cfg)
    return LaplacianResult(posteriors=posteriors, raw_scores=raw)
