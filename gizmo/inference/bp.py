"""
Loopy belief propagation over a GIZMO graph.

Public entry point: :func:`run_bayesian_inference`.

The core is a synchronous message-passing loop over a sparse edge list.
Messages are 3-vectors (one per state), stored as a ``(n_directed_edges, 3)``
numpy array. Each undirected graph edge is represented as two directed
messages ``i→j`` and ``j→i`` that are updated in lockstep using pairwise
potentials built from the edge's ``role``/``edge_type``.

Convergence is declared when the maximum per-node marginal Δ drops below
``tol`` or ``max_iter`` is reached. Damping λ on message updates avoids
oscillation on tight cycles.

Status (April 2026)
-------------------
Retained for reproducibility but does not extract label-specific
signal beyond a permutation null on IDH multi-omics benchmarks (see
``benchmarks/results/METHOD_HISTORY.md``). For new work prefer
:func:`gizmo.inference.laplacian.run_laplacian_inference`, which is
direction-aware (DR=0.889, +6.3σ above null on IDH at α=1.0),
topology-controlled (D^{-1/2} normalization), and ~200× faster.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from gizmo.inference.model import (
    DEFAULT_COUPLINGS, DEFAULT_ASSAY_SIGMAS, Couplings,
    unary_prior, gaussian_obs_likelihood, hard_confirmation_likelihood,
    pairwise_potential, DOWN, NORMAL, UP, SIGN,
)

log = logging.getLogger(__name__)


# Currency metabolite IDs to skip. These are high-degree hubs that act as
# message sinks in BP and wash out real structure. This is a conservative
# short list — extend as needed.
CURRENCY_SUBSTRINGS = (
    "water", "h2o", "atp", "adp", "amp", "nad+", "nadh", "nadp+", "nadph",
    "proton", "oxygen", "carbon dioxide", "phosphate", "diphosphate",
    "ammonia", "ammonium", "coenzyme a", "coa",
)
CURRENCY_CHEBI_IDS = {
    "CHEBI:15377",  # water
    "CHEBI:15378",  # H+
    "CHEBI:30616",  # ATP
    "CHEBI:456216", # ADP
    "CHEBI:456215", # AMP
    "CHEBI:57540",  # NAD+
    "CHEBI:57945",  # NADH
    "CHEBI:58349",  # NADP+
    "CHEBI:57783",  # NADPH
    "CHEBI:15379",  # O2
    "CHEBI:16526",  # CO2
    "CHEBI:43474",  # HPO4 2-
    "CHEBI:33019",  # diphosphate
    "CHEBI:16134",  # NH3
    "CHEBI:28938",  # NH4+
    "CHEBI:57287",  # CoA
}


def _sigma_for(rec, cfg: "BPConfig") -> float:
    """Per-assay σ with fallback to ``cfg.obs_sigma``."""
    atype = (rec.assay_type or "").lower()
    if atype and atype in cfg.assay_sigmas:
        return float(cfg.assay_sigmas[atype])
    return float(cfg.obs_sigma)


def _is_currency(nid: str, attrs: dict) -> bool:
    if attrs.get("is_currency"):
        return True
    name = (attrs.get("name") or "").strip().lower()
    if not name:
        return False
    if any(s == name for s in CURRENCY_SUBSTRINGS):
        return True
    chebi = str(attrs.get("chebi_id") or "").strip()
    if chebi and chebi.upper() in CURRENCY_CHEBI_IDS:
        return True
    return False


# ---------------------------------------------------------------------------
# Config / result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class BPConfig:
    max_iter:        int   = 50
    tol:             float = 1e-3
    damping:         float = 0.5
    normal_bias:     float = 1.5
    obs_sigma:       float = 1.0
    obs_mu:          float = 1.0
    # Per-assay override of obs_sigma keyed by ``EvidenceRecord.assay_type``.
    # Defaults derived from mRNA–protein correlation literature (see
    # ``DEFAULT_ASSAY_SIGMAS`` in gizmo.inference.model). Unknown assay types
    # fall back to ``obs_sigma``.
    assay_sigmas:    dict[str, float] = field(
        default_factory=lambda: dict(DEFAULT_ASSAY_SIGMAS)
    )
    skip_currency:   bool  = True
    # Edge-level conditional-currency skip set. Populate with
    # ``compute_conditional_currency_edges(mg)``. Edges in this set
    # are removed from the message-passing graph — useful for α-KG /
    # SAM / glutamate cofactor hubs that should be cofactors in
    # transamination / α-KG dioxygenase / methyltransferase reactions
    # but features in TCA / glutamate metabolism / methionine cycle.
    skip_edges:      Optional[set] = None
    restrict_to_observed_hops: Optional[int] = 3  # None = whole graph
    couplings:       Couplings = field(default_factory=lambda: Couplings(
        values=dict(DEFAULT_COUPLINGS.values)
    ))
    # Anchor-aware scoring: node IDs whose state is pinned to NORMAL with
    # high confidence before BP runs. Features the analyst flagged as
    # reference compounds go here.
    anchor_nodes:    Optional[set] = None
    anchor_confidence: float = 0.95
    log_every:       int   = 10


@dataclass
class BPResult:
    posteriors: dict[str, np.ndarray]   # node_id → 3-vector
    converged:  bool
    iterations: int
    final_delta: float

    def p_up(self, node_id: str) -> float:
        v = self.posteriors.get(node_id)
        return float(v[UP]) if v is not None else 0.0

    def p_down(self, node_id: str) -> float:
        v = self.posteriors.get(node_id)
        return float(v[DOWN]) if v is not None else 0.0

    def p_perturbed(self, node_id: str) -> float:
        """Total probability the node is not normal."""
        v = self.posteriors.get(node_id)
        return float(v[UP] + v[DOWN]) if v is not None else 0.0

    def signed_score(self, node_id: str) -> float:
        """
        Signed expected state in [−1, +1]:
            E[sgn(X)] = p(UP) − p(DOWN).
        """
        v = self.posteriors.get(node_id)
        if v is None:
            return 0.0
        return float(v[UP] - v[DOWN])


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def run_bayesian_inference(mg, ctx, config: Optional[BPConfig] = None) -> BPResult:
    """
    Run loopy BP on ``mg`` using evidence from ``ctx``.

    Parameters
    ----------
    mg     : GizmoGraph (compartment-collapsing recommended before calling)
    ctx    : SampleContext with EvidenceRecords carrying effect_size + confidence
    config : optional BPConfig; defaults to module defaults
    """
    cfg = config or BPConfig()
    g = mg.graph

    # ----- 1. Choose the node set ---------------------------------------
    obs_ids: set[str] = set()
    for rec in ctx.records(mapped_only=True):
        if rec.node_id and g.has_node(rec.node_id):
            obs_ids.add(rec.node_id)
    if not obs_ids:
        return BPResult({}, True, 0, 0.0)

    if cfg.restrict_to_observed_hops is None:
        candidate = set(g.nodes)
    else:
        # BFS up to k hops out of any observed node (undirected)
        frontier = set(obs_ids)
        candidate = set(frontier)
        for _ in range(cfg.restrict_to_observed_hops):
            next_frontier: set[str] = set()
            for nid in frontier:
                next_frontier.update(g.successors(nid))
                next_frontier.update(g.predecessors(nid))
            next_frontier -= candidate
            if not next_frontier:
                break
            candidate.update(next_frontier)
            frontier = next_frontier

    if cfg.skip_currency:
        candidate = {
            nid for nid in candidate
            if not (g.nodes[nid].get("node_type") == "metabolite"
                    and _is_currency(nid, g.nodes[nid]))
        }

    if not candidate:
        return BPResult({}, True, 0, 0.0)

    # ----- 2. Index nodes, build edge list ------------------------------
    node_ids = sorted(candidate)
    index = {nid: i for i, nid in enumerate(node_ids)}
    n = len(node_ids)

    # Directed edge list (we'll treat each undirected edge as two messages)
    src_idx: list[int] = []
    tgt_idx: list[int] = []
    psi_list: list[np.ndarray] = []   # 3x3 potential for each directed edge

    skip_edges = cfg.skip_edges or set()
    for u, v, attrs in g.edges(data=True):
        if u not in index or v not in index:
            continue
        # Conditional-currency skip: edges where a borderline metabolite
        # (α-KG, Glu, SAM, etc.) acts as a cofactor cosubstrate in this
        # specific reaction's EC class. Skip in either direction.
        if (u, v) in skip_edges or (v, u) in skip_edges:
            continue
        role = (attrs.get("role") or attrs.get("edge_type") or "").lower()
        theta, sign = cfg.couplings.theta_sign(role)
        if theta == 0.0:
            continue
        psi = pairwise_potential(theta, sign)
        i, j = index[u], index[v]
        # Message u→v uses psi(u, v); message v→u uses psi.T
        src_idx.append(i); tgt_idx.append(j); psi_list.append(psi)
        src_idx.append(j); tgt_idx.append(i); psi_list.append(psi.T)

    if not src_idx:
        # No usable edges → everything is prior × observation only.
        posteriors = _compute_unary_only(
            g, node_ids, index, ctx, cfg,
        )
        return BPResult(posteriors, True, 0, 0.0)

    src_arr = np.asarray(src_idx, dtype=np.int64)
    tgt_arr = np.asarray(tgt_idx, dtype=np.int64)
    psi_arr = np.stack(psi_list, axis=0)   # (E, 3, 3)
    n_edges = psi_arr.shape[0]

    # ----- 3. Unary potentials (prior × observation likelihood) --------
    unary = np.empty((n, 3), dtype=np.float64)
    for nid in node_ids:
        attrs = g.nodes[nid]
        node_type = attrs.get("node_type") or "unknown"
        unary[index[nid]] = unary_prior(node_type, normal_bias=cfg.normal_bias)

    anchor_nodes = cfg.anchor_nodes or set()

    # Overlay anchor pins first (strong NORMAL prior). Regular observations
    # applied afterwards can still overcome these if the evidence is
    # overwhelming, but anchors dampen typical noise.
    for nid in anchor_nodes:
        if nid not in index:
            continue
        like = hard_confirmation_likelihood(NORMAL, confidence=cfg.anchor_confidence)
        idx = index[nid]
        unary[idx] = unary[idx] * like
        unary[idx] /= unary[idx].sum()

    # Overlay observations
    for rec in ctx.records(mapped_only=True):
        if rec.node_id not in index:
            continue
        # Skip observations on explicit anchor nodes — they're reference
        # compounds, not signal-bearing measurements.
        if rec.node_id in anchor_nodes:
            continue
        idx = index[rec.node_id]
        ntype = (rec.node_type or g.nodes[rec.node_id].get("node_type") or "").lower()
        if ntype in ("disease", "phenotype"):
            # Confirmation of presence (UP state)
            like = hard_confirmation_likelihood(UP, confidence=0.95)
        else:
            like = gaussian_obs_likelihood(
                rec.effect_size if rec.effect_size is not None else 0.0,
                confidence=rec.confidence or 1.0,
                sigma=_sigma_for(rec, cfg),
                mu_up=cfg.obs_mu,
            )
        unary[idx] = unary[idx] * like
        unary[idx] /= unary[idx].sum()

    # ----- 4. Initialise messages --------------------------------------
    msgs = np.ones((n_edges, 3), dtype=np.float64) / 3.0

    # Precompute, for each node j, the list of incoming directed-edge indices.
    # This lets us accumulate incoming products quickly.
    incoming: list[np.ndarray] = [np.empty(0, dtype=np.int64)] * n
    # Build lists then convert
    _tmp: list[list[int]] = [[] for _ in range(n)]
    for e, j in enumerate(tgt_arr):
        _tmp[j].append(e)
    incoming = [np.asarray(lst, dtype=np.int64) for lst in _tmp]

    # Map src→list of directed-edge indices as well (for skipping reverse)
    # For each directed edge e = (i, j), its reverse is at e^1 (XOR 1) because
    # we inserted pairs consecutively.
    # -> reverse index is e XOR 1 (works because pairs are sequential)

    # ----- 5. Main BP loop ---------------------------------------------
    final_delta = 0.0
    iterations = 0
    converged = False
    prev_marginals = _marginals(unary, msgs, incoming)

    for it in range(cfg.max_iter):
        new_msgs = np.empty_like(msgs)
        # For each directed edge e = (i → j):
        #   outgoing factor at i = unary[i] × Π_{incoming_to_i, excluding rev(e)} msgs
        for e in range(n_edges):
            i = src_arr[e]
            j = tgt_arr[e]
            rev = e ^ 1   # reverse directed edge index
            # Product of all messages coming into i, except the one from j
            inc = incoming[i]
            if inc.size == 0:
                prod = np.ones(3)
            else:
                prod = np.prod(msgs[inc], axis=0)
                # Divide out rev message (guard against zero)
                rev_msg = msgs[rev]
                prod = np.where(rev_msg > 1e-30, prod / rev_msg, prod)
            node_factor = unary[i] * prod
            # Message m_{i→j}(x_j) = Σ_{x_i} ψ(x_i, x_j) × node_factor(x_i)
            m = psi_arr[e].T @ node_factor   # psi[e] rows=x_i, cols=x_j
            s = m.sum()
            if s <= 0:
                m = np.ones(3) / 3.0
            else:
                m = m / s
            new_msgs[e] = m

        # Damping
        msgs = cfg.damping * msgs + (1.0 - cfg.damping) * new_msgs
        # Renormalise
        msgs = msgs / msgs.sum(axis=1, keepdims=True).clip(min=1e-30)

        marginals = _marginals(unary, msgs, incoming)
        delta = float(np.max(np.abs(marginals - prev_marginals)))
        prev_marginals = marginals
        final_delta = delta
        iterations = it + 1

        if cfg.log_every and (it + 1) % cfg.log_every == 0:
            log.info("BP iter %d  max Δmarginal = %.2e", it + 1, delta)

        if delta < cfg.tol:
            converged = True
            break

    # ----- 6. Package posteriors ---------------------------------------
    posteriors = {nid: marginals[index[nid]].copy() for nid in node_ids}
    return BPResult(posteriors, converged, iterations, final_delta)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _marginals(unary, msgs, incoming) -> np.ndarray:
    """Compute node marginals = normalised unary × Π incoming messages."""
    n = unary.shape[0]
    out = unary.copy()
    for i in range(n):
        inc = incoming[i]
        if inc.size:
            out[i] *= np.prod(msgs[inc], axis=0)
    # Normalise
    s = out.sum(axis=1, keepdims=True).clip(min=1e-30)
    return out / s


def _compute_unary_only(g, node_ids, index, ctx, cfg) -> dict[str, np.ndarray]:
    """Fallback when no edges match: return prior × observation."""
    posteriors: dict[str, np.ndarray] = {}
    for nid in node_ids:
        attrs = g.nodes[nid]
        p = unary_prior(attrs.get("node_type") or "unknown",
                        normal_bias=cfg.normal_bias)
        posteriors[nid] = p
    for rec in ctx.records(mapped_only=True):
        if rec.node_id not in index:
            continue
        ntype = (rec.node_type or g.nodes[rec.node_id].get("node_type") or "").lower()
        if ntype in ("disease", "phenotype"):
            like = hard_confirmation_likelihood(UP, confidence=0.95)
        else:
            like = gaussian_obs_likelihood(
                rec.effect_size or 0.0,
                confidence=rec.confidence or 1.0,
                sigma=_sigma_for(rec, cfg),
                mu_up=cfg.obs_mu,
            )
        posteriors[rec.node_id] = posteriors[rec.node_id] * like
        s = posteriors[rec.node_id].sum()
        if s > 0:
            posteriors[rec.node_id] /= s
    return posteriors


__all__ = ["run_bayesian_inference", "BPConfig", "BPResult"]
