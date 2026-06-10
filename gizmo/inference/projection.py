"""GIZMO Paper 1 per-patient projection pipeline (importable).

Wraps the substrate-projection methodology from the GIZMO Paper 1
manuscript into a clean importable API for downstream consumers
(GrAndMA Django service layer, notebooks, CLI).

Pipeline:
  1. Build biochem subgraph from a MultiGraph + compute Laplacian + log_PR
  2. Solve MAP per patient → per-patient state vector F over substrate nodes
  3. β/α decompose: β = OLS projection onto log_PR; α = orthogonal residual
  4. PCA on α → top α-PC directions
  5. Extract signed basins per α-PC: two connected sub-graphs naming both
     poles of the patient axis

Reference: stage3_map_reconstruction.py + loo_pc_stability.py in
benchmarks/, refactored here for production use.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, cg, expm_multiply
from sklearn.decomposition import PCA

log = logging.getLogger(__name__)

# Node types eligible for substrate projection (no pathway / disease / drug)
_PROPAGATION_TYPES = frozenset({"reaction", "metabolite", "gene"})


def _within_patient_zscore(
    data: dict[str, dict[str, float]],
    log_transform: bool = True,
    min_features: int = 10,
) -> dict[str, dict[str, float]]:
    """Within-patient log + z-score preprocessing (Paper 1 v6 §5 diagnostic).

    For each patient, take the patient's positive feature values, optionally
    ``log2(v + 1)``-transform, then z-score using **this patient's own** mean
    and std. Patients with fewer than ``min_features`` positive features are
    excluded.

    Used as a fine-grained subtype/anchor-recovery diagnostic preprocessing
    (Manuscript v6 §5 + Methods §"Within-patient z-score") that surfaces
    strong-driver signal absorbed by per-modality global-std normalization.
    NOT a canonical replacement for cohort-level discrimination preprocessing.

    Parameters
    ----------
    data : dict[patient_id, dict[feature_name, value]]
        Per-patient sparse feature observations.
    log_transform : bool, default True
        If True, apply ``log2(v + 1)`` before z-scoring. Already-log
        modalities (Olink NPX, microarray log-ratios) should pass False.
    min_features : int, default 10
        Minimum positive features per patient to compute z-score. Patients
        below this floor get an empty dict.

    Returns
    -------
    dict[patient_id, dict[feature_name, z_scored_value]]
        Same shape as input; values replaced by per-patient z-scores.
    """
    out = {}
    for sid, feat_dict in data.items():
        vals = np.array([v for v in feat_dict.values() if v > 0],
                         dtype=np.float64)
        if len(vals) < min_features:
            out[sid] = {}
            continue
        xs = np.log2(vals + 1.0) if log_transform else vals
        mu, sd = xs.mean(), xs.std() + 1e-9
        zd = {}
        for f, v in feat_dict.items():
            if v > 0:
                x = np.log2(v + 1.0) if log_transform else v
                zd[f] = float((x - mu) / sd)
        out[sid] = zd
    return out


@dataclass
class SubstrateGeometry:
    """Substrate Laplacian + log_PR direction needed for projection."""
    sub: nx.Graph              # undirected biochem subgraph
    nodes: list[str]           # ordered node IDs
    nid_idx: dict[str, int]    # node_id → row index
    L: sp.csr_matrix           # combinatorial Laplacian (D − A)
    log_pr: np.ndarray         # log10 PageRank per node, the β direction


@dataclass
class ModalitySetup:
    """One per modality (prot / metab / etc.). Built by callers."""
    label: str
    sigma: float
    diffusion_t: float
    feature_cols: list[tuple[str, int]]   # [(feature_name, node_col_idx)]
    data: dict[str, dict[str, float]]      # {patient_id: {feature: zscored_value}}


@dataclass
class Paper1Result:
    """Full Paper 1 output for one cohort."""
    patient_ids: list[str]
    F: np.ndarray                         # n_patients × n_nodes
    beta: np.ndarray                      # n_patients
    alpha_norm: np.ndarray                # n_patients (‖α‖₂)
    alpha_pc_scores: np.ndarray           # n_patients × n_components
    alpha_pc_components: np.ndarray       # n_components × n_nodes
    alpha_pc_explained_variance: np.ndarray  # n_components
    signed_basins: list[dict[str, Any]]   # one per α-PC
    smoothness: np.ndarray                # n_patients (f^T L f)
    diagnostics: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Substrate geometry
# ---------------------------------------------------------------------------

def build_biochem_subgraph(mg, hub_cap: int = 200) -> SubstrateGeometry:
    """Build the substrate subgraph used for Paper 1 projection.

    Drops pathway / disease / drug nodes (they're not propagation targets),
    drops currency-flagged metabolites, and caps hubs at ``hub_cap`` degree.

    Parameters
    ----------
    mg : MultiGraph-like with a ``.graph`` attribute (networkx)
        The full GIZMO graph (Reactome + StringDB + HMDB + KEGG + ...).
    hub_cap : int
        Maximum degree allowed in the subgraph.

    Returns
    -------
    SubstrateGeometry
        Laplacian + log_PR + node ordering needed for projection.
    """
    g = mg.graph
    keep = [
        n for n, a in g.nodes(data=True)
        if a.get("node_type") in _PROPAGATION_TYPES
        and not a.get("is_currency")
        and g.degree(n) <= hub_cap
    ]
    sub = g.subgraph(keep).copy()
    if sub.is_directed():
        sub = sub.to_undirected()
    nodes = sorted(sub.nodes())
    nid_idx = {n: i for i, n in enumerate(nodes)}

    # Combinatorial Laplacian L = D − A (NOT normalized; matches stage3 MAP)
    A = nx.adjacency_matrix(sub, nodelist=nodes).astype(float)
    deg = np.asarray(A.sum(axis=1)).ravel()
    L = (sp.diags(deg) - A).tocsr()

    # log_PR for β direction
    pr = nx.pagerank(sub)
    log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)

    return SubstrateGeometry(sub=sub, nodes=nodes, nid_idx=nid_idx, L=L, log_pr=log_pr)


# ---------------------------------------------------------------------------
# MAP solve per patient
# ---------------------------------------------------------------------------

def solve_map(
    geometry: SubstrateGeometry,
    modality_setups: list[ModalitySetup],
    patient_ids: list[str],
    *,
    lambda_smooth: float | None = None,
    lambda_calibration: float = 0.1,
    ridge_alpha: float = 1e-6,
    cg_rtol: float = 1e-8,
    cg_maxiter: int = 2000,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the per-patient MAP system.

    Per patient i, solve:

        (L + Σ_k (1/σ_k²) P_k^T P_k) f_i = Σ_k (1/σ_k²) P_k^T x_{i,k_smoothed}

    where x_{i,k_smoothed} = expm(−t_k L) x_{i,k} is the modality heat-kernel
    pre-smoothed input. Sparse CG with Jacobi preconditioning.

    Returns
    -------
    F : (n_patients, n_nodes) float32 — per-patient state vectors
    smoothness : (n_patients,) float32 — f_i^T L f_i (graph smoothness)
    """
    n = len(geometry.nodes)
    L = geometry.L

    # λ_smooth: if not provided, derive from σ via formula
    if lambda_smooth is None:
        mean_data_weight = float(np.mean([1.0 / (s.sigma ** 2) for s in modality_setups]))
        # crude lambda_median proxy from substrate spectrum (use 0.05 as fallback)
        lambda_median = 0.05
        lambda_smooth = lambda_calibration * mean_data_weight / (lambda_median * n)

    base_M = lambda_smooth * L + sp.diags(np.full(n, ridge_alpha))
    F = np.zeros((len(patient_ids), n), dtype=np.float32)
    smoothness = np.zeros(len(patient_ids), dtype=np.float32)

    for i, sid in enumerate(patient_ids):
        rhs = np.zeros(n)
        per_patient_anchor = np.zeros(n)
        for s in modality_setups:
            if sid not in s.data:
                continue
            patient_obs = s.data[sid]
            sigma_sq_inv = 1.0 / (s.sigma ** 2)
            x_sparse = np.zeros(n)
            for f, col in s.feature_cols:
                if f in patient_obs:
                    x_sparse[col] += float(patient_obs[f])
                    per_patient_anchor[col] += sigma_sq_inv
            # Modality heat-kernel pre-smoothing
            if s.diffusion_t > 0:
                x_smoothed = expm_multiply(-s.diffusion_t * L, x_sparse)
            else:
                x_smoothed = x_sparse
            rhs += x_smoothed * sigma_sq_inv

        M_i = (base_M + sp.diags(per_patient_anchor)).tocsr()
        M_diag = M_i.diagonal()
        M_inv = 1.0 / np.where(M_diag > 1e-12, M_diag, 1.0)
        Minv_op = LinearOperator(
            M_i.shape, matvec=lambda v, w=M_inv: v * w, dtype=M_i.dtype
        )
        try:
            # scipy <1.12 uses tol=, ≥1.12 uses rtol=. Detect once per call.
            try:
                f_sol, info = cg(M_i, rhs, M=Minv_op, rtol=cg_rtol, maxiter=cg_maxiter)
            except TypeError:
                f_sol, info = cg(M_i, rhs, M=Minv_op, tol=cg_rtol, maxiter=cg_maxiter)
            if info > 0:
                log.warning("MAP CG did not converge for patient %s (%d iters)", sid, info)
        except Exception as exc:
            log.exception("MAP solve failed for patient %s: %s", sid, exc)
            f_sol = np.zeros(n)
        F[i] = f_sol.astype(np.float32)
        smoothness[i] = float(f_sol @ (L @ f_sol))

    return F, smoothness


# ---------------------------------------------------------------------------
# β/α decomposition + α-PCA
# ---------------------------------------------------------------------------

def decompose_beta_alpha(
    F: np.ndarray,
    log_pr: np.ndarray,
    n_components: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, PCA]:
    """Decompose F into β (log_PR projection) + α (orthogonal residual).

    Then run PCA on α.

    Returns
    -------
    beta : (n_patients,) — per-patient projection onto log_PR direction
    alpha_norm : (n_patients,) — ‖α_i‖₂
    alpha_pc_scores : (n_patients, n_components) — α projected onto top PCs
    pca : sklearn PCA object (so callers can inspect components_, explained_variance_ratio_)
    """
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr
    x_mean = x.mean()
    x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    alpha_norm = np.linalg.norm(alpha, axis=1)
    pca = PCA(n_components=min(n_components, alpha.shape[0]), random_state=0)
    alpha_pc_scores = pca.fit_transform(alpha)
    return beta, alpha_norm, alpha_pc_scores, pca


# ---------------------------------------------------------------------------
# Signed-basin extraction
# ---------------------------------------------------------------------------

def extract_signed_basin(
    pc_loadings: np.ndarray,
    geometry: SubstrateGeometry,
    *,
    top_k: int = 15,
    min_mass_fraction: float = 0.5,
) -> dict[str, Any]:
    """Extract two connected sub-graphs (one per sign of the PC) that name
    the poles of the patient axis along this α-PC.

    Each basin = largest connected component within the substrate subgraph
    induced by the top-K most-loaded nodes on that sign.

    Returns dict with:
      - pos_basin : {nodes, sum_loading, size}
      - neg_basin : {nodes, sum_loading, size}
      - top_k : the K used
      - mass_fraction : fraction of total |loadings|^2 captured by selected nodes
    """
    sub = geometry.sub
    nodes = geometry.nodes
    abs_sq = pc_loadings ** 2
    total = abs_sq.sum() + 1e-12

    # Sort separately by sign
    pos_idx = np.where(pc_loadings > 0)[0]
    neg_idx = np.where(pc_loadings < 0)[0]
    pos_top = pos_idx[np.argsort(-pc_loadings[pos_idx])][:top_k]
    neg_top = neg_idx[np.argsort(pc_loadings[neg_idx])][:top_k]

    def _largest_cc(node_ids: list[str]) -> tuple[list[str], int]:
        if not node_ids:
            return [], 0
        induced = sub.subgraph(node_ids).copy()
        if induced.number_of_nodes() == 0:
            return [], 0
        ccs = list(nx.connected_components(induced))
        ccs.sort(key=len, reverse=True)
        largest = list(ccs[0]) if ccs else []
        return largest, len(largest)

    pos_node_ids = [nodes[i] for i in pos_top]
    neg_node_ids = [nodes[i] for i in neg_top]
    pos_cc, pos_size = _largest_cc(pos_node_ids)
    neg_cc, neg_size = _largest_cc(neg_node_ids)

    mass_top = abs_sq[pos_top].sum() + abs_sq[neg_top].sum()

    return {
        "pos_basin": {
            "nodes": pos_cc,
            "sum_loading": float(pc_loadings[pos_top].sum()),
            "size": pos_size,
            "candidates": pos_node_ids,
        },
        "neg_basin": {
            "nodes": neg_cc,
            "sum_loading": float(pc_loadings[neg_top].sum()),
            "size": neg_size,
            "candidates": neg_node_ids,
        },
        "top_k": top_k,
        "mass_fraction_top_2k": float(mass_top / total),
    }


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------

def run_paper1_pipeline(
    geometry: SubstrateGeometry,
    modality_setups: list[ModalitySetup],
    patient_ids: list[str],
    *,
    n_components: int = 5,
    basin_top_k: int = 15,
    lambda_smooth: float | None = None,
) -> Paper1Result:
    """End-to-end Paper 1 projection: MAP solve → β/α → α-PCA → signed basins.

    This is the load-bearing entry point used by the GrAndMA service layer
    (apps/gizmo/services.py:run_gizmo_paper1_projection).
    """
    log.info(
        "Paper 1 projection: %d patients × %d nodes × %d modalities",
        len(patient_ids), len(geometry.nodes), len(modality_setups),
    )

    # 1. MAP solve per patient
    F, smoothness = solve_map(
        geometry, modality_setups, patient_ids, lambda_smooth=lambda_smooth
    )

    # 2. β/α decomposition + α-PCA
    beta, alpha_norm, alpha_pc_scores, pca = decompose_beta_alpha(
        F, geometry.log_pr, n_components=n_components
    )

    # 3. Signed-basin extraction per α-PC
    signed_basins = []
    for k in range(pca.components_.shape[0]):
        basin = extract_signed_basin(
            pca.components_[k], geometry, top_k=basin_top_k
        )
        basin["pc_index"] = k + 1
        basin["explained_variance_ratio"] = float(pca.explained_variance_ratio_[k])
        signed_basins.append(basin)

    diagnostics = {
        "n_patients": len(patient_ids),
        "n_nodes": len(geometry.nodes),
        "n_modalities": len(modality_setups),
        "modalities": [
            {"label": s.label, "sigma": s.sigma, "n_features": len(s.feature_cols)}
            for s in modality_setups
        ],
        "mean_smoothness": float(smoothness.mean()),
        "alpha_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
    }

    return Paper1Result(
        patient_ids=patient_ids,
        F=F,
        beta=beta,
        alpha_norm=alpha_norm,
        alpha_pc_scores=alpha_pc_scores,
        alpha_pc_components=pca.components_,
        alpha_pc_explained_variance=pca.explained_variance_ratio_,
        signed_basins=signed_basins,
        smoothness=smoothness,
        diagnostics=diagnostics,
    )


def result_to_json_summary(
    result: Paper1Result,
    geometry: SubstrateGeometry,
    *,
    top_node_subset: int = 200,
) -> dict[str, Any]:
    """Convert a Paper1Result into JSON-serialisable dict for GizmoRun storage.

    Drops the full F matrix (38k × N is too large for JSONField); keeps
    a top-K substrate-node subset by mean |F| for visualization.
    """
    mean_abs_F = np.abs(result.F).mean(axis=0)
    top_node_idx = np.argsort(-mean_abs_F)[:top_node_subset]
    top_node_ids = [geometry.nodes[i] for i in top_node_idx]
    F_subset = result.F[:, top_node_idx]

    return {
        "patient_ids": list(result.patient_ids),
        "beta_per_patient": result.beta.astype(float).tolist(),
        "alpha_norm_per_patient": result.alpha_norm.astype(float).tolist(),
        "alpha_pc_scores": result.alpha_pc_scores.astype(float).tolist(),
        "alpha_pc_explained_variance": result.alpha_pc_explained_variance.astype(float).tolist(),
        "smoothness_per_patient": result.smoothness.astype(float).tolist(),
        "signed_basins": [
            {
                "pc_index": b["pc_index"],
                "explained_variance_ratio": b["explained_variance_ratio"],
                "mass_fraction_top_2k": b["mass_fraction_top_2k"],
                "pos_basin": {
                    "nodes": list(b["pos_basin"]["nodes"]),
                    "sum_loading": b["pos_basin"]["sum_loading"],
                    "size": b["pos_basin"]["size"],
                },
                "neg_basin": {
                    "nodes": list(b["neg_basin"]["nodes"]),
                    "sum_loading": b["neg_basin"]["sum_loading"],
                    "size": b["neg_basin"]["size"],
                },
            }
            for b in result.signed_basins
        ],
        "top_substrate_nodes": top_node_ids,
        "F_top_subset": F_subset.astype(float).tolist(),
        "diagnostics": result.diagnostics,
    }
