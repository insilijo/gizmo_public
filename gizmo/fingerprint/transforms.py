"""Per-patient multi-resolution fingerprint transform.

Collapses the (n_patients, n_nodes) F matrix into ~30 engineered features
suitable as a supervised-learning input table (see
``discomarker.methods.gizmo_fingerprint``).

Fingerprint columns (per patient):

    beta              — projection onto log_PR direction
    alpha_norm        — ‖α‖₂, mechanism magnitude
    frac_metabolite   — |α| mass fraction on metabolite nodes
    frac_gene         — |α| mass fraction on gene nodes
    frac_reaction     — |α| mass fraction on reaction nodes
    module_0..K-1     — signed α mass on top-K modules (ranked by group ‖α_mod‖_F)
    peak_0..K-1       — signed α at top-K nodes (ranked by mean |α|)

Module / peak provenance is attached to the returned DataFrame's ``attrs``
dict so callers can map column names back to substrate node identities.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

_NODE_TYPE_FRACTION_KEYS = ("metabolite", "gene", "reaction")


def _decompose_beta_alpha(
    F: np.ndarray, log_pr: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """β/α decomposition matching ``gizmo.inference.projection.decompose_beta_alpha``.

    Returns ``(beta, alpha, alpha_norm)``. Unlike the projection version,
    this one returns the full α residual matrix (needed for fingerprint
    features and downstream α-residual clustering).
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
    return beta, alpha, alpha_norm


def multi_resolution_fingerprint(
    F: np.ndarray,
    log_pr: np.ndarray,
    modules: dict[int, list[str]],
    node_types: list[str],
    nodes: list[str],
    *,
    top_k_modules: int = 10,
    top_k_peaks: int = 10,
) -> dict[str, Any]:
    """Engineer a per-patient fingerprint from a F matrix.

    Parameters
    ----------
    F : (n_patients, n_nodes) state vectors. Caller is responsible for
        subsetting F rows to the patients of interest.
    log_pr : (n_nodes,) per-node log10 PageRank from the bundle.
    modules : ``{cluster_idx: [node_id]}``, as returned by
        :func:`gizmo.fingerprint.load_bundle`. Node IDs reference ``nodes``.
    node_types : per-column node type (one of metabolite/gene/reaction/other),
        parallel to F columns.
    nodes : per-column node ID, parallel to F columns. Used to resolve
        module node IDs back to column indices.
    top_k_modules : number of modules to expose as columns.
    top_k_peaks : number of single-node peaks to expose as columns.

    Returns
    -------
    dict with two entries:
        ``fingerprint``     : pd.DataFrame, rows = patients (anonymous int
                              index — caller is responsible for re-indexing),
                              cols = engineered features.
        ``alpha_residual``  : np.ndarray (n_patients, n_nodes), used for
                              α-residual clustering downstream.

    The DataFrame's ``attrs`` carries provenance: ``module_node_ids`` maps
    each ``module_i`` column to its substrate node IDs; ``peak_node_ids``
    maps each ``peak_i`` column to its node ID.
    """
    F = np.asarray(F)
    log_pr = np.asarray(log_pr)
    if F.ndim != 2:
        raise ValueError(f"F must be 2D, got shape {F.shape}")
    n_patients, n_nodes = F.shape
    if log_pr.shape[0] != n_nodes:
        raise ValueError(
            f"log_pr length {log_pr.shape[0]} != F columns {n_nodes}"
        )
    if len(node_types) != n_nodes or len(nodes) != n_nodes:
        raise ValueError(
            f"node_types ({len(node_types)}) / nodes ({len(nodes)}) must "
            f"match F columns ({n_nodes})"
        )
    if n_patients == 0:
        raise ValueError("F has zero patient rows")

    # β / α decomposition
    beta, alpha, alpha_norm = _decompose_beta_alpha(F, log_pr)
    abs_alpha = np.abs(alpha)  # (n_patients, n_nodes)

    columns: dict[str, np.ndarray] = {
        "beta": beta.astype(np.float64),
        "alpha_norm": alpha_norm.astype(np.float64),
    }

    # Modality fractions: per-patient |α| mass split by node_type
    node_types_arr = np.asarray(node_types)
    per_patient_total = abs_alpha.sum(axis=1) + 1e-12
    for key in _NODE_TYPE_FRACTION_KEYS:
        mask = node_types_arr == key
        if mask.any():
            mass = abs_alpha[:, mask].sum(axis=1)
            columns[f"frac_{key}"] = (mass / per_patient_total).astype(np.float64)
        else:
            columns[f"frac_{key}"] = np.zeros(n_patients, dtype=np.float64)

    # Modules: select top-K by group-level ‖α_mod‖_F, then per-patient signed sum
    node_id_to_idx = {nid: i for i, nid in enumerate(nodes)}
    module_keys = sorted(modules.keys())
    module_provenance: list[tuple[int, list[str]]] = []
    module_scores_global: list[tuple[int, float, np.ndarray]] = []
    for cid in module_keys:
        idx = np.array(
            [node_id_to_idx[nid] for nid in modules[cid] if nid in node_id_to_idx],
            dtype=int,
        )
        if idx.size == 0:
            continue
        mod_alpha = alpha[:, idx]
        frob = float(np.linalg.norm(mod_alpha, ord="fro"))
        module_scores_global.append((cid, frob, idx))
    module_scores_global.sort(key=lambda t: -t[1])
    selected_modules = module_scores_global[:top_k_modules]

    for j, (cid, _frob, idx) in enumerate(selected_modules):
        signed_score = alpha[:, idx].sum(axis=1)
        columns[f"module_{j}"] = signed_score.astype(np.float64)
        module_provenance.append((cid, [nodes[i] for i in idx]))

    # Pad with zero columns if fewer than top_k_modules survived
    for j in range(len(selected_modules), top_k_modules):
        columns[f"module_{j}"] = np.zeros(n_patients, dtype=np.float64)
        module_provenance.append((-1, []))

    # Peaks: top-K nodes by mean |α|, exposed as signed α per patient
    mean_abs_alpha = abs_alpha.mean(axis=0)
    k_peaks = min(top_k_peaks, n_nodes)
    peak_node_idx = np.argsort(-mean_abs_alpha)[:k_peaks]
    peak_node_ids: list[str] = []
    for j, node_idx in enumerate(peak_node_idx):
        columns[f"peak_{j}"] = alpha[:, node_idx].astype(np.float64)
        peak_node_ids.append(nodes[node_idx])
    for j in range(k_peaks, top_k_peaks):
        columns[f"peak_{j}"] = np.zeros(n_patients, dtype=np.float64)
        peak_node_ids.append("")

    fp = pd.DataFrame(columns)
    fp.attrs["module_node_ids"] = {
        f"module_{j}": prov[1] for j, prov in enumerate(module_provenance)
    }
    fp.attrs["module_cluster_ids"] = {
        f"module_{j}": prov[0] for j, prov in enumerate(module_provenance)
    }
    fp.attrs["peak_node_ids"] = {
        f"peak_{j}": peak_node_ids[j] for j in range(top_k_peaks)
    }

    return {"fingerprint": fp, "alpha_residual": alpha}
