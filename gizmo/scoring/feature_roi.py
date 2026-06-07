"""Per-feature ROI = predictive contribution / unit cost.

Extracted from benchmarks/diagnostics/per_feature_value_cost_v2.py so other
tools (SQuID-INC panel designer, downstream notebooks) can import the
decomposition + global-table lookup without re-running the cohort sweep.

The benchmark script remains the source of truth for *computing* the global
ROI table. This module owns the math primitives + a thin loader.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def _repo_root() -> Path:
    """Return the GIZMO repo root from this module's location."""
    return Path(__file__).resolve().parent.parent.parent

# Per-feature unit costs ($/sample) — derived from typical commercial panels.
COST_PROT_BROAD = 0.40      # Olink Explore
COST_PROT_TARGETED = 2.60   # Olink 96-plex
COST_RNA_BROAD = 0.027      # RNA-seq WTS
COST_RNA_PANEL = 0.31       # nanoString
COST_METAB_BROAD = 0.80     # Metabolon HD4
COST_METAB_TARGET = 6.00    # targeted LC-MS
COST_NMR_PER_METAB = 5.00


_DEFAULT_GLOBAL_ROI_PATH = _repo_root() / "benchmarks" / "results" / "per_feature_global_roi.tsv"


def decompose_alpha(F: np.ndarray, log_pr: np.ndarray) -> np.ndarray:
    """Project F onto the orthogonal complement of the log-PageRank axis.

    Used by Paper-1 substrate-Laplacian inference to separate β (severity-
    aligned, log_PR direction) from α (orthogonal residual carrying the
    discriminative biology).
    """
    F_norm = np.linalg.norm(F, axis=1, keepdims=True) + 1e-12
    F = F / F_norm
    x = log_pr
    x_mean = x.mean()
    x_var = x.var() + 1e-12
    F_mean = F.mean(axis=1, keepdims=True)
    cov = ((F - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    return F - F_mean - beta[:, None] * (x - x_mean)[None, :]


def compute_per_feature_roi(
    *,
    pc_loadings: np.ndarray,
    nid_idx: dict[str, int],
    feat_to_node: dict[str, str],
    feature_values: dict[str, dict[str, float]],
    patients: Iterable[str],
    cost_per_feature: float,
) -> pd.DataFrame:
    """Compute per-feature contribution + ROI for one cohort × one modality.

    Inputs mirror the per-cohort loop in the benchmark driver:
      - pc_loadings: best discriminating α-PC's loading vector (n_nodes,)
      - nid_idx: node_id → row in pc_loadings
      - feat_to_node: feature_name → substrate node_id
      - feature_values: {patient_id: {feature: value}}
      - patients: ordered patient ids whose values we should pull
      - cost_per_feature: $/sample for one feature on this modality
    """
    patients = list(patients)
    data_lc = {str(p).lower(): v for p, v in feature_values.items()}
    rows = []
    for feat, nid in feat_to_node.items():
        if nid not in nid_idx:
            continue
        loading = abs(float(pc_loadings[nid_idx[nid]]))
        vals = np.array(
            [data_lc[str(p).lower()].get(feat, np.nan)
             for p in patients if str(p).lower() in data_lc],
            dtype=float,
        )
        if vals.size == 0:
            continue
        var = float(np.nanstd(vals)) if (~np.isnan(vals)).any() else 0.0
        if var == 0:
            continue
        contribution = loading * var
        roi = contribution / cost_per_feature if cost_per_feature > 0 else 0.0
        rows.append({
            "feature": str(feat),
            "anchor_node": str(nid),
            "loading_at_node": loading,
            "feature_variance": var,
            "predictive_contribution": float(contribution),
            "cost_per_feature_usd": float(cost_per_feature),
            "roi": float(roi),
        })
    return pd.DataFrame(rows)


@dataclass
class GlobalROITable:
    """Loaded global ROI table indexed for fast lookups."""
    df: pd.DataFrame
    by_anchor: dict[str, dict]

    @classmethod
    def load(cls, path: Path | str | None = None) -> "GlobalROITable | None":
        p = Path(path) if path else _DEFAULT_GLOBAL_ROI_PATH
        if not p.exists():
            return None
        df = pd.read_csv(p, sep="\t")
        by_anchor: dict[str, dict] = {}
        for _, r in df.iterrows():
            by_anchor[str(r["anchor_node"])] = r.to_dict()
        return cls(df=df, by_anchor=by_anchor)

    def lookup(self, anchor_node: str) -> dict | None:
        return self.by_anchor.get(str(anchor_node))

    def filter_by_modality(self, modality: str) -> pd.DataFrame:
        return self.df[self.df["modality"] == modality]
