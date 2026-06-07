"""Augment multi-PC × metadata table with per-PC variance fraction +
Bonferroni-corrected significance (R9 vectors 2 + 4).

Reads:
  - multi_pc_vs_mofa_factors.tsv (per-cohort × metadata × {gizmo, mofa} strength + raw_p)
  - stage3_F_<cohort>.npz (to recompute α-PCA + explained_variance_ratio)
  - mofa_weights_<cohort>.json (factor_explained_variance if available)

Adds columns:
  - gizmo_pc_var_frac   — fraction of α-variance explained by the winning α-PC
  - mofa_factor_var_frac — fraction of MOFA-total-factor-variance for winning factor
  - n_cohort_tests       — N_PC × N_metadata search space per cohort
  - bonf_cohort_thr      — 0.05 / n_cohort_tests
  - bonf_global_thr      — 0.05 / total_tests_panel
  - gizmo_survives_bonf_cohort, gizmo_survives_bonf_global  (bool)
  - mofa_survives_bonf_cohort, mofa_survives_bonf_global    (bool)

Also writes a per-cohort search-space disclosure table.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

REPO = Path("/home/jgardner/GIZMO")
RESULTS = REPO / "benchmarks/results"
UR = RESULTS / "unsupervised"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))


def get_log_pr():
    import networkx as nx
    from gizmo.export.json_export import read_json
    from per_patient_wlsp_v2 import biochem_subgraph
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    sub_dir, nodes, _ = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    return np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)


def cohort_alpha_pc_variance(cohort, log_pr):
    """Return dict pc_idx (1-based) → explained_variance_ratio."""
    F_path = UR / f"stage3_F_{cohort}.npz"
    if not F_path.exists(): return {}
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"].astype(np.float64)
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    pca = PCA(n_components=5, random_state=0)
    pca.fit_transform(alpha)
    return {k + 1: float(pca.explained_variance_ratio_[k]) for k in range(5)}


def cohort_mofa_factor_variance(cohort):
    """Return dict factor_idx (1-based) → explained_variance_ratio (relative to total
    factor variance), or {} if unavailable."""
    w = UR / "mofa_weights" / f"mofa_weights_{cohort}.json"
    if not w.exists(): return {}
    d = json.load(open(w))
    scores = d.get("factor_scores")
    if not scores: return {}
    scores = np.array(scores)
    var = scores.var(axis=0)
    total = var.sum() + 1e-12
    return {k + 1: float(var[k] / total) for k in range(scores.shape[1])}


def main():
    print("Loading log_PR…", flush=True)
    log_pr = get_log_pr()

    df = pd.read_csv(RESULTS / "multi_pc_vs_mofa_factors.tsv", sep="\t")

    # Per-cohort variance dicts
    cohorts = sorted(df["cohort"].unique())
    giz_var = {c: cohort_alpha_pc_variance(c, log_pr) for c in cohorts}
    mof_var = {c: cohort_mofa_factor_variance(c) for c in cohorts}

    # Per-cohort search space
    n_per_cohort = df.groupby("cohort").size()
    n_total = len(df)
    bonf_global = 0.05 / n_total
    df["n_cohort_tests"] = df["cohort"].map(lambda c: int(n_per_cohort[c]))
    df["bonf_cohort_thr"] = 0.05 / df["n_cohort_tests"]
    df["bonf_global_thr"] = bonf_global

    # Variance fractions
    df["gizmo_pc_var_frac"] = df.apply(
        lambda r: giz_var.get(r["cohort"], {}).get(int(r["gizmo_best_pc"]), np.nan)
        if pd.notna(r["gizmo_best_pc"]) else np.nan, axis=1)
    df["mofa_factor_var_frac"] = df.apply(
        lambda r: mof_var.get(r["cohort"], {}).get(int(r["mofa_best_factor"]), np.nan)
        if pd.notna(r["mofa_best_factor"]) else np.nan, axis=1)

    # Bonferroni-survival flags
    df["gizmo_survives_bonf_cohort"] = (df["gizmo_raw_p"] < df["bonf_cohort_thr"]) & (df["gizmo_strength"] >= 0.40)
    df["gizmo_survives_bonf_global"] = (df["gizmo_raw_p"] < df["bonf_global_thr"]) & (df["gizmo_strength"] >= 0.40)
    df["mofa_survives_bonf_cohort"] = (df["mofa_raw_p"] < df["bonf_cohort_thr"]) & (df["mofa_strength"] >= 0.40)
    df["mofa_survives_bonf_global"] = (df["mofa_raw_p"] < df["bonf_global_thr"]) & (df["mofa_strength"] >= 0.40)

    out = RESULTS / "multi_pc_vs_mofa_factors_augmented.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")

    # Per-cohort search-space disclosure
    disclosure = []
    for c in cohorts:
        sub = df[df["cohort"] == c]
        n_t = len(sub)
        bonf_c = 0.05 / n_t
        # Pick the single winning row per method (max strength)
        if not sub.empty:
            giz_winner = sub.loc[sub["gizmo_strength"].idxmax()]
            mof_winner = sub.loc[sub["mofa_strength"].idxmax()]
            disclosure.append({
                "cohort": c, "n_cohort_tests": n_t,
                "bonf_cohort_thr": bonf_c, "bonf_global_thr": bonf_global,
                "gizmo_best_pc": int(giz_winner["gizmo_best_pc"])
                    if pd.notna(giz_winner["gizmo_best_pc"]) else None,
                "gizmo_best_metadata": giz_winner["metadata"],
                "gizmo_strength": giz_winner["gizmo_strength"],
                "gizmo_raw_p": giz_winner["gizmo_raw_p"],
                "gizmo_pc_var_frac": giz_winner["gizmo_pc_var_frac"],
                "gizmo_survives_bonf_cohort": bool(giz_winner["gizmo_survives_bonf_cohort"]),
                "gizmo_survives_bonf_global": bool(giz_winner["gizmo_survives_bonf_global"]),
                "mofa_best_factor": int(mof_winner["mofa_best_factor"])
                    if pd.notna(mof_winner["mofa_best_factor"]) else None,
                "mofa_best_metadata": mof_winner["metadata"],
                "mofa_strength": mof_winner["mofa_strength"],
                "mofa_raw_p": mof_winner["mofa_raw_p"],
                "mofa_factor_var_frac": mof_winner["mofa_factor_var_frac"],
                "mofa_survives_bonf_cohort": bool(mof_winner["mofa_survives_bonf_cohort"]),
                "mofa_survives_bonf_global": bool(mof_winner["mofa_survives_bonf_global"]),
            })
    df_disc = pd.DataFrame(disclosure)
    out2 = RESULTS / "multi_pc_winners_disclosure.tsv"
    df_disc.to_csv(out2, sep="\t", index=False)
    print(f"Wrote {out2}")

    print(f"\n=== Per-cohort search-space disclosure ===")
    print(f"Global Bonferroni threshold: α = 0.05/{n_total} = {bonf_global:.2e}\n")
    cols = ["cohort", "n_cohort_tests", "gizmo_best_pc", "gizmo_pc_var_frac",
            "gizmo_strength", "gizmo_raw_p", "gizmo_survives_bonf_cohort",
            "gizmo_survives_bonf_global", "mofa_best_factor", "mofa_factor_var_frac",
            "mofa_strength", "mofa_raw_p", "mofa_survives_bonf_global"]
    print(df_disc[cols].to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
