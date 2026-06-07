"""Symmetric multi-axis test — both GIZMO and MOFA+ get to use all their axes.

For each cohort × metadata field, compute the discrimination AUC (or |ρ|) of:
  - α-PC1..PC5 (GIZMO, 5 axes)
  - MOFA factor scores 1..K (MOFA+, K factors — typically 5-15)

Pick the best axis per method per (cohort, metadata). Compare pass rates at:
  raw strength ≥ 0.40, per-cohort Bonferroni (α=0.05/n_tests_cohort),
  full-grid Bonferroni (α=0.05/n_tests_total).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

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
    sub_dir, nodes, nid_idx = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    return np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)


def gizmo_pc_scores(cohort, log_pr):
    """Returns (patient_ids, scores_matrix [n_pat × n_pcs])."""
    # Some cohorts are saved with a suffix (_edge_informed / _combined /
    # _node_informed) reflecting which MAP variant was used. Pick the first
    # available, in priority order matching the rest of the codebase.
    F_path = None
    for suffix in ("", "_edge_informed", "_combined", "_node_informed"):
        cand = UR / f"stage3_F_{cohort}{suffix}.npz"
        if cand.exists():
            F_path = cand
            break
    if F_path is None: return None, None
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"].astype(np.float64)
    pids = [str(p).lower() for p in fd["patient_ids"]]
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    pca = PCA(n_components=5, random_state=0)
    return pids, pca.fit_transform(alpha)


def mofa_factor_scores(cohort):
    """Returns (sample_ids_lowered, factor_scores [n_samples × K]) or (None, None)."""
    w_file = UR / "mofa_weights" / f"mofa_weights_{cohort}.json"
    if not w_file.exists(): return None, None
    d = json.load(open(w_file))
    samples = d.get("samples", [])
    scores = d.get("factor_scores")
    if not samples or scores is None: return None, None
    pids = [str(s).lower() for s in samples]
    return pids, np.array(scores)


def compute_strength_and_p(scores_col, vals, metric):
    """Compute (best_strength, raw_p_value) where strength = 2|AUC-0.5| or |ρ|.
    metric='auc' for binary; 'spearman' for continuous."""
    mask = ~pd.isna(vals)
    if mask.sum() < 5:
        return np.nan, np.nan
    s_sub = scores_col[mask]; v_sub = vals[mask]
    if metric == "auc":
        v_bin = np.asarray(v_sub).astype(int)
        if len(np.unique(v_bin)) < 2 or (v_bin == 0).sum() < 3 or (v_bin == 1).sum() < 3:
            return np.nan, np.nan
        try:
            auc = roc_auc_score(v_bin, s_sub)
            auc_strength = 2 * abs(auc - 0.5)
            grp0 = s_sub[v_bin == 0]
            grp1 = s_sub[v_bin == 1]
            _, p = mannwhitneyu(grp0, grp1, alternative="two-sided")
            return float(auc_strength), float(p)
        except Exception:
            return np.nan, np.nan
    if metric == "spearman":
        try:
            r, p = spearmanr(s_sub, v_sub)
            return float(abs(r)), float(p)
        except Exception:
            return np.nan, np.nan
    return np.nan, np.nan


def best_axis_for_metadata(scores_matrix, pids_scores, vals_dict, metric, pids_meta):
    """For each axis (column of scores), compute strength against vals.
    Returns (best_axis_idx, best_strength, best_raw_p)."""
    pid_to_row = {p: i for i, p in enumerate(pids_scores)}
    aligned_idx = [pid_to_row.get(p) for p in pids_meta]
    valid = [i for i, idx in enumerate(aligned_idx) if idx is not None]
    if len(valid) < 5: return None, None, None
    valid_meta_pids = [pids_meta[i] for i in valid]
    valid_scores_rows = np.array([aligned_idx[i] for i in valid])
    vals = np.array([vals_dict[p] for p in valid_meta_pids])
    n_axes = scores_matrix.shape[1]
    best_axis = None; best_strength = -np.inf; best_p = np.nan
    for k in range(n_axes):
        col = scores_matrix[valid_scores_rows, k]
        strength, p = compute_strength_and_p(col, vals, metric)
        if not np.isnan(strength) and strength > best_strength:
            best_strength = strength; best_axis = k; best_p = p
    if best_axis is None: return None, None, None
    return best_axis, float(best_strength), float(best_p)


# (combined multi-axis ridge/logistic removed — best-single-axis is the
# apples-to-apples symmetric test the hostile reviewer asked for)


def main():
    print("Loading substrate hub direction…", flush=True)
    log_pr = get_log_pr()

    # Load full discrimination atlas
    df = pd.read_csv(RESULTS / "axis_metadata_discrimination_extended.tsv", sep="\t")
    # Get per-cohort metadata field list + per-patient values
    # We need to re-derive per-patient values — but the discrimination_extended
    # was computed from labels. For multi-axis testing we need to load per-cohort
    # metadata from the original cohort loaders.
    # For now, work from the discrimination_extended TSV which has the labels
    # we used to compute it for α-PC1..PC5. We can extend MOFA+ by computing
    # factor-scores × metadata directly if we have the cohort metadata loaders.
    #
    # Simpler approach: use the existing α-PC × metadata × strength table for
    # GIZMO. For MOFA+, compute per-factor strength on the SAME metadata set
    # by loading cohort metadata via existing loaders.
    from benchmarks.diagnostics.axis_metadata_extended import (
        load_kmplot_metadata, load_cptac_metadata, load_tcga_idh_glioma_metadata,
        load_idh_glioma_trautwein_metadata, load_gao_ra_metadata,
        load_crohn_metadata, load_filbin_metadata, load_erawijantari_metadata,
        load_hmp2_metadata, load_tcga_luad_metadata, load_su_covid_metadata,
        load_gse_series_metadata,
    )
    # Importing axis_metadata_extended overwrites the module-level REPO with
    # its own Path object. Re-anchor for our own use below.
    REPO_LOCAL = Path("/home/jgardner/GIZMO")

    REPO_PATH = REPO_LOCAL
    loader_map = {
        "KMPLOT_BRCA": load_kmplot_metadata,
        "CPTAC_CCRCC": lambda: load_cptac_metadata("CPTAC_CCRCC"),
        "CPTAC_COAD": lambda: load_cptac_metadata("CPTAC_COAD"),
        "CPTAC_OV": lambda: load_cptac_metadata("CPTAC_OV"),
        "TCGA_IDH_glioma": load_tcga_idh_glioma_metadata,
        "IDH_glioma": load_idh_glioma_trautwein_metadata,
        "Gao_RA": load_gao_ra_metadata,
        "Crohn": load_crohn_metadata,
        "Filbin_COVID": load_filbin_metadata,
        "Erawijantari": load_erawijantari_metadata,
        "HMP2_IBD_CD": load_hmp2_metadata,
        "TCGA_LUAD": load_tcga_luad_metadata,
        "Su_COVID": load_su_covid_metadata,
        "GSE65391_SLE": lambda: load_gse_series_metadata(
            REPO_PATH / "data/cohorts/GSE65391_SLE/GSE65391_series_matrix.txt.gz"),
        "GSE65682_sepsis": lambda: load_gse_series_metadata(
            REPO_PATH / "data/cohorts/GSE65682_sepsis/GSE65682_series_matrix.txt.gz"),
    }

    cohorts = [c for c in loader_map.keys()]
    rows = []
    for cohort in cohorts:
        print(f"\n=== {cohort} ===", flush=True)
        # Load metadata
        try:
            md = loader_map[cohort]()
        except Exception as e:
            print(f"  metadata loader failed: {e}", flush=True); continue
        if md is None or md.empty:
            print(f"  no metadata", flush=True); continue
        id_col = None
        for c in ("patient_id", "sample_id", "id"):
            if c in md.columns: id_col = c; break
        if id_col is None: id_col = md.columns[0]
        md["pid"] = md[id_col].astype(str).str.lower()
        # Strip CPTAC _T/_N
        if cohort.startswith("CPTAC"):
            import re
            md["pid"] = md["pid"].apply(lambda p: re.sub(r"_[tn]$", "", str(p)))
        # GIZMO scores
        giz_pids, giz_scores = gizmo_pc_scores(cohort, log_pr)
        if giz_pids is None:
            print(f"  no GIZMO F", flush=True); continue
        # MOFA factor scores
        mof_pids, mof_scores = mofa_factor_scores(cohort)
        if cohort.startswith("CPTAC"):
            # CPTAC outputs preserve the per-sample _T/_N suffix; the
            # metadata table is per-patient. Strip both score-source ids
            # to the patient root so the metadata join lands.
            import re
            if giz_pids:
                giz_pids = [re.sub(r"_[tn]$", "", p) for p in giz_pids]
            if mof_pids:
                mof_pids = [re.sub(r"_[tn]$", "", p) for p in mof_pids]
        if mof_pids is None:
            print(f"  no MOFA factors", flush=True); continue

        # For each metadata field, find best GIZMO PC + best MOFA factor
        skip_cols = {id_col, "pid"}
        for col in md.columns:
            if col in skip_cols: continue
            vals = pd.to_numeric(md[col], errors="coerce")
            if vals.isna().sum() == len(vals): continue
            # Determine metric: binary if 2 unique values; else spearman
            unique_vals = vals.dropna().unique()
            if len(unique_vals) == 2 and set(unique_vals).issubset({0, 1}):
                metric = "auc"
            else:
                metric = "spearman"
            vals_dict = dict(zip(md["pid"], vals))
            pids_meta = list(md["pid"])

            g_axis, g_strength, g_p = best_axis_for_metadata(
                giz_scores, giz_pids, vals_dict, metric, pids_meta)
            m_axis, m_strength, m_p = best_axis_for_metadata(
                mof_scores, mof_pids, vals_dict, metric, pids_meta)
            if g_strength is None and m_strength is None:
                continue
            rows.append({
                "cohort": cohort, "metadata": col, "metric": metric,
                "n_meta": int((~vals.isna()).sum()),
                "gizmo_best_pc": (g_axis + 1) if g_axis is not None else None,
                "gizmo_strength": g_strength, "gizmo_raw_p": g_p,
                "mofa_best_factor": (m_axis + 1) if m_axis is not None else None,
                "mofa_strength": m_strength, "mofa_raw_p": m_p,
            })
        n_metadata_tested = len([r for r in rows if r["cohort"] == cohort])
        print(f"  {n_metadata_tested} metadata fields tested", flush=True)

    df_out = pd.DataFrame(rows)
    df_out.to_csv(RESULTS / "multi_pc_vs_mofa_factors.tsv", sep="\t", index=False)
    print(f"\nWrote {RESULTS / 'multi_pc_vs_mofa_factors.tsv'}", flush=True)

    # ---------- Pass rate analysis ----------
    n_tests_total = len(df_out)
    bonferroni_global = 0.05 / n_tests_total if n_tests_total else 1.0
    # Per-cohort Bonferroni
    n_per_cohort = df_out.groupby("cohort").size()
    df_out["bonferroni_cohort_thr"] = df_out["cohort"].map(lambda c: 0.05 / n_per_cohort[c])

    # Cohort-level pass: does cohort have AT LEAST ONE (metadata, axis) pair
    # passing the threshold for each method?
    print(f"\nTotal tests: {n_tests_total}, Bonferroni global: {bonferroni_global:.2e}")
    print(f"Per-cohort Bonferroni (median): {n_per_cohort.median():.0f} tests × {0.05/n_per_cohort.median():.2e}")

    summary_rows = []
    for cohort in df_out["cohort"].unique():
        c = df_out[df_out["cohort"] == cohort]
        n_t = len(c)
        # Pass: any test with strength >= 0.40 (raw); MW p < threshold
        bonf_c = 0.05 / n_t
        giz_pass_raw    = (c["gizmo_strength"] >= 0.40).any()
        giz_pass_bonfc  = ((c["gizmo_strength"] >= 0.40) & (c["gizmo_raw_p"] < bonf_c)).any()
        giz_pass_bonfg  = ((c["gizmo_strength"] >= 0.40) & (c["gizmo_raw_p"] < bonferroni_global)).any()
        mof_pass_raw    = (c["mofa_strength"] >= 0.40).any()
        mof_pass_bonfc  = ((c["mofa_strength"] >= 0.40) & (c["mofa_raw_p"] < bonf_c)).any()
        mof_pass_bonfg  = ((c["mofa_strength"] >= 0.40) & (c["mofa_raw_p"] < bonferroni_global)).any()
        summary_rows.append({
            "cohort": cohort, "n_tests": n_t,
            "GIZMO_pass_raw": giz_pass_raw, "GIZMO_pass_bonf_cohort": giz_pass_bonfc,
            "GIZMO_pass_bonf_global": giz_pass_bonfg,
            "MOFA_pass_raw": mof_pass_raw, "MOFA_pass_bonf_cohort": mof_pass_bonfc,
            "MOFA_pass_bonf_global": mof_pass_bonfg,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULTS / "multi_pc_vs_mofa_pass_rates.tsv", sep="\t", index=False)
    print(f"\n== Per-cohort pass summary ==\n{summary.to_string(index=False)}")
    n_c = len(summary)
    print(f"\nAggregate pass rates (cohorts where any metadata × axis passes):")
    for col in ["GIZMO_pass_raw", "GIZMO_pass_bonf_cohort", "GIZMO_pass_bonf_global",
                  "MOFA_pass_raw", "MOFA_pass_bonf_cohort", "MOFA_pass_bonf_global"]:
        print(f"  {col}: {int(summary[col].sum())}/{n_c} ({100*summary[col].mean():.0f}%)")


if __name__ == "__main__":
    main()
