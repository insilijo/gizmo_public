"""Symmetric multi-axis test — substrate-matched MOFA+ variant.

Same test as multi_pc_vs_mofa_factors.py, but MOFA+ uses the substrate-restricted
weights (only features mappable to substrate, the same set GIZMO sees).
Loads from unsupervised/mofa_weights_substrate_matched/.
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
MOFA_DIR = UR / "mofa_weights_substrate_matched"
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
    F_path = None
    for suffix in ("", "_edge_informed", "_combined", "_node_informed"):
        cand = UR / f"stage3_F_{cohort}{suffix}.npz"
        if cand.exists():
            F_path = cand; break
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
    """Load substrate-matched MOFA+ factor scores. Tries both file layouts:
      - mofa_weights_substrate_matched/mofa_weights_<cohort>_substrate.json
        (original 6-cohort sm pipeline)
      - mofa_weights/mofa_weights_<cohort>_sm.json
        (extended streaming-MOFA+ runner, 5 additional cohorts)
    """
    candidates = [
        MOFA_DIR / f"mofa_weights_{cohort}_substrate.json",
        UR / "mofa_weights" / f"mofa_weights_{cohort}_sm.json",
    ]
    w_file = next((p for p in candidates if p.exists()), None)
    if w_file is None: return None, None
    d = json.load(open(w_file))
    samples = d.get("samples", [])
    scores = d.get("factor_scores")
    if not samples or scores is None: return None, None
    pids = [str(s).lower() for s in samples]
    return pids, np.array(scores)


def compute_strength_and_p(scores_col, vals, metric):
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
            grp0 = s_sub[v_bin == 0]; grp1 = s_sub[v_bin == 1]
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


def main():
    print("Loading substrate hub direction…", flush=True)
    log_pr = get_log_pr()

    from benchmarks.diagnostics.axis_metadata_extended import (
        load_kmplot_metadata, load_cptac_metadata, load_tcga_idh_glioma_metadata,
        load_idh_glioma_trautwein_metadata, load_gao_ra_metadata,
        load_crohn_metadata, load_filbin_metadata, load_erawijantari_metadata,
        load_hmp2_metadata, load_tcga_luad_metadata, load_su_covid_metadata,
        load_gse_series_metadata,
    )
    REPO_LOCAL = Path("/home/jgardner/GIZMO")

    # Cohorts where substrate-matched MOFA+ ran. Original 6 + 5 new ones
    # added via the streaming runner (CPTAC trio via subsample MOFA+_sm;
    # SLE / sepsis via IncrementalPCA on substrate-mappable inputs).
    loader_map = {
        "Crohn": load_crohn_metadata,
        "Gao_RA": load_gao_ra_metadata,
        "Su_COVID": load_su_covid_metadata,
        "Erawijantari": load_erawijantari_metadata,
        "IDH_glioma": load_idh_glioma_trautwein_metadata,
        "Filbin_COVID": load_filbin_metadata,
        "CPTAC_CCRCC": lambda: load_cptac_metadata("CPTAC_CCRCC"),
        "CPTAC_COAD":  lambda: load_cptac_metadata("CPTAC_COAD"),
        "CPTAC_OV":    lambda: load_cptac_metadata("CPTAC_OV"),
        "GSE65391_SLE": lambda: load_gse_series_metadata(
            REPO_LOCAL / "data/cohorts/GSE65391_SLE/GSE65391_series_matrix.txt.gz"),
        "GSE65682_sepsis": lambda: load_gse_series_metadata(
            REPO_LOCAL / "data/cohorts/GSE65682_sepsis/GSE65682_series_matrix.txt.gz"),
    }

    rows = []
    for cohort in loader_map.keys():
        print(f"\n=== {cohort} ===", flush=True)
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
        if cohort.startswith("CPTAC"):
            import re
            md["pid"] = md["pid"].apply(lambda p: re.sub(r"_[tn]$", "", str(p)))

        giz_pids, giz_scores = gizmo_pc_scores(cohort, log_pr)
        if giz_pids is None:
            print(f"  no GIZMO F", flush=True); continue
        mof_pids, mof_scores = mofa_factor_scores(cohort)
        if mof_pids is None:
            print(f"  no substrate-matched MOFA factors", flush=True); continue
        if cohort.startswith("CPTAC"):
            import re
            giz_pids = [re.sub(r"_[tn]$", "", p) for p in giz_pids]
            mof_pids = [re.sub(r"_[tn]$", "", p) for p in mof_pids]
        print(f"  GIZMO {len(giz_pids)} pids × {giz_scores.shape[1]} PCs; "
              f"MOFA(substrate) {len(mof_pids)} samples × {mof_scores.shape[1]} factors",
              flush=True)

        skip_cols = {id_col, "pid"}
        n_for_cohort = 0
        for col in md.columns:
            if col in skip_cols: continue
            vals = pd.to_numeric(md[col], errors="coerce")
            if vals.isna().sum() == len(vals): continue
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
            if g_strength is None and m_strength is None: continue
            rows.append({
                "cohort": cohort, "metadata": col, "metric": metric,
                "n_meta": int((~vals.isna()).sum()),
                "gizmo_best_pc": (g_axis + 1) if g_axis is not None else None,
                "gizmo_strength": g_strength, "gizmo_raw_p": g_p,
                "mofa_best_factor": (m_axis + 1) if m_axis is not None else None,
                "mofa_strength": m_strength, "mofa_raw_p": m_p,
            })
            n_for_cohort += 1
        print(f"  {n_for_cohort} metadata fields tested", flush=True)

    df_out = pd.DataFrame(rows)
    df_out.to_csv(RESULTS / "multi_pc_vs_mofa_substrate_matched.tsv", sep="\t", index=False)
    print(f"\nWrote {RESULTS / 'multi_pc_vs_mofa_substrate_matched.tsv'}", flush=True)

    n_tests_total = len(df_out)
    bonferroni_global = 0.05 / n_tests_total if n_tests_total else 1.0
    n_per_cohort = df_out.groupby("cohort").size()

    print(f"\nTotal tests: {n_tests_total}, Bonferroni global: {bonferroni_global:.2e}")

    summary_rows = []
    for cohort in df_out["cohort"].unique():
        c = df_out[df_out["cohort"] == cohort]
        n_t = len(c); bonf_c = 0.05 / n_t
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
    summary.to_csv(RESULTS / "multi_pc_vs_mofa_substrate_matched_pass_rates.tsv",
                   sep="\t", index=False)
    print(f"\n== Per-cohort pass summary (substrate-matched MOFA+) ==")
    print(summary.to_string(index=False))
    n_c = len(summary)
    print(f"\nAggregate pass rates ({n_c} cohorts):")
    for col in ["GIZMO_pass_raw", "GIZMO_pass_bonf_cohort", "GIZMO_pass_bonf_global",
                "MOFA_pass_raw", "MOFA_pass_bonf_cohort", "MOFA_pass_bonf_global"]:
        print(f"  {col}: {int(summary[col].sum())}/{n_c} ({100*summary[col].mean():.0f}%)")


if __name__ == "__main__":
    main()
