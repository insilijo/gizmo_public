"""Leave-one-out α-PC stability test on cohorts with narrow-mass PCs.

For each (cohort, winning PC), drop each patient, recompute β/α decomposition
+ α-PCA, measure |cos(PC_LOO, PC_original)|. Report bottom-5% LOO cosine and
fraction with |cos| < 0.5 (PC-direction-flip threshold).

Outputs:
  loo_pc_stability.tsv  — per (cohort, pc) row with stability stats
  loo_pc_stability_detail.tsv  — per (cohort, pc, dropped_patient) row
"""
from __future__ import annotations

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


def alpha_pca(F_unit, log_pr, n_components=5):
    """Same β/α decomposition + α-PCA as paper."""
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    pca = PCA(n_components=n_components, random_state=0)
    scores = pca.fit_transform(alpha)
    return pca.components_, pca.explained_variance_ratio_, scores


def narrow_mass_stats(loadings):
    """Fraction of |loadings|^2 on top-K nodes for K = 10, 50, 100."""
    sq = loadings ** 2
    total = sq.sum() + 1e-12
    sorted_sq = np.sort(sq)[::-1]
    return {
        "top10_mass": float(sorted_sq[:10].sum() / total),
        "top50_mass": float(sorted_sq[:50].sum() / total),
        "top100_mass": float(sorted_sq[:100].sum() / total),
    }


def cohort_loo(cohort_name, F_path, log_pr, target_pcs=(1, 2, 3, 4, 5)):
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"].astype(np.float64)
    pids = [str(p) for p in fd["patient_ids"]]
    n_pat = F.shape[0]
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)

    components_full, evr_full, _ = alpha_pca(F_unit, log_pr, n_components=max(target_pcs))
    full_mass = {k: narrow_mass_stats(components_full[k - 1]) for k in target_pcs}
    full_evr = {k: float(evr_full[k - 1]) for k in target_pcs}

    print(f"  full α-PC explained variance: " +
          ", ".join(f"PC{k}={full_evr[k]:.3f}" for k in target_pcs), flush=True)
    print(f"  full α-PC top-10 mass: " +
          ", ".join(f"PC{k}={full_mass[k]['top10_mass']:.3f}" for k in target_pcs), flush=True)

    # LOO
    cos_per_pc = {k: [] for k in target_pcs}
    for i in range(n_pat):
        idx = [j for j in range(n_pat) if j != i]
        F_sub = F_unit[idx]
        comps_sub, _, _ = alpha_pca(F_sub, log_pr, n_components=max(target_pcs))
        for k in target_pcs:
            # PC sign ambiguous — take absolute cosine
            v0 = components_full[k - 1]; v1 = comps_sub[k - 1]
            cos = float(np.abs(np.dot(v0, v1) / (np.linalg.norm(v0) * np.linalg.norm(v1) + 1e-12)))
            cos_per_pc[k].append({"dropped_patient": pids[i], "cos_to_original": cos})
        if (i + 1) % 50 == 0 or (i + 1) == n_pat:
            print(f"    {i+1}/{n_pat} patients dropped", flush=True)

    summary_rows = []
    detail_rows = []
    for k in target_pcs:
        cos_vals = np.array([d["cos_to_original"] for d in cos_per_pc[k]])
        summary_rows.append({
            "cohort": cohort_name, "pc": k, "n_patients": n_pat,
            "explained_variance_ratio": full_evr[k],
            "top10_mass": full_mass[k]["top10_mass"],
            "top50_mass": full_mass[k]["top50_mass"],
            "top100_mass": full_mass[k]["top100_mass"],
            "loo_cos_mean": float(cos_vals.mean()),
            "loo_cos_min": float(cos_vals.min()),
            "loo_cos_p05": float(np.percentile(cos_vals, 5)),
            "loo_cos_p25": float(np.percentile(cos_vals, 25)),
            "loo_cos_below_05": int((cos_vals < 0.5).sum()),
            "loo_cos_below_09": int((cos_vals < 0.9).sum()),
            "loo_cos_below_095": int((cos_vals < 0.95).sum()),
            "verdict": ("STABLE (p05≥0.9)" if np.percentile(cos_vals, 5) >= 0.9
                        else "MODERATE (p05≥0.7)" if np.percentile(cos_vals, 5) >= 0.7
                        else "UNSTABLE (p05<0.7)"),
        })
        for d in cos_per_pc[k]:
            detail_rows.append({"cohort": cohort_name, "pc": k,
                                "dropped_patient": d["dropped_patient"],
                                "cos_to_original": d["cos_to_original"]})
    return summary_rows, detail_rows


def main():
    print("Loading log_PR…", flush=True)
    log_pr = get_log_pr()

    # Cohorts where v4 manuscript surfaces a narrow-mass / sub-PC1 winning axis
    cohorts_to_test = [
        ("Filbin_COVID",        UR / "stage3_F_Filbin_COVID.npz"),
        ("TCGA_IDH_glioma",     UR / "stage3_F_TCGA_IDH_glioma.npz"),
        ("KMPLOT_BRCA",         UR / "stage3_F_KMPLOT_BRCA.npz"),
        ("CPTAC_CCRCC",         UR / "stage3_F_CPTAC_CCRCC.npz"),
        ("CPTAC_COAD",          UR / "stage3_F_CPTAC_COAD.npz"),
        ("CPTAC_OV",            UR / "stage3_F_CPTAC_OV.npz"),
        ("GSE65391_SLE",        UR / "stage3_F_GSE65391_SLE.npz"),
        ("Gao_RA",              UR / "stage3_F_Gao_RA.npz"),
        ("Su_COVID",            UR / "stage3_F_Su_COVID.npz"),
        ("Erawijantari",        UR / "stage3_F_Erawijantari.npz"),
    ]

    all_summary = []
    all_detail = []
    for name, path in cohorts_to_test:
        if not path.exists():
            print(f"\n=== {name}: SKIP (no F file) ===", flush=True); continue
        print(f"\n=== {name} ===", flush=True)
        try:
            s, d = cohort_loo(name, path, log_pr)
            all_summary.extend(s); all_detail.extend(d)
        except Exception as e:
            print(f"  ERROR: {e}", flush=True); continue

    df_s = pd.DataFrame(all_summary)
    df_d = pd.DataFrame(all_detail)
    df_s.to_csv(RESULTS / "loo_pc_stability.tsv", sep="\t", index=False)
    df_d.to_csv(RESULTS / "loo_pc_stability_detail.tsv", sep="\t", index=False)
    print(f"\nWrote {RESULTS / 'loo_pc_stability.tsv'}")
    print(f"Wrote {RESULTS / 'loo_pc_stability_detail.tsv'}")
    print(f"\n== Summary ==")
    print(df_s.to_string(index=False))


if __name__ == "__main__":
    main()
