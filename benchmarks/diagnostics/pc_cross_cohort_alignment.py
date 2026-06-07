"""Are α-PC directions actually similar across cohorts (38,148-D unit vectors)?

For each pair of cohorts (c1, c2), compute |cosine(PC_k of c1, PC_k of c2)|
across PC ranks 1..5. High cosine = same substrate direction across cohorts.

Output: benchmarks/results/pc_cross_cohort_alignment.tsv + heatmap PNG.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

REPO = Path("/home/jgardner/GIZMO")
RESULTS = REPO / "benchmarks/results"
UR = RESULTS / "unsupervised"
FIG_DIR = RESULTS / "figures"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))


def load_substrate_log_pr():
    import networkx as nx
    from gizmo.export.json_export import read_json
    from per_patient_wlsp_v2 import biochem_subgraph
    print("Loading substrate…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    sub_dir, nodes, nid_idx = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)
    return log_pr, nodes


def find_F_path(cohort):
    for cand in [UR / f"stage3_F_{cohort}.npz",
                 UR / f"stage3_F_{cohort}_combined.npz",
                 UR / f"stage3_F_{cohort}_edge_informed.npz",
                 UR / f"stage3_F_{cohort}_node_informed.npz"]:
        if cand.exists():
            return cand
    return None


def decompose_unit_norm(F, log_pr):
    F_norm = np.linalg.norm(F, axis=1, keepdims=True) + 1e-12
    F = F / F_norm
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F.mean(axis=1, keepdims=True)
    cov = ((F - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F - F_mean - beta[:, None] * (x - x_mean)[None, :]
    return beta, alpha


def main():
    log_pr, nodes = load_substrate_log_pr()
    cohorts = ["CPTAC_CCRCC", "CPTAC_COAD", "CPTAC_OV",
               "GSE65391_SLE", "GSE65682_sepsis",
               "Crohn", "Su_COVID", "Erawijantari", "Filbin_COVID",
               "GSE89408_RA", "Gao_RA", "HMP2_IBD_CD",
               "IDH_glioma", "TCGA_IDH_glioma", "KMPLOT_BRCA", "TCGA_LUAD"]
    n_pcs = 5

    # Extract per-cohort PC directions
    cohort_pcs = {}  # cohort -> (n_pcs, n_nodes) array
    for cohort in cohorts:
        F_path = find_F_path(cohort)
        if F_path is None:
            continue
        fd = np.load(F_path, allow_pickle=True)
        F = fd["F"].astype(np.float64)
        if F.shape[1] != len(log_pr):
            continue
        _, alpha = decompose_unit_norm(F, log_pr)
        try:
            pca = PCA(n_components=min(n_pcs, alpha.shape[0] - 1), random_state=0)
            pca.fit(alpha)
            cohort_pcs[cohort] = pca.components_  # (n_pcs, 38148)
            print(f"  {cohort}: {pca.components_.shape[0]} PCs extracted", flush=True)
        except Exception as e:
            print(f"  {cohort}: PCA failed: {e}", flush=True)

    cohort_list = sorted(cohort_pcs.keys())
    n_cohorts = len(cohort_list)
    print(f"\nComputing pairwise PC alignment across {n_cohorts} cohorts × {n_pcs} PCs…", flush=True)

    # For each PC rank k, compute n_cohorts × n_cohorts matrix of |cosine|
    # Also: for each pair, find best-matching PC rank in the other cohort
    rows = []
    for k in range(n_pcs):
        mat = np.eye(n_cohorts)
        for i, c1 in enumerate(cohort_list):
            if k >= cohort_pcs[c1].shape[0]: continue
            v1 = cohort_pcs[c1][k]
            v1 = v1 / (np.linalg.norm(v1) + 1e-12)
            for j, c2 in enumerate(cohort_list):
                if i == j: continue
                if k >= cohort_pcs[c2].shape[0]: continue
                v2 = cohort_pcs[c2][k]
                v2 = v2 / (np.linalg.norm(v2) + 1e-12)
                cos_val = float(np.abs(v1 @ v2))
                mat[i, j] = cos_val
                rows.append({"cohort_1": c1, "cohort_2": c2, "pc_rank": k + 1,
                              "cosine": cos_val})
        # Plot per-PC heatmap
        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(mat, cmap="viridis", aspect="auto", vmin=0, vmax=1)
        ax.set_xticks(range(n_cohorts))
        ax.set_xticklabels(cohort_list, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(n_cohorts))
        ax.set_yticklabels(cohort_list, fontsize=7)
        for i in range(n_cohorts):
            for j in range(n_cohorts):
                if i == j: continue
                v = mat[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6, color="white" if v > 0.5 else "black")
        plt.colorbar(im, ax=ax, label="|cosine|")
        ax.set_title(f"Pairwise α-PC{k+1} direction alignment across cohorts\n"
                     f"(|cosine| of 38,148-D unit vectors; 1=identical direction, 0=orthogonal)",
                     fontsize=11)
        plt.tight_layout()
        plt.savefig(FIG_DIR / f"fig_pc_alignment_pc{k+1}.png", dpi=120, bbox_inches="tight")
        plt.close()
        print(f"  PC{k+1}: max off-diagonal |cosine| = {mat[mat<0.99].max():.3f}; "
              f"median off-diagonal = {np.median(mat[mat<0.99]):.3f}", flush=True)

    # Also: best-match PC across ranks
    # For each pair (c1, c2), find best |cosine| between any PC of c1 and any PC of c2
    best_match_rows = []
    for c1 in cohort_list:
        for c2 in cohort_list:
            if c1 == c2: continue
            pcs1 = cohort_pcs[c1]; pcs2 = cohort_pcs[c2]
            for k1 in range(pcs1.shape[0]):
                for k2 in range(pcs2.shape[0]):
                    v1 = pcs1[k1] / (np.linalg.norm(pcs1[k1]) + 1e-12)
                    v2 = pcs2[k2] / (np.linalg.norm(pcs2[k2]) + 1e-12)
                    best_match_rows.append({"c1": c1, "c2": c2,
                                              "pc_c1": k1+1, "pc_c2": k2+1,
                                              "cosine": float(np.abs(v1 @ v2))})

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "pc_cross_cohort_alignment.tsv", sep="\t", index=False)
    bm_df = pd.DataFrame(best_match_rows)
    bm_df.to_csv(RESULTS / "pc_cross_cohort_alignment_all_pairs.tsv", sep="\t", index=False)

    # Print best matches: for each cohort pair, find top match
    print(f"\nTop cross-cohort PC matches (any PC × any PC, cosine > 0.5):", flush=True)
    high_matches = bm_df[bm_df["cosine"] > 0.5].sort_values("cosine", ascending=False)
    print(f"  Total pairs with cosine > 0.5: {len(high_matches)}", flush=True)
    print(high_matches.head(30).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
