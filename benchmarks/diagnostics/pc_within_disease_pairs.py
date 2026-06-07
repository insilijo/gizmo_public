"""Within-disease pair α-PC alignment: do same-disease cohorts converge on shared PCs?

Pairs:
  - Su_COVID vs Filbin_COVID (both COVID; different studies)
  - GSE89408_RA vs Gao_RA (both rheumatoid arthritis)
  - HMP2_IBD_CD vs Crohn (both IBD/CD-spectrum)
  - IDH_glioma vs TCGA_IDH_glioma (positive control — same disease, two cohorts)

For each pair, compute the 5×5 matrix of |cosine| between PC1..PC5 of c1 and
PC1..PC5 of c2 (38,148-D unit vectors), and the across-cohort distribution baseline.
"""
from __future__ import annotations

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
    return log_pr


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


def get_pcs(cohort, log_pr, n_pcs=5):
    F_path = find_F_path(cohort)
    if F_path is None:
        return None
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"].astype(np.float64)
    if F.shape[1] != len(log_pr):
        return None
    _, alpha = decompose_unit_norm(F, log_pr)
    n = min(n_pcs, alpha.shape[0] - 1)
    if n < 1:
        return None
    pca = PCA(n_components=n, random_state=0)
    pca.fit(alpha)
    return pca.components_


def pairwise_cosine(pcs1, pcs2):
    """Return (n1, n2) |cosine| matrix between row-vectors of pcs1 and pcs2."""
    n1, n2 = pcs1.shape[0], pcs2.shape[0]
    A = pcs1 / (np.linalg.norm(pcs1, axis=1, keepdims=True) + 1e-12)
    B = pcs2 / (np.linalg.norm(pcs2, axis=1, keepdims=True) + 1e-12)
    return np.abs(A @ B.T)


def main():
    log_pr = load_substrate_log_pr()
    n_pcs = 5

    pairs = [
        ("Su_COVID", "Filbin_COVID", "COVID"),
        ("GSE89408_RA", "Gao_RA", "RA"),
        ("HMP2_IBD_CD", "Crohn", "IBD"),
        ("IDH_glioma", "TCGA_IDH_glioma", "IDH-glioma (control)"),
    ]

    # All other cohorts to compute null
    all_cohorts = ["CPTAC_CCRCC", "CPTAC_COAD", "CPTAC_OV",
                   "GSE65391_SLE", "GSE65682_sepsis",
                   "Crohn", "Su_COVID", "Erawijantari", "Filbin_COVID",
                   "GSE89408_RA", "Gao_RA", "HMP2_IBD_CD",
                   "IDH_glioma", "TCGA_IDH_glioma", "KMPLOT_BRCA", "TCGA_LUAD"]

    print("Extracting PCs per cohort…", flush=True)
    all_pcs = {}
    for c in all_cohorts:
        p = get_pcs(c, log_pr, n_pcs=n_pcs)
        if p is not None:
            all_pcs[c] = p
            print(f"  {c}: {p.shape[0]} PCs", flush=True)

    # Within-disease pair analysis
    rows = []
    print("\n=== Within-disease-pair PC alignment ===", flush=True)
    for c1, c2, label in pairs:
        if c1 not in all_pcs or c2 not in all_pcs:
            print(f"  [skip] {label}: missing F for {c1 if c1 not in all_pcs else c2}", flush=True)
            continue
        M = pairwise_cosine(all_pcs[c1], all_pcs[c2])
        print(f"\n  {label} ({c1} vs {c2})", flush=True)
        print(f"  PC×PC |cosine| matrix:", flush=True)
        header = "          " + "  ".join([f"{c2}-PC{j+1}" for j in range(M.shape[1])])
        print("  " + header, flush=True)
        for i in range(M.shape[0]):
            row = f"  {c1}-PC{i+1}:  " + "  ".join([f"  {M[i,j]:.3f}   " for j in range(M.shape[1])])
            print(row, flush=True)
        diag = np.diag(M)
        max_per_row = M.max(axis=1)
        print(f"  Diagonal (PC1↔PC1, PC2↔PC2…): {diag.tolist()}", flush=True)
        print(f"  Best-match per c1-PC (max over c2-PCs): {max_per_row.tolist()}", flush=True)
        print(f"  Mean |cosine|: {M.mean():.3f}; Max |cosine|: {M.max():.3f}", flush=True)

        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                rows.append({"pair_label": label, "c1": c1, "c2": c2,
                             "pc1": i+1, "pc2": j+1,
                             "cosine": float(M[i, j]),
                             "kind": "within_disease"})

    # Null distribution: all cross-disease pairs
    print("\n=== Cross-disease (null) baseline ===", flush=True)
    within_disease_set = set()
    for c1, c2, _ in pairs:
        within_disease_set.add(tuple(sorted([c1, c2])))

    null_cosines = []
    cohort_list = sorted(all_pcs.keys())
    for i, c1 in enumerate(cohort_list):
        for c2 in cohort_list[i+1:]:
            key = tuple(sorted([c1, c2]))
            if key in within_disease_set:
                continue
            M = pairwise_cosine(all_pcs[c1], all_pcs[c2])
            for ii in range(M.shape[0]):
                for jj in range(M.shape[1]):
                    null_cosines.append(M[ii, jj])
                    rows.append({"pair_label": "cross_disease_null",
                                 "c1": c1, "c2": c2,
                                 "pc1": ii+1, "pc2": jj+1,
                                 "cosine": float(M[ii, jj]),
                                 "kind": "cross_disease"})
    null_arr = np.array(null_cosines)
    print(f"  N cross-disease cosines: {len(null_arr)}", flush=True)
    print(f"  Mean: {null_arr.mean():.3f}; "
          f"Median: {np.median(null_arr):.3f}; "
          f"95th percentile: {np.percentile(null_arr, 95):.3f}; "
          f"99th percentile: {np.percentile(null_arr, 99):.3f}; "
          f"Max: {null_arr.max():.3f}", flush=True)

    # Z-score within-disease pairs against null
    print("\n=== Z-scored within-disease alignments ===", flush=True)
    null_mu = null_arr.mean()
    null_sd = null_arr.std() + 1e-12
    for c1, c2, label in pairs:
        if c1 not in all_pcs or c2 not in all_pcs:
            continue
        M = pairwise_cosine(all_pcs[c1], all_pcs[c2])
        z_max = (M.max() - null_mu) / null_sd
        # Empirical p
        p_emp = float((null_arr >= M.max()).mean())
        print(f"  {label}: max cosine = {M.max():.3f}; "
              f"Z vs null = {z_max:.2f}; empirical p = {p_emp:.4f}", flush=True)

    # Save table
    df = pd.DataFrame(rows)
    out = RESULTS / "pc_within_disease_pairs.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}", flush=True)

    # Heatmap figure: 2×2 grid of within-disease pair matrices
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for ax, (c1, c2, label) in zip(axes.ravel(), pairs):
        if c1 not in all_pcs or c2 not in all_pcs:
            ax.set_title(f"{label} — missing"); ax.axis("off"); continue
        M = pairwise_cosine(all_pcs[c1], all_pcs[c2])
        im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(M.shape[1]))
        ax.set_xticklabels([f"PC{j+1}" for j in range(M.shape[1])])
        ax.set_yticks(range(M.shape[0]))
        ax.set_yticklabels([f"PC{i+1}" for i in range(M.shape[0])])
        ax.set_xlabel(c2); ax.set_ylabel(c1)
        ax.set_title(f"{label}\n(max |cos| = {M.max():.2f})", fontsize=10)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center",
                        fontsize=8, color="white" if M[i,j] > 0.5 else "black")
        plt.colorbar(im, ax=ax, fraction=0.04)
    fig.suptitle("Within-disease α-PC alignment (does the disease signature replicate across studies?)",
                 fontsize=12)
    plt.tight_layout()
    fig_path = FIG_DIR / "fig_pc_within_disease_pairs.png"
    plt.savefig(fig_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Wrote {fig_path}", flush=True)


if __name__ == "__main__":
    main()
