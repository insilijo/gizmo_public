"""β ablation: PCA directly on F vs F-α-PC (β removed).

If PCA-on-F captures the same α-PC structure, β is unnecessary.
If β-removal changes the principal directions, β provides interpretive separation.
"""
import sys, json
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, '/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS = REPO / 'benchmarks/results/unsupervised/zscored'

mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr

cohorts = ['TCGA_IDH_glioma', 'KMPLOT_BRCA', 'TCGA_LUAD', 'Filbin_COVID']
print(f'{"Cohort":<20} | F-α-PC explanation | Direct-PCA-on-F explanation | Top-7 cos similarity')
print('-' * 100)

for c in cohorts:
    fp = SNAPSHOT / f'stage3_F_{c}_edge_informed.npz'
    if not fp.exists(): continue
    F = np.load(fp, allow_pickle=True)['F']
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)

    # F-α-PC pipeline (β/α decomposition + PCA on α)
    x = log_pr; xm = x.mean(); xv = x.var() + 1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm) * (x - xm)).mean(axis=1, keepdims=True)
    beta = (cov / xv).ravel()
    alpha = F_unit - Fm - beta[:, None] * (x - xm)[None, :]
    pca_a = PCA(n_components=7, random_state=0); pca_a.fit(alpha)

    # Direct PCA on F (no β removal)
    pca_f = PCA(n_components=7, random_state=0); pca_f.fit(F_unit)

    # Compute pairwise cosine between F-α components and PCA-on-F components
    sims = np.zeros((7, 7))
    for i in range(7):
        for j in range(7):
            sims[i, j] = abs(np.dot(pca_a.components_[i], pca_f.components_[j]))
    # Best-match cosine per α-PC
    best_match = sims.max(axis=1)

    # Variance fraction comparison
    f_pcs_evr = pca_f.explained_variance_ratio_
    a_pcs_evr = pca_a.explained_variance_ratio_

    print(f'{c:<20} | α-EVR top-7: {a_pcs_evr.sum()*100:.1f}% | F-EVR top-7: {f_pcs_evr.sum()*100:.1f}% | '
          f'best-match cos: {best_match.mean():.3f} ± {best_match.std():.3f}')
    print(f'  α-PC EVR per PC: {(a_pcs_evr*100).round(1).tolist()}')
    print(f'  F-PC EVR per PC: {(f_pcs_evr*100).round(1).tolist()}')
    print(f'  best-match cos per α-PC: {best_match.round(3).tolist()}')
    print()

# Also: do α-PC and PCA-on-F give the same prognostic discrimination?
# This is the actual "is β necessary?" test
print('\n=== Does β/α split materially change discrimination? ===')
print('Reference: TCGA_IDH F-α-PC4 = 0.762; TCGA_LUAD F-α-PC6 = 0.610 under v5-canonical')
print('TBD: compute matching PCA-on-F best-PC C-index for direct comparison')
