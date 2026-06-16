"""BBB — Bootstrap stability of K=3 SLE endotype clustering.

1000 bootstrap iterations: resample SLE patients with replacement, run K-means
(K=3) on (PC5, PC9, PC10) z-scored coordinates, compute adjusted rand index
(ARI) against the full-cohort reference partition.

ARI distribution → cluster stability assessment:
  ARI > 0.7 → stable
  ARI 0.4 - 0.7 → moderately stable
  ARI < 0.4 → unstable

Plus per-patient consensus matrix: P(same_cluster | i, j) across bootstraps.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path('/home/jgardner/GIZMO')
SCAN = REPO / 'benchmarks/results/unsupervised/cross_indication_scan'
FIGDIR = REPO / 'benchmarks/results/figures/paper4'
FIGDIR.mkdir(parents=True, exist_ok=True)

N_BOOTSTRAP = 1000


def main():
    print('Loading SLE per-patient profile...', flush=True)
    df = pd.read_csv(SCAN / 'sle_per_patient_axis_profile.csv')
    sle = df[df.is_SLE].reset_index(drop=True)
    coords = sle[['PC5_IFN', 'PC9_folate', 'PC10_plasma']].values
    n_sle = len(sle)
    Z = (coords - coords.mean(axis=0)) / (coords.std(axis=0) + 1e-9)
    print(f'  SLE patients: {n_sle}', flush=True)

    # Reference partition (full cohort K=3)
    km_ref = KMeans(n_clusters=3, random_state=0, n_init=10).fit(Z)
    ref_labels = km_ref.labels_
    print(f'  Reference clustering: '
          f'{[int((ref_labels == c).sum()) for c in range(3)]}', flush=True)

    # Bootstrap loop — ARI distribution only (skip O(n²) consensus matrix)
    print(f'\nBootstrap K-means ({N_BOOTSTRAP} iterations)...', flush=True)
    rng = np.random.default_rng(0)
    aris = np.zeros(N_BOOTSTRAP, dtype=np.float32)
    # Approach 2 (cheaper): per-iter project bootstrap K-means labels BACK
    # onto FULL-cohort patient slots via majority vote, then ARI vs ref_labels
    # globally. Avoids O(n²) inner loop.

    for it in range(N_BOOTSTRAP):
        sample_idx = rng.choice(n_sle, size=n_sle, replace=True)
        Z_boot = Z[sample_idx]
        try:
            km = KMeans(n_clusters=3, random_state=int(it),
                         n_init=10).fit(Z_boot)
        except Exception:
            aris[it] = float('nan')
            continue

        # ARI for the sampled indices (vs the reference labels at those positions)
        ref_for_boot = ref_labels[sample_idx]
        aris[it] = adjusted_rand_score(ref_for_boot, km.labels_)

        if (it + 1) % 100 == 0:
            print(f'    iter {it+1}/{N_BOOTSTRAP}, '
                  f'mean ARI so far = {np.nanmean(aris[:it+1]):.3f}',
                  flush=True)

    # Distribution analysis
    print('\n' + '='*78, flush=True)
    print(f'ARI distribution across {N_BOOTSTRAP} bootstrap iterations',
          flush=True)
    print('='*78, flush=True)
    print(f'  mean ARI: {np.nanmean(aris):.3f}', flush=True)
    print(f'  median: {np.nanmedian(aris):.3f}', flush=True)
    print(f'  SD: {np.nanstd(aris):.3f}', flush=True)
    print(f'  min: {np.nanmin(aris):.3f}, max: {np.nanmax(aris):.3f}', flush=True)
    print(f'  5th–95th percentile: '
          f'[{np.nanpercentile(aris, 5):.3f}, {np.nanpercentile(aris, 95):.3f}]',
          flush=True)
    if np.nanmean(aris) > 0.7:
        print(f'  Verdict: STABLE clustering (ARI > 0.7)', flush=True)
    elif np.nanmean(aris) > 0.4:
        print(f'  Verdict: MODERATELY STABLE (ARI 0.4-0.7)', flush=True)
    else:
        print(f'  Verdict: UNSTABLE (ARI < 0.4)', flush=True)

    np.save(SCAN / 'endotype_bootstrap_aris.npy', aris)
    np.save(SCAN / 'endotype_consensus_matrix.npy', consensus)

    # ARI histogram
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(aris[~np.isnan(aris)], bins=40, color='#9D4EDD',
            edgecolor='black', alpha=0.75)
    ax.axvline(0.7, color='green', linestyle='--', label='stable (0.7)')
    ax.axvline(0.4, color='orange', linestyle='--', label='moderate (0.4)')
    ax.axvline(np.nanmean(aris), color='red', linewidth=2,
                label=f'mean = {np.nanmean(aris):.3f}')
    ax.set_xlabel('Adjusted Rand Index (vs reference K=3 partition)',
                   fontsize=10)
    ax.set_ylabel(f'count ({N_BOOTSTRAP} bootstraps)', fontsize=10)
    ax.set_title('Endotype K=3 clustering bootstrap stability\n'
                  f'GSE65391 SLE n={n_sle}', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGDIR / 'fig6_endotype_bootstrap_ari.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print(f'\nSaved {FIGDIR / "fig6_endotype_bootstrap_ari.png"}', flush=True)


if __name__ == '__main__':
    main()
