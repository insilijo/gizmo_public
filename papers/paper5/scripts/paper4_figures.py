"""Paper 4 figure rendering — main figures for the manuscript.

Generates:
  Fig 2: GSE45291 4-group axis projections (HC / RA-DMARD-IR / RA-TNF-IR / SLE)
         — the strongest single replication result, p<10⁻¹²⁰
  Fig 3: SLE per-patient (PC5, PC9, PC10) 3-panel scatter w/ K=3 endotype overlay
  Fig 4: Novel-target validation forest plot — KDM8/PRMT1/CFI/GLO1 validated,
         PDIA3/CTSB failed; substrate-connectivity predictor
  Fig 5: NNLS prescription heatmap — patient × regime weights, ordered by SLEDAI

All figures saved as PNG (300 dpi) to benchmarks/results/figures/paper4/
"""
from __future__ import annotations
import sys, gc
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
CACHE = RESULTS / 'cohort_alpha_cache'
SCAN = RESULTS / 'cross_indication_scan'
FIGDIR = REPO / 'benchmarks/results/figures/paper4'
FIGDIR.mkdir(parents=True, exist_ok=True)

# Brand palette
HC_COLOR = '#3A86FF'          # blue
SLE_COLOR = '#FF006E'          # magenta
RA_DMARD_COLOR = '#FFBE0B'     # gold
RA_TNF_COLOR = '#FB5607'       # orange
SEVERE_COLOR = '#D62828'       # crimson
QUIESCENT_COLOR = '#06A77D'    # teal
MODERATE_COLOR = '#9D4EDD'     # purple


def fig2_gse45291_four_group():
    print('Fig 2: GSE45291 4-group projections...', flush=True)
    F = np.load(SCAN / 'GSE45291_F.npy')
    meta = pd.read_csv(SCAN / 'GSE45291_meta.csv')
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / 'benchmarks'))
    import loocv_longitudinal_audit as audit
    geom = audit.geom
    lpn = (geom.log_pr / (np.linalg.norm(geom.log_pr) + 1e-9)).astype(np.float32)
    a0 = np.load(CACHE / 'Wang_RA_alpha_t0.npy')
    pca = PCA(n_components=15, svd_solver='randomized', random_state=0).fit(a0)
    axes = [pca.components_[k] for k in [4, 8, 9]]
    names = ['PC5 IFN/STAT1', 'PC9 folate', 'PC10 plasma cell']
    beta = F @ lpn
    proj = {n: F @ v - beta * float(lpn @ v) for n, v in zip(names, axes)}

    dist = meta['disease'].astype(str).str.strip()
    groups = [
        ('HC',          dist == 'Control', HC_COLOR),
        ('RA DMARD-IR', dist.str.contains('DMARD-IR'), RA_DMARD_COLOR),
        ('RA TNF-IR',   dist.str.contains('TNF-IR'), RA_TNF_COLOR),
        ('SLE',         dist.str.contains('SLE'), SLE_COLOR),
    ]

    fig, axs = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, name in zip(axs, names):
        data = []; colors = []
        for gname, mask, color in groups:
            data.append(proj[name][mask.values])
            colors.append(color)
        bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                        showfliers=False)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color); patch.set_alpha(0.6)
        for med in bp['medians']:
            med.set_color('black'); med.set_linewidth(1.5)
        ax.set_xticklabels([g[0] for g in groups], rotation=20, ha='right',
                            fontsize=9)
        ax.set_ylabel(f'projection on {name}', fontsize=10)
        ax.axhline(0, color='gray', lw=0.5, linestyle='--', alpha=0.6)
        ax.grid(axis='y', alpha=0.3)
    fig.suptitle('GSE45291 (n=805) — 4-group axis projections.\n'
                  'PC9 KW p<10⁻¹²⁰, PC10 KW p<10⁻¹²⁴. SLE and RA on '
                  'opposite sides of HC.', fontsize=11)
    plt.tight_layout()
    plt.savefig(FIGDIR / 'fig2_gse45291_four_group.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print(f'  Saved {FIGDIR / "fig2_gse45291_four_group.png"}', flush=True)


def fig3_sle_per_patient_endotypes():
    print('Fig 3: SLE per-patient axis profile + endotypes...', flush=True)
    df = pd.read_csv(SCAN / 'sle_per_patient_axis_profile.csv')

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    pairs = [('PC5_IFN', 'PC10_plasma'),
             ('PC5_IFN', 'PC9_folate'),
             ('PC9_folate', 'PC10_plasma')]
    sle_mask = df.is_SLE
    sle_only = df.loc[sle_mask, ['PC5_IFN', 'PC9_folate', 'PC10_plasma']].values
    Zn = (sle_only - sle_only.mean(axis=0)) / (sle_only.std(axis=0) + 1e-9)
    km = KMeans(n_clusters=3, random_state=0, n_init=10).fit(Zn)
    # Label clusters by severity (PC5 highest = severe)
    centers = km.cluster_centers_
    order = np.argsort(-centers[:, 0])  # PC5 high to low
    label_map = {old: new for new, old in enumerate(order)}
    labels_remapped = np.array([label_map[c] for c in km.labels_])
    cluster_colors = [SEVERE_COLOR, MODERATE_COLOR, QUIESCENT_COLOR]
    cluster_names = ['Severe (IFN+plasma+folate)',
                      'Moderate', 'HC-like (quiescent)']

    for ax, (a, b) in zip(axs, pairs):
        # HC
        ax.scatter(df.loc[~sle_mask, a], df.loc[~sle_mask, b],
                    c=HC_COLOR, s=14, alpha=0.5, label=f'HC (n={(~sle_mask).sum()})',
                    edgecolors='none')
        # SLE clusters
        sle_df = df.loc[sle_mask].reset_index(drop=True)
        for c in range(3):
            mask = labels_remapped == c
            ax.scatter(sle_df.loc[mask, a], sle_df.loc[mask, b],
                        c=cluster_colors[c], s=12, alpha=0.45,
                        label=f'{cluster_names[c]} (n={mask.sum()})',
                        edgecolors='none')
        ax.set_xlabel(a, fontsize=10); ax.set_ylabel(b, fontsize=10)
        ax.axhline(0, color='gray', lw=0.5, linestyle='--', alpha=0.4)
        ax.axvline(0, color='gray', lw=0.5, linestyle='--', alpha=0.4)
        ax.grid(alpha=0.25)
        if ax == axs[0]:
            ax.legend(fontsize=7, loc='lower right')
    fig.suptitle('GSE65391 SLE per-patient axis-engagement profile (n=996). '
                  'K=3 endotypes match clinical disease activity (χ²p=1.4×10⁻³).',
                  fontsize=11)
    plt.tight_layout()
    plt.savefig(FIGDIR / 'fig3_sle_per_patient_endotypes.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print(f'  Saved {FIGDIR / "fig3_sle_per_patient_endotypes.png"}', flush=True)


def fig4_novel_target_validation():
    print('Fig 4: novel-target validation forest plot...', flush=True)
    # Hand-coded from the GSE49454 CYP-stratified validation results
    targets = [
        # (symbol, axis, z, neighbors_co_loading, log_mwu_p, validated)
        ('KDM8',   'PC10',  28.5,  2,  np.log10(0.008), True),
        ('PRMT1',  'PC5',   28.5, 26,  np.log10(0.013), True),
        ('CFI',    'PC5',   30.9, 10,  np.log10(0.023), True),
        ('GLO1',   'PC10',  33.1,  1,  np.log10(0.031), True),
        ('SAMHD1', 'PC10',  29.9,  1,  np.log10(0.062), False),
        ('CTSB',   'PC10',  30.7,  1,  np.log10(0.123), False),
        ('PDIA3',  'PC10',  43.0,  1,  np.log10(0.595), False),
    ]
    df = pd.DataFrame(targets, columns=[
        'symbol', 'axis', 'z_score', 'n_coload_neighbors',
        'log10_p', 'validated'
    ])
    df = df.sort_values('log10_p')   # most significant on top

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))

    # Left: forest plot of -log10(p)
    ax = axs[0]
    colors = ['#06A77D' if v else '#D62828' for v in df.validated]
    y = np.arange(len(df))
    ax.barh(y, -df.log10_p, color=colors, edgecolor='black', linewidth=0.5)
    ax.axvline(-np.log10(0.05), color='black', lw=1, linestyle='--',
               label='p=0.05 threshold')
    ax.set_yticks(y)
    ax.set_yticklabels([f'{r.symbol} ({r.axis})' for _, r in df.iterrows()],
                       fontsize=10)
    ax.set_xlabel('-log₁₀(MWU p) in GSE49454 CYP+ vs CYP- SLE', fontsize=10)
    ax.set_title('Independent validation of novel SLE targets', fontsize=10)
    ax.legend(loc='lower right', fontsize=8)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.25)

    # Right: z vs n_co-loading neighbors (substrate-connectivity predictor)
    ax = axs[1]
    for _, r in df.iterrows():
        color = '#06A77D' if r.validated else '#D62828'
        ax.scatter(r.n_coload_neighbors, r.z_score, s=160, c=color,
                    edgecolors='black', linewidths=0.8)
        ax.annotate(r.symbol, (r.n_coload_neighbors, r.z_score),
                    xytext=(7, 3), textcoords='offset points', fontsize=9)
    ax.set_xlabel('# of 1-hop substrate neighbors with |z|>3 on same axis',
                   fontsize=10)
    ax.set_ylabel('self z-score on axis', fontsize=10)
    ax.set_title('Substrate-connectivity predicts validation', fontsize=10)
    # Legend
    import matplotlib.patches as mpatches
    h1 = mpatches.Patch(color='#06A77D', label='Validated (p<0.05)')
    h2 = mpatches.Patch(color='#D62828', label='Failed')
    ax.legend(handles=[h1, h2], loc='upper right', fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle('Novel SLE target validation in GSE49454 (n=157)\n'
                 '4/7 candidates validate; substrate-connectivity discriminates',
                 fontsize=11)
    plt.tight_layout()
    plt.savefig(FIGDIR / 'fig4_novel_target_validation.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print(f'  Saved {FIGDIR / "fig4_novel_target_validation.png"}', flush=True)


def fig5_nnls_prescription_heatmap():
    print('Fig 5: NNLS prescription heatmap...', flush=True)
    p = SCAN / 'sle_signaware_recommendations.csv'
    if not p.exists():
        print(f'  {p} missing; skip', flush=True); return
    df = pd.read_csv(p)
    # Take signaware weights only
    w_cols = [c for c in df.columns if c.startswith('w_signaware_')]
    if not w_cols:
        print(f'  no signaware weight columns; skip', flush=True); return
    # Sort by SLEDAI if available
    if 'sledai' in df.columns:
        df_sorted = df.sort_values('sledai', ascending=False, na_position='last')
    else:
        df_sorted = df
    W = df_sorted[w_cols].values
    fig, ax = plt.subplots(figsize=(8, 10))
    im = ax.imshow(np.log1p(W), aspect='auto', cmap='Reds')
    ax.set_xticks(np.arange(len(w_cols)))
    ax.set_xticklabels([c.replace('w_signaware_', '') for c in w_cols],
                       rotation=45, ha='right', fontsize=9)
    ax.set_ylabel(f'patients (n={len(df_sorted)}) — sorted by SLEDAI ↓', fontsize=9)
    ax.set_title('Sign-aware NNLS prescription weights (log1p)', fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.04, label='log1p(weight)')
    plt.tight_layout()
    plt.savefig(FIGDIR / 'fig5_nnls_prescription_heatmap.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print(f'  Saved {FIGDIR / "fig5_nnls_prescription_heatmap.png"}', flush=True)


def main():
    print(f'Output directory: {FIGDIR}', flush=True)
    try: fig2_gse45291_four_group()
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f'  Fig 2 failed: {e}', flush=True)
    try: fig3_sle_per_patient_endotypes()
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f'  Fig 3 failed: {e}', flush=True)
    try: fig4_novel_target_validation()
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f'  Fig 4 failed: {e}', flush=True)
    try: fig5_nnls_prescription_heatmap()
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f'  Fig 5 failed: {e}', flush=True)


if __name__ == '__main__':
    main()
