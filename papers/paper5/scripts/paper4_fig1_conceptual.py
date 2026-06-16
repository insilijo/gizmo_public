"""DDD — Figure 1: conceptual schematic of the coordinate transformation chain.

5-panel schematic illustrating:
  A. Raw expression (uninterpretable, no biology prior)
  B. F = substrate-aware Laplacian smoothing
  C. β/α orthogonal decomposition (universal burden + multi-modal complement)
  D. α → PC eigenmodes via spectral PCA
  E. Per-patient projection onto canonical axes

Pure matplotlib schematic; no real data. Saves to:
  benchmarks/results/figures/paper4/fig1_conceptual_schematic.png
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO = Path('/home/jgardner/GIZMO')
FIGDIR = REPO / 'benchmarks/results/figures/paper4'
FIGDIR.mkdir(parents=True, exist_ok=True)

# Brand colors
HC = '#3A86FF'
SLE = '#FF006E'
SUBSTRATE_BG = '#F0F4F8'
BOX_BG = '#FFFFFF'
ACCENT = '#D62828'


def main():
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(3, 5, hspace=0.55, wspace=0.4)

    # ====== TOP ROW: Coordinate transformation chain ======
    # A: Raw expression
    ax = fig.add_subplot(gs[0, 0])
    np.random.seed(0)
    expr = np.random.normal(0, 1, (12, 12))
    im = ax.imshow(expr, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
    ax.set_title('A. Raw expression\n(genes × samples)', fontsize=10,
                  fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(0.5, -0.18, 'no biology prior\nuninterpretable factors',
             transform=ax.transAxes, ha='center', va='top',
             fontsize=8, style='italic', color='#555')

    # B: Substrate graph + smoothing
    ax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor(SUBSTRATE_BG)
    # Node-link layout
    rng = np.random.default_rng(2)
    n_nodes = 18
    pos = rng.uniform(0.1, 0.9, (n_nodes, 2))
    # Edges: nearest neighbors
    edges = []
    for i in range(n_nodes):
        d = np.linalg.norm(pos - pos[i], axis=1)
        for j in np.argsort(d)[1:3]:
            if (j, i) not in edges and (i, j) not in edges:
                edges.append((i, j))
    for (i, j) in edges:
        ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
                color='#888', lw=1.0, alpha=0.7, zorder=1)
    sig = np.zeros(n_nodes)
    sig[2] = 1.0; sig[5] = 0.7
    for _ in range(2):
        new_sig = sig.copy()
        for k in range(n_nodes):
            neighbors = [j for (i, j) in edges if i == k] + \
                        [i for (i, j) in edges if j == k]
            if neighbors:
                new_sig[k] = 0.5 * sig[k] + 0.5 * np.mean(
                    [sig[n] for n in neighbors])
        sig = new_sig
    ax.scatter(pos[:, 0], pos[:, 1], c=sig, cmap='Reds', s=200,
                edgecolors='black', linewidths=0.6, zorder=2,
                vmin=0, vmax=1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title('B. F = solve_map\n(graph-smoothed signal)', fontsize=10,
                  fontweight='bold')
    ax.text(0.5, -0.18, 'substrate ≡ biology prior\ncrosstalk via Laplacian',
             transform=ax.transAxes, ha='center', va='top',
             fontsize=8, style='italic', color='#555')

    # C: β/α decomposition
    ax = fig.add_subplot(gs[0, 2])
    # Draw F vector decomposed into β (along log_PR axis) + α (orthogonal)
    F_vec = np.array([0.7, 0.5])
    lpn = np.array([1.0, 0.0])
    beta = F_vec @ lpn
    alpha = F_vec - beta * lpn
    ax.arrow(0, 0, F_vec[0], F_vec[1], head_width=0.04, head_length=0.04,
              fc='black', ec='black', length_includes_head=True, lw=2)
    ax.text(F_vec[0] + 0.04, F_vec[1] + 0.02, 'F', fontsize=14,
             fontweight='bold')
    ax.arrow(0, 0, beta, 0, head_width=0.04, head_length=0.04,
              fc=HC, ec=HC, length_includes_head=True, lw=2)
    ax.text(beta/2, -0.08, 'β · log_PR', fontsize=10, color=HC, ha='center',
             fontweight='bold')
    ax.arrow(beta, 0, 0, alpha[1], head_width=0.04, head_length=0.04,
              fc=SLE, ec=SLE, length_includes_head=True, lw=2)
    ax.text(beta + 0.05, alpha[1]/2, 'α', fontsize=14, color=SLE,
             fontweight='bold')
    ax.set_xlim(-0.1, 0.95); ax.set_ylim(-0.2, 0.7)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect('equal')
    ax.set_title('C. F = β·log_PR + α\n(orthogonal decomposition)',
                  fontsize=10, fontweight='bold')
    ax.text(0.5, -0.32, 'β: universal burden\nα: multi-modal complement',
             transform=ax.transAxes, ha='center', va='top',
             fontsize=8, style='italic', color='#555')

    # D: α → PCs
    ax = fig.add_subplot(gs[0, 3])
    # Multi-mode amplitude bars (representing PC1..PC5 components)
    pcs = ['PC1', 'PC2', 'PC3', 'PC5\nIFN', 'PC9\nfolate', 'PC10\nplasma']
    amps = [3.5, 1.5, 1.0, 0.95, 0.65, 0.72]
    colors = ['#999', '#999', '#999', SLE, '#FFB84D', '#6A0DAD']
    ax.bar(range(len(pcs)), amps, color=colors, edgecolor='black',
            linewidth=0.6)
    ax.set_xticks(range(len(pcs)))
    ax.set_xticklabels(pcs, fontsize=8)
    ax.set_ylabel('|amplitude|', fontsize=9)
    ax.set_title('D. PCA on α\n(spectral eigenmodes)', fontsize=10,
                  fontweight='bold')
    ax.text(0.5, -0.42, 'mid-low PCs = interpretable\nbiological modules',
             transform=ax.transAxes, ha='center', va='top',
             fontsize=8, style='italic', color='#555')

    # E: per-patient projection 3D-ish scatter
    ax = fig.add_subplot(gs[0, 4])
    rng2 = np.random.default_rng(3)
    # HC cluster
    hc_pts = rng2.normal(0, 0.4, (20, 2))
    ax.scatter(hc_pts[:, 0], hc_pts[:, 1], c=HC, s=40, alpha=0.7,
                label='HC', edgecolors='none')
    # SLE cluster shifted
    sle_pts = rng2.normal(1.5, 0.5, (40, 2))
    ax.scatter(sle_pts[:, 0], sle_pts[:, 1], c=SLE, s=40, alpha=0.7,
                label='SLE', edgecolors='none')
    ax.set_xlabel('projection on PC5 (IFN)', fontsize=9)
    ax.set_ylabel('projection on PC10 (plasma)', fontsize=9)
    ax.set_title('E. Per-patient profile\n(disease coords)',
                  fontsize=10, fontweight='bold')
    ax.axhline(0, color='gray', lw=0.5, linestyle='--', alpha=0.4)
    ax.axvline(0, color='gray', lw=0.5, linestyle='--', alpha=0.4)
    ax.legend(fontsize=8, loc='upper left')

    # ====== MIDDLE ROW: arrows + label ======
    ax = fig.add_subplot(gs[1, :])
    ax.axis('off')
    # Chain label
    ax.annotate('', xy=(0.96, 0.6), xytext=(0.05, 0.6),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='black', lw=2.5))
    label_text = ('Coordinate transformation chain — biology-aware basis '
                  'where modes approximately decouple (PDE-style)')
    ax.text(0.5, 0.85, label_text, ha='center', va='center',
             fontsize=11, fontweight='bold',
             transform=ax.transAxes)

    # ====== BOTTOM ROW: applications ======
    ax = fig.add_subplot(gs[2, 0])
    # Endotype K-means
    rng3 = np.random.default_rng(4)
    cluster_colors = ['#D62828', '#9D4EDD', '#06A77D']
    cluster_names = ['Severe', 'Moderate', 'Quiescent']
    cluster_centers = [(1.8, -1.5), (0.5, 0), (-1.2, 1.2)]
    for c, cn, color in zip(cluster_centers, cluster_names, cluster_colors):
        pts = rng3.normal(c, 0.45, (35, 2))
        ax.scatter(pts[:, 0], pts[:, 1], c=color, alpha=0.6, s=30,
                    edgecolors='none', label=cn)
    ax.set_xlabel('PC5 (IFN) →', fontsize=9)
    ax.set_ylabel('PC10 (plasma) →', fontsize=9)
    ax.set_title('Endotype stratification\n(K=3 SLE)', fontsize=10,
                  fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')

    # Drug-target output
    ax = fig.add_subplot(gs[2, 1])
    targets = ['STAT1', 'JAK1', 'PSMB9', 'CD20', 'DHFR', 'PRMT1*', 'KDM8*']
    statuses = ['approved', 'approved', 'approved', 'approved', 'approved',
                 'novel', 'novel']
    z_scores = [40, 30, 28, 24, 32, 28, 28]
    colors = ['#06A77D' if s == 'approved' else '#D62828' for s in statuses]
    ax.barh(range(len(targets)), z_scores, color=colors,
             edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(targets)))
    ax.set_yticklabels(targets, fontsize=9)
    ax.set_xlabel('axis |z|', fontsize=9)
    ax.invert_yaxis()
    ax.set_title('Drug-target output\n(approved + novel)', fontsize=10,
                  fontweight='bold')
    h_approved = mpatches.Patch(color='#06A77D', label='Approved drug')
    h_novel = mpatches.Patch(color='#D62828', label='Novel candidate')
    ax.legend(handles=[h_approved, h_novel], fontsize=7, loc='lower right')

    # NNLS prescription
    ax = fig.add_subplot(gs[2, 2])
    regimes = ['MTX', 'JAKi', 'anti-CD20', 'bortezomib']
    weights = [0.6, 0.45, 0.3, 0.2]
    bars = ax.barh(regimes, weights, color='#FF006E',
                    edgecolor='black', linewidth=0.5)
    ax.set_xlabel('NNLS weight', fontsize=9)
    ax.invert_yaxis()
    ax.set_title('Combination prescription\n(per-patient NNLS)',
                  fontsize=10, fontweight='bold')
    ax.set_xlim(0, 0.75)

    # Cross-disease comparison
    ax = fig.add_subplot(gs[2, 3])
    groups = ['HC', 'RA\nDMARD-IR', 'RA\nTNF-IR', 'SLE']
    pc10_means = [-3.85, +4.06, +2.34, -5.92]
    group_colors = [HC, '#FFBE0B', '#FB5607', SLE]
    bars = ax.bar(range(len(groups)), pc10_means, color=group_colors,
                   edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, fontsize=8)
    ax.set_ylabel('PC10 projection mean', fontsize=9)
    ax.set_title('Cross-disease separation\n(GSE45291 p<10⁻¹²⁰)',
                  fontsize=10, fontweight='bold')
    ax.axhline(0, color='gray', lw=0.5, linestyle='--', alpha=0.5)

    # Sign-aware contraindication
    ax = fig.add_subplot(gs[2, 4])
    drugs = ['anti-TNF', 'Filbin', 'Wang_RA', 'MTX', 'bariatric']
    contra_pct = [96, 95, 98, 10, 10]
    colors = ['#D62828' if p > 50 else '#06A77D' for p in contra_pct]
    ax.barh(drugs, contra_pct, color=colors, edgecolor='black',
             linewidth=0.5)
    ax.set_xlabel('% SLE patients\ncontraindicated', fontsize=8)
    ax.invert_yaxis()
    ax.axvline(50, color='gray', lw=0.5, linestyle='--', alpha=0.5)
    ax.set_title('Sign-aware filter\n(anti-TNF auto-flagged)',
                  fontsize=10, fontweight='bold')

    plt.suptitle('Substrate-aware coordinate system for multi-axis disease '
                 'biology', fontsize=13, fontweight='bold', y=0.998)
    plt.savefig(FIGDIR / 'fig1_conceptual_schematic.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print(f'Saved {FIGDIR / "fig1_conceptual_schematic.png"}')


if __name__ == '__main__':
    main()
