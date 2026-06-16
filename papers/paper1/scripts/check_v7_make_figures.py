"""Generate v7 manuscript figures from deposited results."""
import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

REPO = Path('/home/jgardner/GIZMO')
ZS = REPO / 'benchmarks/results/unsupervised/zscored'
FIG_DIR = REPO / 'figures_v7'
FIG_DIR.mkdir(exist_ok=True)

# === FIGURE 2: Phase 3b interpretability heatmap ===
print('Building Figure 2...', flush=True)
d = json.load(open(ZS / 'v7_interpretability_eval_v4.json'))
cells = d['cells']
cohorts = sorted({c['cohort'] for c in cells})
# Build cohort × (pc, sign) matrix; PCs labeled PC1+, PC1-, PC2+, ..., PC5-
labels = []
for pc in range(1, 6):
    labels.append(f'PC{pc}+')
    labels.append(f'PC{pc}-')
auroc_M = np.full((len(cohorts), len(labels)), np.nan)
p_M = np.full((len(cohorts), len(labels)), np.nan)
for c in cells:
    if c['pc'] > 5: continue
    i = cohorts.index(c['cohort'])
    j = labels.index(f'PC{c["pc"]}{c["sign"]}')
    auroc_M[i, j] = c['auroc']
    p_M[i, j] = c['p_empirical']

fig, ax = plt.subplots(figsize=(11, 6))
im = ax.imshow(auroc_M, aspect='auto', cmap='RdBu_r', vmin=0.3, vmax=0.7)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=0)
ax.set_yticks(range(len(cohorts)))
ax.set_yticklabels(cohorts)
plt.colorbar(im, ax=ax, label='AUROC (truth-gene vs |loading|)')
# Annotate cells with p<0.05 stars + AUROC value
for i in range(len(cohorts)):
    for j in range(len(labels)):
        if np.isnan(auroc_M[i, j]): continue
        v = auroc_M[i, j]; p = p_M[i, j]
        color = 'white' if v > 0.65 or v < 0.35 else 'black'
        annot = f'{v:.2f}'
        if p < 0.01: annot += '**'
        elif p < 0.05: annot += '*'
        ax.text(j, i, annot, ha='center', va='center', fontsize=7, color=color)
ax.set_title('Figure 2: Quantitative interpretability evaluation\n'
              'Cohort × α-PC × sign AUROC (truth-gene rank in loading vector)\n'
              '* p<0.05 degree-preserving null  ** p<0.01     24/110 cells significant; 10/11 cohorts pass',
              fontsize=10)
ax.set_xlabel('α-PC × sign')
ax.set_ylabel('Cohort')
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure2_interpretability.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure2_interpretability.png')


# === FIGURE 4: Per-basin survival decomposition ===
print('Building Figure 4...', flush=True)
df_surv = pd.read_csv(ZS / 'v7_basin_survival_broad.tsv', sep='\t')
# Add Filbin manually + TCGA_IDH
filbin = pd.DataFrame([
    ('Filbin_COVID', 'PC1-', 'Proteoglycan/MMP', 19, 7, 0.776),
    ('Filbin_COVID', 'PC4+', 'ECM/HMOX1/IL10', 54, 19, 0.768),
    ('Filbin_COVID', 'PC7-', 'Proteoglycan', 25, 5, 0.764),
    ('Filbin_COVID', 'PC5-', 'Cytokine storm', 93, 24, 0.730),
    ('Filbin_COVID', 'PC2-', 'Cadherin/ECM', 175, 21, 0.727),
    ('Filbin_COVID', 'PC3+', 'Innate immune', 94, 9, 0.608),
    ('Filbin_COVID', 'PC7+', 'Mixed signal', 135, 42, 0.605),
    ('Filbin_COVID', 'PC3-', 'RTK signaling', 98, 20, 0.561),
    ('Filbin_COVID', 'PC6+', 'IFN signaling', 97, 17, 0.546),
    ('Filbin_COVID', 'PC2+', 'Coagulation', 27, 9, 0.492),
    ('Filbin_COVID', 'PC1+', 'Apoptotic sig', 170, 26, 0.482),
    ('Filbin_COVID', 'PC6-', 'Mixed MMP', 42, 12, 0.422),
], columns=['cohort','pc','category','n_genes','n_obs','cindex'])
tcga_idh = pd.DataFrame([
    ('TCGA_IDH', 'PC2-', 'Ciliary', 7, 7, 0.590),
    ('TCGA_IDH', 'PC3+', 'OXPHOS', 16, 3, 0.543),
    ('TCGA_IDH', 'PC4-', 'DNA damage', 3, 3, 0.526),
    ('TCGA_IDH', 'PC1+', 'Neuronal', 3, 3, 0.517),
], columns=['cohort','pc','category','n_genes','n_obs','cindex'])
df_all = pd.concat([df_surv, filbin, tcga_idh], ignore_index=True)

# 4-panel figure
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
baseline_F = {'KMPLOT_BRCA': 0.496, 'TCGA_LUAD': 0.580, 'TCGA_IDH': 0.824,
                'Filbin_COVID': 0.787}
baseline_PCA = {'KMPLOT_BRCA': 0.596, 'TCGA_LUAD': 0.599, 'TCGA_IDH': 0.821,
                  'Filbin_COVID': 0.795}
metric = {'KMPLOT_BRCA': 'C-index', 'TCGA_LUAD': 'C-index',
            'TCGA_IDH': 'C-index', 'Filbin_COVID': 'AUC'}

cohort_panels = [
    ('TCGA_LUAD', axes[0, 0]),
    ('Filbin_COVID', axes[0, 1]),
    ('TCGA_IDH', axes[1, 0]),
    ('KMPLOT_BRCA', axes[1, 1]),
]
for cohort, ax in cohort_panels:
    sub = df_all[df_all['cohort'] == cohort].sort_values('cindex', ascending=True)
    if len(sub) == 0:
        ax.text(0.5, 0.5, f'{cohort}\n(no per-basin data)',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        continue
    y = np.arange(len(sub))
    colors = ['#1f77b4' if c > 0.55 else '#aaaaaa' for c in sub['cindex']]
    bars = ax.barh(y, sub['cindex'], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels([f'{r["pc"]}: {r["category"]} ({int(r["n_obs"])}g)'
                         for _, r in sub.iterrows()], fontsize=8)
    ax.axvline(0.5, color='gray', linestyle=':', linewidth=0.8, label='Chance (0.50)')
    ax.axvline(baseline_F[cohort], color='red', linestyle='--', linewidth=1,
               label=f'F-7 ensemble ({baseline_F[cohort]:.2f})')
    ax.axvline(baseline_PCA[cohort], color='blue', linestyle='--', linewidth=1,
               label=f'PCA-on-input ({baseline_PCA[cohort]:.2f})')
    ax.set_xlabel(metric[cohort])
    ax.set_title(f'{cohort} per-basin {metric[cohort]}')
    ax.legend(loc='lower right', fontsize=7)
    ax.set_xlim(0.35, 0.85)
fig.suptitle('Figure 4: Per-basin patient activation scores → survival/mortality discrimination',
             fontsize=12)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure4_basin_survival.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure4_basin_survival.png')


# === FIGURE 5: F-basin vs PCA orthogonal gene panels ===
print('Building Figure 5...', flush=True)
overlap_data = [
    {'cohort': 'TCGA_LUAD', 'F_basin': 'T cell/MHC (PC4+, 68g)',
     'F_C': 0.608, 'PCA_PC': 'PCA-PC1: SFTPC/SCGB1A1 (alveolar)', 'PCA_C': 0.577,
     'overlap': 2, 'F_size': 68, 'PCA_size': 50},
    {'cohort': 'TCGA_IDH', 'F_basin': 'Ciliary DNAAF (PC2-, 7g)',
     'F_C': 0.590, 'PCA_PC': 'PCA-PC1: CLIC1/S100A11 (motility)', 'PCA_C': 0.624,
     'overlap': 0, 'F_size': 7, 'PCA_size': 50},
    {'cohort': 'Filbin_COVID', 'F_basin': 'ECM degradation (PC1-, 18g)',
     'F_C': 0.776, 'PCA_PC': 'PCA-PC1: NTproBNP/FGF23 (organ damage)', 'PCA_C': 0.717,
     'overlap': 0, 'F_size': 18, 'PCA_size': 50},
]
fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
for ax, od in zip(axes, overlap_data):
    # Bar chart: F-only, overlap, PCA-only
    f_only = od['F_size'] - od['overlap']
    pca_only = od['PCA_size'] - od['overlap']
    bars = ['F-basin only', 'Shared', 'PCA-PC only']
    vals = [f_only, od['overlap'], pca_only]
    colors = ['#1f77b4', '#888888', '#ff7f0e']
    ax.bar(bars, vals, color=colors)
    for i, v in enumerate(vals):
        ax.text(i, v + 1, str(v), ha='center', fontsize=10)
    ax.set_ylabel('# genes')
    ax.set_title(f'{od["cohort"]}\n'
                  f'F-basin: {od["F_basin"]}, score={od["F_C"]:.3f}\n'
                  f'PCA-PC top-50: score={od["PCA_C"]:.3f}',
                  fontsize=9)
    ax.tick_params(axis='x', labelrotation=20)
fig.suptitle('Figure 5: F-basins and PCA-best-PC read orthogonal prognostic biology\n'
              '(0-2 shared genes across three cohorts; comparable discrimination)',
              fontsize=11)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure5_f_pca_orthogonal.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure5_f_pca_orthogonal.png')

# Reuse Figure 3 (basin activation matrix - already exists)
import shutil
src = ZS / 'v7_basin_activation_matrix.png'
dst = FIG_DIR / 'figure3_activation_matrix.png'
shutil.copy(src, dst)
print(f'  Copied {dst}')

print('\nAll figures generated in', FIG_DIR)
"""Two additional v7 figures: per-α-PC F vs PCA ablation + preprocessing scope."""
import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path('/home/jgardner/GIZMO')
ZS = REPO / 'benchmarks/results/unsupervised/zscored'
FIG_DIR = REPO / 'figures_v7'

# === FIGURE 6: Per-α-PC F vs PCA ablation across 4 cohorts ===
print('Building Figure 6...', flush=True)
# F-α-PC under v5-canonical (Phase 5a v5)
f_v5 = json.load(open(ZS / 'v7_v5_F_per_pc_ablation.json'))
# PCA on input per-PC
pca = json.load(open(ZS / 'v7_pca_per_pc_ablation.json'))

cohorts = ['KMPLOT_BRCA', 'TCGA_LUAD', 'TCGA_IDH_glioma', 'Filbin_COVID']
labels = ['KMPLOT_BRCA\n(BRCA OS)', 'TCGA_LUAD\n(LUAD OS)',
            'TCGA_IDH_glioma\n(Glioma OS)', 'Filbin_COVID\n(28d Mortality)']
fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), sharey=True)
n_pcs = 7
x = np.arange(1, n_pcs + 1)
for k, (cohort, label) in enumerate(zip(cohorts, labels)):
    ax = axes[k]
    f_scores = [f_v5[cohort].get(f'α-PC{i}', {}).get('mean') if isinstance(f_v5[cohort], dict) and isinstance(f_v5[cohort].get(f'α-PC{i}'), dict) else None for i in range(1, n_pcs + 1)]
    pca_scores = [pca[cohort].get(f'PC{i}', {}).get('mean') if isinstance(pca[cohort].get(f'PC{i}'), dict) else None for i in range(1, n_pcs + 1)]
    f_y = [v if v is not None else np.nan for v in f_scores]
    pca_y = [v if v is not None else np.nan for v in pca_scores]
    ax.plot(x, f_y, 'o-', label='F-α-PC (v5-canonical)', color='#d62728', markersize=8)
    ax.plot(x, pca_y, 's-', label='PCA-on-input', color='#1f77b4', markersize=7)
    ax.axhline(0.5, color='gray', linestyle=':', linewidth=0.7, label='Chance')
    # Mark best F and best PCA
    best_f_i = np.nanargmax(f_y); best_p_i = np.nanargmax(pca_y)
    ax.annotate(f'F best: PC{best_f_i+1}\n{f_y[best_f_i]:.3f}',
                xy=(best_f_i+1, f_y[best_f_i]), xytext=(best_f_i+1, f_y[best_f_i]+0.04),
                fontsize=8, ha='center', color='#d62728',
                arrowprops=dict(arrowstyle='->', color='#d62728', alpha=0.7))
    ax.annotate(f'PCA best: PC{best_p_i+1}\n{pca_y[best_p_i]:.3f}',
                xy=(best_p_i+1, pca_y[best_p_i]), xytext=(best_p_i+1, pca_y[best_p_i]-0.06),
                fontsize=8, ha='center', color='#1f77b4',
                arrowprops=dict(arrowstyle='->', color='#1f77b4', alpha=0.7))
    ax.set_xlabel('PC index')
    ax.set_xticks(x)
    ax.set_title(label, fontsize=10)
    if k == 0:
        ax.set_ylabel('Discrimination (Cox C-index / AUC)')
        ax.legend(loc='lower left', fontsize=8)
    ax.set_ylim(0.35, 0.85)
    # Annotate gap for TCGA_IDH (the big win)
    if cohort == 'TCGA_IDH_glioma':
        gap = f_y[best_f_i] - pca_y[best_p_i]
        ax.text(0.5, 0.97, f'Δ best-PC: +{gap:.3f}',
                transform=ax.transAxes, ha='center', va='top',
                fontsize=11, fontweight='bold', color='darkgreen',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#e0f5e0', edgecolor='darkgreen'))
fig.suptitle('Figure 6: Per-α-PC F vs PCA-on-input ablation (matched preprocessing)\n'
              'F-α-PCs match or exceed PCA-PCs at the best-single-PC level in all 4 cohorts; '
              'TCGA_IDH +0.138 = substrate-projection on IDH-mut catalysis axis',
              fontsize=10)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure6_f_vs_pca_per_pc.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure6_f_vs_pca_per_pc.png')


# === FIGURE 7: Preprocessing scope - z-score vs canonical ===
print('Building Figure 7...', flush=True)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: 2HG anchor rank shift on TCGA_IDH
ax = axes[0]
labels_p = ['Canonical\n(per-mod log + global-std)',
              'Within-patient z-score\n(§5 diagnostic)']
ranks = [6355, 26]
colors = ['#1f77b4', '#d62728']
bars = ax.bar(labels_p, ranks, color=colors)
ax.set_ylabel('Rank of mito-2HG anchor (R-ALL-879997)\n[out of 38,148 substrate nodes; lower = better]')
ax.set_yscale('log')
for b, v in zip(bars, ranks):
    ax.text(b.get_x() + b.get_width()/2, v * 1.3, str(v),
              ha='center', fontsize=11, fontweight='bold')
ax.axhline(38148 * 0.05, color='gray', linestyle='--', alpha=0.6, label='Top 5% cutoff (rank 1907)')
ax.set_title('A) TCGA_IDH-glioma: mito-2HG anchor recovery\nis preprocessing-conditional',
             fontsize=10)
ax.legend(fontsize=8, loc='upper right')

# Panel B: F-feature C-index under each preprocessing
ax = axes[1]
cohorts_b = ['KMPLOT_BRCA', 'TCGA_LUAD', 'TCGA_IDH']
zscore_vals = [0.451, 0.580, 0.790]
canon_vals = [0.496, 0.555, 0.824]  # v5-canonical
pca_input_vals = [0.596, 0.616, 0.821]  # PCA-substrate-matched

x = np.arange(len(cohorts_b))
w = 0.27
ax.bar(x - w, zscore_vals, w, label='F (within-patient z-score)', color='#ff7f0e')
ax.bar(x, canon_vals, w, label='F (canonical: per-mod-std)', color='#d62728')
ax.bar(x + w, pca_input_vals, w, label='PCA-on-substrate-matched input', color='#1f77b4')
ax.axhline(0.5, color='gray', linestyle=':', linewidth=0.7)
ax.set_xticks(x)
ax.set_xticklabels(cohorts_b)
ax.set_ylabel('Cox C-index (5-fold CV)')
ax.set_title('B) Survival discrimination: canonical preprocessing\nis the right choice for magnitude-driven phenotypes',
             fontsize=10)
ax.legend(fontsize=8, loc='lower right')
ax.set_ylim(0.4, 0.9)

fig.suptitle('Figure 7: Scope conditions — preprocessing choice is task-specific',
             fontsize=11)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure7_preprocessing_scope.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure7_preprocessing_scope.png')

print('\nDone — figures_v7 now has 6 generated figures (figure1 schematic still TBD).')
"""Network sub-graph visualizations + patient embeddings for v7."""
import sys
sys.path.insert(0, '/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
import json
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS = REPO / 'benchmarks/results/unsupervised/zscored'
FIG_DIR = REPO / 'figures_v7'

# Load substrate
mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr
sub_node_set = set(geom.nodes)
G = nx.Graph()
G.add_nodes_from(geom.nodes)
for u, v in mg.graph.edges():
    if u in sub_node_set and v in sub_node_set:
        G.add_edge(u, v)


def get_basin_nodes(cohort, pc_idx_0based, sign, thresh_q=0.95):
    """Return nodes in the basin (largest CC of top-5% loading same-sign)."""
    npz = np.load(SNAPSHOT / f'stage3_F_{cohort}_edge_informed.npz', allow_pickle=True)
    F = npz['F']
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm = x.mean(); xv = x.var() + 1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm) * (x - xm)).mean(axis=1, keepdims=True)
    beta = (cov / xv).ravel()
    alpha = F_unit - Fm - beta[:, None] * (x - xm)[None, :]
    pca = PCA(n_components=7, random_state=0); pca.fit(alpha)
    loadings = pca.components_[pc_idx_0based]
    abs_l = np.abs(loadings)
    threshold = np.quantile(abs_l, thresh_q)
    sign_int = 1 if sign == '+' else -1
    mask = (np.sign(loadings) == sign_int) & (abs_l >= threshold)
    keep = [geom.nodes[i] for i in range(len(loadings)) if mask[i]]
    H = G.subgraph(keep)
    ccs = sorted(nx.connected_components(H), key=len, reverse=True)
    if not ccs: return set(), loadings
    return ccs[0], loadings


def node_label(nid):
    a = mg.graph.nodes.get(nid, {})
    if a.get('node_type') == 'gene':
        sym = a.get('symbol') or a.get('name') or nid.replace('symbol:', '')
        return sym
    return None  # reactions and metabolites get no label (just colored)


def node_type(nid):
    return mg.graph.nodes.get(nid, {}).get('node_type', '?')


# === FIGURE 8: Network visualization of conserved OXPHOS basin (TCGA_IDH PC3+ and KMPLOT PC1-) ===
print('Building Figure 8...', flush=True)

# TCGA_IDH PC3+ (mt-OXPHOS) and KMPLOT PC1- (mt-OXPHOS)
idh_basin, _ = get_basin_nodes('TCGA_IDH_glioma', 2, '+')  # PC3 = index 2
kmplot_basin, _ = get_basin_nodes('KMPLOT_BRCA', 0, '-')

# Show side-by-side network sub-graphs
fig, axes = plt.subplots(1, 2, figsize=(15, 8))
NODE_COLORS = {'gene': '#1f77b4', 'reaction': '#ff7f0e', 'metabolite': '#2ca02c'}

for ax, basin, title in [(axes[0], idh_basin, 'TCGA_IDH α-PC3+ (16 gene, 28 reaction nodes)\nMitochondrial OXPHOS basin'),
                            (axes[1], kmplot_basin, 'KMPLOT_BRCA α-PC1- (19 gene, 17 reaction nodes)\nMitochondrial OXPHOS basin')]:
    if not basin:
        ax.text(0.5, 0.5, 'empty basin', ha='center', va='center')
        continue
    H = G.subgraph(basin)
    pos = nx.spring_layout(H, k=1.5/np.sqrt(len(basin)), iterations=80, seed=42)
    # Draw edges
    nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.4, edge_color='gray', width=0.8)
    # Draw nodes by type
    for ntype in ['reaction', 'metabolite', 'gene']:
        nodes = [n for n in basin if node_type(n) == ntype]
        if not nodes: continue
        nx.draw_networkx_nodes(H, pos, nodelist=nodes,
                                  node_color=NODE_COLORS[ntype],
                                  node_size=250 if ntype == 'gene' else 80,
                                  ax=ax, alpha=0.85,
                                  label=ntype)
    # Label only genes
    labels = {n: node_label(n) for n in basin if node_label(n)}
    nx.draw_networkx_labels(H, pos, labels, font_size=8, ax=ax)
    ax.set_title(title, fontsize=11)
    ax.set_axis_off()
    legend_elem = [Patch(facecolor=NODE_COLORS[t], label=t.title()) for t in ['gene', 'reaction', 'metabolite']]
    ax.legend(handles=legend_elem, loc='upper right', fontsize=8)

fig.suptitle('Figure 8: Substrate sub-graph (basin) of conserved mitochondrial OXPHOS biology\n'
              'Same biochemical neighborhood surfaces on different α-PC indices in unrelated cohorts',
              fontsize=11)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure8_oxphos_basin_network.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure8_oxphos_basin_network.png')


# === FIGURE 9: Filbin ECM-degradation basin network (the headline F-vs-PCA winning basin) ===
print('Building Figure 9...', flush=True)
filbin_basin, _ = get_basin_nodes('Filbin_COVID', 0, '-')  # PC1 = index 0
fig, ax = plt.subplots(figsize=(11, 9))
if filbin_basin:
    H = G.subgraph(filbin_basin)
    pos = nx.spring_layout(H, k=1.8/np.sqrt(len(filbin_basin)), iterations=80, seed=42)
    nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.4, edge_color='gray', width=0.7)
    for ntype in ['reaction', 'metabolite', 'gene']:
        nodes = [n for n in filbin_basin if node_type(n) == ntype]
        if not nodes: continue
        nx.draw_networkx_nodes(H, pos, nodelist=nodes,
                                  node_color=NODE_COLORS[ntype],
                                  node_size=350 if ntype == 'gene' else 100,
                                  ax=ax, alpha=0.85, label=ntype)
    labels = {n: node_label(n) for n in filbin_basin if node_label(n)}
    nx.draw_networkx_labels(H, pos, labels, font_size=9, ax=ax)
    ax.set_axis_off()
    legend_elem = [Patch(facecolor=NODE_COLORS[t], label=t.title()) for t in ['gene', 'reaction', 'metabolite']]
    ax.legend(handles=legend_elem, loc='upper right', fontsize=9)
ax.set_title('Figure 9: Filbin_COVID α-PC1- basin — ECM degradation by neutrophil enzymes\n'
              f'19 gene + reaction nodes; AUC 0.776 for 28d mortality (vs PCA 0.717, Δ +0.06)\n'
              'Proteoglycan core (BGN/DCN/ACAN/VCAN/LUM/FMOD/KERA) + cathepsins (CTSL/CTSK) + MMP (MMP13/20)',
              fontsize=10)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure9_filbin_ecm_basin_network.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure9_filbin_ecm_basin_network.png')


# === FIGURE 1: Substrate construction schematic + per-patient pipeline ===
print('Building Figure 1...', flush=True)
fig, axes = plt.subplots(1, 4, figsize=(18, 5))

# Panel A: substrate sub-graph (small visual)
ax = axes[0]
# Generate a small representative subgraph of the substrate for visualization
sample_nodes = list(geom.nodes)[:80]
H_sample = G.subgraph(sample_nodes)
ccs = sorted(nx.connected_components(H_sample), key=len, reverse=True)
if ccs:
    H_show = G.subgraph(ccs[0] | set(list(ccs[1] if len(ccs) > 1 else set())))
    pos = nx.spring_layout(H_show, seed=42)
    nx.draw_networkx_edges(H_show, pos, ax=ax, alpha=0.4, edge_color='gray', width=0.5)
    for ntype in ['reaction', 'metabolite', 'gene']:
        nodes = [n for n in H_show.nodes if node_type(n) == ntype]
        if not nodes: continue
        nx.draw_networkx_nodes(H_show, pos, nodelist=nodes,
                                  node_color=NODE_COLORS[ntype],
                                  node_size=120 if ntype == 'gene' else 60,
                                  ax=ax, alpha=0.85)
ax.set_title('A) Substrate (38,148 nodes)\nReactome + StringDB + HMDB + KEGG\nCC-BY 4.0',
             fontsize=10)
ax.set_axis_off()

# Panel B: MAP solve
ax = axes[1]
ax.text(0.5, 0.85, 'B) Per-patient MAP solve', ha='center', va='top',
        fontsize=11, fontweight='bold')
ax.text(0.5, 0.65,
        r'$\mathcal{L}(F) = (x - A_{obs}F)^T \Sigma^{-1} (x - A_{obs}F)$' + '\n'
        r'$+ \lambda F^T L_{signed} F + \rho \|F\|^2$',
        ha='center', va='center', fontsize=10)
ax.text(0.5, 0.40, 'Strictly convex\n→ unique F per patient', ha='center', va='center',
        fontsize=10, color='darkgreen', style='italic')
ax.text(0.5, 0.15, 'F: 38,148-D substrate-coordinate\nstate vector', ha='center',
        va='center', fontsize=9, color='gray')
ax.set_axis_off()

# Panel C: β/α decomposition diagram
ax = axes[2]
# Show conceptual 2D: log-PR axis + α-residual
ax.set_xlim(-1, 1.5); ax.set_ylim(-1, 1.5)
ax.arrow(0, 0, 1.2, 0, head_width=0.05, head_length=0.05, fc='#1f77b4', ec='#1f77b4', linewidth=2)
ax.text(1.25, -0.1, r'$\log PR$', fontsize=11, color='#1f77b4', fontweight='bold')
ax.text(0.7, -0.15, r'$\beta$', fontsize=12, color='#1f77b4')
# alpha residual
ax.arrow(0.9, 0, 0, 0.9, head_width=0.05, head_length=0.05, fc='#d62728', ec='#d62728', linewidth=2)
ax.text(0.95, 0.5, r'$\alpha$', fontsize=12, color='#d62728')
# F vector
ax.arrow(0, 0, 0.9, 0.9, head_width=0.06, head_length=0.06, fc='#2ca02c', ec='#2ca02c', linewidth=2.5)
ax.text(0.45, 0.55, r'$F_p$', fontsize=13, color='#2ca02c', fontweight='bold')
ax.set_title(r'C) $\beta$/$\alpha$ decomposition' + '\n'
              r'$F_p = \beta_p \cdot \log PR + \alpha_p$',
              fontsize=10)
ax.set_axis_off()

# Panel D: signed-basin extraction
ax = axes[3]
# Show a tiny basin
small_basin = [n for n in geom.nodes if node_type(n) == 'gene'][:8]
# Add a few connecting reaction nodes
H_basin = G.subgraph(small_basin)
if H_basin.number_of_edges() == 0:
    # Pick a more connected starting node
    for try_node in geom.nodes[:200]:
        nbrs = list(G.neighbors(try_node))
        if len(nbrs) >= 8:
            small_basin = [try_node] + nbrs[:8]
            H_basin = G.subgraph(small_basin)
            break

if H_basin.number_of_nodes() > 0:
    pos = nx.spring_layout(H_basin, seed=42)
    # Color by sign (alternate for illustration)
    pos_signs = ['+' if i % 2 == 0 else '-' for i in range(len(H_basin.nodes))]
    pos_color = '#d62728'; neg_color = '#1f77b4'
    nx.draw_networkx_edges(H_basin, pos, ax=ax, alpha=0.5, edge_color='gray', width=1.0)
    for i, n in enumerate(H_basin.nodes):
        color = pos_color if pos_signs[i] == '+' else neg_color
        nx.draw_networkx_nodes(H_basin, pos, nodelist=[n], node_color=color,
                                  node_size=180, ax=ax, alpha=0.85)
ax.set_title('D) Signed-basin extraction\nLargest connected sub-graph of\nsame-sign top-5% loading nodes',
             fontsize=10)
ax.set_axis_off()

fig.suptitle('Figure 1: GIZMO pipeline — substrate construction, MAP projection, β/α decomposition, basin extraction',
             fontsize=11)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure1_pipeline_schematic.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure1_pipeline_schematic.png')

print('\nDone — figures_v7 now has 9 figures (heatmaps, bar charts, line plots, networks, schematic).')
"""Cross-cohort chord diagrams + patient F-space embeddings for v7."""
import sys
sys.path.insert(0, '/home/jgardner/GIZMO')
import json
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.path import Path as MPLPath
from matplotlib.patches import PathPatch
import matplotlib.cm as cm

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS = REPO / 'benchmarks/results/unsupervised/zscored'
FIG_DIR = REPO / 'figures_v7'

# === FIGURE 10: Chord diagram of cross-cohort basin conservation ===
print('Building Figure 10 — chord diagram of basin conservation...', flush=True)
d = json.load(open(ZS / 'v7_basin_conservation.json'))

cohorts = ['CorEvitas_RA','GSE89408_RA','Gao_RA','Crohn','HMP2_IBD_CD',
            'Filbin_COVID','Su_COVID','KMPLOT_BRCA','TCGA_LUAD',
            'IDH_glioma','TCGA_IDH_glioma']
cohort_labels = {
    'CorEvitas_RA':'RA/CorEvitas','GSE89408_RA':'RA/GSE89408','Gao_RA':'RA/Gao',
    'Crohn':'IBD/Crohn','HMP2_IBD_CD':'IBD/HMP2',
    'Filbin_COVID':'COVID/Filbin','Su_COVID':'COVID/Su',
    'KMPLOT_BRCA':'BRCA/KMPLOT','TCGA_LUAD':'LUAD/TCGA',
    'IDH_glioma':'Glioma/Trautwein','TCGA_IDH_glioma':'Glioma/TCGA',
}
# Disease group colors
disease_colors = {
    'CorEvitas_RA':'#e377c2','GSE89408_RA':'#e377c2','Gao_RA':'#e377c2',
    'Crohn':'#ff7f0e','HMP2_IBD_CD':'#ff7f0e',
    'Filbin_COVID':'#1f77b4','Su_COVID':'#1f77b4',
    'KMPLOT_BRCA':'#2ca02c','TCGA_LUAD':'#d62728',
    'IDH_glioma':'#9467bd','TCGA_IDH_glioma':'#9467bd',
}

def categorize(genes):
    gs = set(genes)
    if any(g.startswith('RPL') or g.startswith('RPS') for g in gs): return 'Ribosome'
    if any(g.startswith('MT-') or g.startswith('NDUF') for g in gs): return 'OXPHOS'
    if any(g in {'B2M','CD3D','CD3E','CD3G','CD247','HLA-A','HLA-B','HLA-C'} for g in gs): return 'T cell/MHC'
    if any(g in {'ATM','BRCA1','BARD1','BLM','RAD51D','FANCD2','FANCM','FAAP100'} for g in gs): return 'DNA damage'
    if any(g.startswith('COL') or g in {'ACAN','BGN','DCN','VCAN','BMP1'} for g in gs): return 'Collagen/ECM'
    if any(g.startswith('IFT') or g.startswith('DNAAF') for g in gs): return 'Cilia'
    if any(g in {'C3','C4A','C4B','C5','C6','C7','C8A','C9'} for g in gs): return 'Complement'
    return 'Other'

CATEGORY_COLORS = {
    'Ribosome': '#1f77b4', 'OXPHOS': '#d62728', 'T cell/MHC': '#ff7f0e',
    'DNA damage': '#9467bd', 'Collagen/ECM': '#2ca02c', 'Cilia': '#8c564b',
    'Complement': '#e377c2', 'Other': '#7f7f7f',
}

# Position cohorts on circle
n = len(cohorts)
angles = np.linspace(np.pi/2, np.pi/2 - 2*np.pi, n, endpoint=False)
positions = {c: (np.cos(a), np.sin(a)) for c, a in zip(cohorts, angles)}

fig, ax = plt.subplots(figsize=(11, 11))
ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal'); ax.set_axis_off()

# Draw cohort markers and labels
for c in cohorts:
    px, py = positions[c]
    ax.scatter([px], [py], s=400, c=disease_colors[c], edgecolors='black', zorder=5)
    # Label outside the circle
    a = np.arctan2(py, px)
    lx = 1.25 * np.cos(a); ly = 1.25 * np.sin(a)
    ha = 'left' if lx > 0 else 'right'
    ax.text(lx, ly, cohort_labels[c], fontsize=10, ha=ha, va='center',
            fontweight='bold')

# Draw conserved basin pairs as curved chords
for pair in d['conserved_pairs']:
    if pair['jaccard'] < 0.30: continue
    cA = pair['cohort_A']; cB = pair['cohort_B']
    if cA not in positions or cB not in positions: continue
    cat = categorize(pair['shared_genes'])
    color = CATEGORY_COLORS.get(cat, '#7f7f7f')
    pA = positions[cA]; pB = positions[cB]
    # Bezier curve through origin (for chord effect)
    mid = (0, 0)
    verts = [(pA[0]*0.95, pA[1]*0.95), mid, (pB[0]*0.95, pB[1]*0.95)]
    codes = [MPLPath.MOVETO, MPLPath.CURVE3, MPLPath.CURVE3]
    path = MPLPath(verts, codes)
    width = pair['jaccard'] * 5
    patch = PathPatch(path, facecolor='none', edgecolor=color,
                       linewidth=width, alpha=0.4)
    ax.add_patch(patch)

# Legend
legend_handles = [Patch(facecolor=c, label=cat) for cat, c in CATEGORY_COLORS.items()
                    if cat != 'Other']
legend_handles += [Patch(facecolor='#7f7f7f', label='Other biology')]
ax.legend(handles=legend_handles, loc='center', fontsize=10,
            title='Conserved basin\nbiology category',
            title_fontsize=11)

ax.set_title('Figure 10: Cross-cohort basin conservation as a substrate-fixed activation vocabulary\n'
              '36 cross-cohort basin pairs at Jaccard ≥ 0.30 (≥ 5 shared genes)\n'
              'Line width ∝ Jaccard similarity; color = biology category',
              fontsize=11)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure10_chord_conservation.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure10_chord_conservation.png')


# === FIGURE 11: Patient F-space 2D projection across cohorts ===
print('Building Figure 11 — patient F-space embeddings...', flush=True)
all_F = []
all_pids = []
all_cohorts = []
for c in cohorts:
    fp = SNAPSHOT / f'stage3_F_{c}_edge_informed.npz'
    if not fp.exists():
        continue
    npz = np.load(fp, allow_pickle=True)
    F = npz['F']
    pids = [str(p) for p in npz['patient_ids']]
    all_F.append(F)
    all_pids.extend(pids)
    all_cohorts.extend([c] * len(pids))
F_combined = np.vstack(all_F)
print(f'  Combined F: {F_combined.shape}')

# Quick L2 normalize per patient
F_norm = F_combined / (np.linalg.norm(F_combined, axis=1, keepdims=True) + 1e-12)
# PCA-2 for visualization
pca = PCA(n_components=2, random_state=42)
emb = pca.fit_transform(F_norm)
print(f'  Top 2 PCs EVR: {pca.explained_variance_ratio_}')

fig, ax = plt.subplots(figsize=(11, 9))
for c in cohorts:
    mask = np.array(all_cohorts) == c
    ax.scatter(emb[mask, 0], emb[mask, 1], s=20, c=disease_colors[c],
                label=cohort_labels[c], alpha=0.6, edgecolors='none')
ax.set_xlabel(f'F-space PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)')
ax.set_ylabel(f'F-space PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
ax.set_title('Figure 11: All patients projected into shared F-space coordinates\n'
              f'n = {len(all_pids)} patients from 11 cohorts; top-2 PC visualization\n'
              'Each cohort occupies a distinct sub-region (substrate-coordinate disease separation)',
              fontsize=10)
ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
ax.grid(alpha=0.2)
plt.tight_layout()
plt.savefig(FIG_DIR / 'figure11_patient_f_space.png', dpi=150)
plt.close()
print(f'  Wrote {FIG_DIR}/figure11_patient_f_space.png')

print('\nDone — figures_v7 now has chord diagrams + patient embeddings.')
