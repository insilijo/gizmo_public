"""v7 Phase 5b: cross-cohort basin conservation across 11 cohorts.

For each cohort under v5-canonical preprocessing:
  1. Load F + decompose α + PCA top-7
  2. Extract basin (largest connected same-sign top-5% loading sub-graph) per α-PC × sign
  3. Get gene members of each basin

Then identify conserved biochemical neighborhoods by computing pairwise basin-gene
Jaccard overlap across cohorts.
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
import numpy as np
import networkx as nx
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS_DIR = REPO / 'benchmarks/results/unsupervised/zscored'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph

COHORTS = [
    'CorEvitas_RA', 'Crohn', 'Filbin_COVID', 'GSE89408_RA',
    'Gao_RA', 'HMP2_IBD_CD', 'IDH_glioma', 'KMPLOT_BRCA',
    'Su_COVID', 'TCGA_IDH_glioma', 'TCGA_LUAD',
]

print('Loading substrate...', flush=True)
mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr
sub_node_set = set(geom.nodes)

G = nx.Graph()
G.add_nodes_from(geom.nodes)
for u, v in mg.graph.edges():
    if u in sub_node_set and v in sub_node_set:
        G.add_edge(u, v)
print(f'Substrate: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges', flush=True)

node_to_idx = {n: i for i, n in enumerate(geom.nodes)}
def node_sym(nid):
    attrs = mg.graph.nodes.get(nid, {})
    if attrs.get('node_type') == 'gene':
        return attrs.get('symbol') or attrs.get('name') or nid.replace('symbol:', '')
    return None


def extract_basin(loadings, sign, thresh_q=0.95, min_size=5):
    abs_load = np.abs(loadings)
    threshold = np.quantile(abs_load, thresh_q)
    sign_match = (np.sign(loadings) == sign) & (abs_load >= threshold)
    keep = [geom.nodes[i] for i in range(len(loadings)) if sign_match[i]]
    H = G.subgraph(keep)
    if H.number_of_nodes() == 0: return None
    ccs = sorted(nx.connected_components(H), key=len, reverse=True)
    if not ccs or len(ccs[0]) < min_size: return None
    return ccs[0]


def cohort_basins(cohort):
    """Return list of (cohort, pc_idx, sign, basin_genes) for all basins."""
    # Find v5-canonical F
    paths = [
        SNAPSHOT / f'stage3_F_{cohort}_edge_informed.npz',
        SNAPSHOT / f'stage3_F_{cohort}.npz',
    ]
    f_path = next((p for p in paths if p.exists()), None)
    if f_path is None:
        print(f'  {cohort}: no F found', flush=True)
        return []
    npz = np.load(f_path, allow_pickle=True)
    F = npz['F']
    if F.shape[1] != len(geom.nodes):
        print(f'  {cohort}: F shape {F.shape} mismatch (substrate {len(geom.nodes)})',
              flush=True)
        return []

    # Decompose
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm = x.mean(); xv = x.var() + 1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm) * (x - xm)).mean(axis=1, keepdims=True)
    beta = (cov / xv).ravel()
    alpha = F_unit - Fm - beta[:, None] * (x - xm)[None, :]
    pca = PCA(n_components=7, random_state=0)
    pca.fit(alpha)

    basins = []
    for pc_k in range(7):
        for sign_label, sign in [('+', 1), ('-', -1)]:
            basin = extract_basin(pca.components_[pc_k], sign)
            if basin is None: continue
            genes = []
            for nid in basin:
                sym = node_sym(nid)
                if sym: genes.append(sym)
            basins.append({
                'cohort': cohort,
                'pc': pc_k + 1, 'sign': sign_label,
                'basin_size': len(basin),
                'n_genes': len(genes),
                'genes': sorted(genes),
            })
    return basins


print('\n=== Extracting basins across cohorts ===\n', flush=True)
all_basins = []
for c in COHORTS:
    print(f'{c}...', flush=True)
    cb = cohort_basins(c)
    print(f'  {len(cb)} basins', flush=True)
    all_basins.extend(cb)
print(f'\nTotal basins: {len(all_basins)}', flush=True)

# Build basin gene sets
basin_genes = []
for b in all_basins:
    if b['n_genes'] >= 3:
        basin_genes.append((b['cohort'], f'PC{b["pc"]}{b["sign"]}',
                             set(b['genes'])))
print(f'Basins with ≥3 gene members: {len(basin_genes)}', flush=True)

# Pairwise Jaccard between cross-cohort basins
print('\n=== Cross-cohort basin conservation (Jaccard ≥ 0.30, n_genes ≥ 5) ===\n',
      flush=True)
conserved = []
for i in range(len(basin_genes)):
    cA, lA, gA = basin_genes[i]
    for j in range(i+1, len(basin_genes)):
        cB, lB, gB = basin_genes[j]
        if cA == cB: continue
        inter = gA & gB
        union = gA | gB
        jac = len(inter) / len(union) if union else 0
        if jac >= 0.30 and len(inter) >= 5:
            conserved.append({
                'cohort_A': cA, 'basin_A': lA, 'cohort_B': cB, 'basin_B': lB,
                'n_genes_A': len(gA), 'n_genes_B': len(gB),
                'n_shared': len(inter), 'jaccard': jac,
                'shared_genes': sorted(inter),
            })

# Sort by jaccard descending, show
conserved_sorted = sorted(conserved, key=lambda x: -x['jaccard'])
print(f'Found {len(conserved)} cross-cohort basin pairs at Jaccard ≥ 0.30 + ≥5 shared genes\n',
      flush=True)
print(f'{"Cohort A":<18} {"PC":<6} {"Cohort B":<18} {"PC":<6} '
      f'{"|A|":<5} {"|B|":<5} {"∩":<4} {"Jac":<6} sample shared',
      flush=True)
print('-' * 110, flush=True)
for c in conserved_sorted[:50]:
    sample = ','.join(c['shared_genes'][:5])
    print(f'{c["cohort_A"]:<18} {c["basin_A"]:<6} {c["cohort_B"]:<18} '
          f'{c["basin_B"]:<6} {c["n_genes_A"]:<5} {c["n_genes_B"]:<5} '
          f'{c["n_shared"]:<4} {c["jaccard"]:.3f}  {sample}', flush=True)

out = ZS_DIR / 'v7_basin_conservation.json'
out.write_text(json.dumps({
    'cohorts': COHORTS,
    'all_basins': all_basins,
    'conserved_pairs': conserved_sorted,
    'compute_seconds': time.time(),  # rough
}, indent=2))
print(f'\nWrote {out}', flush=True)
