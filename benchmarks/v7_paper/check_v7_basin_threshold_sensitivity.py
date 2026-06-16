"""Basin 5% threshold sensitivity: re-extract at {1, 2, 5, 10, 15}% quantiles
and report cross-cohort Jaccard conservation rate at each."""
import sys, json
import numpy as np
import networkx as nx
from pathlib import Path
sys.path.insert(0, '/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'

mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr
sub_node_set = set(geom.nodes)
G = nx.Graph()
G.add_nodes_from(geom.nodes)
for u, v in mg.graph.edges():
    if u in sub_node_set and v in sub_node_set:
        G.add_edge(u, v)

cohorts = ['CorEvitas_RA','GSE89408_RA','Gao_RA','Crohn','HMP2_IBD_CD',
            'Filbin_COVID','Su_COVID','KMPLOT_BRCA','TCGA_LUAD',
            'IDH_glioma','TCGA_IDH_glioma']

def extract_basins(F, threshold_quantile):
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm = x.mean(); xv = x.var() + 1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm) * (x - xm)).mean(axis=1, keepdims=True)
    beta = (cov / xv).ravel()
    alpha = F_unit - Fm - beta[:, None] * (x - xm)[None, :]
    pca = PCA(n_components=7, random_state=0); pca.fit(alpha)
    basins = []
    for pc_k in range(7):
        comp = pca.components_[pc_k]
        abs_load = np.abs(comp)
        threshold = np.quantile(abs_load, 1 - threshold_quantile)
        for sign_label, sign in [('+', 1), ('-', -1)]:
            mask = (np.sign(comp) == sign) & (abs_load >= threshold)
            keep = [geom.nodes[i] for i in range(len(comp)) if mask[i]]
            H = G.subgraph(keep)
            ccs = sorted(nx.connected_components(H), key=len, reverse=True)
            if not ccs or len(ccs[0]) < 5: continue
            genes = set()
            for nid in ccs[0]:
                attrs = mg.graph.nodes.get(nid, {})
                if attrs.get('node_type') == 'gene':
                    sym = attrs.get('symbol') or attrs.get('name')
                    if sym: genes.add(sym)
            basins.append((pc_k+1, sign_label, genes))
    return basins

print(f'{"Threshold":<10} | {"Total basins":<14} | {"Cross-cohort pairs at Jac≥0.30":<35}')
print('-' * 70)
for q in [0.01, 0.02, 0.05, 0.10, 0.15]:
    all_basins_at_q = []
    for c in cohorts:
        for variant in ['_edge_informed', '_node_informed', '']:
            fp = SNAPSHOT / f'stage3_F_{c}{variant}.npz'
            if fp.exists():
                F = np.load(fp, allow_pickle=True)['F']
                if F.shape[1] == len(geom.nodes):
                    for pc, sign, genes in extract_basins(F, q):
                        if len(genes) >= 3:
                            all_basins_at_q.append((c, pc, sign, genes))
                    break
    # Compute pairwise Jaccard ≥ 0.30
    n_pairs = 0
    for i in range(len(all_basins_at_q)):
        cA, pA, sA, gA = all_basins_at_q[i]
        for j in range(i+1, len(all_basins_at_q)):
            cB, pB, sB, gB = all_basins_at_q[j]
            if cA == cB: continue
            inter = gA & gB
            jac = len(inter) / max(1, len(gA | gB))
            if jac >= 0.30 and len(inter) >= 5:
                n_pairs += 1
    print(f'{q*100:.0f}%        | {len(all_basins_at_q):<14} | {n_pairs:<35}')
