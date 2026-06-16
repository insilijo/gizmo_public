"""Re-run pathway-null for HMP2 only (was silently skipped — plain .npz, no edge_informed)."""
import sys, json, gzip
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
sys.path.insert(0, '/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS = REPO / 'benchmarks/results/unsupervised/zscored'

mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr

# Alias map + Reactome-degree
gene_reaction_degree = defaultdict(int)
for u, v in mg.graph.edges():
    u_type = mg.graph.nodes.get(u, {}).get('node_type', '')
    v_type = mg.graph.nodes.get(v, {}).get('node_type', '')
    if u_type == 'gene' and v_type == 'reaction':
        sym = mg.graph.nodes[u].get('symbol') or mg.graph.nodes[u].get('name') or u.replace('symbol:', '')
        if sym: gene_reaction_degree[sym] += 1
    elif v_type == 'gene' and u_type == 'reaction':
        sym = mg.graph.nodes[v].get('symbol') or mg.graph.nodes[v].get('name') or v.replace('symbol:', '')
        if sym: gene_reaction_degree[sym] += 1

sub_genes = set()
sym_to_node = {}
for nid in geom.nodes:
    attrs = mg.graph.nodes.get(nid, {})
    if attrs.get('node_type') == 'gene':
        sym = attrs.get('symbol') or attrs.get('name') or nid.replace('symbol:', '')
        if sym:
            sub_genes.add(sym); sym_to_node[sym] = nid

gene_rxn_degs = np.array([gene_reaction_degree.get(g, 0) for g in sub_genes])
sub_gene_list = sorted(sub_genes)
gene_to_decile = {}
for g in sub_genes:
    d = gene_reaction_degree.get(g, 0)
    decile = min(9, int(np.searchsorted(np.quantile(gene_rxn_degs, np.linspace(0, 1, 11)), d, side='right') - 1))
    gene_to_decile[g] = decile
decile_to_genes = defaultdict(list)
for g, dec in gene_to_decile.items():
    decile_to_genes[dec].append(g)

cur = pd.read_csv(REPO / 'data/curation/v7_cohort_key_genes.tsv', sep='\t')
truth_genes = [g for g in cur[cur['cohort']=='HMP2_IBD_CD']['key_gene_symbol'] if g in sub_genes]
print(f'HMP2 truth genes in substrate: {truth_genes}')

# Load HMP2 F (plain version)
fp = SNAPSHOT / 'stage3_F_HMP2_IBD_CD.npz'
F = np.load(fp, allow_pickle=True)['F']
print(f'HMP2 F shape: {F.shape}, substrate size: {len(geom.nodes)}')
if F.shape[1] != len(geom.nodes):
    print('SHAPE MISMATCH; cannot run')
    sys.exit()

F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
x = log_pr; xm = x.mean(); xv = x.var() + 1e-12
Fm = F_unit.mean(axis=1, keepdims=True)
cov = ((F_unit - Fm) * (x - xm)).mean(axis=1, keepdims=True)
beta = (cov / xv).ravel()
alpha = F_unit - Fm - beta[:, None] * (x - xm)[None, :]
pca = PCA(n_components=7, random_state=0); pca.fit(alpha)

truth_deciles = [gene_to_decile[g] for g in truth_genes]
decile_counts = np.bincount(truth_deciles, minlength=10)
decile_pools = {d: [g for g in decile_to_genes[d] if g not in truth_genes] for d in range(10)}

rng = np.random.default_rng(42)
results = []
for pc_k in range(5):
    loadings = pca.components_[pc_k]
    gene_to_score = {}
    for g in sub_genes:
        nid = sym_to_node[g]
        if nid in geom.nid_idx:
            gene_to_score[g] = loadings[geom.nid_idx[nid]]
    for sign_label, sign_mult in [('+', 1), ('-', -1)]:
        score_arr = np.array([sign_mult * gene_to_score.get(g, 0) for g in sub_gene_list])
        label_arr = np.array([1 if g in truth_genes else 0 for g in sub_gene_list])
        obs_auroc = roc_auc_score(label_arr, score_arr)
        null_aurocs = []
        for _ in range(1000):
            null_set = []
            for d in range(10):
                pool = decile_pools[d]
                n_d = decile_counts[d]
                if n_d == 0 or len(pool) < n_d: continue
                null_set.extend(rng.choice(pool, size=n_d, replace=False))
            null_label_arr = np.array([1 if g in null_set else 0 for g in sub_gene_list])
            if sum(null_label_arr) < 2: continue
            null_aurocs.append(roc_auc_score(null_label_arr, score_arr))
        null_aurocs = np.array(null_aurocs)
        p_emp = ((null_aurocs >= obs_auroc).sum() + 1) / (len(null_aurocs) + 1)
        results.append({'cohort':'HMP2_IBD_CD','pc':pc_k+1,'sign':sign_label,
                         'truth_n':len(truth_genes),'obs_auroc':obs_auroc,
                         'null_median':float(np.median(null_aurocs)),
                         'p_pathway_matched':p_emp})

df_new = pd.DataFrame(results)
print('\nHMP2 results:')
print(df_new.to_string(index=False))
n_p05 = (df_new['p_pathway_matched']<0.05).sum()
print(f'\nHMP2 passes: {n_p05}/10 cells at p<0.05')

# Merge into existing TSV
existing = pd.read_csv(ZS / 'v7_pathway_matched_null.tsv', sep='\t')
combined = pd.concat([existing, df_new], ignore_index=True)
combined.to_csv(ZS / 'v7_pathway_matched_null.tsv', sep='\t', index=False)
print(f'\nUpdated {ZS}/v7_pathway_matched_null.tsv with HMP2 cells')

# Final summary
n_cohorts_passing = combined.groupby('cohort').apply(lambda g: (g['p_pathway_matched']<0.05).any()).sum()
print(f'\nFINAL: {n_cohorts_passing}/11 cohorts pass pathway-null at p<0.05')
