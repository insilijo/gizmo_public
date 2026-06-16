"""Pathway-membership-matched null for interpretability test.

Critique: degree-matched null doesn't control for annotation density. Reactome-curated
genes (the "key genes" in curation) are over-represented at substrate nodes that
participate in many Reactome reactions. A stricter null samples from the same
Reactome pathway-membership distribution as the curated set.

For each cohort × α-PC × sign cell:
  1. Curated key gene set: count their Reactome reaction participation degree
  2. Null sampling: sample N genes from substrate matched to this Reactome-reaction-degree
  3. Recompute AUROC, repeat 1000x, empirical p-value
"""
import sys, json, time, gzip
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, '/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS = REPO / 'benchmarks/results/unsupervised/zscored'
CUR = REPO / 'data/curation/v7_cohort_key_genes.tsv'

print('Loading substrate...', flush=True)
mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr

# Build per-gene Reactome reaction-participation degree
print('Computing per-gene Reactome reaction degree...', flush=True)
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
print(f'  {len(gene_reaction_degree)} genes with ≥1 Reactome reaction', flush=True)

# Substrate gene symbols
sub_genes = set()
sym_to_node = {}
for nid in geom.nodes:
    attrs = mg.graph.nodes.get(nid, {})
    if attrs.get('node_type') == 'gene':
        sym = attrs.get('symbol') or attrs.get('name') or nid.replace('symbol:', '')
        if sym:
            sub_genes.add(sym)
            sym_to_node[sym] = nid
print(f'  {len(sub_genes)} substrate gene symbols', flush=True)

# Stratify substrate genes by Reactome-degree decile
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

# Load curated key genes
print('Loading curated key genes...', flush=True)
cur_df = pd.read_csv(CUR, sep='\t')
cohort_to_truth = defaultdict(list)
for _, r in cur_df.iterrows():
    cohort_to_truth[r['cohort']].append(r['key_gene_symbol'])

# For each cohort, decompose F-v5-canonical → α-PCA → 7 components
def auroc(scores, labels):
    """AUROC of binary labels vs continuous scores."""
    from sklearn.metrics import roc_auc_score
    if len(set(labels)) < 2: return np.nan
    return roc_auc_score(labels, scores)


print('Computing pathway-matched null per cohort × α-PC × sign...', flush=True)
COHORTS = ['CorEvitas_RA','GSE89408_RA','Gao_RA','Crohn','HMP2_IBD_CD',
            'Filbin_COVID','Su_COVID','KMPLOT_BRCA','TCGA_LUAD',
            'IDH_glioma','TCGA_IDH_glioma']
rng = np.random.default_rng(42)
results = []

for cohort in COHORTS:
    truth_genes = [g for g in cohort_to_truth[cohort] if g in sub_genes]
    if len(truth_genes) < 2:
        continue
    # Compute truth's decile distribution
    truth_deciles = [gene_to_decile[g] for g in truth_genes]
    decile_counts = np.bincount(truth_deciles, minlength=10)
    # Sample pool: from each decile, the available genes excluding truth
    decile_pools = {d: [g for g in decile_to_genes[d] if g not in truth_genes]
                     for d in range(10)}

    # Load F + decompose
    fp = SNAPSHOT / f'stage3_F_{cohort}_edge_informed.npz'
    if not fp.exists(): continue
    F = np.load(fp, allow_pickle=True)['F']
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm = x.mean(); xv = x.var() + 1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm) * (x - xm)).mean(axis=1, keepdims=True)
    beta = (cov / xv).ravel()
    alpha = F_unit - Fm - beta[:, None] * (x - xm)[None, :]
    pca = PCA(n_components=7, random_state=0); pca.fit(alpha)

    # Build gene → node index map for THIS cohort's geometry
    truth_node_idx = []
    for g in truth_genes:
        nid = sym_to_node.get(g)
        if nid in geom.nid_idx:
            truth_node_idx.append(geom.nid_idx[nid])

    # For each PC × sign, compute observed AUROC + null distribution
    for pc_k in range(5):
        loadings = pca.components_[pc_k]
        # Build gene-only score vector (substrate gene nodes only)
        gene_to_score = {}
        for g in sub_genes:
            nid = sym_to_node[g]
            if nid in geom.nid_idx:
                gene_to_score[g] = loadings[geom.nid_idx[nid]]

        for sign_label, sign_mult in [('+', 1), ('-', -1)]:
            # Test signed score on gene labels (truth = 1 if in curated)
            score_arr = np.array([sign_mult * gene_to_score.get(g, 0) for g in sub_gene_list])
            label_arr = np.array([1 if g in truth_genes else 0 for g in sub_gene_list])
            obs_auroc = auroc(score_arr, label_arr)

            # Null: sample pathway-matched random gene set, recompute AUROC
            null_aurocs = []
            for _ in range(1000):
                null_set = []
                for d in range(10):
                    pool = decile_pools[d]
                    n_d = decile_counts[d]
                    if n_d == 0: continue
                    if len(pool) < n_d: continue
                    null_set.extend(rng.choice(pool, size=n_d, replace=False))
                null_label_arr = np.array([1 if g in null_set else 0 for g in sub_gene_list])
                if sum(null_label_arr) < 2:
                    continue
                null_aurocs.append(auroc(score_arr, null_label_arr))
            null_aurocs = np.array(null_aurocs)
            n_better = (null_aurocs >= obs_auroc).sum()
            p_emp = (n_better + 1) / (len(null_aurocs) + 1)
            results.append({
                'cohort': cohort, 'pc': pc_k+1, 'sign': sign_label,
                'truth_n': len(truth_genes), 'obs_auroc': obs_auroc,
                'null_median': float(np.median(null_aurocs)),
                'p_pathway_matched': p_emp,
            })

# Save
df = pd.DataFrame(results)
df.to_csv(ZS / 'v7_pathway_matched_null.tsv', sep='\t', index=False)
print(f'Wrote {ZS}/v7_pathway_matched_null.tsv  ({len(df)} cells)')

# Summary
n_p05 = (df['p_pathway_matched'] < 0.05).sum()
n_cohorts_passing = df.groupby('cohort').apply(lambda g: (g['p_pathway_matched'] < 0.05).any()).sum()
print(f'\nCells significant at p_pathway-matched < 0.05: {n_p05} / {len(df)}')
print(f'Cohorts passing at least one cell: {n_cohorts_passing} / 11')
print(f'\nDistribution by cohort:')
print(df.groupby('cohort').apply(lambda g: f"{(g['p_pathway_matched']<0.05).sum()} of {len(g)} cells significant").to_string())
