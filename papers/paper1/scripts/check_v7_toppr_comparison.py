"""TopPR (static PageRank, no data) vs F-α-PC on v7 substrate.

For each cohort × curated truth gene set, rank substrate genes by:
  (a) TopPR = static substrate PageRank (no data input)
  (b) Best F-α-PC |loading|
And compute fold-enrichment of truth genes at top-K = {50, 200, 500}.
"""
import sys, json
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0,'/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS = REPO / 'benchmarks/results/unsupervised/zscored'

mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr

# Substrate gene nodes
sub_genes = []
gene_to_idx = {}
for nid in geom.nodes:
    attrs = mg.graph.nodes.get(nid, {})
    if attrs.get('node_type') == 'gene':
        sym = attrs.get('symbol') or attrs.get('name') or nid.replace('symbol:', '')
        if sym:
            sub_genes.append(sym)
            gene_to_idx[sym] = geom.nid_idx[nid]
sub_gene_set = set(sub_genes)
n_genes = len(sub_genes)
print(f'Substrate gene nodes: {n_genes}', flush=True)

# Static PageRank = log_pr (already substrate-coordinate)
# TopPR ranking: by descending substrate PageRank, restricted to gene nodes
pr_values = np.array([log_pr[gene_to_idx[g]] for g in sub_genes])
pr_rank = np.argsort(-pr_values)
pr_ranked_genes = [sub_genes[i] for i in pr_rank]

# Curated truth
cur = pd.read_csv(REPO / 'data/curation/v7_cohort_key_genes.tsv', sep='\t')

# Cohorts
COHORTS = ['CorEvitas_RA','GSE89408_RA','Gao_RA','Crohn','HMP2_IBD_CD',
            'Filbin_COVID','Su_COVID','KMPLOT_BRCA','TCGA_LUAD',
            'IDH_glioma','TCGA_IDH_glioma']

# For each cohort: compute F-α-PC top-K vs TopPR top-K vs random
results = []
for cohort in COHORTS:
    truth = [g for g in cur[cur['cohort']==cohort]['key_gene_symbol'] if g in sub_gene_set]
    if len(truth) < 2: continue

    # Load F
    fp = None
    for v in ['_edge_informed', '_node_informed', '']:
        p = SNAPSHOT / f'stage3_F_{cohort}{v}.npz'
        if p.exists(): fp = p; break
    if fp is None: continue
    F = np.load(fp, allow_pickle=True)['F']
    if F.shape[1] != len(geom.nodes): continue

    # F-α-PC: get top-K from best α-PC × sign
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm=x.mean(); xv=x.var()+1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm)*(x-xm)).mean(axis=1, keepdims=True)
    beta = (cov/xv).ravel()
    alpha = F_unit - Fm - beta[:,None]*(x-xm)[None,:]
    pca = PCA(n_components=7, random_state=0); pca.fit(alpha)

    # Best PC × sign by truth-recovery (oracle — for the comparison)
    best_alpha_recall = {50: 0, 200: 0, 500: 0}
    for pc_k in range(7):
        loadings = pca.components_[pc_k]
        gene_scores = np.array([loadings[gene_to_idx[g]] for g in sub_genes])
        for sign_mult in [1, -1]:
            signed = sign_mult * gene_scores
            rank_order = np.argsort(-signed)
            ranked_genes = [sub_genes[i] for i in rank_order]
            for K in [50, 200, 500]:
                topK = set(ranked_genes[:K])
                recall = sum(1 for t in truth if t in topK)
                if recall > best_alpha_recall[K]:
                    best_alpha_recall[K] = recall

    for K in [50, 200, 500]:
        topK_pr = set(pr_ranked_genes[:K])
        pr_recall = sum(1 for t in truth if t in topK_pr)
        # Expected null random recovery: K/n_genes * len(truth)
        null_recall = K / n_genes * len(truth)

        fold_alpha = best_alpha_recall[K] / max(null_recall, 0.01)
        fold_pr = pr_recall / max(null_recall, 0.01)

        results.append({
            'cohort': cohort, 'K': K, 'truth_n': len(truth),
            'F_alpha_recall': best_alpha_recall[K],
            'TopPR_recall': pr_recall,
            'null_expected': round(null_recall, 2),
            'F_alpha_fold': round(fold_alpha, 2),
            'TopPR_fold': round(fold_pr, 2),
        })

df = pd.DataFrame(results)
print('\n=== TopPR vs F-α-PC truth-gene recovery (best α-PC × sign per cohort) ===')
print(df.to_string(index=False))
df.to_csv(ZS / 'v7_toppr_comparison.tsv', sep='\t', index=False)
print(f'\nWrote {ZS}/v7_toppr_comparison.tsv')

# Summary: at each K, what's the median fold-enrichment for each method?
print('\n=== Summary by K ===')
for K in [50, 200, 500]:
    sub = df[df['K']==K]
    print(f'K={K}: median F-α-PC fold = {sub["F_alpha_fold"].median():.2f},  '
          f'median TopPR fold = {sub["TopPR_fold"].median():.2f}')
