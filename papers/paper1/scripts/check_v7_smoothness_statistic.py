"""Compute smoothness statistic s_cohort across 11 cohorts.

s_cohort = cos(v_data_top1_SVD, v_low_eigen_substrate_Laplacian)
where:
  v_data = top-1 right singular vector of cohort's substrate-mapped input matrix
  v_low_eigen = first non-trivial Laplacian eigenvector of substrate restricted to observed nodes

Validation: does s predict whether canonical preprocessing works?
"""
import sys, json
import numpy as np
import pandas as pd
import openpyxl, gzip
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
sys.path.insert(0, '/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import TruncatedSVD

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
ZS = REPO / 'benchmarks/results/unsupervised/zscored'

print('Loading substrate...', flush=True)
mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)

# Build undirected adjacency for Laplacian
print('Building substrate adjacency...', flush=True)
N = len(geom.nodes)
node_idx = {n: i for i, n in enumerate(geom.nodes)}
rows, cols = [], []
for u, v in mg.graph.edges():
    if u in node_idx and v in node_idx:
        i, j = node_idx[u], node_idx[v]
        rows.append(i); cols.append(j)
        rows.append(j); cols.append(i)
data = np.ones(len(rows))
A = csr_matrix((data, (rows, cols)), shape=(N, N))
# Symmetric normalized Laplacian L = I - D^{-1/2} A D^{-1/2}
deg = np.array(A.sum(axis=1)).flatten()
deg_safe = np.where(deg > 0, 1.0/np.sqrt(deg), 0)
D_inv_sqrt = csr_matrix((deg_safe, (np.arange(N), np.arange(N))), shape=(N, N))
L = csr_matrix(np.eye(N)) - D_inv_sqrt @ A @ D_inv_sqrt
print(f'  Substrate {N} nodes, computing low Laplacian eigenvector...', flush=True)
# First non-trivial eigenvector = 2nd smallest eigenvalue
# Use eigsh with which='SM' (smallest magnitude)
vals, vecs = eigsh(L, k=3, which='SM')
v_low_eigen = vecs[:, 1]  # Skip trivial constant
v_low_eigen = v_low_eigen / np.linalg.norm(v_low_eigen)
print(f'  Done. Eigenvalue: {vals[1]:.6f}', flush=True)

# For each cohort, compute s_cohort
COHORTS = ['CorEvitas_RA','GSE89408_RA','Gao_RA','Crohn','HMP2_IBD_CD',
            'Filbin_COVID','Su_COVID','KMPLOT_BRCA','TCGA_LUAD',
            'IDH_glioma','TCGA_IDH_glioma']
results = []
for cohort in COHORTS:
    # Load F
    fp_e = SNAPSHOT / f'stage3_F_{cohort}_edge_informed.npz'
    fp_p = SNAPSHOT / f'stage3_F_{cohort}.npz'
    fp = fp_e if fp_e.exists() else fp_p
    if not fp.exists():
        print(f'  {cohort}: no F file', flush=True)
        continue
    F = np.load(fp, allow_pickle=True)['F']
    if F.shape[1] != N:
        print(f'  {cohort}: shape mismatch', flush=True)
        continue

    # F top-1 right singular vector
    svd = TruncatedSVD(n_components=1, random_state=0)
    svd.fit(F)
    v_data = svd.components_[0]
    v_data = v_data / np.linalg.norm(v_data)

    # s = |cos(v_data, v_low_eigen)|
    s = abs(float(np.dot(v_data, v_low_eigen)))
    results.append({'cohort': cohort, 's_cohort': s, 'F_shape': str(F.shape)})
    print(f'  {cohort:<20} s = {s:.3f}  F: {F.shape}', flush=True)

df = pd.DataFrame(results)
df.to_csv(ZS / 'v7_smoothness_statistic.tsv', sep='\t', index=False)
print(f'\nWrote {ZS}/v7_smoothness_statistic.tsv')

# Look at relationship: s vs degree-null pass / pathway-null pass / Gao failure
print('\n=== Validation: does s predict canonical preprocessing success? ===')
# Cohorts that pass degree-null: 10/11 (all except Gao_RA)
# Cohorts that fail under pathway null: KMPLOT, CorEvitas_RA, GSE89408_RA, Gao_RA
# Gao_RA AUROC 0.406 (below chance) is the catastrophic failure
gao_s = df[df['cohort']=='Gao_RA']['s_cohort'].iloc[0]
print(f'Gao_RA s = {gao_s:.3f} (the cohort that fails canonical at AUROC 0.406)')
print(f'TCGA_IDH s = {df[df["cohort"]=="TCGA_IDH_glioma"]["s_cohort"].iloc[0]:.3f}')
print(f'  (the cohort that succeeds under z-score after canonical mis-aligns 2HG)')
print()
print('s distribution:')
print(df.sort_values('s_cohort'))
