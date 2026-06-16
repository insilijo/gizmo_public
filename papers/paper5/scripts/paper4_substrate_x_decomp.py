"""Substrate-aware × decomposition-method 2D comparison.

For each combination of:
  axis 1: substrate-awareness (raw expression vs substrate-aware α)
  axis 2: decomposition method (PCA, NMF, FastICA)

Apply to GSE65391 SLE; report SLE-vs-HC discrimination (best PC AUC) +
biological interpretability (top-30 gene overlap with curated SLE pathways).

Tests whether substrate-aware decomposition consistently outperforms non-
substrate-aware regardless of the decomposition method choice (PCA-only is
not load-bearing; substrate-awareness is).
"""
from __future__ import annotations
import sys, gc
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, NMF, FastICA, FactorAnalysis
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
CACHE = RESULTS / 'cohort_alpha_cache'
SCAN = RESULTS / 'cross_indication_scan'
COHORT_SLE = REPO / 'data/cohorts/GSE65391_SLE'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


ISG_GENES = {'STAT1','STAT2','IRF1','IRF7','IRF9','ISG15','ISG20','IFI6',
             'IFI27','IFI35','IFI44','IFI44L','IFIT1','IFIT2','IFIT3','IFITM1',
             'IFITM2','IFITM3','MX1','MX2','OAS1','OAS2','OAS3','OASL',
             'USP18','RSAD2','DDX60','DDX58','IFIH1','PARP9','PARP12','PARP14',
             'EPSTI1','CMPK2','XAF1','SAMD9','SAMD9L','LY6E','LAMP3','BST2',
             'HERC5','HERC6','CXCL10','CXCL11'}
FOLATE_GENES = {'DHFR','TYMS','MTHFR','MTR','MTHFD1','MTHFD2','MTHFD1L','SHMT1',
                'SHMT2','FPGS','GGH','GART','ATIC','AICDA','CTH','BHMT','CBS',
                'MTAP','MAT2A','AHCY','GLDC','GCSH','AMT','DLD'}
PLASMA_GENES = {'PSMB8','PSMB9','PSMB10','PSME1','PSME2','IGHG1','IGHG2','IGHG3',
                'IGHG4','IGHA1','IGHA2','IGHM','IGHD','IGHE','IGKC','IGLC1',
                'MS4A1','CD19','CD27','CD38','TNFSF13B','XBP1','PRDM1','IRF4',
                'PDIA3','PDIA4','PDIA6','RPS6','RPL5','RPL10A','RPL18A','RPL27A',
                'EEF1A1','EEF1G'}
COMPLEMENT_GENES = {'C1QA','C1QB','C1QC','C1S','C1R','C2','C3','C4A','C4B','C5',
                    'CFB','CFD','CFH','CFI','CFP','C3AR1','C5AR1','CR1','CR2',
                    'CD55','CD46','CD59','MASP1','MASP2','MBL2'}
ALL_BIOLOGY = (ISG_GENES | FOLATE_GENES | PLASMA_GENES | COMPLEMENT_GENES)


def parse_sle_raw_expression(geom):
    from sle_random_pc_null import parse_series_matrix, parse_platform_probes
    probes, expr, sample_cols, meta = parse_series_matrix(
        COHORT_SLE / 'GSE65391_series_matrix.txt.gz')
    probe_to_sym = parse_platform_probes(COHORT_SLE / 'GSE65391_family.soft.gz')
    probe_to_node = {p: f'symbol:{probe_to_sym[p]}' for p in probes
                      if p in probe_to_sym
                      and f'symbol:{probe_to_sym[p]}' in geom.nid_idx}
    unique_nodes = sorted(set(probe_to_node.values()))
    probe_idx = {p: i for i, p in enumerate(probes)}
    node_to_probes = {nid: [] for nid in unique_nodes}
    for p, nid in probe_to_node.items():
        node_to_probes[nid].append(probe_idx[p])
    X = np.zeros((expr.shape[1], len(unique_nodes)), dtype=np.float32)
    for j, nid in enumerate(unique_nodes):
        idxs = node_to_probes[nid]
        if idxs: X[:, j] = expr[idxs, :].mean(axis=0)
    is_sle = np.array([(m.get('disease state') == 'SLE') for m in meta])
    valid = ~np.isnan(X).any(axis=1)
    X = X[valid]; is_sle = is_sle[valid]
    return X, is_sle, unique_nodes


def score_components(components, X, y, node_id_resolver, n_top=30):
    """For each component, AUC + interpretability. Returns dicts."""
    aucs = []
    interps = []
    n_comp = components.shape[0]
    for k in range(n_comp):
        v = components[k]
        proj = X @ v
        try:
            a = roc_auc_score(y.astype(int), proj)
            a = float(max(a, 1 - a))
        except Exception:
            a = 0.5
        aucs.append(a)
        # Top-30 gene features
        top_idx = np.argsort(-np.abs(v))[:n_top]
        syms = set()
        for j in top_idx:
            s = node_id_resolver(j)
            if s: syms.add(s)
        if syms:
            interps.append(len(syms & ALL_BIOLOGY) / len(syms))
        else:
            interps.append(0.0)
    best_k = int(np.argmax(aucs))
    return {
        'best_AUC': float(max(aucs)),
        'mean_AUC': float(np.mean(aucs)),
        'best_PC': best_k + 1,
        'best_interp': float(interps[best_k]),
        'mean_interp': float(np.mean(interps)),
    }


def main():
    print('Loading substrate...', flush=True)
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom; mg = audit_mod.mg
    lpn = (geom.log_pr / (np.linalg.norm(geom.log_pr) + 1e-9)).astype(np.float32)

    print('Parsing raw SLE expression...', flush=True)
    X_sle, is_sle, unique_nodes = parse_sle_raw_expression(geom)
    print(f'  Shape: {X_sle.shape}', flush=True)
    mu = X_sle.mean(axis=0); sd = X_sle.std(axis=0) + 1e-9
    Xz = (X_sle - mu) / sd

    # Resolver for raw expression (gene-only)
    def resolve_raw(j):
        nid = unique_nodes[j]
        return nid.replace('symbol:', '').strip() if nid.startswith('symbol:') else None

    # Substrate-aware α for SLE (from cached F)
    F_sle = np.load(SCAN / 'SLE_F.npy')
    beta = F_sle @ lpn
    # α residual = F - outer(β, lpn) — too large to store; instead fit on Wang RA α
    # for decomposition derivation (matches the framework approach for axis library).

    # Wang RA α (the framework's "source cohort" for axis derivation)
    a0 = np.load(CACHE / 'Wang_RA_alpha_t0.npy')
    print(f'  Wang α: {a0.shape}', flush=True)

    # Substrate-aware projections need gene-only top contributors
    # Pre-build gene-node mask
    is_gene_node = np.array([
        bool(geom.nodes[j].startswith('symbol:')
             or mg.graph.nodes.get(geom.nodes[j], {}).get('gene_symbols'))
        for j in range(len(geom.nodes))
    ], dtype=bool)
    gene_indices = np.where(is_gene_node)[0]

    def resolve_substrate(j):
        nid = geom.nodes[j]
        attrs = mg.graph.nodes.get(nid, {})
        syms = attrs.get('gene_symbols', []) or []
        if syms: return syms[0]
        if nid.startswith('symbol:'): return nid.replace('symbol:', '').strip()
        return None

    def score_substrate_method(name, components, X_proj, y, n_top=30):
        """Substrate-aware components project SLE F via α @ v = F @ v - β·lpn·v.
        Top contributors restricted to gene nodes for fair comparison."""
        aucs, interps = [], []
        for k in range(components.shape[0]):
            v = components[k]
            proj = F_sle @ v - beta * float(lpn @ v)
            try:
                a = roc_auc_score(y.astype(int), proj)
                a = float(max(a, 1 - a))
            except Exception:
                a = 0.5
            aucs.append(a)
            # Top-30 GENE contributors only
            gene_v = v[gene_indices]
            top_gene_local = np.argsort(-np.abs(gene_v))[:n_top]
            syms = set()
            for j_local in top_gene_local:
                j = gene_indices[j_local]
                s = resolve_substrate(j)
                if s: syms.add(s)
            if syms:
                interps.append(len(syms & ALL_BIOLOGY) / len(syms))
            else:
                interps.append(0.0)
        best_k = int(np.argmax(aucs))
        return {
            'best_AUC': float(max(aucs)),
            'mean_AUC': float(np.mean(aucs)),
            'best_PC': best_k + 1,
            'best_interp': float(interps[best_k]),
            'mean_interp': float(np.mean(interps)),
        }

    N_COMP = 15
    results = []

    print('\n' + '='*100, flush=True)
    print('2D comparison — substrate-awareness × decomposition method', flush=True)
    print('='*100, flush=True)

    # Row 1: PCA
    print('\nPCA on raw expression...', flush=True)
    pca_raw = PCA(n_components=N_COMP, svd_solver='randomized',
                   random_state=0).fit(Xz)
    r1 = score_components(pca_raw.components_, Xz, is_sle, resolve_raw)
    r1.update({'substrate_aware': False, 'method': 'PCA'})
    results.append(r1)

    print('PCA on substrate-aware α (Wang RA derivation)...', flush=True)
    pca_sub = PCA(n_components=N_COMP, svd_solver='randomized',
                   random_state=0).fit(a0)
    r2 = score_substrate_method('PCA-on-α', pca_sub.components_, F_sle, is_sle)
    r2.update({'substrate_aware': True, 'method': 'PCA'})
    results.append(r2)

    # Row 2: NMF
    print('NMF on raw expression (shifted positive)...', flush=True)
    X_pos = Xz - Xz.min(axis=0) + 0.01
    nmf_raw = NMF(n_components=N_COMP, random_state=0, max_iter=300,
                   init='nndsvd').fit(X_pos)
    r3 = score_components(nmf_raw.components_, Xz, is_sle, resolve_raw)
    r3.update({'substrate_aware': False, 'method': 'NMF'})
    results.append(r3)

    print('NMF on substrate-aware α (Wang, shifted positive)...', flush=True)
    a0_pos = a0 - a0.min(axis=0) + 0.01
    nmf_sub = NMF(n_components=N_COMP, random_state=0, max_iter=300,
                   init='nndsvd').fit(a0_pos)
    r4 = score_substrate_method('NMF-on-α', nmf_sub.components_, F_sle, is_sle)
    r4.update({'substrate_aware': True, 'method': 'NMF'})
    results.append(r4)

    # Row 3: FastICA
    print('FastICA on raw expression...', flush=True)
    ica_raw = FastICA(n_components=N_COMP, random_state=0, max_iter=300,
                       whiten='unit-variance').fit(Xz)
    r5 = score_components(ica_raw.components_, Xz, is_sle, resolve_raw)
    r5.update({'substrate_aware': False, 'method': 'FastICA'})
    results.append(r5)

    print('FastICA on substrate-aware α...', flush=True)
    ica_sub = FastICA(n_components=N_COMP, random_state=0, max_iter=300,
                       whiten='unit-variance').fit(a0)
    r6 = score_substrate_method('FastICA-on-α', ica_sub.components_, F_sle, is_sle)
    r6.update({'substrate_aware': True, 'method': 'FastICA'})
    results.append(r6)

    # Report 2D table
    print('\n' + '='*100, flush=True)
    print('2D comparison: substrate-awareness × decomposition method', flush=True)
    print('='*100, flush=True)
    print(f'  {"method":<10}{"substrate-aware?":<18}{"Best AUC":>10}'
          f'{"Mean AUC":>10}{"Best interp":>14}{"Mean interp":>14}', flush=True)
    for r in results:
        sub = 'YES' if r['substrate_aware'] else 'no'
        print(f'  {r["method"]:<10}{sub:<18}{r["best_AUC"]:>10.3f}'
              f'{r["mean_AUC"]:>10.3f}{r["best_interp"]:>14.3f}'
              f'{r["mean_interp"]:>14.3f}', flush=True)

    # Per-method delta (substrate-aware vs raw)
    print('\nSubstrate-awareness improvement (substrate-aware − raw) per method:',
          flush=True)
    print(f'  {"method":<10}{"ΔBest AUC":>12}{"ΔMean AUC":>12}'
          f'{"ΔBest interp":>16}', flush=True)
    methods = ['PCA', 'NMF', 'FastICA']
    for m in methods:
        raw = [r for r in results if r['method'] == m
               and not r['substrate_aware']][0]
        sub = [r for r in results if r['method'] == m
               and r['substrate_aware']][0]
        d_auc = sub['best_AUC'] - raw['best_AUC']
        d_mean = sub['mean_AUC'] - raw['mean_AUC']
        d_interp = sub['best_interp'] - raw['best_interp']
        print(f'  {m:<10}{d_auc:>+12.3f}{d_mean:>+12.3f}{d_interp:>+16.3f}',
              flush=True)

    df = pd.DataFrame(results)
    df.to_csv(SCAN / 'substrate_x_decomp_2D.csv', index=False)
    print(f'\n  Saved: {SCAN / "substrate_x_decomp_2D.csv"}', flush=True)


if __name__ == '__main__':
    main()
