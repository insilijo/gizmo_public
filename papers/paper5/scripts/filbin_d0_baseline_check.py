"""Filbin Day 0 baseline β + α-PC AUCs (COVID+ vs COVID−).

Quick standalone check to confirm whether β degenerates in Olink-only data
OR whether β is genuinely the disease-state axis that doesn't shift longitudinally.
"""
from __future__ import annotations
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

REPO = Path('/home/jgardner/GIZMO')
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import (
    build_biochem_subgraph, solve_map, decompose_beta_alpha, ModalitySetup,
)
from sklearn.metrics import roc_auc_score


CACHE = Path.home() / '.cache' / 'filbin_mgh_covid'


def main():
    print('Filbin Day 0 baseline (COVID+ vs COVID−) β + α-PC check')
    print('=' * 70)

    t2a = pd.read_excel(CACHE / 'Suppl_T2_Olink_Assays_NPX.xlsx',
                         sheet_name='2A-Olink-Assay', header=1)
    oid_to_sym = dict(zip(t2a['OlinkID'].astype(str), t2a['Assay'].astype(str)))

    clin = pd.read_excel(CACHE / 'Clinical_Metadata.xlsx', sheet_name='Subject-level metadata')
    pid_covid = dict(zip(clin['Public ID'].astype(int).astype(str),
                          clin['COVID'].astype(int)))

    ol = pd.read_excel(CACHE / 'Olink_Proteomics.xlsx')
    ol = ol[ol.Day == 0].copy()
    ol['pid'] = ol['Public ID'].astype(str).str.replace(r'_D0$', '', regex=True)

    data = {}
    y = {}
    for _, row in ol.iterrows():
        pid = row['pid']
        if pid not in pid_covid: continue
        d = {}
        for col, v in row.items():
            if not isinstance(col, str) or not col.startswith('OID'): continue
            sym = oid_to_sym.get(col)
            if not sym or sym == 'nan': continue
            try: vf = float(v)
            except (TypeError, ValueError): continue
            if pd.isna(vf): continue
            d[sym] = max(d.get(sym, -1e9), vf)
        if d:
            data[pid] = d
            y[pid] = 'active' if pid_covid[pid] == 1 else 'control'

    sids = sorted(data.keys())
    print(f'n={len(sids)}  labels={Counter(y[s] for s in sids)}')

    mg = read_json(str(REPO / 'data/processed/human_full_rhea_full/graph.json'))
    geom = build_biochem_subgraph(mg, hub_cap=500)
    print(f'Substrate: {len(geom.nodes)} nodes')

    all_genes = sorted({k for sid in sids for k in data[sid].keys()})
    node_ids = {k: f'symbol:{k}' for k in all_genes if f'symbol:{k}' in geom.nid_idx}
    kept = [k for k in all_genes if k in node_ids]
    print(f'Mapped genes: {len(kept)}/{len(all_genes)}')

    X = np.zeros((len(sids), len(kept)))
    for i, sid in enumerate(sids):
        for j, k in enumerate(kept):
            X[i, j] = data[sid].get(k, 0.0)
    mu = X.mean(axis=0); sd = X.std(axis=0) + 1e-9
    Xz = (X - mu) / sd
    feat_cols = [(f'feat_{k}', geom.nid_idx[node_ids[kept[k]]]) for k in range(len(kept))]
    pdata = {sid: {f'feat_{k}': float(Xz[i, k]) for k in range(len(kept))}
             for i, sid in enumerate(sids)}
    main = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=pdata)

    print('Solving MAP...')
    F, _ = solve_map(geom, [main], sids)
    beta, _, alpha_pc, pca = decompose_beta_alpha(F, geom.log_pr, n_components=5)

    is_active = np.array([1 if y[s] == 'active' else 0 for s in sids])
    def auc(s): a = roc_auc_score(is_active, s); return max(a, 1 - a)

    beta_auc = auc(beta)
    alpha_aucs = [auc(alpha_pc[:, k]) for k in range(5)]
    print(f'\nβ AUC (COVID+ vs COVID-): {beta_auc:.3f}')
    print(f'α-PC AUCs: {[f"{a:.3f}" for a in alpha_aucs]}')
    print(f'EV: {[f"{v:.2%}" for v in pca.explained_variance_ratio_]}')

    # Group means
    print(f'\nβ means: active={beta[is_active.astype(bool)].mean():+.3f}  '
          f'control={beta[~is_active.astype(bool)].mean():+.3f}')


if __name__ == '__main__':
    main()
