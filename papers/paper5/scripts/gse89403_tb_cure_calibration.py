"""GSE89403 TB treatment cure calibration.

Tests whether β shifts under CURATIVE 6-month antibiotic regimen for TB.
TB patients followed at DX (diagnosis) → Day 7 → Week 4 → Week 24 (post-cure).
Treatment outcome labels: Definite Cure / Probable / Possible / Not Cured.

Critical second test (after HCV DAA) of the "drugs don't move β" claim:
  - DX TB vs Healthy Controls → does β work here?
  - DX → Week 24 in Definite Cure patients → does β shift toward healthy?

The 6-month antibiotic regimen is a different kind of cure (eradicates the
pathogen by killing it, vs DAA which inhibits viral replication). If β stays
flat HERE too, the "β = immune-metabolic set point that persists past cure"
claim generalizes to TB.
"""
from __future__ import annotations
import sys, json, gzip
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
COHORT_DIR = REPO / 'data/cohorts/GSE89403_TB_treat'

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import (
    build_biochem_subgraph, solve_map, decompose_beta_alpha, ModalitySetup,
)
from sklearn.metrics import roc_auc_score


def parse_series_matrix(path):
    """Return dict {sample_code: {disease, treatmentresult, time, subject}}."""
    out = {}
    chars = []
    titles = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_characteristics_ch1'):
                vals = [s.strip('"') for s in line.rstrip().split('\t')[1:]]
                chars.append(vals)
            elif line.startswith('!Sample_title'):
                titles = [s.strip('"').split('-')[0] for s in line.rstrip().split('\t')[1:]]
    if not chars: return out
    # chars[0]=tissue, [1]=sample_code, [2]=subject, [3]=disease, [4]=treatmentresult, [5]=time
    for i in range(len(titles)):
        sample_code = chars[1][i].replace('sample_code: ', '')
        # Take the canonical record per sample_code (first occurrence, since technical reps share metadata)
        if sample_code in out: continue
        out[sample_code] = {
            'subject': chars[2][i].replace('subject: ', ''),
            'disease': chars[3][i].replace('disease state: ', ''),
            'treatment_result': chars[4][i].replace('treatmentresult: ', ''),
            'time': chars[5][i].replace('time: ', ''),
        }
    return out


def main():
    print('=' * 80)
    print('GSE89403 TB cure paired calibration')
    print('=' * 80)

    annots = parse_series_matrix(COHORT_DIR / 'GSE89403_series_matrix.txt.gz')
    print(f'Annotated samples: {len(annots)}')
    print(f'Treatment results: {Counter(a["treatment_result"] for a in annots.values())}')
    print(f'Disease states: {Counter(a["disease"] for a in annots.values())}')

    # Load expression
    print('\nLoading expression matrix (454 sample columns, ENSG rows)...')
    df = pd.read_csv(COHORT_DIR / 'GSE89403_log2ExpGeneNames_AllSamples.csv.gz',
                      compression='gzip', index_col=0)
    print(f'Expression: {df.shape}')
    # Drop 'symbol' column if present
    if 'symbol' in df.columns:
        gene_symbols = df['symbol']
        df = df.drop(columns=['symbol'])
        print(f'After dropping symbol col: {df.shape}')

    # Map ENSG → substrate
    mg = read_json(str(REPO / 'data/processed/human_full_rhea_full/graph.json'))
    geom = build_biochem_subgraph(mg, hub_cap=500)
    ensg_ids = list(df.index)
    mapped_ensg = {e: f'ENSG:{e}' for e in ensg_ids if f'ENSG:{e}' in geom.nid_idx}
    print(f'ENSG mapped to substrate: {len(mapped_ensg)}/{len(ensg_ids)} = '
          f'{100*len(mapped_ensg)/len(ensg_ids):.1f}%')
    kept = [e for e in ensg_ids if e in mapped_ensg]
    feat_node_ids = [mapped_ensg[e] for e in kept]
    expr = df.loc[kept].T  # samples × genes
    print(f'Final: {expr.shape}')

    # Build sample → annotation lookup
    sample_cols = list(expr.index)
    print(f'Total samples in expression matrix: {len(sample_cols)}')
    annot_match = {s: annots.get(s) for s in sample_cols}
    matched = sum(1 for v in annot_match.values() if v is not None)
    print(f'Matched annotations: {matched}/{len(sample_cols)}')

    # Filter: TB subjects with Definite Cure + paired DX/Week24
    tb_dc = {s: a for s, a in annot_match.items()
             if a and a['disease'] == 'TB Subjects' and a['treatment_result'] == 'Definite Cure'}
    # Find patients with both DX and week_24
    by_subj = {}
    for sid, a in tb_dc.items():
        by_subj.setdefault(a['subject'], {})[a['time']] = sid
    paired_subjs = [s for s, tps in by_subj.items() if 'DX' in tps and 'week_24' in tps]
    print(f'\nDefinite Cure TB patients with DX + Week 24 paired: {len(paired_subjs)}')

    # Healthy controls
    hd_sids = [s for s, a in annot_match.items() if a and a['disease'] == 'Healthy Controls']
    print(f'Healthy Controls: {len(hd_sids)}')

    # Baseline cohort: DX TB (Definite Cure) vs Healthy
    dx_sids = [by_subj[s]['DX'] for s in paired_subjs]
    cohort = dx_sids + hd_sids
    print(f'\n[Baseline] TB DX (Definite Cure n={len(dx_sids)}) vs Healthy (n={len(hd_sids)})')
    X = expr.loc[cohort].values
    mu = X.mean(axis=0); sd = X.std(axis=0) + 1e-9
    Xz = (X - mu) / sd
    feat_cols = [(f'feat_{k}', geom.nid_idx[feat_node_ids[k]]) for k in range(len(feat_node_ids))]
    data = {sid: {f'feat_{k}': float(Xz[i, k]) for k in range(len(feat_node_ids))}
            for i, sid in enumerate(cohort)}
    main = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=data)
    print(f'  Solving MAP on {len(cohort)} samples × {len(feat_node_ids)} features...')
    F_base, _ = solve_map(geom, [main], cohort)
    beta, _, alpha_pc, pca_base = decompose_beta_alpha(F_base, geom.log_pr, n_components=5)
    is_tb = np.array([1 if a.get('disease', '') == 'TB Subjects' else 0
                       for a in [annot_match[s] for s in cohort]])
    def auc(s): a = roc_auc_score(is_tb, s); return max(a, 1-a)
    beta_auc = auc(beta)
    a_aucs = [auc(alpha_pc[:, k]) for k in range(5)]
    print(f'  β AUC (TB DX vs Healthy): {beta_auc:.3f}')
    print(f'  α-PC AUCs: {[f"{a:.3f}" for a in a_aucs]}')

    # Paired DX → Week 24
    print(f'\n[Paired DX→Week24] {len(paired_subjs)} patients (all Definite Cure)')
    dx_paired = [by_subj[s]['DX'] for s in paired_subjs]
    w24_paired = [by_subj[s]['week_24'] for s in paired_subjs]
    cohort = dx_paired + w24_paired
    X = expr.loc[cohort].values
    mu = X.mean(axis=0); sd = X.std(axis=0) + 1e-9
    Xz = (X - mu) / sd
    data = {sid: {f'feat_{k}': float(Xz[i, k]) for k in range(len(feat_node_ids))}
            for i, sid in enumerate(cohort)}
    main = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=data)
    print(f'  Solving MAP on {len(cohort)} samples × {len(feat_node_ids)} features...')
    F, _ = solve_map(geom, [main], cohort)
    beta_p, _, alpha_pc_p, pca = decompose_beta_alpha(F, geom.log_pr, n_components=5)
    n = len(paired_subjs)
    beta_pre = beta_p[:n]; beta_post = beta_p[n:]
    dbeta = beta_post - beta_pre
    alpha_pre = alpha_pc_p[:n]; alpha_post = alpha_pc_p[n:]
    dalpha = alpha_post - alpha_pre

    F_pre = F[:n]; F_post = F[n:]
    dF = F_post - F_pre
    dF_norm = np.linalg.norm(dF, axis=1)
    dF_node_mean = dF.mean(axis=0)

    print(f'  ‖ΔF‖ per patient: mean={dF_norm.mean():.2f}, σ={dF_norm.std():.2f}')
    print(f'  β: DX={beta_pre.mean():+.3f} → Wk24={beta_post.mean():+.3f}, '
          f'Δβ={dbeta.mean():+.4f}, σ_Δβ={dbeta.std():.4f}')
    for k in range(5):
        v = pca.explained_variance_ratio_[k]
        print(f'  α-PC{k+1}: DX={alpha_pre[:,k].mean():+.3f} → Wk24={alpha_post[:,k].mean():+.3f}  '
              f'Δmean={dalpha[:,k].mean():+.3f}, σ={dalpha[:,k].std():.3f}  EV={v:.2%}')

    # Top shifted nodes
    top_idx = np.argsort(-np.abs(dF_node_mean))[:15]
    print(f'  Top shifted nodes (DX → Week 24 cure):')
    for j in top_idx:
        nid = geom.nodes[j]
        a = mg.graph.nodes.get(nid, {})
        name = (a.get('name', '') or nid)[:50]
        print(f'    {name:<52} Δ={dF_node_mean[j]:+.3f}')

    out = {
        'cohort': 'GSE89403_TB_treatment',
        'baseline': {
            'n_TB_DX': len(dx_sids), 'n_HD': len(hd_sids),
            'n_features': len(feat_node_ids),
            'beta_AUC': float(beta_auc),
            'alpha_PC_AUCs': [float(a) for a in a_aucs],
        },
        'paired_DX_Wk24': {
            'n_paired': len(paired_subjs),
            'mean_dF_norm': float(dF_norm.mean()),
            'std_dF_norm': float(dF_norm.std()),
            'mean_dbeta': float(dbeta.mean()),
            'std_dbeta': float(dbeta.std()),
            'beta_DX_mean': float(beta_pre.mean()),
            'beta_Wk24_mean': float(beta_post.mean()),
            'alpha_PC_dmean': [float(dalpha[:,k].mean()) for k in range(5)],
            'alpha_PC_dstd':  [float(dalpha[:,k].std())  for k in range(5)],
        }
    }
    out_json = RESULTS / 'gse89403_tb_cure_calibration.json'
    def to_native(o):
        if isinstance(o, dict): return {k: to_native(v) for k, v in o.items()}
        if isinstance(o, list): return [to_native(v) for v in o]
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        return o
    with open(out_json, 'w') as fh:
        json.dump(to_native(out), fh, indent=2)
    print(f'\nSaved: {out_json}')


if __name__ == '__main__':
    main()
