"""TB progressor replication: project GSE94438 + GSE107994 baseline samples
onto Filbin-PC5 TNF axis; test cross-sectional discrimination.

This is the external-replication move for the cross-cohort axis-projection paper.
Tests a different question than the original GSE89403 cure outcome: does TNF
axis engagement at LTBI baseline predict who progresses to active TB?
"""
from __future__ import annotations
import sys, gzip, json, gc
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_1samp

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


def parse_series_matrix(path):
    """Return list of per-sample dicts indexed by sample order in series matrix."""
    titles = []; chars_rows = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles = [s.strip('"') for s in line.rstrip().split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'):
                chars_rows.append([s.strip('"') for s in line.rstrip().split('\t')[1:]])
    samples = []
    for i, title in enumerate(titles):
        d = {'title': title}
        for row in chars_rows:
            if i < len(row):
                if ':' in row[i]:
                    key, val = row[i].split(':', 1)
                    d[key.strip()] = val.strip()
        samples.append(d)
    return samples


def solve_F(geom, X_arr, feat_node_ids, sample_ids):
    from gizmo.inference.projection import solve_map, ModalitySetup
    mu = X_arr.mean(axis=0); sd = X_arr.std(axis=0) + 1e-9
    Xz = (X_arr - mu) / sd
    feat_cols = [(f'feat_{k}', geom.nid_idx[feat_node_ids[k]]) for k in range(len(feat_node_ids))]
    pdata = {sid: {f'feat_{k}': float(Xz[i, k]) for k in range(Xz.shape[1])}
             for i, sid in enumerate(sample_ids)}
    setup = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=pdata)
    F, _ = solve_map(geom, [setup], sample_ids)
    return F.astype(np.float32)


def load_gse94438(geom):
    """Load GSE94438 baseline (T0) samples; return X, feat_node_ids, metadata list."""
    print('\n[GSE94438] loading...', flush=True)
    series_path = REPO / 'data/cohorts/GSE94438/GSE94438_series_matrix.txt.gz'
    expr_path = REPO / 'data/cohorts/GSE94438/GSE94438_rawCounts_GeneNames_AllSamples.csv.gz'
    samples = parse_series_matrix(series_path)
    print(f'  Total samples in series matrix: {len(samples)}', flush=True)
    # Per-sample metadata: code → group / time.from.exposure / time.to.tb / site
    for s in samples:
        s['expr_col'] = f'X{s.get("code", "?")}'
    print(f'  Group: {Counter(s.get("group") for s in samples)}', flush=True)
    print(f'  time.from.exposure: {Counter(s.get("time.from.exposure.months") for s in samples)}',
          flush=True)
    # Filter to BASELINE samples (time.from.exposure == 0) with known group
    baseline_samples = [s for s in samples
                         if s.get('time.from.exposure.months') == '0'
                         and s.get('group') in ('case (TB)', 'Control')]
    print(f'  Baseline samples with case/control label: {len(baseline_samples)}', flush=True)
    print(f'    cases: {sum(1 for s in baseline_samples if s["group"] == "case (TB)")}',
          flush=True)
    print(f'    controls: {sum(1 for s in baseline_samples if s["group"] == "Control")}',
          flush=True)

    # Load expression matrix
    df = pd.read_csv(expr_path, index_col=0)
    if 'symbol' in df.columns:
        symbol_col = df['symbol']
        df = df.drop(columns=['symbol'])
    else:
        symbol_col = None
    print(f'  Expression: {df.shape}', flush=True)

    # log2 CPM normalize (since raw counts)
    lib_sizes = df.sum(axis=0)
    cpm = df.divide(lib_sizes, axis=1) * 1e6
    logcpm = np.log2(cpm + 1).astype(np.float32)
    del df, cpm; gc.collect()

    # Map ENSG → substrate
    ensg_ids = list(logcpm.index)
    feat_node_ids = [f'ENSG:{e}' for e in ensg_ids if f'ENSG:{e}' in geom.nid_idx]
    ensg_keep = [n.replace('ENSG:', '') for n in feat_node_ids]
    print(f'  Mapped features: {len(feat_node_ids)}', flush=True)

    # Build sample list — drop missing
    use_meta = []; use_cols = []
    for s in baseline_samples:
        col = s['expr_col']
        if col in logcpm.columns:
            use_meta.append(s); use_cols.append(col)
    print(f'  Samples with expression: {len(use_cols)}', flush=True)

    X = logcpm.loc[ensg_keep, use_cols].values.T.astype(float)
    del logcpm; gc.collect()
    return X, feat_node_ids, use_cols, use_meta


def load_gse107994(geom):
    """Load GSE107994 baseline samples."""
    print('\n[GSE107994 Leicester] loading...', flush=True)
    series_path = REPO / 'data/cohorts/GSE107994/GSE107994_series_matrix.txt.gz'
    expr_path = REPO / 'data/cohorts/GSE107994/GSE107994_edgeR_normalized_Leicester_with_progressor_longitudinal.xlsx'
    samples = parse_series_matrix(series_path)
    print(f'  Total samples: {len(samples)}', flush=True)
    print(f'  Group: {Counter(s.get("group") for s in samples)}', flush=True)
    print(f'  Timepoint: {Counter(s.get("timepoint_months") for s in samples)}', flush=True)

    # Filter to Baseline timepoint with LTBI groups (progressor vs non-progressor)
    baseline_samples = [s for s in samples
                         if s.get('timepoint_months') == 'Baseline']
    print(f'  Baseline samples: {len(baseline_samples)}', flush=True)
    print(f'    by group: {Counter(s.get("group") for s in baseline_samples)}', flush=True)

    # Load expression — sample columns are "Leicester_with_progressor_longitudinal_SampleN"
    df = pd.read_excel(expr_path, sheet_name=0)
    print(f'  Expression: {df.shape}', flush=True)
    # First 3 cols: Genes / Gene_name / Gene_biotype
    meta_cols = ['Genes', 'Gene_name', 'Gene_biotype']
    sample_cols = [c for c in df.columns if c not in meta_cols]
    # Build sample-N → expression column index
    df = df.set_index('Genes')
    df = df[sample_cols].astype(np.float32)
    print(f'  Sample columns count: {len(sample_cols)}', flush=True)

    # Map sample by NAME (series matrix has gaps in Sample numbering)
    title_to_meta = {s['title']: s for s in samples}
    use_meta = []; use_col_idx = []
    for col in df.columns:
        if col in title_to_meta and title_to_meta[col] in baseline_samples:
            use_meta.append(title_to_meta[col]); use_col_idx.append(col)
    print(f'  Mapped baseline samples (by name): {len(use_col_idx)}', flush=True)
    print(f'  Baseline group breakdown after mapping: '
          f'{Counter(s.get("group") for s in use_meta)}', flush=True)

    # Map ENSG → substrate
    ensg_ids = list(df.index)
    feat_node_ids = [f'ENSG:{e}' for e in ensg_ids if f'ENSG:{e}' in geom.nid_idx]
    ensg_keep = [n.replace('ENSG:', '') for n in feat_node_ids]
    print(f'  Mapped features: {len(feat_node_ids)}', flush=True)

    X = df.loc[ensg_keep, use_col_idx].values.T.astype(float)
    del df; gc.collect()
    return X, feat_node_ids, use_col_idx, use_meta


def project_and_test(name, X, feat_node_ids, sample_ids, meta, geom, lpn, filbin_pc5,
                      progressor_pred, label_field):
    """MAP solve, β/α decomp, project onto Filbin-PC5, test discrimination."""
    print(f'\n  [{name}] Solving MAP on {X.shape[0]} samples × {X.shape[1]} features...',
          flush=True)
    F = solve_F(geom, X, feat_node_ids, sample_ids)
    beta = F @ lpn
    alpha = (F - np.outer(beta, lpn)).astype(np.float32)
    del F; gc.collect()

    # Project
    proj_pc5 = alpha @ filbin_pc5
    del alpha; gc.collect()

    # Labels
    labels = np.array([1 if progressor_pred(m) else 0 for m in meta])
    print(f'  Progressor count: {labels.sum()}/{len(labels)}', flush=True)
    if labels.sum() < 2 or (1 - labels).sum() < 2:
        print(f'  ⚠ Too few samples in one class; skipping discrimination test', flush=True)
        return None

    # Stats
    mwu = mannwhitneyu(proj_pc5[labels == 1], proj_pc5[labels == 0])
    p_grp = float(mwu.pvalue)
    t, p_eng = ttest_1samp(proj_pc5, 0.0)
    p_eng = float(p_eng)
    proj_pos = float(proj_pc5[labels == 1].mean())
    proj_neg = float(proj_pc5[labels == 0].mean())
    print(f'  PC5 mean: progressors={proj_pos:+.3f}  non-progressors={proj_neg:+.3f}',
          flush=True)
    print(f'  MWU p_grp = {p_grp:.4f}', flush=True)
    print(f'  Cohort-wide PC5 one-sample t: p_engage = {p_eng:.4f}', flush=True)

    # β discrimination too
    mwu_b = mannwhitneyu(beta[labels == 1], beta[labels == 0])
    p_grp_b = float(mwu_b.pvalue)
    print(f'  β mean: progressors={beta[labels==1].mean():+.3f}  '
          f'non-progressors={beta[labels==0].mean():+.3f}, MWU p_β = {p_grp_b:.4f}', flush=True)

    # Random-axis null (300 trials)
    rng = np.random.default_rng(0)
    # Need alpha back for random-axis null — recompute
    F = solve_F(geom, X, feat_node_ids, sample_ids)
    beta2 = F @ lpn
    alpha2 = (F - np.outer(beta2, lpn)).astype(np.float32)
    del F; gc.collect()
    n_random = 300
    rand_pgrp = np.zeros(n_random)
    for s in range(n_random):
        v = rng.standard_normal(filbin_pc5.shape[0]).astype(np.float32)
        v /= (np.linalg.norm(v) + 1e-9)
        proj_v = alpha2 @ v
        rand_pgrp[s] = float(mannwhitneyu(proj_v[labels == 1], proj_v[labels == 0]).pvalue)
    p_random_below = float((rand_pgrp <= p_grp).sum() / n_random)
    p_random_below_05 = float((rand_pgrp < 0.05).sum() / n_random)
    del alpha2; gc.collect()
    print(f'  Random-axis null (n=300): {p_random_below*n_random:.0f}/{n_random} ≤ obs p; '
          f'{p_random_below_05*100:.1f}% achieve p<.05 (null calibration check)', flush=True)

    return {
        'cohort': name,
        'n_samples': int(X.shape[0]),
        'n_progressor': int(labels.sum()),
        'n_non_progressor': int((1 - labels).sum()),
        'p_grp_pc5': p_grp,
        'p_engage_pc5': p_eng,
        'mean_pc5_progressor': proj_pos,
        'mean_pc5_non_progressor': proj_neg,
        'p_grp_beta': p_grp_b,
        'random_axis_rate_below_obs': p_random_below,
        'random_axis_rate_below_05': p_random_below_05,
    }


def main():
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom
    lpn = (geom.log_pr / (np.linalg.norm(geom.log_pr) + 1e-9)).astype(np.float32)
    filbin_pc5 = np.load(RESULTS / 'filbin_pc5_eigvec.npy').astype(np.float32)
    print(f'Loaded Filbin-PC5: {filbin_pc5.shape}', flush=True)

    out = {}

    # ===== GSE94438 (Zak GC6-74) =====
    print('=' * 78, flush=True)
    print('Cohort 1: GSE94438 (Zak GC6-74) — multi-country progressor surveillance', flush=True)
    print('=' * 78, flush=True)
    X1, fnids1, sids1, meta1 = load_gse94438(geom)
    if len(meta1) >= 10:
        r = project_and_test('GSE94438', X1, fnids1, sids1, meta1, geom, lpn, filbin_pc5,
                              progressor_pred=lambda m: m.get('group') == 'case (TB)',
                              label_field='group')
        if r: out['GSE94438'] = r
    del X1, fnids1, sids1, meta1; gc.collect()

    # ===== GSE107994 (Leicester) =====
    print('\n' + '=' * 78, flush=True)
    print('Cohort 2: GSE107994 (Leicester) — UK contact progressor surveillance', flush=True)
    print('=' * 78, flush=True)
    X2, fnids2, sids2, meta2 = load_gse107994(geom)
    # Two subtests: LTBI progressor vs LTBI non-progressor (the headline test)
    # AND Active vs all non-active (for cross-sectional validation)
    if len(meta2) >= 10:
        # Test 1: LTBI Progressor vs LTBI (the headline)
        # Filter to LTBI samples only
        ltbi_idx = [i for i, m in enumerate(meta2) if m.get('group') in ('LTBI_Progressor', 'LTBI')]
        if len(ltbi_idx) >= 10:
            X2_ltbi = X2[ltbi_idx]
            sids2_ltbi = [sids2[i] for i in ltbi_idx]
            meta2_ltbi = [meta2[i] for i in ltbi_idx]
            print('\n  --- Sub-test 1: LTBI Progressor vs LTBI Non-progressor ---', flush=True)
            r1 = project_and_test('GSE107994_LTBI_progression', X2_ltbi, fnids2, sids2_ltbi,
                                    meta2_ltbi, geom, lpn, filbin_pc5,
                                    progressor_pred=lambda m: m.get('group') == 'LTBI_Progressor',
                                    label_field='group')
            if r1: out['GSE107994_LTBI_progression'] = r1

        # Test 2: Active TB vs all others (cross-sectional disease state)
        print('\n  --- Sub-test 2: Active TB vs LTBI+Control (full baseline) ---', flush=True)
        r2 = project_and_test('GSE107994_active_vs_all', X2, fnids2, sids2, meta2, geom, lpn,
                               filbin_pc5,
                               progressor_pred=lambda m: m.get('group') == 'Active_TB',
                               label_field='group')
        if r2: out['GSE107994_active_vs_all'] = r2

    out_path = RESULTS / 'tb_progressor_replication.json'
    with open(out_path, 'w') as fh: json.dump(out, fh, indent=2)
    print(f'\nSaved: {out_path}', flush=True)

    # Summary
    print('\n' + '=' * 78, flush=True)
    print('REPLICATION SUMMARY', flush=True)
    print('=' * 78, flush=True)
    for k, v in out.items():
        print(f'\n  {k}:', flush=True)
        print(f'    n_progressor / n_non_progressor: {v["n_progressor"]} / {v["n_non_progressor"]}',
              flush=True)
        print(f'    Filbin-PC5 p_grp: {v["p_grp_pc5"]:.4f}', flush=True)
        print(f'    β p_grp: {v["p_grp_beta"]:.4f}', flush=True)
        print(f'    Random-axis null rate: {v["random_axis_rate_below_obs"]:.4f}',
              flush=True)


if __name__ == '__main__':
    main()
