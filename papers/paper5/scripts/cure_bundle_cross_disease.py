"""Test cure-axis bundle in additional cohorts:

(A) Wang RA trajectory (n=38 paired DX→followup MTX response): does the bundle
    coordinate in RA, or is the coordination TB-specific?
(B) Filbin COVID trajectory (n=40 paired D0→D3): does the bundle coordinate
    in acute viral disease?
(C) GSE94438 progressors vs non-progressors: does bundle direction differ
    between LTBI patients who progress to active TB vs those who don't?

Bundle = {Filbin-PC3 mitochondrial, Filbin-PC5 TNF/integrin, Wang-PC5 STAT1 IFN}
Originally bootstrap-validated in GSE89403 TB DX→Wk24 (n=76):
  PC3↔PC5: ρ=+0.64; PC3↔WP5: ρ=+0.75; PC5↔WP5: ρ=+0.54
Replicated in GSE94438 multi-timepoint (n=69):
  PC3↔PC5: ρ=+0.60; PC3↔WP5: ρ=+0.73; PC5↔WP5: ρ=+0.77
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from collections import Counter

import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, ttest_1samp

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
CACHE = RESULTS / 'cohort_alpha_cache'


def project_dpc(alpha_t0, alpha_t1, v):
    return (alpha_t1 @ v - alpha_t0 @ v).astype(np.float64)


def bundle_correlations(dpc_fp3, dpc_fp5, dpc_wp5, n_boot=1000, seed=0):
    pairs = [
        ('ΔFilbin-PC3', 'ΔFilbin-PC5', dpc_fp3, dpc_fp5, +0.639),
        ('ΔFilbin-PC3', 'ΔWang-PC5',   dpc_fp3, dpc_wp5, +0.749),
        ('ΔFilbin-PC5', 'ΔWang-PC5',   dpc_fp5, dpc_wp5, +0.537),
    ]
    n = len(dpc_fp3)
    rng = np.random.default_rng(seed)
    out = {}
    print(f'  {"pair":<32}{"obs ρ":>9}{"95% CI":>22}{"GSE89403 ρ":>14}{"verdict":>14}',
          flush=True)
    for name_a, name_b, x, y, orig_rho in pairs:
        rho, _ = spearmanr(x, y)
        boots = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            r, _ = spearmanr(x[idx], y[idx])
            if not np.isnan(r): boots.append(r)
        if len(boots) < 100:
            ci = (float('nan'), float('nan'))
            verdict = '?'
        else:
            ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))
            same_sign = np.sign(rho) == np.sign(orig_rho)
            ci_nonzero = (ci[0] > 0) or (ci[1] < 0)
            ci_overlap = (orig_rho >= ci[0]) and (orig_rho <= ci[1])
            if ci_overlap and ci_nonzero and same_sign: verdict = '✓ strong'
            elif same_sign and ci_nonzero:               verdict = '~ direction'
            elif same_sign:                              verdict = '? underpowered'
            else:                                         verdict = '✗ different'
        ci_str = f'[{ci[0]:+.2f}, {ci[1]:+.2f}]'
        print(f'  {name_a} ↔ {name_b:<14}{rho:>+9.3f}{ci_str:>22}{orig_rho:>+14.3f}{verdict:>14}',
              flush=True)
        out[f'{name_a}__{name_b}'] = {
            'obs_rho': float(rho), 'CI95': list(ci),
            'gse89403_rho': float(orig_rho), 'verdict': verdict,
        }
    return out


def main():
    lpn = np.load(CACHE / 'lpn.npy').astype(np.float32)
    filbin_pc3 = np.load(RESULTS / 'filbin_pc3_eigvec.npy').astype(np.float32)
    filbin_pc5 = np.load(RESULTS / 'filbin_pc5_eigvec.npy').astype(np.float32)
    wang_pc5 = np.load(RESULTS / 'wang_pc5_eigvec.npy').astype(np.float32)

    out = {}

    # ============================================================
    # (A) Wang RA trajectory bundle test
    # ============================================================
    print('='*80, flush=True)
    print('(A) Wang RA paired DX→followup (n=38, MTX response)', flush=True)
    print('='*80, flush=True)
    wang_t0 = np.load(CACHE / 'Wang_RA_alpha_t0.npy')
    wang_t1 = np.load(CACHE / 'Wang_RA_alpha_t1.npy')
    wang_labels = np.load(CACHE / 'Wang_RA_labels.npy')
    print(f'  α matrix: {wang_t0.shape}; labels: Response={wang_labels.sum()}, '
          f'No-Resp={(1-wang_labels).sum()}', flush=True)
    d_fp3_w = project_dpc(wang_t0, wang_t1, filbin_pc3)
    d_fp5_w = project_dpc(wang_t0, wang_t1, filbin_pc5)
    d_wp5_w = project_dpc(wang_t0, wang_t1, wang_pc5)
    print('  Bundle pairwise correlations in Wang RA trajectory:', flush=True)
    out['Wang_RA'] = bundle_correlations(d_fp3_w, d_fp5_w, d_wp5_w)

    # Test bundle direction discrimination of outcome
    print('\n  Bundle direction → Wang outcome discrimination:', flush=True)
    for name, dpc in [('ΔFilbin-PC3', d_fp3_w), ('ΔFilbin-PC5', d_fp5_w), ('ΔWang-PC5', d_wp5_w)]:
        mwu = mannwhitneyu(dpc[wang_labels == 1], dpc[wang_labels == 0])
        rmean = dpc[wang_labels == 1].mean(); nrmean = dpc[wang_labels == 0].mean()
        print(f'    {name:<18}  R={rmean:+.2f}  NR={nrmean:+.2f}  MWU p={mwu.pvalue:.4f}',
              flush=True)

    # ============================================================
    # (B) Filbin COVID trajectory bundle test
    # ============================================================
    print('\n' + '='*80, flush=True)
    print('(B) Filbin COVID paired D0→D3 (n=40, Improved vs Worsened)', flush=True)
    print('='*80, flush=True)
    filbin_t0 = np.load(CACHE / 'Filbin_alpha_t0.npy')
    filbin_t1 = np.load(CACHE / 'Filbin_alpha_t1.npy')
    filbin_labels = np.load(CACHE / 'Filbin_labels.npy')
    print(f'  α matrix: {filbin_t0.shape}; labels: Improved={filbin_labels.sum()}, '
          f'Worsened={(1-filbin_labels).sum()}', flush=True)
    d_fp3_f = project_dpc(filbin_t0, filbin_t1, filbin_pc3)
    d_fp5_f = project_dpc(filbin_t0, filbin_t1, filbin_pc5)
    d_wp5_f = project_dpc(filbin_t0, filbin_t1, wang_pc5)
    print('  Bundle pairwise correlations in Filbin COVID trajectory:', flush=True)
    out['Filbin'] = bundle_correlations(d_fp3_f, d_fp5_f, d_wp5_f)

    print('\n  Bundle direction → Filbin outcome discrimination:', flush=True)
    for name, dpc in [('ΔFilbin-PC3', d_fp3_f), ('ΔFilbin-PC5', d_fp5_f), ('ΔWang-PC5', d_wp5_f)]:
        mwu = mannwhitneyu(dpc[filbin_labels == 1], dpc[filbin_labels == 0])
        imean = dpc[filbin_labels == 1].mean(); wmean = dpc[filbin_labels == 0].mean()
        print(f'    {name:<18}  Imp={imean:+.2f}  Wors={wmean:+.2f}  MWU p={mwu.pvalue:.4f}',
              flush=True)

    # ============================================================
    # (C) GSE94438 progressor vs non-progressor bundle direction
    # ============================================================
    print('\n' + '='*80, flush=True)
    print('(C) GSE94438 multi-timepoint subset: progressor vs non-progressor bundle direction',
          flush=True)
    print('='*80, flush=True)
    # Reload from the saved bundle replication json (has subject-level group info)
    # Or recompute — easier to just rebuild the projections from cached data
    # Actually we don't have GSE94438 α in cache (different cohort), so re-derive from script run
    # For simplicity, re-run the GSE94438 pipeline subset just for projection
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / 'benchmarks'))
    import gzip, pandas as pd
    from collections import defaultdict
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom
    series_path = REPO / 'data/cohorts/GSE94438/GSE94438_series_matrix.txt.gz'
    expr_path = REPO / 'data/cohorts/GSE94438/GSE94438_rawCounts_GeneNames_AllSamples.csv.gz'

    titles = []; chars_rows = []
    with gzip.open(series_path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles = [s.strip('"') for s in line.rstrip().split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'):
                chars_rows.append([s.strip('"') for s in line.rstrip().split('\t')[1:]])
    samples = []
    for i, t in enumerate(titles):
        d = {'title': t}
        for row in chars_rows:
            if i < len(row) and ':' in row[i]:
                k, v = row[i].split(':', 1); d[k.strip()] = v.strip()
        samples.append(d)
    for s in samples:
        s['expr_col'] = f"X{s.get('code', '?')}"
        try: s['t'] = float(s.get('time.from.exposure.months', 'NA'))
        except ValueError: s['t'] = float('nan')

    by_subj = defaultdict(list)
    for s in samples:
        if s.get('subjectid') and not np.isnan(s['t']): by_subj[s['subjectid']].append(s)
    pairs = []
    for sid, ss in by_subj.items():
        by_t = {s['t']: s for s in ss}
        if 0.0 in by_t:
            for t1m in (6.0, 18.0):
                if t1m in by_t:
                    pairs.append((sid, by_t[0.0], by_t[t1m], t1m))

    df = pd.read_csv(expr_path, index_col=0)
    if 'symbol' in df.columns: df = df.drop(columns=['symbol'])
    lib_sizes = df.sum(axis=0); cpm = df.divide(lib_sizes, axis=1) * 1e6
    logcpm = np.log2(cpm + 1).astype(np.float32)
    ensg_ids = list(logcpm.index)
    feat_node_ids = [f'ENSG:{e}' for e in ensg_ids if f'ENSG:{e}' in geom.nid_idx]
    ensg_keep = [n.replace('ENSG:', '') for n in feat_node_ids]
    all_cols = []; pair_indices = []
    for sid, t0_s, t1_s, t1m in pairs:
        if t0_s['expr_col'] in logcpm.columns and t1_s['expr_col'] in logcpm.columns:
            i0 = len(all_cols); all_cols.append(t0_s['expr_col'])
            i1 = len(all_cols); all_cols.append(t1_s['expr_col'])
            pair_indices.append((i0, i1, sid, t0_s.get('group', 'NA'), t1m))
    print(f'  GSE94438 pairs: {len(pair_indices)}', flush=True)

    from gizmo.inference.projection import solve_map, ModalitySetup
    X = logcpm.loc[ensg_keep, all_cols].values.T.astype(float)
    mu = X.mean(axis=0); sd = X.std(axis=0) + 1e-9
    Xz = (X - mu) / sd
    feat_cols = [(f'feat_{k}', geom.nid_idx[feat_node_ids[k]]) for k in range(len(feat_node_ids))]
    pdata = {sid: {f'feat_{k}': float(Xz[i, k]) for k in range(Xz.shape[1])}
             for i, sid in enumerate(all_cols)}
    setup = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=pdata)
    F, _ = solve_map(geom, [setup], all_cols)
    F = F.astype(np.float32)
    beta = F @ lpn
    alpha = (F - np.outer(beta, lpn)).astype(np.float32)
    del F; import gc; gc.collect()

    proj_fp3 = alpha @ filbin_pc3
    proj_fp5 = alpha @ filbin_pc5
    proj_wp5 = alpha @ wang_pc5
    del alpha; gc.collect()

    d_fp3_g = np.array([proj_fp3[i1] - proj_fp3[i0] for (i0, i1, _, _, _) in pair_indices])
    d_fp5_g = np.array([proj_fp5[i1] - proj_fp5[i0] for (i0, i1, _, _, _) in pair_indices])
    d_wp5_g = np.array([proj_wp5[i1] - proj_wp5[i0] for (i0, i1, _, _, _) in pair_indices])
    groups = np.array([g for (_, _, _, g, _) in pair_indices])
    is_progressor = (groups == 'case (TB)').astype(int)
    print(f'  Progressors: {is_progressor.sum()}, Non-progressors: {(1-is_progressor).sum()}',
          flush=True)

    print('\n  Bundle direction → GSE94438 progressor discrimination:', flush=True)
    for name, dpc in [('ΔFilbin-PC3', d_fp3_g), ('ΔFilbin-PC5', d_fp5_g), ('ΔWang-PC5', d_wp5_g)]:
        mwu = mannwhitneyu(dpc[is_progressor == 1], dpc[is_progressor == 0])
        pmean = dpc[is_progressor == 1].mean(); nmean = dpc[is_progressor == 0].mean()
        print(f'    {name:<18}  Progressor={pmean:+.2f}  Non-Prog={nmean:+.2f}  '
              f'MWU p={mwu.pvalue:.4f}', flush=True)

    out['GSE94438_progressor_test'] = {
        'n_progressor': int(is_progressor.sum()),
        'n_non_progressor': int((1 - is_progressor).sum()),
        'per_axis': {
            'Filbin-PC3': {'progressor_mean': float(d_fp3_g[is_progressor==1].mean()),
                            'non_progressor_mean': float(d_fp3_g[is_progressor==0].mean()),
                            'p_grp': float(mannwhitneyu(d_fp3_g[is_progressor==1], d_fp3_g[is_progressor==0]).pvalue)},
            'Filbin-PC5': {'progressor_mean': float(d_fp5_g[is_progressor==1].mean()),
                            'non_progressor_mean': float(d_fp5_g[is_progressor==0].mean()),
                            'p_grp': float(mannwhitneyu(d_fp5_g[is_progressor==1], d_fp5_g[is_progressor==0]).pvalue)},
            'Wang-PC5': {'progressor_mean': float(d_wp5_g[is_progressor==1].mean()),
                          'non_progressor_mean': float(d_wp5_g[is_progressor==0].mean()),
                          'p_grp': float(mannwhitneyu(d_wp5_g[is_progressor==1], d_wp5_g[is_progressor==0]).pvalue)},
        },
    }

    out_path = RESULTS / 'cure_bundle_cross_disease.json'
    with open(out_path, 'w') as fh: json.dump(out, fh, indent=2)
    print(f'\nSaved: {out_path}', flush=True)


if __name__ == '__main__':
    main()
