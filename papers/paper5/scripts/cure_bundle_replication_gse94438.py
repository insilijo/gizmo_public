"""Replicate cure-axis bundle correlation structure in GSE94438 multi-timepoint subset.

Original finding (GSE89403 DX→Wk24, n=76):
  ΔFilbin-PC3 ↔ ΔFilbin-PC5  ρ=+0.64 [+0.48, +0.76]
  ΔFilbin-PC3 ↔ ΔWang-PC5    ρ=+0.75 [+0.61, +0.84]
  ΔFilbin-PC5 ↔ ΔWang-PC5    ρ=+0.54 [+0.35, +0.70]

GSE94438 has paired samples at months 0/6/18 from exposure for some subjects.
We compute ΔPC_T = projection(later) - projection(month_0) per subject per axis,
then test pairwise Spearman correlation of the bundle axes' trajectory signals.

Caveat: GSE94438 is progression monitoring, not cure trajectory. The replication
question is whether the three orthogonal axes coordinate in *any* TB-related
biology trajectory, or only in cure trajectories specifically.
"""
from __future__ import annotations
import sys, gzip, json, gc
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


def parse_series_matrix(path):
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
            if i < len(row) and ':' in row[i]:
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


def main():
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom
    lpn = (geom.log_pr / (np.linalg.norm(geom.log_pr) + 1e-9)).astype(np.float32)
    filbin_pc3 = np.load(RESULTS / 'filbin_pc3_eigvec.npy').astype(np.float32)
    filbin_pc5 = np.load(RESULTS / 'filbin_pc5_eigvec.npy').astype(np.float32)
    wang_pc5 = np.load(RESULTS / 'wang_pc5_eigvec.npy').astype(np.float32)
    print(f'Loaded 3 bundle eigvecs ({filbin_pc3.shape})', flush=True)

    # ---- Identify multi-timepoint subjects in GSE94438 ----
    print('\nParsing GSE94438 metadata...', flush=True)
    series_path = REPO / 'data/cohorts/GSE94438/GSE94438_series_matrix.txt.gz'
    samples = parse_series_matrix(series_path)
    for s in samples:
        s['expr_col'] = f"X{s.get('code', '?')}"
        try: s['t'] = float(s.get('time.from.exposure.months', 'NA'))
        except ValueError: s['t'] = float('nan')

    # Group by subject
    by_subject = defaultdict(list)
    for s in samples:
        if s.get('subjectid') and not np.isnan(s['t']):
            by_subject[s['subjectid']].append(s)

    # Find subjects with month_0 AND a later timepoint
    pairs = []  # list of (subject, t0_sample, t1_sample, t1_months)
    for sid, subj_samples in by_subject.items():
        by_t = {s['t']: s for s in subj_samples}
        if 0.0 in by_t:
            t0 = by_t[0.0]
            for t1m in (6.0, 18.0):
                if t1m in by_t:
                    pairs.append((sid, t0, by_t[t1m], t1m))
    print(f'  Subjects with paired month_0 + later timepoint: {len(pairs)}', flush=True)
    print(f'  Pair breakdown: {Counter(p[3] for p in pairs)}', flush=True)
    print(f'  Outcome breakdown (case vs control): '
          f'{Counter(p[1].get("group") for p in pairs)}', flush=True)

    if len(pairs) < 20:
        print(f'  ⚠ Only {len(pairs)} paired subjects — limited power', flush=True)

    # ---- Load expression and project ----
    print('\nLoading GSE94438 expression...', flush=True)
    expr_path = REPO / 'data/cohorts/GSE94438/GSE94438_rawCounts_GeneNames_AllSamples.csv.gz'
    df = pd.read_csv(expr_path, index_col=0)
    if 'symbol' in df.columns: df = df.drop(columns=['symbol'])
    lib_sizes = df.sum(axis=0)
    cpm = df.divide(lib_sizes, axis=1) * 1e6
    logcpm = np.log2(cpm + 1).astype(np.float32)
    del df, cpm; gc.collect()
    ensg_ids = list(logcpm.index)
    feat_node_ids = [f'ENSG:{e}' for e in ensg_ids if f'ENSG:{e}' in geom.nid_idx]
    ensg_keep = [n.replace('ENSG:', '') for n in feat_node_ids]
    print(f'  Mapped features: {len(feat_node_ids)}', flush=True)

    # Collect all sample columns we need
    all_cols = []
    pair_indices = []  # (t0_idx, t1_idx, subject, group, t1m)
    for sid, t0_s, t1_s, t1m in pairs:
        if t0_s['expr_col'] in logcpm.columns and t1_s['expr_col'] in logcpm.columns:
            idx0 = len(all_cols)
            all_cols.append(t0_s['expr_col'])
            idx1 = len(all_cols)
            all_cols.append(t1_s['expr_col'])
            pair_indices.append((idx0, idx1, sid, t0_s.get('group', 'NA'), t1m))
    print(f'  Pairs with expression available: {len(pair_indices)}', flush=True)

    X = logcpm.loc[ensg_keep, all_cols].values.T.astype(float)
    del logcpm; gc.collect()

    # Pooled z-score across all samples in this analysis (consistent with other cohorts)
    print(f'\nSolving MAP on {X.shape[0]} samples × {X.shape[1]} features...', flush=True)
    F = solve_F(geom, X, feat_node_ids, all_cols)
    beta = F @ lpn
    alpha = (F - np.outer(beta, lpn)).astype(np.float32)
    del F; gc.collect()
    print(f'  α matrix: {alpha.shape}', flush=True)

    # Project onto each bundle axis
    proj_fp3 = alpha @ filbin_pc3
    proj_fp5 = alpha @ filbin_pc5
    proj_wp5 = alpha @ wang_pc5
    del alpha; gc.collect()

    # Compute ΔPC per pair
    dpc_fp3 = np.array([proj_fp3[i1] - proj_fp3[i0] for (i0, i1, _, _, _) in pair_indices])
    dpc_fp5 = np.array([proj_fp5[i1] - proj_fp5[i0] for (i0, i1, _, _, _) in pair_indices])
    dpc_wp5 = np.array([proj_wp5[i1] - proj_wp5[i0] for (i0, i1, _, _, _) in pair_indices])
    group_label = np.array([g for (_, _, _, g, _) in pair_indices])
    t1m_label = np.array([t for (_, _, _, _, t) in pair_indices])

    n = len(dpc_fp3)
    print(f'\nBundle replication test (n={n}):', flush=True)
    print(f'  Group breakdown: {Counter(group_label)}', flush=True)
    print(f'  Timepoint breakdown: {Counter(t1m_label.astype(int))}', flush=True)

    # ---- Pairwise observed Spearman + bootstrap CIs ----
    pair_specs = [
        ('ΔFilbin-PC3', 'ΔFilbin-PC5', dpc_fp3, dpc_fp5, +0.639),
        ('ΔFilbin-PC3', 'ΔWang-PC5',   dpc_fp3, dpc_wp5, +0.749),
        ('ΔFilbin-PC5', 'ΔWang-PC5',   dpc_fp5, dpc_wp5, +0.537),
    ]
    print(f'\n  Pairwise Spearman ρ (observed) + 95% bootstrap CIs (1000 resamples):',
          flush=True)
    print(f'  {"pair":<30}{"obs ρ":>9}{"95% CI":>22}{"GSE89403 ρ":>14}{"replicates?":>14}',
          flush=True)
    rng = np.random.default_rng(0)
    n_boot = 1000
    bootstrap_results = {}
    for name_a, name_b, x, y, orig_rho in pair_specs:
        rho, _ = spearmanr(x, y)
        boot_rhos = []
        for b in range(n_boot):
            idx = rng.integers(0, n, n)
            r, _ = spearmanr(x[idx], y[idx])
            if not np.isnan(r): boot_rhos.append(r)
        if len(boot_rhos) < 100:
            ci_lo, ci_hi = float('nan'), float('nan')
            verdict = '?'
        else:
            ci_lo = float(np.percentile(boot_rhos, 2.5))
            ci_hi = float(np.percentile(boot_rhos, 97.5))
            # Replication: CI overlaps original ρ AND CI excludes 0 (same direction)
            ci_overlap = (orig_rho >= ci_lo) and (orig_rho <= ci_hi)
            ci_nonzero = (ci_lo > 0) or (ci_hi < 0)
            same_sign = np.sign(rho) == np.sign(orig_rho)
            if ci_overlap and ci_nonzero and same_sign: verdict = '✓ strong'
            elif same_sign and ci_nonzero:              verdict = '~ direction'
            elif same_sign:                              verdict = '? underpowered'
            else:                                         verdict = '✗ different'
        ci_str = f'[{ci_lo:+.2f}, {ci_hi:+.2f}]'
        print(f'  {name_a} ↔ {name_b:<14}{rho:>+9.3f}{ci_str:>22}{orig_rho:>+14.3f}{verdict:>14}',
              flush=True)
        bootstrap_results[f'{name_a}__{name_b}'] = {
            'obs_rho': float(rho), 'CI95': [ci_lo, ci_hi],
            'gse89403_rho': float(orig_rho), 'verdict': verdict,
        }

    out = {
        'cohort': 'GSE94438',
        'n_pairs': int(n),
        'timepoint_distribution': {str(k): int(v) for k, v in Counter(t1m_label.astype(int)).items()},
        'group_distribution': {str(k): int(v) for k, v in Counter(group_label).items()},
        'pairwise': bootstrap_results,
        'gse89403_reference': {
            'Filbin-PC3__Filbin-PC5': 0.639,
            'Filbin-PC3__Wang-PC5':   0.749,
            'Filbin-PC5__Wang-PC5':   0.537,
        },
    }
    out_path = RESULTS / 'cure_bundle_replication_gse94438.json'
    with open(out_path, 'w') as fh: json.dump(out, fh, indent=2)
    print(f'\nSaved: {out_path}', flush=True)


if __name__ == '__main__':
    main()
