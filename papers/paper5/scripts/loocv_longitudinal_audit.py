"""LOOCV-corrected trajectory α-PC AUCs across paired longitudinal cohorts.

For each paired cohort, solve MAP jointly on (T0 + T_later) samples (2N rows),
compute β/α decomposition. Trajectory deltas per patient:
  Δβ_i = β_i(T_later) − β_i(T0)
  Δα-PC_k,i = α-PC_k,i(T_later) − α-PC_k,i(T0)
AUC tested against binary trajectory label (Responder/Non-Responder,
Improved/Worsened, Cure/Not-Cure).

LOOCV: hold out patient i's BOTH timepoints, refit PCA on (2N−2) rows,
sign-align components to the naive full-data fit, project i's two timepoints,
compute held-out Δα-PC. Aggregate held-out Δα-PCs across N folds, compute AUC.

β is graph-fixed (log_PR projection) so no CV is needed for Δβ AUC.

Saves _scores so DeLong CI + label-shuffle permutation null can be run
post-hoc by loocv_perm_and_ci_longitudinal.py.
"""
from __future__ import annotations
import sys, json, gzip, re
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import (
    build_biochem_subgraph, solve_map, ModalitySetup,
)

print('Loading substrate...')
mg = read_json(str(REPO / 'data/processed/human_full_rhea_full/graph.json'))
geom = build_biochem_subgraph(mg, hub_cap=500)
print(f'  {len(geom.nodes)} nodes')


def safe_auc(y, s):
    a = roc_auc_score(y, s); return max(a, 1 - a)


# ---------- cohort loaders: return (F, trajectory_label, n_pairs) -----------
# F is (2N, n_nodes) stacked: rows [0:N] = T0 in patient order, [N:2N] = T_later
# trajectory_label is (N,) binary

def cohort_wang_ra():
    """Wang 2025 RA MTX paired baseline (A) → followup (B), binary Resp vs No-Resp."""
    print('\n[Wang_RA_trajectory] loading...')
    WANG_DIR = REPO / 'data/cohorts/Wang_RA_MTX/hesy1191569605-rheumatoid-arthritis-0c94b6d'
    df = pd.read_csv(WANG_DIR / 'Figure6/csv/RA_DATAKNN1.csv')
    PROT = list(df.columns[15:])
    # Group1 ∈ {Health, RA_A, RA_B, RA_C}; pid via Sample column
    df['pid'] = df['Sample'].astype(str).str[1:].str.lstrip('0')
    df['tp'] = df['Sample'].astype(str).str[0]  # A / B / C
    by_pid = df.groupby('pid').agg({'tp': list, 'Drugs Response': 'first'})
    # Require A and B paired with non-null Drugs Response
    paired_pids = [p for p, r in by_pid.iterrows()
                   if 'A' in r['tp'] and 'B' in r['tp']
                   and r['Drugs Response'] in ('Response', 'No Response')]
    print(f'  Paired A+B patients: {len(paired_pids)}')
    print(f'  Response: {Counter(by_pid.loc[p, "Drugs Response"] for p in paired_pids)}')

    feat = sorted(set(PROT))
    feat_node_ids = [f'symbol:{g}' for g in feat if f'symbol:{g}' in geom.nid_idx]
    print(f'  Mapped features: {len(feat_node_ids)}/{len(feat)}')
    feat = [n.replace('symbol:', '') for n in feat_node_ids]

    X_t0 = np.zeros((len(paired_pids), len(feat)))
    X_t1 = np.zeros((len(paired_pids), len(feat)))
    for i, pid in enumerate(paired_pids):
        sub = df[df.pid == pid]
        rowA = sub[sub.tp == 'A'].iloc[0]
        rowB = sub[sub.tp == 'B'].iloc[0]
        for j, g in enumerate(feat):
            X_t0[i, j] = float(rowA[g]) if pd.notna(rowA[g]) else 0.0
            X_t1[i, j] = float(rowB[g]) if pd.notna(rowB[g]) else 0.0
    y = np.array([1 if by_pid.loc[p, 'Drugs Response'] == 'Response' else 0
                  for p in paired_pids])
    return _solve_paired(X_t0, X_t1, feat_node_ids, paired_pids, y)


def cohort_filbin_d0d3():
    """Filbin COVID D0→D3 paired, binary Improved vs Worsened."""
    print('\n[Filbin_D0D3_trajectory] loading...')
    CACHE = Path.home() / '.cache' / 'filbin_mgh_covid'
    t2a = pd.read_excel(CACHE / 'Suppl_T2_Olink_Assays_NPX.xlsx', sheet_name='2A-Olink-Assay', header=1)
    oid_to_sym = dict(zip(t2a['OlinkID'].astype(str), t2a['Assay'].astype(str)))
    clin = pd.read_excel(CACHE / 'Clinical_Metadata.xlsx', sheet_name='Subject-level metadata')
    clin['core_pid'] = clin['Public ID'].astype(str)
    pid_covid = dict(zip(clin['core_pid'], clin['COVID'].astype(int)))
    pid_acuity = {str(r['core_pid']): {0: r.get('Acuity 0'), 3: r.get('Acuity 3')}
                  for _, r in clin.iterrows()}

    ol = pd.read_excel(CACHE / 'Olink_Proteomics.xlsx')
    ol['core_pid'] = ol['Public ID'].astype(str).str.replace(r'_D\d+$', '', regex=True).str.replace(r'_E$', '', regex=True)
    by_day = {0: {}, 3: {}}
    for _, row in ol.iterrows():
        try: d = int(row['Day'])
        except Exception: continue
        if d not in (0, 3): continue
        pid = str(row['core_pid'])
        if pid_covid.get(pid) != 1: continue
        feat_d = {}
        for col, v in row.items():
            if not isinstance(col, str) or not col.startswith('OID'): continue
            sym = oid_to_sym.get(col)
            if not sym or sym == 'nan': continue
            try: vf = float(v)
            except (TypeError, ValueError): continue
            if pd.isna(vf): continue
            feat_d[sym] = max(feat_d.get(sym, -1e9), vf)
        if feat_d:
            by_day[d][pid] = feat_d

    paired = sorted(set(by_day[0].keys()) & set(by_day[3].keys()))
    # Improved (Acuity↑) vs Worsened (Acuity↓); Stable dropped
    use = []
    y = []
    for pid in paired:
        a0 = pid_acuity[pid][0]; a3 = pid_acuity[pid][3]
        if pd.isna(a0) or pd.isna(a3): continue
        d = int(a3) - int(a0)
        if d > 0: use.append(pid); y.append(1)  # Improved
        elif d < 0: use.append(pid); y.append(0)  # Worsened
    print(f'  Improved vs Worsened paired: {len(use)} ({Counter(y)})')

    all_genes = sorted({k for pid in use for d in (0, 3) for k in by_day[d][pid].keys()})
    feat_node_ids = [f'symbol:{g}' for g in all_genes if f'symbol:{g}' in geom.nid_idx]
    feat = [n.replace('symbol:', '') for n in feat_node_ids]
    print(f'  Mapped features: {len(feat_node_ids)}/{len(all_genes)}')

    X_t0 = np.zeros((len(use), len(feat)))
    X_t1 = np.zeros((len(use), len(feat)))
    for i, pid in enumerate(use):
        for j, g in enumerate(feat):
            X_t0[i, j] = by_day[0][pid].get(g, 0.0)
            X_t1[i, j] = by_day[3][pid].get(g, 0.0)
    return _solve_paired(X_t0, X_t1, feat_node_ids, use, np.array(y))


def cohort_tb_dx_wk24():
    """GSE89403 TB DX→Wk24 paired, binary Definite Cure vs (Probable|Possible|Not Cured)."""
    print('\n[TB_DX_Wk24_trajectory] loading...')
    COHORT = REPO / 'data/cohorts/GSE89403_TB_treat'
    annots = {}
    titles = []; chars = []
    with gzip.open(COHORT / 'GSE89403_series_matrix.txt.gz', 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles = [s.strip('"').split('-')[0] for s in line.rstrip().split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'):
                chars.append([s.strip('"') for s in line.rstrip().split('\t')[1:]])
    for i in range(len(titles)):
        sc = chars[1][i].replace('sample_code: ', '')
        subj = chars[2][i].replace('subject: ', '')
        if sc in annots: continue
        annots[sc] = {
            'subject': subj,
            'disease': chars[3][i].replace('disease state: ', ''),
            'tr': chars[4][i].replace('treatmentresult: ', ''),
            'time': chars[5][i].replace('time: ', ''),
        }
    df = pd.read_csv(COHORT / 'GSE89403_log2ExpGeneNames_AllSamples.csv.gz',
                     compression='gzip', index_col=0)
    if 'symbol' in df.columns: df = df.drop(columns=['symbol'])
    # Find subjects with TB + DX + week_24 paired AND treatment_result populated
    by_subj = {}
    for sid, a in annots.items():
        if a['disease'] != 'TB Subjects': continue
        if sid not in df.columns: continue
        by_subj.setdefault(a['subject'], {})[a['time']] = (sid, a['tr'])
    paired = [(s, tps['DX'][0], tps['week_24'][0], tps['DX'][1])
              for s, tps in by_subj.items()
              if 'DX' in tps and 'week_24' in tps and tps['DX'][1] in
              ('Definite Cure', 'Probable Cure', 'Possible Cure', 'Not Cured')]
    paired = [(s, dx, wk, tr) for s, dx, wk, tr in paired
              if tr in ('Definite Cure', 'Not Cured')]  # binary contrast
    print(f'  TB DX→Wk24 paired (Cure vs Not Cured): {len(paired)}')
    print(f'  Outcome: {Counter(tr for _, _, _, tr in paired)}')

    feat_node_ids = [f'ENSG:{e}' for e in df.index if f'ENSG:{e}' in geom.nid_idx]
    ensg_keep = [n.replace('ENSG:', '') for n in feat_node_ids]
    print(f'  Mapped features: {len(feat_node_ids)}')
    X_t0 = df.loc[ensg_keep, [dx for _, dx, _, _ in paired]].values.T.astype(float)
    X_t1 = df.loc[ensg_keep, [wk for _, _, wk, _ in paired]].values.T.astype(float)
    y = np.array([1 if tr == 'Definite Cure' else 0 for _, _, _, tr in paired])
    subj_ids = [s for s, _, _, _ in paired]
    return _solve_paired(X_t0, X_t1, feat_node_ids, subj_ids, y)


# ---------------------- shared MAP solver helpers --------------------------

def _solve_paired(X_t0, X_t1, feat_node_ids, patient_ids, labels):
    """Stack T0+T1, pooled z-score, solve MAP on 2N samples.

    Returns dict with F_t0, F_t1, naive Δβ, naive Δα-PC scores,
    LOOCV Δα-PC scores, labels.
    """
    N = X_t0.shape[0]
    X_pool = np.vstack([X_t0, X_t1])
    mu = X_pool.mean(axis=0); sd = X_pool.std(axis=0) + 1e-9
    X_t0z = (X_t0 - mu) / sd; X_t1z = (X_t1 - mu) / sd
    feat_cols = [(f'feat_{k}', geom.nid_idx[feat_node_ids[k]])
                 for k in range(len(feat_node_ids))]
    t0_sids = [f'{p}_T0' for p in patient_ids]
    t1_sids = [f'{p}_T1' for p in patient_ids]
    all_sids = t0_sids + t1_sids
    pdata = {}
    for i, sid in enumerate(t0_sids):
        pdata[sid] = {f'feat_{k}': float(X_t0z[i, k]) for k in range(X_t0z.shape[1])}
    for i, sid in enumerate(t1_sids):
        pdata[sid] = {f'feat_{k}': float(X_t1z[i, k]) for k in range(X_t1z.shape[1])}
    setup = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=pdata)
    print(f'  Solving MAP on 2N={len(all_sids)} samples × {len(feat_node_ids)} features...')
    F, _ = solve_map(geom, [setup], all_sids)
    F_t0 = F[:N]; F_t1 = F[N:]

    log_pr = geom.log_pr
    log_pr_norm = log_pr / (np.linalg.norm(log_pr) + 1e-9)
    # β and α at each timepoint
    beta_t0 = F_t0 @ log_pr_norm; beta_t1 = F_t1 @ log_pr_norm
    alpha_t0 = F_t0 - np.outer(beta_t0, log_pr_norm)
    alpha_t1 = F_t1 - np.outer(beta_t1, log_pr_norm)

    # Naive: PCA on the full 2N-row α matrix
    alpha_all = np.vstack([alpha_t0, alpha_t1])
    pca = PCA(n_components=5, svd_solver='randomized', random_state=0).fit(alpha_all)
    naive_components = pca.components_.copy()
    naive_alpha_t0 = pca.transform(alpha_t0)
    naive_alpha_t1 = pca.transform(alpha_t1)
    naive_dalpha = naive_alpha_t1 - naive_alpha_t0
    naive_aucs = [safe_auc(labels, naive_dalpha[:, k]) for k in range(5)]
    naive_ev = pca.explained_variance_ratio_.tolist()

    # LOOCV — hold out patient i (both rows), refit PCA on (2N-2) rows
    loocv_dalpha = np.zeros((N, 5))
    for i in range(N):
        mask = np.ones(2 * N, dtype=bool)
        mask[i] = False; mask[N + i] = False  # drop both timepoints
        pca_i = PCA(n_components=5, svd_solver='randomized', random_state=0).fit(alpha_all[mask])
        for k in range(5):
            if np.dot(pca_i.components_[k], naive_components[k]) < 0:
                pca_i.components_[k] *= -1
        s_t0 = pca_i.transform(alpha_t0[i:i+1])[0]
        s_t1 = pca_i.transform(alpha_t1[i:i+1])[0]
        loocv_dalpha[i] = s_t1 - s_t0
    loocv_aucs = [safe_auc(labels, loocv_dalpha[:, k]) for k in range(5)]

    # Δβ (graph-fixed)
    dbeta = beta_t1 - beta_t0
    dbeta_auc = safe_auc(labels, dbeta)

    return {
        'n_pairs': int(N),
        'n_features': int(len(feat_node_ids)),
        'labels_pos': int(labels.sum()),
        'labels_neg': int((1 - labels).sum()),
        'dbeta_AUC': float(dbeta_auc),
        'naive_dalpha_PC_AUCs': [float(a) for a in naive_aucs],
        'loocv_dalpha_PC_AUCs': [float(a) for a in loocv_aucs],
        'bias_per_PC': [float(naive_aucs[k] - loocv_aucs[k]) for k in range(5)],
        'explained_variance': naive_ev,
        '_scores': {
            'dbeta': dbeta.tolist(),
            'loocv_dalpha_pc': loocv_dalpha.tolist(),
            'labels': labels.tolist(),
            'patient_ids': list(patient_ids),
        },
    }


def main():
    cohorts = {
        'Wang_RA_traj': cohort_wang_ra,
        'Filbin_D0D3_traj': cohort_filbin_d0d3,
        'TB_DX_Wk24_traj': cohort_tb_dx_wk24,
    }
    out_json = RESULTS / 'loocv_longitudinal_audit.json'
    if out_json.exists():
        with open(out_json) as fh: results = json.load(fh)
        print(f'Resuming, already have: {list(results.keys())}')
    else:
        results = {}
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for name, fn in cohorts.items():
        if only and name != only: continue
        if name in results:
            print(f'\n[{name}] skipped — already done.')
            continue
        try:
            out = fn()
            results[name] = out
            print(f'  Δβ AUC:              {out["dbeta_AUC"]:.3f}')
            print(f'  naive Δα-PC AUCs:    {[f"{a:.3f}" for a in out["naive_dalpha_PC_AUCs"]]}')
            print(f'  LOOCV Δα-PC AUCs:    {[f"{a:.3f}" for a in out["loocv_dalpha_PC_AUCs"]]}')
            print(f'  bias (naive-LOOCV):  {[f"{b:+.3f}" for b in out["bias_per_PC"]]}')
            with open(out_json, 'w') as fh: json.dump(results, fh, indent=2)
            print(f'  saved partial to {out_json}')
        except Exception as e:
            print(f'  FAILED: {e}')
            import traceback; traceback.print_exc()
    print(f'\nFinal save: {out_json}')


if __name__ == '__main__':
    main()
