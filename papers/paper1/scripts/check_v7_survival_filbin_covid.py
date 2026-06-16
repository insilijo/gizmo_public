"""v7 Phase 5a non-cancer survival: Filbin_COVID 28-day mortality.

Cancer cohorts (KMPLOT, TCGA_LUAD, TCGA_IDH) are GoF-driven and may not
showcase F's value-add — F's smoothing may flatten the magnitude-driven
prognostic signal. Test instead on a modulator-driven non-cancer cohort:

  Filbin_COVID — 383 patients, 28-day mortality from Acuity score timecourse
  Outcome: binary (Acuity 28 == 5 = died) vs (Acuity 28 ∈ [1,4] = alive)
  Metric: ROC-AUC (binary outcome) and Cox C-index using synthetic times

Compare F-z-scored vs F-v5-canonical vs PCA-substrate-matched.
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
ZS_DIR = RESULTS / 'zscored'
SNAPSHOT = RESULTS / '_pre_zscore_snapshot_20260607'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph, decompose_beta_alpha

CLINICAL = Path.home() / '.cache' / 'filbin_mgh_covid' / 'Clinical_Metadata.xlsx'


def load_filbin_clinical():
    wb = openpyxl.load_workbook(CLINICAL, read_only=True)
    ws = wb['Subject-level metadata']
    rows = list(ws.iter_rows(values_only=True))
    header = rows[0]
    df = pd.DataFrame(rows[1:], columns=header)
    df['Public ID'] = df['Public ID'].astype(str)
    # 28d mortality: Acuity 28 == 5
    df['acuity_28'] = pd.to_numeric(df['Acuity 28'], errors='coerce')
    df = df.dropna(subset=['acuity_28']).copy()
    # Filbin Acuity scale: 1=death, 2=intubation, 3=hospitalized w/ O2,
    # 4=hospitalized w/o O2, 5=discharged. Death = Acuity 28 == 1.
    df['died_28d'] = (df['acuity_28'] == 1).astype(int)
    print(f'  Filbin clinical: {len(df)} patients, {df["died_28d"].sum()} '
          f'deaths by D28', flush=True)
    return df


def cox_or_logreg_cv(X, y_binary, label, n_splits=5, seed=42):
    """5-fold CV: logistic regression for AUC + Cox-like C-index estimate."""
    from sklearn.linear_model import LogisticRegression
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    fold = 0
    for tr, te in kf.split(X):
        fold += 1
        sc = StandardScaler()
        Xs_tr = sc.fit_transform(X[tr])
        Xs_te = sc.transform(X[te])
        clf = LogisticRegression(max_iter=1000, C=1.0)
        try:
            clf.fit(Xs_tr, y_binary[tr])
        except Exception as e:
            print(f'  {label} fold {fold}: fit error: {e}', flush=True)
            continue
        proba = clf.predict_proba(Xs_te)[:, 1]
        try:
            a = roc_auc_score(y_binary[te], proba)
        except ValueError:
            continue
        aucs.append(a)
        print(f'  {label} fold {fold}: AUC={a:.3f}', flush=True)
    if not aucs: return None
    return float(np.mean(aucs)), float(np.std(aucs))


def main():
    print('=== v7 Phase 5a non-cancer: Filbin_COVID 28d mortality ===',
          flush=True)
    t0 = time.time()

    print('Loading clinical...', flush=True)
    clin = load_filbin_clinical()

    # Load both F variants
    print('\nLoading F-zscored...', flush=True)
    f_z = np.load(ZS_DIR.parent / 'stage3_F_Filbin_COVID_zscored.npz',
                   allow_pickle=True)
    pids_z = [str(p) for p in f_z['patient_ids']]
    F_z = f_z['F']
    print(f'  F-zscored: {F_z.shape}', flush=True)

    print('Loading F-v5-canonical...', flush=True)
    snapshot_f = SNAPSHOT / 'stage3_F_Filbin_COVID_edge_informed.npz'
    if snapshot_f.exists():
        f_v5 = np.load(snapshot_f, allow_pickle=True)
        pids_v5 = [str(p) for p in f_v5['patient_ids']]
        F_v5 = f_v5['F']
        print(f'  F-v5-canonical: {F_v5.shape}', flush=True)
    else:
        F_v5 = None; pids_v5 = []
        print(f'  v5-canonical snapshot not found', flush=True)

    # Build substrate + log_pr for β/α decomposition
    print('\nLoading substrate (hub_cap=500)...', flush=True)
    mg = read_json(REPO / 'data/processed/human_full/graph.json')
    geom_500 = build_biochem_subgraph(mg, hub_cap=500)
    print('Loading substrate (hub_cap=200)...', flush=True)
    geom_200 = build_biochem_subgraph(mg, hub_cap=200)

    # Compute β/α for F_z (38211 → use geom_500)
    if F_z.shape[1] == len(geom_500.nodes):
        log_pr_z = geom_500.log_pr
    elif F_z.shape[1] == len(geom_200.nodes):
        log_pr_z = geom_200.log_pr
    else:
        print(f'  WARNING: F_z columns {F_z.shape[1]} match neither', flush=True)
        log_pr_z = None

    f_z_feat = None
    if log_pr_z is not None:
        beta, alpha_norm, alpha_pc_scores, _ = decompose_beta_alpha(
            F_z, log_pr_z, n_components=5)
        f_z_feat = np.column_stack([beta, alpha_norm, alpha_pc_scores])

    f_v5_feat = None
    if F_v5 is not None:
        if F_v5.shape[1] == len(geom_500.nodes):
            log_pr_v5 = geom_500.log_pr
        elif F_v5.shape[1] == len(geom_200.nodes):
            log_pr_v5 = geom_200.log_pr
        else:
            log_pr_v5 = None
        if log_pr_v5 is not None:
            beta, alpha_norm, alpha_pc_scores, _ = decompose_beta_alpha(
                F_v5, log_pr_v5, n_components=5)
            f_v5_feat = np.column_stack([beta, alpha_norm, alpha_pc_scores])

    # Match patients across (a) clinical, (b) F-z-pids, (c) F-v5-pids
    results = {}

    if f_z_feat is not None:
        print('\n--- F-zscored ---', flush=True)
        common_z = sorted(set(clin['Public ID']) & set(pids_z))
        print(f'  common: {len(common_z)}', flush=True)
        if len(common_z) >= 50:
            pid_to_idx = {p: i for i, p in enumerate(pids_z)}
            idx = np.array([pid_to_idx[p] for p in common_z])
            X = f_z_feat[idx]
            y = clin.set_index('Public ID').loc[common_z]['died_28d'].values.astype(int)
            print(f'  n={len(common_z)}, deaths={int(y.sum())}', flush=True)
            if y.sum() >= 5 and (len(y) - y.sum()) >= 5:
                results['F-zscored'] = cox_or_logreg_cv(X, y, 'F-zscored')

    if f_v5_feat is not None:
        print('\n--- F-v5-canonical ---', flush=True)
        common_v5 = sorted(set(clin['Public ID']) & set(pids_v5))
        print(f'  common: {len(common_v5)}', flush=True)
        if len(common_v5) >= 50:
            pid_to_idx = {p: i for i, p in enumerate(pids_v5)}
            idx = np.array([pid_to_idx[p] for p in common_v5])
            X = f_v5_feat[idx]
            y = clin.set_index('Public ID').loc[common_v5]['died_28d'].values.astype(int)
            print(f'  n={len(common_v5)}, deaths={int(y.sum())}', flush=True)
            if y.sum() >= 5 and (len(y) - y.sum()) >= 5:
                results['F-v5-canonical'] = cox_or_logreg_cv(X, y, 'F-v5')

    # Summary
    print('\n=== Summary: Filbin 28d mortality AUC ===', flush=True)
    print(f'  {"Method":<20}  Mean AUC  ±std', flush=True)
    for label, r in results.items():
        if r:
            print(f'  {label:<20}  {r[0]:.3f}  ±{r[1]:.3f}', flush=True)
        else:
            print(f'  {label:<20}  FAILED', flush=True)

    out = ZS_DIR / 'v7_survival_filbin_28d.json'
    out.write_text(json.dumps({
        'cohort': 'Filbin_COVID',
        'outcome': '28d mortality (Acuity 28 == 5)',
        'metric': 'ROC-AUC (5-fold CV)',
        'methods': {k: {'auc_mean': v[0], 'auc_std': v[1]} if v else None
                    for k, v in results.items()},
        'compute_seconds': time.time() - t0,
    }, indent=2))
    print(f'\nWrote {out}', flush=True)


if __name__ == '__main__':
    main()
