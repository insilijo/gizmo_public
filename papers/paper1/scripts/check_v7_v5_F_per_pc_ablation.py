"""v7: per-α-PC ablation on v5-canonical F across all 4 cohorts.

Companion to check_v7_pca_per_pc_ablation.py — same survival/AUC test
on each individual α-PC, but using the v5-canonical F (from the
_pre_zscore_snapshot_20260607/) instead of the mixed-preprocessing F.
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl
from sklearn.linear_model import LogisticRegression
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


def cv_metric(X, y_or_time_event, label, kind='auc', n_splits=5, seed=42):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    scores = []
    for tr, te in kf.split(X):
        if X.ndim == 1:
            Xt = X.reshape(-1, 1)
        else:
            Xt = X
        sc = StandardScaler()
        Xs_tr = sc.fit_transform(Xt[tr]); Xs_te = sc.transform(Xt[te])
        if kind == 'auc':
            y = y_or_time_event
            clf = LogisticRegression(max_iter=1000, C=1.0)
            try:
                clf.fit(Xs_tr, y[tr])
                proba = clf.predict_proba(Xs_te)[:, 1]
                scores.append(roc_auc_score(y[te], proba))
            except (ValueError, Exception):
                continue
        elif kind == 'cox':
            from lifelines import CoxPHFitter
            from lifelines.utils import concordance_index
            time_arr, event_arr = y_or_time_event
            df = pd.DataFrame(Xs_tr,
                                columns=[f'f{i}' for i in range(Xt.shape[1])])
            df['time'] = time_arr[tr]; df['event'] = event_arr[tr]
            cph = CoxPHFitter(penalizer=0.01)
            try:
                cph.fit(df, duration_col='time', event_col='event',
                        show_progress=False)
            except Exception:
                continue
            df_te = pd.DataFrame(Xs_te,
                                   columns=[f'f{i}' for i in range(Xt.shape[1])])
            rs = cph.predict_partial_hazard(df_te).values
            scores.append(concordance_index(time_arr[te], -rs, event_arr[te]))
    if not scores: return None
    return float(np.mean(scores)), float(np.std(scores))


def load_v5_alpha_pcs(cohort, hub_cap=200):
    """Load v5-canonical F + decompose β/α/α-PCs."""
    f_path = SNAPSHOT / f'stage3_F_{cohort}_edge_informed.npz'
    if not f_path.exists():
        f_path = SNAPSHOT / f'stage3_F_{cohort}.npz'
    if not f_path.exists():
        return None, None, None
    npz = np.load(f_path, allow_pickle=True)
    F = npz['F']
    pids = [str(p) for p in npz['patient_ids']]
    if 'tcga' in cohort.lower():
        pids = [p.lower() for p in pids]
    mg = read_json(REPO / 'data/processed/human_full/graph.json')
    geom = build_biochem_subgraph(mg, hub_cap=hub_cap)
    if F.shape[1] != len(geom.nodes):
        for hc in (500, 200):
            geom = build_biochem_subgraph(mg, hub_cap=hc)
            if F.shape[1] == len(geom.nodes):
                break
        else:
            return None, None, None
    beta, alpha_norm, alpha_pc_scores, _ = decompose_beta_alpha(
        F, geom.log_pr, n_components=7)
    return beta, alpha_norm, alpha_pc_scores, pids


def per_pc_ablation(cohort_name, alpha_pcs, y_or_te, kind='auc'):
    print(f'\n=== {cohort_name} v5-canonical α-PC ablation ({kind.upper()}) ===',
          flush=True)
    print(f'  {"α-PC":<10}  Mean    ±std', flush=True)
    res = {}
    for k in range(alpha_pcs.shape[1]):
        score = alpha_pcs[:, k]
        r = cv_metric(score, y_or_te, f'α-PC{k+1}', kind=kind)
        if r:
            res[f'α-PC{k+1}'] = {'mean': r[0], 'std': r[1]}
            print(f'  α-PC{k+1:<7}  {r[0]:.3f}  ±{r[1]:.3f}', flush=True)
        else:
            res[f'α-PC{k+1}'] = None
            print(f'  α-PC{k+1:<7}  FAILED', flush=True)
    return res


def load_kmplot_surv():
    KMPLOT_SURV = Path('/home/jgardner/gitlab-old/'
                        'd2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/'
                        'KMPLOT_BRCA_EXPRESSION/DATA/KMPLOT_BRCA_SURVIVAL.txt')
    surv = pd.read_csv(KMPLOT_SURV, sep='\t')
    surv['AffyID'] = surv['AffyID'].astype(str)
    surv = surv.dropna(subset=['Death_event (1=death)', 'Death_time'])
    surv['event'] = surv['Death_event (1=death)'].astype(int)
    surv['time'] = pd.to_numeric(surv['Death_time'], errors='coerce')
    surv = surv.dropna(subset=['time'])
    surv = surv[surv['time'] > 0]
    return surv.rename(columns={'AffyID': 'patient_id'})


def load_luad_surv():
    CLIN = Path.home() / '.cache' / 'tcga_luad' / 'gdac.broadinstitute.org_LUAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0' / 'LUAD.clin.merged.picked.txt'
    cdf = pd.read_csv(CLIN, sep='\t', header=None, low_memory=False)
    attrs = cdf.iloc[:, 0].astype(str).tolist()
    pids = [str(p).strip().lower() for p in cdf.iloc[0, 1:].tolist()]
    vital = [str(v).strip() for v in cdf.iloc[attrs.index('vital_status'), 1:].tolist()]
    dtd = [str(v).strip() for v in cdf.iloc[attrs.index('days_to_death'), 1:].tolist()]
    dtf = [str(v).strip() for v in cdf.iloc[attrs.index('days_to_last_followup'), 1:].tolist()]
    rows = []
    for i, pid in enumerate(pids):
        try: v = int(float(vital[i]))
        except (ValueError, IndexError): continue
        if v == 1:
            try: t = float(dtd[i]); ev = 1
            except (ValueError, TypeError): continue
        elif v == 0:
            try: t = float(dtf[i]); ev = 0
            except (ValueError, TypeError): continue
        else: continue
        if t > 0 and np.isfinite(t):
            rows.append({'patient_id': pid, 'time': t, 'event': ev})
    return pd.DataFrame(rows)


def load_idh_surv():
    CLIN = Path.home() / '.cache' / 'tcga_idh' / 'lgggbm_tcga_pub_clinical.tsv'
    clin = pd.read_csv(CLIN, sep='\t')
    clin = clin.dropna(subset=['OS_MONTHS', 'OS_STATUS'])
    clin['event'] = clin['OS_STATUS'].astype(str).str.startswith('1').astype(int)
    clin['time'] = pd.to_numeric(clin['OS_MONTHS'], errors='coerce')
    clin = clin.dropna(subset=['time'])
    clin = clin[clin['time'] > 0]
    clin['patient_id'] = clin['patient_id'].astype(str).str.lower()
    return clin


def load_filbin_outcome():
    wb = openpyxl.load_workbook('/home/jgardner/.cache/filbin_mgh_covid/Clinical_Metadata.xlsx',
                                  read_only=True)
    ws = wb['Subject-level metadata']
    rows = list(ws.iter_rows(values_only=True))
    df = pd.DataFrame(rows[1:], columns=rows[0])
    df['patient_id'] = df['Public ID'].astype(str)
    df['acuity_28'] = pd.to_numeric(df['Acuity 28'], errors='coerce')
    df = df.dropna(subset=['acuity_28']).copy()
    df['died_28d'] = (df['acuity_28'] == 1).astype(int)
    return df[['patient_id', 'died_28d']]


def main():
    print('=== v7 v5-canonical F per-α-PC ablation across 4 cohorts ===',
          flush=True)
    t0 = time.time()

    out = {}
    # KMPLOT (Cox)
    beta, an, apc, pids = load_v5_alpha_pcs('KMPLOT_BRCA')
    if apc is not None:
        surv = load_kmplot_surv()
        common = sorted(set(surv['patient_id']) & set(pids))
        pid_to_idx = {p: i for i, p in enumerate(pids)}
        idx = np.array([pid_to_idx[p] for p in common])
        s = surv.set_index('patient_id').loc[common]
        print(f'KMPLOT: n={len(common)}, deaths={int(s["event"].sum())}',
              flush=True)
        out['KMPLOT_BRCA'] = per_pc_ablation(
            'KMPLOT_BRCA', apc[idx], (s['time'].values, s['event'].values.astype(int)),
            kind='cox')

    # TCGA_LUAD (Cox)
    beta, an, apc, pids = load_v5_alpha_pcs('TCGA_LUAD')
    if apc is not None:
        surv = load_luad_surv()
        common = sorted(set(surv['patient_id']) & set(pids))
        pid_to_idx = {p: i for i, p in enumerate(pids)}
        idx = np.array([pid_to_idx[p] for p in common])
        s = surv.set_index('patient_id').loc[common]
        print(f'TCGA_LUAD: n={len(common)}, deaths={int(s["event"].sum())}',
              flush=True)
        out['TCGA_LUAD'] = per_pc_ablation(
            'TCGA_LUAD', apc[idx], (s['time'].values, s['event'].values.astype(int)),
            kind='cox')

    # TCGA_IDH (Cox)
    beta, an, apc, pids = load_v5_alpha_pcs('TCGA_IDH_glioma')
    if apc is not None:
        surv = load_idh_surv()
        common = sorted(set(surv['patient_id']) & set(pids))
        pid_to_idx = {p: i for i, p in enumerate(pids)}
        idx = np.array([pid_to_idx[p] for p in common])
        s = surv.set_index('patient_id').loc[common]
        print(f'TCGA_IDH: n={len(common)}, deaths={int(s["event"].sum())}',
              flush=True)
        out['TCGA_IDH_glioma'] = per_pc_ablation(
            'TCGA_IDH_glioma', apc[idx],
            (s['time'].values, s['event'].values.astype(int)), kind='cox')

    # Filbin (AUC)
    beta, an, apc, pids = load_v5_alpha_pcs('Filbin_COVID')
    if apc is not None:
        out_clin = load_filbin_outcome()
        common = sorted(set(out_clin['patient_id']) & set(pids))
        pid_to_idx = {p: i for i, p in enumerate(pids)}
        idx = np.array([pid_to_idx[p] for p in common])
        s = out_clin.set_index('patient_id').loc[common]
        print(f'Filbin: n={len(common)}, deaths={int(s["died_28d"].sum())}',
              flush=True)
        out['Filbin_COVID'] = per_pc_ablation(
            'Filbin_COVID', apc[idx], s['died_28d'].values.astype(int),
            kind='auc')

    p = ZS_DIR / 'v7_v5_F_per_pc_ablation.json'
    p.write_text(json.dumps({**out, 'compute_seconds': time.time() - t0},
                              indent=2))
    print(f'\nWrote {p}', flush=True)


if __name__ == '__main__':
    main()
