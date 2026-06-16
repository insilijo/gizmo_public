"""v7 Phase 5a v3: v5-canonical F survival on all 3 cohorts (KMPLOT + LUAD + IDH).

Confirms whether the TCGA_IDH rescue (v5-canonical 0.824 vs z-scored 0.790)
generalizes to KMPLOT_BRCA and TCGA_LUAD.
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
ZS_DIR = RESULTS / 'zscored'
SNAPSHOT = RESULTS / '_pre_zscore_snapshot_20260607'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph, decompose_beta_alpha


def cox_cindex_cv(X, time_arr, event_arr, label, n_splits=5, seed=42):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    cidx = []
    fold = 0
    for tr, te in kf.split(X):
        fold += 1
        sc = StandardScaler()
        Xs_tr = sc.fit_transform(X[tr]); Xs_te = sc.transform(X[te])
        df = pd.DataFrame(Xs_tr, columns=[f'f{i}' for i in range(X.shape[1])])
        df['time'] = time_arr[tr]; df['event'] = event_arr[tr]
        cph = CoxPHFitter(penalizer=0.01)
        try:
            cph.fit(df, duration_col='time', event_col='event', show_progress=False)
        except Exception as e:
            print(f'  {label} fold {fold}: fit error: {e}', flush=True)
            continue
        df_te = pd.DataFrame(Xs_te, columns=[f'f{i}' for i in range(X.shape[1])])
        rs = cph.predict_partial_hazard(df_te).values
        ci = concordance_index(time_arr[te], -rs, event_arr[te])
        cidx.append(ci)
        print(f'  {label} fold {fold}: C={ci:.3f}', flush=True)
    if not cidx: return None
    return float(np.mean(cidx)), float(np.std(cidx))


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
    surv = surv[surv['time'] > 0].copy()
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
    clin = clin[clin['time'] > 0].copy()
    clin['patient_id'] = clin['patient_id'].astype(str).str.lower()
    return clin


def process_cohort(name, F_file, surv_loader, mg, hub_cap=200):
    print(f'\n=== {name} (v5-canonical preprocessing) ===', flush=True)
    f_path = SNAPSHOT / F_file
    if not f_path.exists():
        print(f'  SKIP: {f_path} missing', flush=True)
        return None
    npz = np.load(f_path, allow_pickle=True)
    F = npz['F']
    # Case-preserve for KMPLOT (GSM IDs); lowercase for TCGA (case-insensitive barcodes)
    if 'tcga' in name.lower():
        pids = [str(p).lower() for p in npz['patient_ids']]
    else:
        pids = [str(p) for p in npz['patient_ids']]
    print(f'  F: {F.shape}', flush=True)

    geom = build_biochem_subgraph(mg, hub_cap=hub_cap)
    if F.shape[1] != len(geom.nodes):
        for hc in (500, 200):
            geom = build_biochem_subgraph(mg, hub_cap=hc)
            if F.shape[1] == len(geom.nodes):
                print(f'  matched at hub_cap={hc}', flush=True)
                break
        else:
            print(f'  WARNING: F has {F.shape[1]} cols, no hub_cap match', flush=True)
            return None
    log_pr = geom.log_pr

    beta, alpha_norm, alpha_pc_scores, pca = decompose_beta_alpha(
        F, log_pr, n_components=5)
    print(f'  EVR top-5: {pca.explained_variance_ratio_.round(3)}', flush=True)
    f_feat = np.column_stack([beta, alpha_norm, alpha_pc_scores])

    surv = surv_loader()
    pid_to_idx = {p: i for i, p in enumerate(pids)}
    common = sorted(set(surv['patient_id']) & set(pid_to_idx))
    print(f'  common F ∩ surv: {len(common)}', flush=True)
    if len(common) < 80:
        print(f'  SKIP: only {len(common)}', flush=True)
        return None

    idx = np.array([pid_to_idx[p] for p in common])
    sub = f_feat[idx]
    s = surv.set_index('patient_id').loc[common]
    time_arr = s['time'].values.astype(float)
    event_arr = s['event'].values.astype(int)
    print(f'  n={len(common)}, deaths={int(event_arr.sum())}', flush=True)
    res = cox_cindex_cv(sub, time_arr, event_arr, f'  F-v5-{name}')
    return res


def main():
    print('=== v7 Phase 5a v3: v5-canonical F survival on 3 cohorts ===',
          flush=True)
    t0 = time.time()

    print('Loading substrate...', flush=True)
    mg = read_json(REPO / 'data/processed/human_full/graph.json')

    cohorts = [
        ('KMPLOT_BRCA', 'stage3_F_KMPLOT_BRCA_edge_informed.npz', load_kmplot_surv),
        ('TCGA_LUAD', 'stage3_F_TCGA_LUAD_edge_informed.npz', load_luad_surv),
        ('TCGA_IDH_glioma', 'stage3_F_TCGA_IDH_glioma_edge_informed.npz', load_idh_surv),
    ]
    results = {}
    for name, fn, sl in cohorts:
        results[name] = process_cohort(name, fn, sl, mg, hub_cap=200)

    print('\n=== Summary: F survival under different preprocessings ===', flush=True)
    PRIOR = {
        'KMPLOT_BRCA': {'zscored': 0.451, 'sm_pca': 0.589, 'all_pca': 0.596},
        'TCGA_LUAD':   {'zscored': 0.580, 'sm_pca': 0.616, 'all_pca': 0.600},
        'TCGA_IDH_glioma': {'zscored': 0.790, 'sm_pca': 0.819, 'all_pca': 0.821},
    }
    print(f'  {"Cohort":<18} {"F-v5-canonical":>15} {"F-zscored":>11} '
          f'{"SM-PCA":>8} {"All-PCA":>8}', flush=True)
    for name in PRIOR:
        r = results.get(name)
        v5_str = f'{r[0]:.3f}±{r[1]:.3f}' if r else '   FAILED   '
        p = PRIOR[name]
        print(f'  {name:<18} {v5_str:>15} {p["zscored"]:>11.3f} '
              f'{p["sm_pca"]:>8.3f} {p["all_pca"]:>8.3f}', flush=True)

    out = ZS_DIR / 'v7_v5_canonical_F_all3.json'
    out_data = {}
    for name in PRIOR:
        r = results.get(name)
        out_data[name] = {
            'F_v5_canonical_mean': r[0] if r else None,
            'F_v5_canonical_std': r[1] if r else None,
            'F_zscored_prior': PRIOR[name]['zscored'],
            'PCA_substrate_matched_prior': PRIOR[name]['sm_pca'],
            'PCA_all_genes_prior': PRIOR[name]['all_pca'],
        }
    out_data['compute_seconds'] = time.time() - t0
    out.write_text(json.dumps(out_data, indent=2))
    print(f'\nWrote {out}', flush=True)


if __name__ == '__main__':
    main()
