"""Discrimination ablation: F-features (β + α-PC1..5 + ||α||_2) vs PCA-on-F top-7.

If F-features C-index ≈ PCA-on-F C-index, β/α adds no discrimination value over direct PCA.
If F-features > PCA-on-F, β/α provides extra useful features.
"""
import sys, json
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, '/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'
mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr

# Loaders (re-using minimal from prior scripts)
def load_kmplot():
    KMPLOT_SURV = Path('/home/jgardner/gitlab-old/d2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/KMPLOT_BRCA_EXPRESSION/DATA/KMPLOT_BRCA_SURVIVAL.txt')
    surv = pd.read_csv(KMPLOT_SURV, sep='\t')
    surv['AffyID'] = surv['AffyID'].astype(str)
    surv = surv.dropna(subset=['Death_event (1=death)','Death_time']).copy()
    surv['event'] = surv['Death_event (1=death)'].astype(int)
    surv['time'] = pd.to_numeric(surv['Death_time'], errors='coerce')
    surv = surv[(surv['time'].notna()) & (surv['time']>0)]
    return surv.set_index('AffyID')[['time','event']]

def load_luad():
    CLIN = Path.home()/'.cache'/'tcga_luad'/'gdac.broadinstitute.org_LUAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0'/'LUAD.clin.merged.picked.txt'
    cdf = pd.read_csv(CLIN, sep='\t', header=None, low_memory=False)
    attrs = cdf.iloc[:,0].astype(str).tolist()
    pids = [str(p).strip().lower() for p in cdf.iloc[0,1:].tolist()]
    vital = [str(v).strip() for v in cdf.iloc[attrs.index('vital_status'),1:].tolist()]
    dtd = [str(v).strip() for v in cdf.iloc[attrs.index('days_to_death'),1:].tolist()]
    dtf = [str(v).strip() for v in cdf.iloc[attrs.index('days_to_last_followup'),1:].tolist()]
    rows = []
    for i,pid in enumerate(pids):
        try: v = int(float(vital[i]))
        except (ValueError,IndexError): continue
        if v==1:
            try: t=float(dtd[i]); ev=1
            except (ValueError,TypeError): continue
        elif v==0:
            try: t=float(dtf[i]); ev=0
            except (ValueError,TypeError): continue
        else: continue
        if t>0 and np.isfinite(t):
            rows.append({'patient_id':pid,'time':t,'event':ev})
    return pd.DataFrame(rows).set_index('patient_id')

def load_idh():
    CLIN = Path.home()/'.cache'/'tcga_idh'/'lgggbm_tcga_pub_clinical.tsv'
    c = pd.read_csv(CLIN,sep='\t').dropna(subset=['OS_MONTHS','OS_STATUS']).copy()
    c['event'] = c['OS_STATUS'].astype(str).str.startswith('1').astype(int)
    c['time'] = pd.to_numeric(c['OS_MONTHS'],errors='coerce')
    c = c[(c['time'].notna())&(c['time']>0)]
    c['patient_id'] = c['patient_id'].astype(str).str.lower()
    return c.set_index('patient_id')[['time','event']]

def cox_cv(X, time_arr, event_arr, n_splits=5, seed=42):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    if X.ndim==1: X = X.reshape(-1,1)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    scores = []
    for tr, te in kf.split(X):
        sc=StandardScaler(); Xs_tr=sc.fit_transform(X[tr]); Xs_te=sc.transform(X[te])
        df = pd.DataFrame(Xs_tr,columns=[f'f{i}' for i in range(X.shape[1])])
        df['time']=time_arr[tr]; df['event']=event_arr[tr]
        cph = CoxPHFitter(penalizer=0.01)
        try: cph.fit(df,duration_col='time',event_col='event',show_progress=False)
        except: continue
        df_te = pd.DataFrame(Xs_te,columns=[f'f{i}' for i in range(X.shape[1])])
        rs = cph.predict_partial_hazard(df_te).values
        scores.append(concordance_index(time_arr[te],-rs,event_arr[te]))
    return (float(np.mean(scores)), float(np.std(scores))) if scores else None


cohort_loaders = {
    'KMPLOT_BRCA': (load_kmplot, False),
    'TCGA_LUAD': (load_luad, True),
    'TCGA_IDH_glioma': (load_idh, True),
}

print(f'{"Cohort":<20} {"F-features 7":<15} {"PCA-on-F top-7":<15} {"Δ":<8}')
print('-'*65)
for c, (loader, tcga_lower) in cohort_loaders.items():
    fp = SNAPSHOT / f'stage3_F_{c}_edge_informed.npz'
    if not fp.exists(): continue
    npz = np.load(fp, allow_pickle=True)
    F = npz['F']
    pids = [str(p).lower() if tcga_lower else str(p) for p in npz['patient_ids']]

    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm=x.mean(); xv=x.var()+1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm)*(x-xm)).mean(axis=1, keepdims=True)
    beta = (cov/xv).ravel()
    alpha = F_unit - Fm - beta[:,None]*(x-xm)[None,:]
    alpha_norm = np.linalg.norm(alpha, axis=1)
    pca_a = PCA(n_components=5, random_state=0)
    alpha_pcs = pca_a.fit_transform(alpha)
    f_feat = np.column_stack([beta, alpha_norm, alpha_pcs])  # 7 features

    pca_f = PCA(n_components=7, random_state=0)
    f_pcs = pca_f.fit_transform(F_unit)  # 7 features

    surv = loader()
    common = sorted(set(pids) & set(surv.index))
    if len(common) < 50: continue
    pid_idx = {p:i for i,p in enumerate(pids)}
    idx = np.array([pid_idx[p] for p in common])
    t_arr = surv.loc[common]['time'].values
    e_arr = surv.loc[common]['event'].values.astype(int)

    r_ff = cox_cv(f_feat[idx], t_arr, e_arr)
    r_fpca = cox_cv(f_pcs[idx], t_arr, e_arr)
    print(f'{c:<20} {r_ff[0]:.3f}±{r_ff[1]:.3f}  {r_fpca[0]:.3f}±{r_fpca[1]:.3f}  '
          f'{r_ff[0]-r_fpca[0]:+.3f}')
