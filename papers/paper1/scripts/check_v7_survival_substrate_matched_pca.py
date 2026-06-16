"""v7 Phase 5a refinement: substrate-matched PCA baseline.

The full PCA-on-raw uses ALL genes in the expression matrix; F only sees
substrate-mappable genes. A fair apples-to-apples test restricts PCA's
input to the same substrate-mappable subset that F uses.

For each cohort:
  - Identify substrate gene symbols (gene nodes in mg.graph.json)
  - Restrict expression matrix to substrate-mappable genes only
  - PCA-top-7 on the restricted matrix
  - 5-fold Cox CV → C-index
  - Compare to (a) full-raw PCA, (b) F-features

If substrate-matched-PCA drops toward F-features, F's underperformance was
due to input-universe asymmetry (not the substrate projection itself).
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
ZS_DIR = RESULTS / 'zscored'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json


def get_substrate_gene_symbols():
    """Return set of substrate-mappable gene symbols."""
    print('Loading substrate to extract gene symbols...', flush=True)
    mg = read_json(REPO / 'data/processed/human_full/graph.json')
    symbols = set()
    for n, a in mg.graph.nodes(data=True):
        if a.get('node_type') == 'gene':
            sym = a.get('symbol') or a.get('name') or n.replace('symbol:', '')
            if sym:
                symbols.add(sym)
    print(f'  {len(symbols)} unique gene symbols in substrate', flush=True)
    return symbols


def cox_cindex_cv(X, time_arr, event_arr, label, n_splits=5, seed=42):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    cidx = []
    for tr, te in kf.split(X):
        sc = StandardScaler()
        Xs_tr = sc.fit_transform(X[tr]); Xs_te = sc.transform(X[te])
        df = pd.DataFrame(Xs_tr, columns=[f'f{i}' for i in range(X.shape[1])])
        df['time'] = time_arr[tr]; df['event'] = event_arr[tr]
        cph = CoxPHFitter(penalizer=0.01)
        try:
            cph.fit(df, duration_col='time', event_col='event', show_progress=False)
        except Exception:
            continue
        rs = cph.predict_partial_hazard(
            pd.DataFrame(Xs_te, columns=[f'f{i}' for i in range(X.shape[1])])).values
        cidx.append(concordance_index(time_arr[te], -rs, event_arr[te]))
    if not cidx: return None
    return float(np.mean(cidx)), float(np.std(cidx))


def load_kmplot_setup():
    KMPLOT_SURV = Path('/home/jgardner/gitlab-old/'
                        'd2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/'
                        'KMPLOT_BRCA_EXPRESSION/DATA/KMPLOT_BRCA_SURVIVAL.txt')
    KMPLOT_XP = Path('/home/jgardner/gitlab-old/'
                      'd2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/'
                      'KMPLOT_BRCA_EXPRESSION/KMPLOT_BRCA_XP_NORMALIZED_CLEANED.tsv')
    surv = pd.read_csv(KMPLOT_SURV, sep='\t')
    surv['AffyID'] = surv['AffyID'].astype(str)
    surv = surv.dropna(subset=['Death_event (1=death)', 'Death_time'])
    surv['event'] = surv['Death_event (1=death)'].astype(int)
    surv['time'] = pd.to_numeric(surv['Death_time'], errors='coerce')
    surv = surv.dropna(subset=['time'])
    surv = surv[surv['time'] > 0].copy()
    ba = pd.read_csv(ZS_DIR / 'KMPLOT_BRCA' / 'beta_alpha.tsv', sep='\t')
    ba['patient_id'] = ba['patient_id'].astype(str)
    xp = pd.read_csv(KMPLOT_XP, sep='\t')
    xp = xp.rename(columns={xp.columns[0]: 'sample_id'})
    xp['sample_id'] = xp['sample_id'].astype(str)
    xp = xp.set_index('sample_id')
    common = sorted(set(surv['AffyID']) & set(ba['patient_id']) & set(xp.index))
    return (surv.set_index('AffyID').loc[common].reset_index(),
            ba.set_index('patient_id').loc[common].reset_index(),
            xp.loc[common])


def load_tcga_luad_setup():
    CLIN = Path.home() / '.cache' / 'tcga_luad' / 'gdac.broadinstitute.org_LUAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0' / 'LUAD.clin.merged.picked.txt'
    EXPR = Path.home() / '.cache' / 'tcga_luad' / 'gdac.broadinstitute.org_LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0' / 'LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt'
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
    surv = pd.DataFrame(rows)
    expr = pd.read_csv(EXPR, sep='\t', skiprows=[1])
    expr = expr.set_index(expr.columns[0])
    expr.columns = [c.lower()[:12] if len(c) >= 12 else c.lower()
                    for c in expr.columns]
    # Broad TCGA gene IDs use "SYMBOL|ENTREZ" format; strip the entrez suffix
    expr.index = [str(i).split('|')[0] for i in expr.index]
    expr = expr.groupby(level=0).mean()  # collapse duplicate symbols
    expr_t = expr.T.groupby(level=0).mean()
    ba = pd.read_csv(ZS_DIR / 'TCGA_LUAD' / 'beta_alpha.tsv', sep='\t')
    ba['patient_id'] = ba['patient_id'].astype(str)
    common = sorted(set(surv['patient_id']) & set(ba['patient_id']) & set(expr_t.index))
    return (surv.set_index('patient_id').loc[common].reset_index(),
            ba.set_index('patient_id').loc[common].reset_index(),
            expr_t.loc[common])


def load_tcga_idh_setup():
    CLIN = Path.home() / '.cache' / 'tcga_idh' / 'lgggbm_tcga_pub_clinical.tsv'
    EXPR = Path('/home/jgardner/gitlab-old/'
                 'd2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/'
                 'TCGA_BRAIN_EXPRESSION/TCGA_GBM_and_LGG_PREPROCESSED_RNASEQ_EXPRESSION.tsv')
    clin = pd.read_csv(CLIN, sep='\t')
    clin = clin.dropna(subset=['OS_MONTHS', 'OS_STATUS']).copy()
    clin['event'] = clin['OS_STATUS'].astype(str).str.startswith('1').astype(int)
    clin['time'] = pd.to_numeric(clin['OS_MONTHS'], errors='coerce')
    clin = clin.dropna(subset=['time'])
    clin = clin[clin['time'] > 0].copy()
    clin['patient_id'] = clin['patient_id'].astype(str).str.lower()
    expr = pd.read_csv(EXPR, sep='\t', index_col=0)
    if expr.shape[0] > expr.shape[1]:
        expr.columns = [c.lower().replace('.', '-') for c in expr.columns]
        expr_t = expr.T
    else:
        expr.index = [i.lower().replace('.', '-') for i in expr.index]
        expr_t = expr
    expr_t.index = [s[:12] if len(s) >= 12 else s for s in expr_t.index]
    expr_t = expr_t.groupby(expr_t.index).mean()
    ba = pd.read_csv(ZS_DIR / 'TCGA_IDH_glioma' / 'beta_alpha.tsv', sep='\t')
    ba['patient_id'] = ba['patient_id'].astype(str).str.lower()
    common = sorted(set(clin['patient_id']) & set(ba['patient_id']) & set(expr_t.index))
    return (clin.set_index('patient_id').loc[common].reset_index(),
            ba.set_index('patient_id').loc[common].reset_index(),
            expr_t.loc[common])


def evaluate_cohort(name, surv, ba, expr, substrate_syms, log_transform=False):
    print(f'\n=== {name} (n={len(surv)}, deaths={int(surv["event"].sum())}) ===',
          flush=True)
    time_arr = surv['time'].values.astype(float)
    event_arr = surv['event'].values.astype(int)
    # Clean expression: NaN → column mean
    expr_arr = expr.values.astype(float)
    col_mean = np.nanmean(expr_arr, axis=0)
    inds = np.where(np.isnan(expr_arr))
    expr_arr[inds] = np.take(col_mean, inds[1])
    if log_transform:
        expr_arr = np.log2(np.maximum(expr_arr, 0.0) + 1.0)

    # All-gene PCA (existing baseline)
    pca_full = PCA(n_components=7, random_state=42)
    X_PCA_full = pca_full.fit_transform(expr_arr)
    r_full = cox_cindex_cv(X_PCA_full, time_arr, event_arr, 'PCA-all-genes')

    # Substrate-matched PCA
    expr_cols = list(expr.columns)
    matched_idx = [i for i, c in enumerate(expr_cols) if c in substrate_syms]
    print(f'  expression genes: {len(expr_cols)}, '
          f'substrate-mappable: {len(matched_idx)}', flush=True)
    expr_matched = expr_arr[:, matched_idx]
    pca_sm = PCA(n_components=7, random_state=42)
    X_PCA_sm = pca_sm.fit_transform(expr_matched)
    r_sm = cox_cindex_cv(X_PCA_sm, time_arr, event_arr, 'PCA-substrate-matched')

    # F-features (all 7)
    f_cols = ['beta', 'alpha_norm', 'alpha_pc1', 'alpha_pc2', 'alpha_pc3',
              'alpha_pc4', 'alpha_pc5']
    X_F = ba[f_cols].values.astype(float)
    r_F = cox_cindex_cv(X_F, time_arr, event_arr, 'F-features')

    out = {}
    print(f'  {"Method":<32}  Mean C  ±std', flush=True)
    for label, r in [('F-features (β + α-PC1..5 + ‖α‖₂)', r_F),
                      ('PCA-substrate-matched (top-7)', r_sm),
                      ('PCA-all-genes (top-7)', r_full)]:
        if r:
            print(f'  {label:<32}  {r[0]:.3f}  ±{r[1]:.3f}', flush=True)
            out[label] = {'cindex_mean': r[0], 'cindex_std': r[1]}
        else:
            out[label] = None
    out['n_substrate_mapped_genes'] = len(matched_idx)
    out['n_total_genes'] = len(expr_cols)
    return out


def main():
    print('=== v7 Phase 5a substrate-matched PCA test ===', flush=True)
    t0 = time.time()
    sub_syms = get_substrate_gene_symbols()

    print('\nLoading KMPLOT...', flush=True)
    s, b, e = load_kmplot_setup()
    r_kmplot = evaluate_cohort('KMPLOT_BRCA', s, b, e, sub_syms,
                                  log_transform=False)

    print('\nLoading TCGA_LUAD...', flush=True)
    s, b, e = load_tcga_luad_setup()
    r_luad = evaluate_cohort('TCGA_LUAD', s, b, e, sub_syms,
                                log_transform=True)

    print('\nLoading TCGA_IDH_glioma...', flush=True)
    s, b, e = load_tcga_idh_setup()
    r_idh = evaluate_cohort('TCGA_IDH_glioma', s, b, e, sub_syms,
                               log_transform=False)

    out = ZS_DIR / 'v7_survival_substrate_matched.json'
    out.write_text(json.dumps({
        'KMPLOT_BRCA': r_kmplot, 'TCGA_LUAD': r_luad,
        'TCGA_IDH_glioma': r_idh, 'compute_seconds': time.time() - t0,
    }, indent=2))
    print(f'\nWrote {out}', flush=True)


if __name__ == '__main__':
    main()
