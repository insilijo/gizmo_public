"""v7: test F-α-PCs as a feature-selection mechanism (raw PCA on F-selected genes).

For each cohort:
  1. Sanity check: PCA on v5-canonical-rescaled input == PCA on raw log input
     (true if cohort is single-modality, because global-std is uniform scaling).
  2. Feature-selection test: take the best F-α-PC, extract top-K most heavily
     loaded gene nodes, restrict raw expression to those genes, do PCA on
     the restricted matrix, compute Cox/AUC discrimination.

If "PCA on F-selected genes" achieves discrimination close to F-α-PC score,
then F's contribution is identifying the right gene set (feature selection);
the substrate-smoothing is the *mechanism* by which F finds those genes but
not the *source* of the discrimination.

If "PCA on F-selected genes" stays close to "PCA on all genes," then F's
contribution is the smoothing-transformed signal (not just selection).
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl
from sklearn.decomposition import PCA
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


def cv_metric(X, y_or_te, kind='auc', n_splits=5, seed=42):
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
            y = y_or_te
            clf = LogisticRegression(max_iter=1000, C=1.0)
            try:
                clf.fit(Xs_tr, y[tr])
                proba = clf.predict_proba(Xs_te)[:, 1]
                scores.append(roc_auc_score(y[te], proba))
            except Exception:
                continue
        elif kind == 'cox':
            from lifelines import CoxPHFitter
            from lifelines.utils import concordance_index
            time_arr, event_arr = y_or_te
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


def load_v5_F_with_alpha(cohort):
    """Load v5-canonical F + decompose, return F, pids, beta, alpha_norm, alpha_pcs, pca_components."""
    f_path = SNAPSHOT / f'stage3_F_{cohort}_edge_informed.npz'
    if not f_path.exists():
        f_path = SNAPSHOT / f'stage3_F_{cohort}.npz'
    npz = np.load(f_path, allow_pickle=True)
    F = npz['F']
    pids = [str(p) for p in npz['patient_ids']]
    if 'tcga' in cohort.lower():
        pids = [p.lower() for p in pids]
    mg = read_json(REPO / 'data/processed/human_full/graph.json')
    for hc in (200, 500):
        geom = build_biochem_subgraph(mg, hub_cap=hc)
        if F.shape[1] == len(geom.nodes):
            break
    else:
        return None
    log_pr = geom.log_pr
    # Recompute α/PCA explicitly to get components
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm = x.mean(); xv = x.var() + 1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm) * (x - xm)).mean(axis=1, keepdims=True)
    beta = (cov / xv).ravel()
    alpha = F_unit - Fm - beta[:, None] * (x - xm)[None, :]
    pca_alpha = PCA(n_components=min(7, alpha.shape[0]), random_state=0)
    alpha_pcs = pca_alpha.fit_transform(alpha)
    alpha_components = pca_alpha.components_  # (n_pcs, n_nodes)

    # Map node IDs to gene symbols
    node_sym = {}
    for nid in geom.nodes:
        attrs = mg.graph.nodes.get(nid, {})
        if attrs.get('node_type') == 'gene':
            sym = attrs.get('symbol') or attrs.get('name') or nid.replace('symbol:', '')
            node_sym[nid] = sym

    return F, pids, beta, alpha, alpha_pcs, alpha_components, geom.nodes, node_sym, mg


def top_loaded_genes_for_pc(alpha_components, node_ids, node_sym, pc_idx,
                             top_k=50, restrict_to=None):
    """Return top-K gene symbols by |loading| on a given α-PC.

    If restrict_to is a set of gene symbols, only those genes are considered
    (used to restrict to OBSERVED genes in the expression matrix, filtering
    out propagation-only substrate nodes that have no input data).
    """
    loadings = alpha_components[pc_idx]
    abs_load = np.abs(loadings)
    order = np.argsort(-abs_load)
    out = []
    for i in order:
        nid = node_ids[i]
        sym = node_sym.get(nid)
        if not sym:
            continue
        if restrict_to is not None and sym not in restrict_to:
            continue
        out.append((sym, float(loadings[i])))
        if len(out) >= top_k:
            break
    return out


def load_kmplot_expr_surv():
    KMPLOT_XP = Path('/home/jgardner/gitlab-old/'
                      'd2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/'
                      'KMPLOT_BRCA_EXPRESSION/KMPLOT_BRCA_XP_NORMALIZED_CLEANED.tsv')
    KMPLOT_SURV = Path('/home/jgardner/gitlab-old/'
                        'd2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/'
                        'KMPLOT_BRCA_EXPRESSION/DATA/KMPLOT_BRCA_SURVIVAL.txt')
    xp = pd.read_csv(KMPLOT_XP, sep='\t')
    xp = xp.rename(columns={xp.columns[0]: 'sample_id'})
    xp['sample_id'] = xp['sample_id'].astype(str)
    xp = xp.set_index('sample_id')
    surv = pd.read_csv(KMPLOT_SURV, sep='\t')
    surv['AffyID'] = surv['AffyID'].astype(str)
    surv = surv.dropna(subset=['Death_event (1=death)', 'Death_time'])
    surv['event'] = surv['Death_event (1=death)'].astype(int)
    surv['time'] = pd.to_numeric(surv['Death_time'], errors='coerce')
    surv = surv.dropna(subset=['time'])
    surv = surv[surv['time'] > 0]
    return xp, surv.rename(columns={'AffyID': 'patient_id'})


def load_tcga_idh_expr_surv():
    CLIN = Path.home() / '.cache' / 'tcga_idh' / 'lgggbm_tcga_pub_clinical.tsv'
    EXPR = Path('/home/jgardner/gitlab-old/'
                 'd2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/'
                 'TCGA_BRAIN_EXPRESSION/TCGA_GBM_and_LGG_PREPROCESSED_RNASEQ_EXPRESSION.tsv')
    clin = pd.read_csv(CLIN, sep='\t')
    clin = clin.dropna(subset=['OS_MONTHS', 'OS_STATUS'])
    clin['event'] = clin['OS_STATUS'].astype(str).str.startswith('1').astype(int)
    clin['time'] = pd.to_numeric(clin['OS_MONTHS'], errors='coerce')
    clin = clin.dropna(subset=['time'])
    clin = clin[clin['time'] > 0]
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
    return expr_t, clin


def test_cohort(cohort_name, best_pc_1based, kind, top_k=50):
    print(f'\n=== {cohort_name} F-α-PC{best_pc_1based} as feature selector (top {top_k}) ===',
          flush=True)
    loaded = load_v5_F_with_alpha(cohort_name)
    if loaded is None:
        print(f'  load failed', flush=True)
        return None
    F, pids, beta, alpha, alpha_pcs, alpha_comps, node_ids, node_sym, mg = loaded

    # Load expression + survival FIRST so we can restrict top-loaded-gene
    # ranking to the OBSERVED gene set (filtering out propagation-only nodes)
    if cohort_name == 'KMPLOT_BRCA':
        xp, surv = load_kmplot_expr_surv()
        kind = 'cox'
        clin_id = 'patient_id'
        time_col = 'time'; event_col = 'event'
        log_transform = False
    elif cohort_name == 'TCGA_IDH_glioma':
        xp, surv = load_tcga_idh_expr_surv()
        kind = 'cox'
        clin_id = 'patient_id'
        time_col = 'time'; event_col = 'event'
        log_transform = False
    else:
        print(f'  {cohort_name} not implemented in test_cohort yet', flush=True)
        return None

    # Build set of observed gene symbols from expression matrix
    expr_genes_set = set(xp.columns)
    # Handle TCGA SYMBOL|ENTREZ format
    expr_genes_set |= {str(c).split('|')[0] for c in xp.columns if '|' in str(c)}

    pc_idx = best_pc_1based - 1
    top_genes = top_loaded_genes_for_pc(alpha_comps, node_ids, node_sym,
                                          pc_idx, top_k=top_k,
                                          restrict_to=expr_genes_set)
    top_symbols = [g for g, _ in top_genes]
    print(f'  Top-{top_k} OBSERVED F-α-PC{best_pc_1based} loadings (sample 10): '
          f'{[g for g, _ in top_genes[:10]]}', flush=True)

    # Intersect patients
    common = sorted(set(surv[clin_id]) & set(xp.index) & set(pids))
    print(f'  common (F ∩ expr ∩ surv): {len(common)}', flush=True)
    if len(common) < 50:
        print(f'  too few patients', flush=True)
        return None

    s = surv.set_index(clin_id).loc[common]
    xp_common = xp.loc[common]

    # Strip TCGA SYMBOL|ENTREZ if needed
    if any('|' in str(c) for c in xp_common.columns[:5]):
        xp_common.columns = [str(c).split('|')[0] for c in xp_common.columns]
        xp_common = xp_common.groupby(level=0, axis=1).mean()
    X_full = xp_common.values.astype(float)
    if log_transform:
        X_full = np.log2(np.maximum(X_full, 0.0) + 1.0)
    col_mean = np.nanmean(X_full, axis=0)
    inds = np.where(np.isnan(X_full))
    X_full[inds] = np.take(col_mean, inds[1])

    # Restrict to F-α-PC top loaded genes (those present in expression)
    expr_genes = list(xp_common.columns)
    expr_gene_set = set(expr_genes)
    selected = [s for s in top_symbols if s in expr_gene_set]
    selected_idx = [expr_genes.index(s) for s in selected]
    print(f'  F-selected genes in expression matrix: '
          f'{len(selected)}/{len(top_symbols)}', flush=True)
    if len(selected) < 5:
        print(f'  too few selected genes; sample: {selected}', flush=True)
        return None

    X_sel = X_full[:, selected_idx]
    print(f'  X_selected: {X_sel.shape}', flush=True)

    # 5-fold CV on PCA-top-1 of F-selected genes
    pca_sel = PCA(n_components=min(7, X_sel.shape[1]), random_state=42)
    X_sel_pcs = pca_sel.fit_transform(X_sel)
    print(f'  EVR top-7 on F-selected: {pca_sel.explained_variance_ratio_.round(3)}',
          flush=True)

    if kind == 'cox':
        te = (s[time_col].values, s[event_col].values.astype(int))
    else:
        te = s[event_col].values.astype(int)

    # Test each PC individually + all 7 together
    print(f'  {"Subset":<22}  Mean  ±std', flush=True)
    for k in range(min(7, X_sel.shape[1])):
        r = cv_metric(X_sel_pcs[:, k], te, kind=kind)
        if r:
            print(f'  PC{k+1}-of-F-selected     {r[0]:.3f}  ±{r[1]:.3f}',
                  flush=True)
    r_all = cv_metric(X_sel_pcs, te, kind=kind)
    if r_all:
        print(f'  All-7-of-F-selected     {r_all[0]:.3f}  ±{r_all[1]:.3f}',
              flush=True)

    # Best single PC + corresponding F-α-PC value for context
    pca_pcs_aucs = []
    for k in range(min(7, X_sel.shape[1])):
        r = cv_metric(X_sel_pcs[:, k], te, kind=kind)
        if r: pca_pcs_aucs.append(r[0])
    best_pca_sel = max(pca_pcs_aucs) if pca_pcs_aucs else None

    return {
        'cohort': cohort_name,
        'F_alpha_pc_used_for_selection': best_pc_1based,
        'n_selected_in_expr': len(selected),
        'top_selected_symbols': selected[:20],
        'pca_on_F_selected_best_PC': best_pca_sel,
        'pca_on_F_selected_all_7': r_all[0] if r_all else None,
    }


def main():
    print('=== v7 F-as-feature-selection test ===', flush=True)
    t0 = time.time()

    # TCGA_IDH: best F-α-PC = α-PC4 (C=0.762 under v5-canonical)
    res_idh = test_cohort('TCGA_IDH_glioma', best_pc_1based=4,
                            kind='cox', top_k=50)
    # KMPLOT: best F-α-PC = α-PC3 (C=0.603 under v5-canonical)
    res_kmplot = test_cohort('KMPLOT_BRCA', best_pc_1based=3,
                                kind='cox', top_k=50)

    out_path = ZS_DIR / 'v7_f_as_feature_selection.json'
    out_path.write_text(json.dumps({
        'TCGA_IDH_glioma': res_idh,
        'KMPLOT_BRCA': res_kmplot,
        'compute_seconds': time.time() - t0,
    }, indent=2))
    print(f'\nWrote {out_path}', flush=True)


if __name__ == '__main__':
    main()
