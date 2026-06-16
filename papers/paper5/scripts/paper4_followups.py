"""Paper 4 follow-ups (items 1, 3, 4):

(1) Cross-cohort α-PC axis correspondence. α-PCs from each cohort are vectors
    in substrate-node space (R^86826). Compute |cos(PC_a_k, PC_b_l)| across
    cohorts; high cosine means the same biological axis is recovered in both.

(3) Per-cluster mechanism naming. For Wang RA's variance-normalized 5D
    k=3 clustering, project each cluster centroid back into substrate-node
    space (via the α-PC components) and report top-loading substrate nodes.
    The 6-patient "PC2-deficit non-responder" cluster gets a named mechanism.

(4) Clinical-baseline benchmark. For Wang RA, train logistic regression with
    LOOCV using:
      (a) clinical features only (DAS28-CRP, VAS, TJC, SJC, CRP, Age,
          anti-CCP, Gender)
      (b) α-PC_T0 only (5 features)
      (c) clinical + α-PC_T0 (additive)
    Compare LOOCV AUC. (c) > (a) means α-PC is additive on clinical baseline.

Re-solves MAP per cohort to obtain α-PC components (not saved in prior runs).
"""
from __future__ import annotations
import sys, json, gzip
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph, solve_map, ModalitySetup
from loocv_prestrat_audit import (
    cohort_wang_ra_baseline_for_response, cohort_filbin_d0_predicts_outcome,
    cohort_tb_dx_predicts_cure, geom, mg,
)


def safe_auc(y, s):
    a = roc_auc_score(y, s); return max(a, 1 - a)


def get_alpha_components(F, log_pr, n_components=5):
    lpn = log_pr / (np.linalg.norm(log_pr) + 1e-9)
    beta = F @ lpn
    alpha = F - np.outer(beta, lpn)
    pca = PCA(n_components=n_components, svd_solver='randomized', random_state=0).fit(alpha)
    return pca.components_, pca.transform(alpha), beta, pca.explained_variance_ratio_


def display_name(node_id):
    a = mg.graph.nodes.get(node_id, {})
    nm = a.get('name', '')
    if nm: return nm[:32]
    return node_id.replace('symbol:', '').replace('ENSG:', '')[:32]


# ---------------- Item 1: cross-cohort α-PC axis cosine ----------------
def cross_cohort_axis_cosine(components_by_cohort):
    """For each pair of cohorts, compute |cos(PC_a_k, PC_b_l)| over substrate.

    Components live on the FULL substrate (R^n_nodes) so cosines are directly
    comparable. Reports top matching PC pair per cohort pair.
    """
    print('\n' + '=' * 70)
    print('Item 1: Cross-cohort α-PC axis correspondence')
    print('=' * 70)
    results = {}
    names = list(components_by_cohort.keys())
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            a_name = names[i]; b_name = names[j]
            A = components_by_cohort[a_name]; B = components_by_cohort[b_name]
            # cosine matrix (5 x 5)
            cos_mat = np.zeros((A.shape[0], B.shape[0]))
            for k in range(A.shape[0]):
                for l in range(B.shape[0]):
                    na = np.linalg.norm(A[k]); nb = np.linalg.norm(B[l])
                    if na == 0 or nb == 0: continue
                    cos_mat[k, l] = abs(np.dot(A[k], B[l]) / (na * nb))
            best_per_row = cos_mat.max(axis=1)
            best_idx_per_row = cos_mat.argmax(axis=1)
            print(f'\n  {a_name} vs {b_name}:')
            for k in range(A.shape[0]):
                print(f'    PC{k+1}_{a_name[:6]} ↔ PC{best_idx_per_row[k]+1}_{b_name[:6]}: '
                      f'|cos| = {best_per_row[k]:.3f}')
            results[f'{a_name}__{b_name}'] = {
                'cos_matrix': cos_mat.tolist(),
                'best_match_per_a_PC': [(int(best_idx_per_row[k])+1, float(best_per_row[k]))
                                          for k in range(A.shape[0])],
            }
    return results


# -------- Item 3: per-cluster mechanism naming (Wang RA) ----------
def cluster_mechanism(cohort, alpha_scores, components, labels, feat_node_ids,
                       n_clusters=3, top_nodes=12):
    """Cluster in variance-normalized 5D α-space; project centroids back
    into substrate-node space; name top loading nodes per cluster.
    """
    print('\n' + '=' * 70)
    print(f'Item 3: Cluster mechanism naming — {cohort}')
    print('=' * 70)
    A = alpha_scores  # (N, 5)
    Az = (A - A.mean(axis=0)) / (A.std(axis=0) + 1e-9)
    km = KMeans(n_clusters=n_clusters, random_state=0, n_init=20).fit(Az)
    centroids_z = km.cluster_centers_  # (n_clusters, 5) in Z-space

    # Un-Z each centroid coord to recover original-units α-PC centroid
    centroids = centroids_z * (A.std(axis=0) + 1e-9) + A.mean(axis=0)

    # Project each centroid into substrate-node space:
    # substrate_centroid_c = sum_k centroid[c, k] * components[k]   (n_nodes,)
    substrate_centroids = centroids @ components  # (n_clusters, n_nodes)

    out = {'cluster_assignment': km.labels_.tolist(),
            'cluster_sizes': [int((km.labels_ == c).sum()) for c in range(n_clusters)],
            'centroids_alpha_pc_space': centroids.tolist(),
            'frac_outcome1_per_cluster': [],
            'top_nodes_per_cluster': []}

    for c in range(n_clusters):
        mask = km.labels_ == c
        frac1 = float(labels[mask].sum() / mask.sum()) if mask.sum() > 0 else 0
        out['frac_outcome1_per_cluster'].append(frac1)
        sc = substrate_centroids[c]
        top_idx = np.argsort(-np.abs(sc))[:top_nodes]
        # Map node index → node_id via geom.nodes
        top_nodes_info = []
        for ni in top_idx:
            nid = geom.nodes[ni]
            top_nodes_info.append({
                'node_id': nid, 'name': display_name(nid),
                'centroid_loading': float(sc[ni]),
            })
        n_size = int(mask.sum())
        dominant_pc = int(np.argmax(np.abs(centroids[c]))) + 1
        dom_sign = '+' if centroids[c, dominant_pc-1] > 0 else '−'
        print(f'\n  Cluster {c} (n={n_size}, outcome=1 fraction={frac1:.2f}):')
        print(f'    Dominant α-PC axis: {dom_sign}PC{dominant_pc}')
        print(f'    Top substrate nodes (signed loading on centroid):')
        for tn in top_nodes_info[:8]:
            print(f'      {tn["name"]:<35} {tn["centroid_loading"]:+.4f}')
        out['top_nodes_per_cluster'].append({
            'cluster_id': c, 'n_patients': n_size,
            'frac_outcome1': frac1,
            'dominant_alpha_PC': dominant_pc,
            'dominant_sign': dom_sign,
            'top_nodes': top_nodes_info,
        })
    return out


# ---- Item 4: clinical-baseline benchmark (Wang RA) ----
def clinical_baseline_benchmark(alpha_scores, labels):
    """LOO-CV logistic regression: clinical / α-only / clinical+α.

    Returns AUC for each variant + relative gain.
    """
    print('\n' + '=' * 70)
    print('Item 4: Clinical-baseline benchmark — Wang RA')
    print('=' * 70)
    WANG = REPO / 'data/cohorts/Wang_RA_MTX/hesy1191569605-rheumatoid-arthritis-0c94b6d'
    df = pd.read_csv(WANG / 'Figure6/csv/RA_DATAKNN1.csv')
    df['pid'] = df['Sample'].astype(str).str[1:].str.lstrip('0')
    df['tp'] = df['Sample'].astype(str).str[0]
    by_pid = df.groupby('pid').agg({
        'tp': list, 'Drugs Response': 'first', 'DAS28-CRP': 'first',
        'VAS': 'first', 'TJC': 'first', 'SJC': 'first', 'CRP': 'first',
        'Age': 'first', 'Gender': 'first',
        'Anti-citrullinated peptide antibodies': 'first', 'Class': 'first',
    })
    paired_pids = [p for p, r in by_pid.iterrows()
                   if 'A' in r['tp'] and 'B' in r['tp']
                   and r['Drugs Response'] in ('Response', 'No Response')]
    print(f'  Paired patients: {len(paired_pids)}')

    # Build clinical matrix (BASELINE row A)
    clin_cols = ['DAS28-CRP', 'VAS', 'TJC', 'SJC', 'CRP', 'Age',
                  'Anti-citrullinated peptide antibodies']
    # Categorical: Gender (M/F), Class (RA category code)
    rowsA = {}
    for pid in paired_pids:
        rA = df[(df.pid == pid) & (df.tp == 'A')].iloc[0]
        d = {}
        for c in clin_cols:
            try: d[c] = float(rA[c]) if pd.notna(rA[c]) else np.nan
            except Exception: d[c] = np.nan
        d['Gender_M'] = 1.0 if str(rA['Gender']).upper().startswith('M') else 0.0
        try: d['Class_num'] = float(rA['Class']) if pd.notna(rA['Class']) else np.nan
        except Exception: d['Class_num'] = np.nan
        rowsA[pid] = d
    clin_features = clin_cols + ['Gender_M', 'Class_num']
    X_clin = np.array([[rowsA[p][c] for c in clin_features] for p in paired_pids])
    # Drop all-NaN columns, then mean-impute remaining per column
    keep_cols = []
    for j, c in enumerate(clin_features):
        col = X_clin[:, j]
        if np.all(np.isnan(col)): print(f'  drop all-NaN: {c}'); continue
        keep_cols.append(j)
    X_clin = X_clin[:, keep_cols]
    clin_features = [clin_features[j] for j in keep_cols]
    for j in range(X_clin.shape[1]):
        col = X_clin[:, j]; mu = np.nanmean(col)
        if np.isnan(mu): mu = 0.0
        X_clin[np.isnan(col), j] = mu
    y = np.array([1 if by_pid.loc[p, 'Drugs Response'] == 'Response' else 0
                  for p in paired_pids])
    print(f'  Clinical features: {clin_features}')
    print(f'  Outcome: {Counter(y)}')

    # α scores are already aligned to the audit's patient order; sanity-check N
    assert len(alpha_scores) == len(paired_pids), \
        f'α-PC scores n={len(alpha_scores)} ≠ paired n={len(paired_pids)}'
    X_alpha = alpha_scores  # (N, 5)
    X_combined = np.hstack([X_clin, X_alpha])

    def loo_auc(X, y, label):
        loo = LeaveOneOut(); pred = np.zeros(len(y))
        for tr, te in loo.split(X):
            sc = StandardScaler().fit(X[tr])
            Xtr = sc.transform(X[tr]); Xte = sc.transform(X[te])
            model = LogisticRegression(max_iter=2000, C=1.0, random_state=0).fit(Xtr, y[tr])
            pred[te] = model.predict_proba(Xte)[:, 1]
        a = safe_auc(y, pred)
        print(f'  LOO-CV AUC ({label}): {a:.3f}')
        return float(a)

    auc_clin = loo_auc(X_clin, y, 'clinical only')
    auc_alpha = loo_auc(X_alpha, y, 'α-PC_T0 only (5)')
    auc_both = loo_auc(X_combined, y, 'clinical + α-PC_T0')
    return {
        'clinical_features': clin_features,
        'auc_clinical_only': auc_clin,
        'auc_alpha_only': auc_alpha,
        'auc_clinical_plus_alpha': auc_both,
        'incremental_AUC_from_alpha': auc_both - auc_clin,
    }


def main():
    cohorts = {
        'Wang_RA': cohort_wang_ra_baseline_for_response,
        'Filbin_D0': cohort_filbin_d0_predicts_outcome,
        'TB_DX': cohort_tb_dx_predicts_cure,
    }
    components_by_cohort = {}
    cohort_data = {}
    for name, fn in cohorts.items():
        F, log_pr, labels = fn()
        comps, scores, beta, ev = get_alpha_components(F, log_pr)
        components_by_cohort[name] = comps
        cohort_data[name] = {'scores': scores, 'beta': beta, 'labels': labels,
                              'explained_variance': ev.tolist()}
        print(f'  [{name}] EV: {[f"{v:.2%}" for v in ev]}')

    out = {}

    # Item 1
    out['item1_cross_cohort_axis_cosine'] = cross_cohort_axis_cosine(components_by_cohort)

    # Item 3 (Wang RA cluster mechanism)
    # Use LOOCV scores from prior prestrat audit for clustering consistency
    cd = cohort_data['Wang_RA']
    prestrat_json = RESULTS / 'loocv_prestrat_audit.json'
    with open(prestrat_json) as fh: prestrat = json.load(fh)
    loocv_scores = np.array(prestrat['Wang_RA_T0_predicts_response']['_scores']['loocv_alpha_pc'])
    print(f'  (Item 3 uses LOOCV α-PC scores from prestrat audit, n={len(loocv_scores)})')
    out['item3_wang_ra_cluster_mechanism'] = cluster_mechanism(
        'Wang_RA', loocv_scores, components_by_cohort['Wang_RA'], cd['labels'],
        feat_node_ids=None, n_clusters=3, top_nodes=12)

    # Item 4
    out['item4_wang_ra_clinical_baseline'] = clinical_baseline_benchmark(
        cd['scores'], cd['labels'])

    out_json = RESULTS / 'paper4_followups.json'
    with open(out_json, 'w') as fh: json.dump(out, fh, indent=2)
    print(f'\nSaved: {out_json}')


if __name__ == '__main__':
    main()
