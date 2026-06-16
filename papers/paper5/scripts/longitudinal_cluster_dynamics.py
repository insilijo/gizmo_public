"""Longitudinal cluster dynamics: do α-PC clusters change over time?

For paired cohorts (Wang RA A→B, Filbin D0→D3):
  1. Solve MAP on pooled (T0+T1) samples
  2. Compute α-PC basis on pooled α
  3. Cluster T0 patients in variance-normalized α-PC2..5 space (k=3)
  4. Project T1 patient scores through the SAME basis + cluster centroids
  5. Assign T1 cluster = nearest centroid in variance-normalized space
  6. Transition matrix P(cluster_T1 | cluster_T0)
  7. Test: are transitions different for responders vs non-responders?

Key question: do non-responder cluster patients STAY in their cluster
(phenotype lock-in), while responder cluster patients SHIFT? That would be
direct evidence the framework captures treatment-induced phenotype change
at the subtype level.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from collections import Counter

import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.stats import fisher_exact

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


def varnorm(A, mu=None, sd=None):
    if mu is None: mu = A.mean(axis=0)
    if sd is None: sd = A.std(axis=0) + 1e-9
    return (A - mu) / sd, mu, sd


def assign_cluster_by_centroid(X_z, centroids):
    """Assign each row of X_z to nearest centroid in Euclidean distance."""
    dists = np.linalg.norm(X_z[:, None, :] - centroids[None, :, :], axis=2)
    return np.argmin(dists, axis=1)


def analyze_paired(name, F_t0, F_t1, log_pr, labels, label_names, n_components=5):
    """Cluster on T0, project T1, transition analysis."""
    print(f'\n{"="*70}\n[{name}] n_pairs={F_t0.shape[0]}\n{"="*70}')
    n = F_t0.shape[0]
    lpn = log_pr / (np.linalg.norm(log_pr) + 1e-9)

    # Pooled β/α basis
    F_all = np.vstack([F_t0, F_t1])
    beta_all = F_all @ lpn
    alpha_all = F_all - np.outer(beta_all, lpn)
    pca = PCA(n_components=n_components, svd_solver='randomized', random_state=0).fit(alpha_all)
    scores_all = pca.transform(alpha_all)  # (2N, 5)

    alpha_pc_t0 = scores_all[:n, 1:]   # PC2..5 only
    alpha_pc_t1 = scores_all[n:, 1:]

    # Variance-normalize using T0 statistics
    z_t0, mu, sd = varnorm(alpha_pc_t0)
    z_t1, _, _   = varnorm(alpha_pc_t1, mu=mu, sd=sd)

    # Cluster T0
    km = KMeans(n_clusters=3, random_state=0, n_init=20).fit(z_t0)
    cluster_t0 = km.labels_
    centroids = km.cluster_centers_
    print(f'\n  T0 cluster sizes: {Counter(cluster_t0)}')

    # Assign T1 via nearest centroid in same variance-normalized space
    cluster_t1 = assign_cluster_by_centroid(z_t1, centroids)
    print(f'  T1 cluster sizes: {Counter(cluster_t1)}')

    # Transition matrix
    print(f'\n  Transition matrix P(cluster_T1 | cluster_T0):')
    trans_matrix = np.zeros((3, 3), dtype=int)
    for i in range(n):
        trans_matrix[cluster_t0[i], cluster_t1[i]] += 1
    print(f'    T0→ /  T1→0   T1→1   T1→2   total  stay-prob')
    for c0 in range(3):
        row_total = trans_matrix[c0].sum()
        stay = trans_matrix[c0, c0] / row_total if row_total > 0 else 0
        print(f'    C{c0}        {trans_matrix[c0, 0]:>4}   {trans_matrix[c0, 1]:>4}   '
              f'{trans_matrix[c0, 2]:>4}   {row_total:>5}   {stay:.2f}')

    # Per-label transition: do labels=1 transition differently than labels=0?
    print(f'\n  By {label_names} outcome:')
    for lbl_val, lbl_name in [(1, label_names[1]), (0, label_names[0])]:
        mask = labels == lbl_val
        n_l = int(mask.sum())
        if n_l == 0: continue
        same_cluster = sum(1 for i in range(n) if labels[i] == lbl_val and cluster_t0[i] == cluster_t1[i])
        print(f'    {lbl_name} (n={n_l}): {same_cluster}/{n_l} '
              f'({100*same_cluster/n_l:.0f}%) stay in same cluster')

    # Stability test: are stay-rates different by label? Fisher 2×2 (stayed vs moved) × (label_pos vs neg)
    stayed_pos = sum(1 for i in range(n) if labels[i] == 1 and cluster_t0[i] == cluster_t1[i])
    moved_pos = sum(1 for i in range(n) if labels[i] == 1 and cluster_t0[i] != cluster_t1[i])
    stayed_neg = sum(1 for i in range(n) if labels[i] == 0 and cluster_t0[i] == cluster_t1[i])
    moved_neg = sum(1 for i in range(n) if labels[i] == 0 and cluster_t0[i] != cluster_t1[i])
    print(f'\n  2x2 stay-vs-move × label:')
    print(f'              stay   move')
    print(f'    pos       {stayed_pos:>4}   {moved_pos:>4}')
    print(f'    neg       {stayed_neg:>4}   {moved_neg:>4}')
    try:
        _, p = fisher_exact([[stayed_pos, moved_pos], [stayed_neg, moved_neg]])
        print(f'    Fisher p (label vs stay/move): {p:.4f}')
    except Exception as e:
        p = float('nan')
        print(f'    Fisher p: nan ({e})')

    # Dominant PC per cluster for naming
    print(f'\n  Dominant α-PC per cluster centroid:')
    for c in range(3):
        dom = int(np.argmax(np.abs(centroids[c]))) + 2  # +2 because we dropped PC1
        sign = '+' if centroids[c, dom-2] > 0 else '−'
        print(f'    C{c}: {sign}PC{dom}')

    return {
        'cluster_t0': cluster_t0.tolist(),
        'cluster_t1': cluster_t1.tolist(),
        'transition_matrix': trans_matrix.tolist(),
        'fisher_p_label_vs_movement': float(p),
        'centroids': centroids.tolist(),
        'label_names': label_names,
    }


def main():
    out = {}

    # === Wang RA paired A→B ===
    # Import audit module FIRST (loads substrate at module level) — don't reload here
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom; mg = audit_mod.mg
    print('\nLoading Wang RA paired...')
    captured = {}
    original = audit_mod._solve_paired
    def capture(X_t0, X_t1, feat_node_ids, patient_ids, labels):
        captured['X_t0'] = X_t0; captured['X_t1'] = X_t1
        captured['feat_node_ids'] = feat_node_ids
        captured['patient_ids'] = patient_ids; captured['labels'] = labels
        return original(X_t0, X_t1, feat_node_ids, patient_ids, labels)
    audit_mod._solve_paired = capture
    try:
        audit_mod.cohort_wang_ra()
    finally:
        audit_mod._solve_paired = original
    X_t0 = captured['X_t0']; X_t1 = captured['X_t1']
    feat_nids = captured['feat_node_ids']
    pids = captured['patient_ids']
    labels = captured['labels']
    # Solve MAP on pooled to get F_t0, F_t1 (reusing audit_mod's geom)
    from gizmo.inference.projection import solve_map, ModalitySetup
    N = X_t0.shape[0]
    X_pool = np.vstack([X_t0, X_t1])
    mu = X_pool.mean(axis=0); sd = X_pool.std(axis=0) + 1e-9
    X_t0z = (X_t0 - mu) / sd; X_t1z = (X_t1 - mu) / sd
    feat_cols = [(f'feat_{k}', geom.nid_idx[feat_nids[k]]) for k in range(len(feat_nids))]
    t0_sids = [f'{p}_T0' for p in pids]; t1_sids = [f'{p}_T1' for p in pids]
    pdata = {sid: {f'feat_{k}': float(X_t0z[i, k]) for k in range(X_t0z.shape[1])}
             for i, sid in enumerate(t0_sids)}
    pdata.update({sid: {f'feat_{k}': float(X_t1z[i, k]) for k in range(X_t1z.shape[1])}
                  for i, sid in enumerate(t1_sids)})
    setup = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=pdata)
    print('Solving MAP on Wang RA pooled (2N)...')
    F, _ = solve_map(geom, [setup], t0_sids + t1_sids)
    F_t0 = F[:N]; F_t1 = F[N:]
    out['Wang_RA'] = analyze_paired('Wang_RA (A→B, MTX)', F_t0, F_t1,
                                     geom.log_pr, labels,
                                     label_names=['No-Response', 'Response'])

    # === Filbin D0→D3 ===
    print('\nLoading Filbin D0→D3...')
    captured.clear()
    audit_mod._solve_paired = capture
    try:
        audit_mod.cohort_filbin_d0d3()
    finally:
        audit_mod._solve_paired = original
    X_t0 = captured['X_t0']; X_t1 = captured['X_t1']
    feat_nids = captured['feat_node_ids']
    pids = captured['patient_ids']
    labels = captured['labels']
    N = X_t0.shape[0]
    X_pool = np.vstack([X_t0, X_t1])
    mu = X_pool.mean(axis=0); sd = X_pool.std(axis=0) + 1e-9
    X_t0z = (X_t0 - mu) / sd; X_t1z = (X_t1 - mu) / sd
    feat_cols = [(f'feat_{k}', geom.nid_idx[feat_nids[k]]) for k in range(len(feat_nids))]
    t0_sids = [f'{p}_T0' for p in pids]; t1_sids = [f'{p}_T1' for p in pids]
    pdata = {sid: {f'feat_{k}': float(X_t0z[i, k]) for k in range(X_t0z.shape[1])}
             for i, sid in enumerate(t0_sids)}
    pdata.update({sid: {f'feat_{k}': float(X_t1z[i, k]) for k in range(X_t1z.shape[1])}
                  for i, sid in enumerate(t1_sids)})
    setup = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                          feature_cols=feat_cols, data=pdata)
    print('Solving MAP on Filbin pooled (2N)...')
    F, _ = solve_map(geom, [setup], t0_sids + t1_sids)
    F_t0 = F[:N]; F_t1 = F[N:]
    out['Filbin_D0D3'] = analyze_paired('Filbin (D0→D3, Acuity)', F_t0, F_t1,
                                          geom.log_pr, labels,
                                          label_names=['Worsened', 'Improved'])

    out_path = RESULTS / 'longitudinal_cluster_dynamics.json'
    with open(out_path, 'w') as fh: json.dump(out, fh, indent=2)
    print(f'\nSaved: {out_path}')


if __name__ == '__main__':
    main()
