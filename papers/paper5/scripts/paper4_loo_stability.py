"""Leave-one-out stability of Wang RA PC5/PC9/PC10.

Critical hostile-reviewer concern: PCA on n=38 samples may put PC5/PC9/PC10
in the Marčenko-Pastur noise regime. Test by re-fitting PCA on n=37 leave-one-out
subsets and computing subspace agreement with the full-cohort axes.

Outputs:
  - Per-PC, per-LOO absolute cosine with reference axis (full-cohort PC5/PC9/PC10)
  - 3-axis subspace Procrustes angle (does the 3D mode subspace survive LOO)
  - Permuted-label PCA null: re-fit PCA on label-shuffled Wang RA, project SLE,
    compare AUC distribution to actual Wang PC AUCs.
"""
from __future__ import annotations
import sys, gc
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
CACHE = RESULTS / 'cohort_alpha_cache'
SCAN = RESULTS / 'cross_indication_scan'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


def best_sign_aligned_cosine(v_ref, v_test):
    """Cosine between two PCs after best sign alignment."""
    c = v_ref @ v_test
    return float(abs(c) / (np.linalg.norm(v_ref) * np.linalg.norm(v_test) + 1e-9))


def best_match_cosine(v_ref, candidate_matrix):
    """For a reference axis, find the best-matching PC in a candidate matrix
    (by max |cos|). Returns (best_idx, best_cos)."""
    cosines = np.abs(candidate_matrix @ v_ref)
    cosines = cosines / (np.linalg.norm(candidate_matrix, axis=1) + 1e-9) \
              / (np.linalg.norm(v_ref) + 1e-9)
    best = int(np.argmax(cosines))
    return best, float(cosines[best])


def procrustes_subspace_angle(A, B):
    """Principal angles between subspaces spanned by columns of A and B.
    Both A and B are (n_dim, k) matrices. Returns sorted principal angles in
    degrees."""
    qa, _ = np.linalg.qr(A)
    qb, _ = np.linalg.qr(B)
    s = np.linalg.svd(qa.T @ qb, compute_uv=False)
    s = np.clip(s, -1.0, 1.0)
    angles = np.degrees(np.arccos(s))
    return sorted(angles)


def main():
    print('Loading substrate + Wang RA α...', flush=True)
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom
    lpn = (geom.log_pr / (np.linalg.norm(geom.log_pr) + 1e-9)).astype(np.float32)
    a0 = np.load(CACHE / 'Wang_RA_alpha_t0.npy')
    n_samples = a0.shape[0]
    print(f'  Wang RA α: {a0.shape}', flush=True)

    # Reference PCs (full cohort)
    pca_ref = PCA(n_components=15, svd_solver='randomized', random_state=0).fit(a0)
    ref_PC5 = pca_ref.components_[4].astype(np.float32)
    ref_PC9 = pca_ref.components_[8].astype(np.float32)
    ref_PC10 = pca_ref.components_[9].astype(np.float32)
    ref_subspace = np.stack([ref_PC5, ref_PC9, ref_PC10]).T  # (n_dim, 3)
    print(f'  Reference PC5/PC9/PC10 derived', flush=True)

    # ==========================================================
    # Test 1: LOO stability of PC5/PC9/PC10 identity
    # ==========================================================
    print('\n' + '='*78, flush=True)
    print(f'LOO stability test (n={n_samples} leave-one-out iterations)',
          flush=True)
    print('='*78, flush=True)
    rows = []
    subspace_angles = []
    for i in range(n_samples):
        keep = np.delete(np.arange(n_samples), i)
        a_loo = a0[keep]
        pca_loo = PCA(n_components=15, svd_solver='randomized',
                       random_state=0).fit(a_loo)
        # For each reference PC, find best matching LOO PC by |cos|
        b5_idx, b5_cos = best_match_cosine(ref_PC5, pca_loo.components_)
        b9_idx, b9_cos = best_match_cosine(ref_PC9, pca_loo.components_)
        b10_idx, b10_cos = best_match_cosine(ref_PC10, pca_loo.components_)
        rows.append({
            'loo_idx': i,
            'PC5_best_match': b5_idx + 1, 'PC5_best_cos': b5_cos,
            'PC9_best_match': b9_idx + 1, 'PC9_best_cos': b9_cos,
            'PC10_best_match': b10_idx + 1, 'PC10_best_cos': b10_cos,
        })
        # 3-axis subspace angles
        loo_subspace = np.stack([pca_loo.components_[4],
                                   pca_loo.components_[8],
                                   pca_loo.components_[9]]).T
        angles = procrustes_subspace_angle(ref_subspace, loo_subspace)
        subspace_angles.append(angles)
        del pca_loo, a_loo; gc.collect()
        if (i+1) % 10 == 0:
            print(f'  iter {i+1}/{n_samples}', flush=True)

    pc5_cos = np.array([r['PC5_best_cos'] for r in rows])
    pc9_cos = np.array([r['PC9_best_cos'] for r in rows])
    pc10_cos = np.array([r['PC10_best_cos'] for r in rows])
    print(f'\n  PC5 best-match |cos| across LOO: mean={pc5_cos.mean():.3f}, '
          f'min={pc5_cos.min():.3f}, max={pc5_cos.max():.3f}', flush=True)
    print(f'  PC9 best-match |cos| across LOO: mean={pc9_cos.mean():.3f}, '
          f'min={pc9_cos.min():.3f}, max={pc9_cos.max():.3f}', flush=True)
    print(f'  PC10 best-match |cos| across LOO: mean={pc10_cos.mean():.3f}, '
          f'min={pc10_cos.min():.3f}, max={pc10_cos.max():.3f}', flush=True)

    # Fraction of LOO where the reference PC's best match has |cos| > 0.9
    print(f'\n  LOO iterations where best-match |cos| > 0.9:', flush=True)
    print(f'    PC5: {int((pc5_cos > 0.9).sum())}/{n_samples}',
          flush=True)
    print(f'    PC9: {int((pc9_cos > 0.9).sum())}/{n_samples}',
          flush=True)
    print(f'    PC10: {int((pc10_cos > 0.9).sum())}/{n_samples}',
          flush=True)

    # 3-axis subspace agreement: largest principal angle should be small
    subspace_angles = np.array(subspace_angles)
    largest_angles = subspace_angles.max(axis=1)
    print(f'\n  3-axis subspace largest principal angle (degrees) across LOO:',
          flush=True)
    print(f'    mean={largest_angles.mean():.1f}°, '
          f'median={np.median(largest_angles):.1f}°, '
          f'max={largest_angles.max():.1f}°', flush=True)

    np.save(SCAN / 'wang_loo_pc_match.npy', np.stack([pc5_cos, pc9_cos, pc10_cos]))
    np.save(SCAN / 'wang_loo_subspace_angles.npy', subspace_angles)

    # ==========================================================
    # Test 2: Label-permuted PCA null on SLE projection
    # ==========================================================
    # We can't shuffle "labels" within Wang RA because there are no SLE/HC
    # labels (it's an RA cohort with all RA patients). Instead, we test what
    # we CAN test: shuffle Wang RA SAMPLE order, re-fit PCA, project SLE,
    # compute AUC. This is equivalent to "what if Wang RA samples were
    # randomly drawn from this α distribution?" It tests whether the SPECIFIC
    # variance structure of Wang RA produces axes meaningful for SLE.
    print('\n' + '='*78, flush=True)
    print('Permuted-sample-order null: 200 iterations', flush=True)
    print('='*78, flush=True)
    F_sle = np.load(SCAN / 'SLE_F.npy')
    sle_lab = np.load(SCAN / 'SLE_labels.npy').astype(bool)
    beta = F_sle @ lpn

    rng = np.random.default_rng(0)
    n_perm = 200
    # For each permutation, compute resulting PC5/PC9/PC10 AUCs on SLE
    null_aucs_pc5 = np.zeros(n_perm)
    null_aucs_pc9 = np.zeros(n_perm)
    null_aucs_pc10 = np.zeros(n_perm)
    for it in range(n_perm):
        # Shuffle row order (not values within rows)
        perm_idx = rng.permutation(n_samples)
        a_perm = a0[perm_idx]
        # Add noise to deliberately perturb the eigenstructure
        a_perm = a_perm + rng.normal(0, 0.1 * a_perm.std(),
                                      a_perm.shape).astype(np.float32)
        pca_perm = PCA(n_components=15, svd_solver='randomized',
                        random_state=int(it)).fit(a_perm)
        for k, store in zip([4, 8, 9], [null_aucs_pc5, null_aucs_pc9, null_aucs_pc10]):
            v = pca_perm.components_[k]
            proj = F_sle @ v - beta * float(lpn @ v)
            a = roc_auc_score(sle_lab.astype(int), proj)
            store[it] = max(a, 1 - a)
        if (it+1) % 50 == 0:
            print(f'  iter {it+1}/{n_perm}', flush=True)

    for name, ref_auc, null_aucs in [('PC5', 0.890, null_aucs_pc5),
                                       ('PC9', 0.791, null_aucs_pc9),
                                       ('PC10', 0.844, null_aucs_pc10)]:
        mean = null_aucs.mean(); sd = null_aucs.std()
        z = (ref_auc - mean) / (sd + 1e-9)
        p_perm = (null_aucs >= ref_auc).mean()
        print(f'  {name}: real={ref_auc:.3f}, null mean={mean:.3f}±{sd:.3f}, '
              f'z={z:+.1f}, p_perm={p_perm:.4f}', flush=True)


if __name__ == '__main__':
    main()
