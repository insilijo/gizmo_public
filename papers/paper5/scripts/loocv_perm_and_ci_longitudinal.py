"""Permutation null + DeLong CI for longitudinal trajectory α-PC AUCs.

Reads `loocv_longitudinal_audit.json` (produced by loocv_longitudinal_audit.py)
which contains _scores per cohort: Δβ per patient, LOOCV Δα-PC per patient,
trajectory label per patient.

Computes:
  - DeLong 95% CI on each AUC (Δβ + 5 Δα-PCs)
  - 1000-permutation label-shuffle null:
      • Δβ: AUC vs shuffled labels
      • Δα: max-over-5-PCs AUC under each shuffle (multiplicity-corrected)
  - Two-sided p-value
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
IN_JSON = REPO / 'benchmarks/results/unsupervised/loocv_longitudinal_audit.json'
OUT_JSON = REPO / 'benchmarks/results/unsupervised/loocv_longitudinal_audit_with_ci.json'


def safe_auc(y, s):
    a = roc_auc_score(y, s); return max(a, 1 - a)


def delong_ci(y_true, y_score, alpha=0.05):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    raw_auc = roc_auc_score(y_true, y_score)
    if raw_auc < 0.5:
        y_score = -y_score; raw_auc = 1 - raw_auc
    pos = y_score[y_true == 1]; neg = y_score[y_true == 0]
    m, n = len(pos), len(neg)
    if m < 2 or n < 2: return float(raw_auc), float('nan'), float('nan')

    def midrank(x):
        order = np.argsort(x); ranks = np.empty(len(x)); i = 0
        while i < len(x):
            j = i
            while j < len(x) - 1 and x[order[j]] == x[order[j+1]]: j += 1
            ranks[order[i:j+1]] = 0.5 * (i + j) + 1; i = j + 1
        return ranks

    Tx = midrank(pos); Ty = midrank(neg); Txy = midrank(np.concatenate([pos, neg]))
    V10 = (Txy[:m] - Tx) / n; V01 = 1 - (Txy[m:] - Ty) / m
    sx = np.var(V10, ddof=1) / m; sy = np.var(V01, ddof=1) / n
    se = np.sqrt(sx + sy)
    if se == 0 or not np.isfinite(se): return float(raw_auc), float('nan'), float('nan')
    z = 1.959963984540054
    return float(raw_auc), float(max(0.0, raw_auc - z * se)), float(min(1.0, raw_auc + z * se))


def permutation_pvals(dbeta, loocv_dalpha, labels, n_perm=1000, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    labels = np.asarray(labels); n = len(labels)
    obs_b = safe_auc(labels, dbeta)
    obs_a = [safe_auc(labels, loocv_dalpha[:, k]) for k in range(loocv_dalpha.shape[1])]
    obs_a_max = max(obs_a)

    n_b = 0; n_a_max = 0; n_a_pc = [0] * loocv_dalpha.shape[1]
    valid = 0
    for _ in range(n_perm):
        perm = rng.permutation(n); yp = labels[perm]
        if len(set(yp)) < 2: continue
        valid += 1
        if safe_auc(yp, dbeta) >= obs_b: n_b += 1
        per_pc = [safe_auc(yp, loocv_dalpha[:, k]) for k in range(loocv_dalpha.shape[1])]
        if max(per_pc) >= obs_a_max: n_a_max += 1
        for k, v in enumerate(per_pc):
            if v >= obs_a[k]: n_a_pc[k] += 1
    return {
        'dbeta_pval': (n_b + 1) / (valid + 1),
        'dalpha_max_pval': (n_a_max + 1) / (valid + 1),
        'dalpha_per_pc_pvals': [(c + 1) / (valid + 1) for c in n_a_pc],
        'n_perm': valid,
    }


def main():
    with open(IN_JSON) as fh: results = json.load(fh)
    enriched = {}
    for name, rec in results.items():
        scores = rec.get('_scores')
        if scores is None:
            print(f'\n[{name}] no _scores — skipping')
            enriched[name] = rec; continue
        labels = np.array(scores['labels'])
        dbeta = np.array(scores['dbeta'])
        loocv_dalpha = np.array(scores['loocv_dalpha_pc'])
        print(f'\n[{name}] n_pairs={len(labels)}, pos={int(labels.sum())}, neg={int((1-labels).sum())}')

        b_auc, b_lo, b_hi = delong_ci(labels, dbeta)
        a_cis = [delong_ci(labels, loocv_dalpha[:, k]) for k in range(loocv_dalpha.shape[1])]
        print(f'  Δβ AUC = {b_auc:.3f} [{b_lo:.3f}, {b_hi:.3f}]')
        for k, (a, lo, hi) in enumerate(a_cis):
            print(f'  Δα-PC{k+1} AUC = {a:.3f} [{lo:.3f}, {hi:.3f}]')

        perm = permutation_pvals(dbeta, loocv_dalpha, labels, n_perm=1000)
        print(f'  Δβ p = {perm["dbeta_pval"]:.4f}')
        print(f'  Δα-max-PC p (multipl-corr) = {perm["dalpha_max_pval"]:.4f}')
        print(f'  Δα per-PC p = {[f"{p:.4f}" for p in perm["dalpha_per_pc_pvals"]]}')

        rec_out = dict(rec)
        rec_out['delong_dbeta'] = {'auc': b_auc, 'lo': b_lo, 'hi': b_hi}
        rec_out['delong_dalpha'] = [{'auc': a, 'lo': lo, 'hi': hi} for a, lo, hi in a_cis]
        rec_out['perm'] = perm
        enriched[name] = rec_out

    with open(OUT_JSON, 'w') as fh: json.dump(enriched, fh, indent=2)
    print(f'\nSaved: {OUT_JSON}')


if __name__ == '__main__':
    main()
