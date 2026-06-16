"""Δα vector direction test — isolates treatment-effect biology.

Going beyond endpoint comparison to compare the SHIFT VECTORS themselves.

For each paired cohort:
  1. Per-patient Δα = (α_T1) − (α_T0) in 4D α-PC2..5 space
  2. Per-patient Δα in SUBSTRATE space (via cohort's PC components)
  3. Mean Δα vector per outcome group (responder vs non-responder)
  4. Angle between Δα_pos and Δα_neg (within-cohort outcome divergence)
  5. Cross-cohort angle between Δα vectors of same outcome class
     (does treatment biology transfer?)

The substrate-space cosines are the proper "Δα direction" measurement
since both vectors live in R^86826 regardless of which cohort they came from.

Cohorts: Wang RA, Filbin D0D3, TB DX→Wk24.
"""
from __future__ import annotations
import sys, json
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


def safe_cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    if na == 0 or nb == 0: return float('nan')
    return float(np.dot(a, b) / (na * nb))


def angle_deg(a, b):
    c = safe_cos(a, b)
    if np.isnan(c): return float('nan')
    c = max(-1.0, min(1.0, c))
    return float(np.degrees(np.arccos(c)))


def analyze(name, F_t0, F_t1, log_pr, labels, label_names):
    n = F_t0.shape[0]
    lpn = log_pr / (np.linalg.norm(log_pr) + 1e-9)
    F_all = np.vstack([F_t0, F_t1])
    beta_all = F_all @ lpn
    alpha_all = F_all - np.outer(beta_all, lpn)
    pca = PCA(n_components=5, svd_solver='randomized', random_state=0).fit(alpha_all)
    components = pca.components_   # (5, n_nodes)
    scores_all = pca.transform(alpha_all)

    alpha_t0 = scores_all[:n, 1:]
    alpha_t1 = scores_all[n:, 1:]
    delta_alpha_pc = alpha_t1 - alpha_t0   # (n, 4) in α-PC2..5 space

    # Project Δα-PC back to substrate-node space
    # delta_alpha_substrate = delta_alpha_pc @ components[1:]   # (n, n_nodes)
    delta_alpha_substrate = delta_alpha_pc @ components[1:]
    # Δβ per patient
    dbeta = beta_all[n:] - beta_all[:n]

    print(f'\n{"="*70}\n[{name}] n_pairs={n}\n{"="*70}')
    print(f'  Labels: 1={label_names[1]} (n={int(labels.sum())}), '
          f'0={label_names[0]} (n={int((1-labels).sum())})')

    # Mean Δα vector per outcome (α-PC space)
    pos_mask = labels == 1; neg_mask = labels == 0
    mean_d_pos_pc = delta_alpha_pc[pos_mask].mean(axis=0)
    mean_d_neg_pc = delta_alpha_pc[neg_mask].mean(axis=0)
    print(f'\n  Mean Δα-PC2..5 vectors:')
    print(f'    {label_names[1]}: {mean_d_pos_pc}, ‖‖={np.linalg.norm(mean_d_pos_pc):.3f}')
    print(f'    {label_names[0]}: {mean_d_neg_pc}, ‖‖={np.linalg.norm(mean_d_neg_pc):.3f}')

    cos_pc = safe_cos(mean_d_pos_pc, mean_d_neg_pc)
    print(f'  Angle in α-PC2..5 space: cos={cos_pc:.3f}, angle={angle_deg(mean_d_pos_pc, mean_d_neg_pc):.1f}°')

    # Mean Δα vector per outcome (SUBSTRATE space) — the fair comparison
    mean_d_pos_sub = delta_alpha_substrate[pos_mask].mean(axis=0)
    mean_d_neg_sub = delta_alpha_substrate[neg_mask].mean(axis=0)
    cos_sub = safe_cos(mean_d_pos_sub, mean_d_neg_sub)
    print(f'  Angle in SUBSTRATE space: cos={cos_sub:.3f}, angle={angle_deg(mean_d_pos_sub, mean_d_neg_sub):.1f}°')

    # Mean Δβ per outcome
    print(f'\n  Mean Δβ:')
    print(f'    {label_names[1]}: {dbeta[pos_mask].mean():+.3f}')
    print(f'    {label_names[0]}: {dbeta[neg_mask].mean():+.3f}')

    return {
        'cohort': name,
        'n_pairs': n,
        'label_names': label_names,
        'mean_dalpha_pc_pos': mean_d_pos_pc.tolist(),
        'mean_dalpha_pc_neg': mean_d_neg_pc.tolist(),
        'mean_dalpha_substrate_pos': mean_d_pos_sub.tolist(),
        'mean_dalpha_substrate_neg': mean_d_neg_sub.tolist(),
        'mean_dbeta_pos': float(dbeta[pos_mask].mean()),
        'mean_dbeta_neg': float(dbeta[neg_mask].mean()),
        'cosine_dalpha_pc_pos_vs_neg': float(cos_pc),
        'angle_dalpha_pc_pos_vs_neg_deg': float(angle_deg(mean_d_pos_pc, mean_d_neg_pc)),
        'cosine_dalpha_substrate_pos_vs_neg': float(cos_sub),
        'angle_dalpha_substrate_pos_vs_neg_deg': float(angle_deg(mean_d_pos_sub, mean_d_neg_sub)),
        # Save substrate Δα for cross-cohort comparison
        '_substrate_d_pos': mean_d_pos_sub.tolist(),
        '_substrate_d_neg': mean_d_neg_sub.tolist(),
    }


def main():
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom
    from gizmo.inference.projection import solve_map, ModalitySetup

    def get_F_pair(loader):
        captured = {}
        original = audit_mod._solve_paired
        def capture(X_t0, X_t1, feat_node_ids, patient_ids, labels):
            captured.update({'X_t0': X_t0, 'X_t1': X_t1,
                             'feat_node_ids': feat_node_ids,
                             'patient_ids': patient_ids, 'labels': labels})
            return original(X_t0, X_t1, feat_node_ids, patient_ids, labels)
        audit_mod._solve_paired = capture
        try: loader()
        finally: audit_mod._solve_paired = original
        X_t0 = captured['X_t0']; X_t1 = captured['X_t1']
        feat_nids = captured['feat_node_ids']
        pids = captured['patient_ids']; labels = captured['labels']
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
        F, _ = solve_map(geom, [setup], t0_sids + t1_sids)
        return F[:N], F[N:], labels

    cohort_specs = [
        ('Wang_RA_MTX',        audit_mod.cohort_wang_ra,    ['No-Response', 'Response']),
        ('Filbin_D0D3_COVID',  audit_mod.cohort_filbin_d0d3, ['Worsened', 'Improved']),
        ('TB_DX_Wk24_cure',    audit_mod.cohort_tb_dx_wk24,  ['Not-Cured', 'Definite-Cure']),
    ]
    results = []
    for name, loader, lbl in cohort_specs:
        print(f'\n>>> {name}: loading + MAP solve...')
        F_t0, F_t1, labels = get_F_pair(loader)
        results.append(analyze(name, F_t0, F_t1, geom.log_pr, labels, lbl))

    # Cross-cohort comparison
    print(f'\n{"="*78}\nCross-cohort Δα-substrate direction transfer\n{"="*78}')
    print(f'\n  Responder-direction transfer:')
    for i in range(len(results)):
        for j in range(i+1, len(results)):
            a = np.array(results[i]['_substrate_d_pos'])
            b = np.array(results[j]['_substrate_d_pos'])
            cos = safe_cos(a, b); ang = angle_deg(a, b)
            print(f'    {results[i]["cohort"]:<22} ↔ {results[j]["cohort"]:<22}  '
                  f'cos={cos:+.3f}  angle={ang:.1f}°')

    print(f'\n  Non-responder-direction transfer:')
    for i in range(len(results)):
        for j in range(i+1, len(results)):
            a = np.array(results[i]['_substrate_d_neg'])
            b = np.array(results[j]['_substrate_d_neg'])
            cos = safe_cos(a, b); ang = angle_deg(a, b)
            print(f'    {results[i]["cohort"]:<22} ↔ {results[j]["cohort"]:<22}  '
                  f'cos={cos:+.3f}  angle={ang:.1f}°')

    print(f'\n  Cross-class within-cohort (pos vs neg Δα angle, summary):')
    for r in results:
        print(f'    {r["cohort"]:<22}  '
              f'substrate angle = {r["angle_dalpha_substrate_pos_vs_neg_deg"]:5.1f}°  '
              f'(cos = {r["cosine_dalpha_substrate_pos_vs_neg"]:+.3f})')

    # Strip large fields before save
    summary = []
    for r in results:
        s = {k: v for k, v in r.items() if not k.startswith('_substrate')}
        summary.append(s)
    out_path = RESULTS / 'delta_alpha_direction_test.json'
    with open(out_path, 'w') as fh: json.dump(summary, fh, indent=2)
    print(f'\nSaved: {out_path}')


if __name__ == '__main__':
    main()
