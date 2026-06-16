"""Step 2: load TB only, load saved Filbin-PC5 eigvec, run B/D/E."""
import sys, json, gc
from pathlib import Path
import numpy as np
from scipy.stats import mannwhitneyu

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


def get_F_pair(loader, audit_mod, geom):
    from gizmo.inference.projection import solve_map, ModalitySetup
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
    return F[:N].astype(np.float32), F[N:].astype(np.float32), labels, pids


def main():
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom
    lpn = (geom.log_pr / (np.linalg.norm(geom.log_pr) + 1e-9)).astype(np.float32)

    filbin_pc5 = np.load(RESULTS / 'filbin_pc5_eigvec.npy').astype(np.float32)
    print(f'Loaded Filbin-PC5: shape={filbin_pc5.shape}', flush=True)

    F_t0, F_t1, labels, pids = get_F_pair(audit_mod.cohort_tb_dx_wk24, audit_mod, geom)
    beta_t0 = F_t0 @ lpn; beta_t1 = F_t1 @ lpn
    alpha_t0 = (F_t0 - np.outer(beta_t0, lpn)).astype(np.float32)
    alpha_t1 = (F_t1 - np.outer(beta_t1, lpn)).astype(np.float32)
    del F_t0, F_t1, beta_t0, beta_t1; gc.collect()
    print(f'TB α matrices: {alpha_t0.shape}', flush=True)

    proj_t0 = alpha_t0 @ filbin_pc5
    proj_t1 = alpha_t1 @ filbin_pc5
    dpc = (proj_t1 - proj_t0).astype(np.float64)
    not_cured = np.where(labels == 0)[0]

    # (D)
    print('\n=== (D) Baseline imbalance ===', flush=True)
    mwu_base = mannwhitneyu(proj_t0[labels==1], proj_t0[labels==0])
    mwu_traj = mannwhitneyu(dpc[labels==1], dpc[labels==0])
    print(f'  PC5_T0 baseline MWU: p={mwu_base.pvalue:.4f}', flush=True)
    print(f'  ΔPC5 trajectory MWU: p={mwu_traj.pvalue:.4f}', flush=True)
    print(f'  PC5_T0: Cure={proj_t0[labels==1].mean():+.2f}, NotCured={proj_t0[labels==0].mean():+.2f}',
          flush=True)
    print(f'  ΔPC5:   Cure={dpc[labels==1].mean():+.2f}, NotCured={dpc[labels==0].mean():+.2f}',
          flush=True)
    if mwu_base.pvalue < 0.10:
        print(f'  ⚠ Baseline imbalance — partially regression-to-mean', flush=True)
    else:
        print(f'  ✓ Baseline balanced — trajectory effect is genuine', flush=True)

    # (B)
    print('\n=== (B) Drop-one TB Not-Cured jackknife ===', flush=True)
    print(f'  Not-Cured IDs: {[pids[i] for i in not_cured]}', flush=True)
    jack_ps = []
    for drop_idx in not_cured:
        keep = np.ones(len(labels), dtype=bool); keep[drop_idx] = False
        if (labels[keep]==0).sum() < 2: continue
        m = mannwhitneyu(dpc[keep][labels[keep]==1], dpc[keep][labels[keep]==0])
        jack_ps.append(float(m.pvalue))
        mark = '⚠' if m.pvalue > 0.05 else ' '
        print(f'  drop {pids[drop_idx]:<6}: p={m.pvalue:.4f} {mark}', flush=True)
    print(f'  jackknife range: [{min(jack_ps):.4f}, {max(jack_ps):.4f}]', flush=True)
    if max(jack_ps) < 0.05:
        print(f'  ✓ All drop-one tests survive at p<0.05', flush=True)
    elif max(jack_ps) < 0.10:
        print(f'  ~ Survives at p<0.10', flush=True)
    else:
        print(f'  ⚠ Some leave-outs push p above 0.05', flush=True)

    # (E) Random-axis null, incremental
    print('\n=== (E) Random-axis null (K=300) ===', flush=True)
    rng = np.random.default_rng(0)
    n_random = 300
    obs_p = float(mwu_traj.pvalue)
    print(f'  Observed p_grp = {obs_p:.4f}', flush=True)
    rand_p_grp = np.zeros(n_random)
    for s in range(n_random):
        v = rng.standard_normal(filbin_pc5.shape[0]).astype(np.float32)
        v /= (np.linalg.norm(v) + 1e-9)
        dpc_v = alpha_t1 @ v - alpha_t0 @ v
        rand_p_grp[s] = float(mannwhitneyu(dpc_v[labels==1], dpc_v[labels==0]).pvalue)
        if (s+1) % 50 == 0:
            print(f'    {s+1}/{n_random} done', flush=True)
    n_below = int((rand_p_grp <= obs_p).sum())
    n_below_05 = int((rand_p_grp < 0.05).sum())
    print(f'\n  Random p_grp ≤ {obs_p:.4f}: {n_below}/{n_random} = {n_below/n_random:.4f}', flush=True)
    print(f'  Random p_grp < 0.05: {n_below_05}/{n_random} = {n_below_05/n_random:.4f}', flush=True)
    for q in [0.05, 0.10, 0.50, 0.90]:
        print(f'    {int(q*100)}%-tile: {float(np.quantile(rand_p_grp, q)):.4f}', flush=True)

    out = {
        'baseline_p': float(mwu_base.pvalue),
        'trajectory_p': float(mwu_traj.pvalue),
        'jackknife_ps': jack_ps,
        'jackknife_max_p': float(max(jack_ps)),
        'jackknife_min_p': float(min(jack_ps)),
        'random_n': n_random,
        'random_rate_below_obs': float(n_below / n_random),
        'random_rate_below_05': float(n_below_05 / n_random),
        'observed_p': obs_p,
    }
    out_path = RESULTS / 'cross_cohort_axis_stress_test_B.json'
    with open(out_path, 'w') as fh: json.dump(out, fh, indent=2)
    print(f'\nSaved: {out_path}', flush=True)


if __name__ == '__main__':
    main()
