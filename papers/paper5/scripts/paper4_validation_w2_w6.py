"""W2 + W6 — Cure-trajectory generalization across paired cohorts + anti-IL6 ↔ COVID test.

W2: For each paired cohort (Wang RA, Filbin COVID, TB cure, HCV DAA, bariatric),
    project mean Δα onto Wang-PC5/9/10. Direction match against HC-recovery
    expectation (PC5↓, PC9↑, PC10↑) tests whether the coordinate system
    consistently captures recovery direction across disease classes.

W6: Anti-IL6 (GSE93272 tocilizumab) axis projected onto Filbin COVID Δα.
    Tocilizumab is FDA-approved for severe COVID-19; the framework should
    detect axis alignment if it correctly predicts this cross-indication.
"""
from __future__ import annotations
import sys, gc
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
CACHE = RESULTS / 'cohort_alpha_cache'
SCAN = RESULTS / 'cross_indication_scan'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))


def main():
    print('Loading substrate + Wang PCs...', flush=True)
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom
    lpn = (geom.log_pr / (np.linalg.norm(geom.log_pr) + 1e-9)).astype(np.float32)
    a0 = np.load(CACHE / 'Wang_RA_alpha_t0.npy')
    pca_w = PCA(n_components=15, svd_solver='randomized', random_state=0).fit(a0)
    AXES = {
        'PC5_IFN':     pca_w.components_[4].astype(np.float32),
        'PC9_folate':  pca_w.components_[8].astype(np.float32),
        'PC10_plasma': pca_w.components_[9].astype(np.float32),
    }
    del a0; gc.collect()

    # ==========================================================
    # W2 — Cure trajectory across all paired cohorts
    # ==========================================================
    print('\n' + '='*80, flush=True)
    print('W2 — Cure-trajectory direction across all paired cohorts',
          flush=True)
    print('='*80, flush=True)
    print(f'  Expected HC-recovery direction: PC5 ↓ (lower IFN), '
          f'PC9 ↑ (folate norm), PC10 ↑ (plasma norm)\n', flush=True)

    paired = [
        ('Wang_RA',      CACHE / 'Wang_RA_alpha_t0.npy',
                         CACHE / 'Wang_RA_alpha_t1.npy',
                         'RA dx→fu (mixed treatment)'),
        ('Filbin_COVID', CACHE / 'Filbin_alpha_t0.npy',
                         CACHE / 'Filbin_alpha_t1.npy',
                         'Acute COVID D0→D3'),
        ('TB_cure',      CACHE / 'TB_alpha_t0.npy',
                         CACHE / 'TB_alpha_t1.npy',
                         'TB pre→post antibiotic'),
        ('HCV_DAA',      CACHE / 'HCV_alpha_t0.npy',
                         CACHE / 'HCV_alpha_t1.npy',
                         'HCV pre→post DAA'),
    ]
    rows = []
    print(f'  {"cohort":<16}{"context":<32}{"ΔPC5 mean":>13}'
          f'{"ΔPC9 mean":>13}{"ΔPC10 mean":>13}{"direction":>14}', flush=True)
    for key, p0, p1, label in paired:
        if not (p0.exists() and p1.exists()):
            print(f'  {key}: α missing; skip', flush=True); continue
        a0_ = np.load(p0); a1_ = np.load(p1)
        n = min(a0_.shape[0], a1_.shape[0])
        delta = (a1_[:n] - a0_[:n]).astype(np.float32)
        d5 = float((delta @ AXES['PC5_IFN']).mean())
        d9 = float((delta @ AXES['PC9_folate']).mean())
        d10 = float((delta @ AXES['PC10_plasma']).mean())
        # Direction match: PC5↓ AND PC9↑ AND PC10↑
        n_correct = sum([d5 < 0, d9 > 0, d10 > 0])
        direction_str = f'{n_correct}/3 correct'
        if n_correct == 3:
            direction_str += ' ✓'
        elif n_correct == 0:
            direction_str += ' ✗ all wrong'
        rows.append({'cohort': key, 'label': label, 'n_paired': n,
                     'dPC5': d5, 'dPC9': d9, 'dPC10': d10,
                     'directions_correct': n_correct})
        print(f'  {key:<16}{label[:31]:<32}{d5:>+13.2f}{d9:>+13.2f}'
              f'{d10:>+13.2f}{direction_str:>14}', flush=True)
        del a0_, a1_, delta; gc.collect()

    # Bariatric (pre/post within RYGB cohort, not paired but mean diff)
    bari = np.load(SCAN / 'bariatric_GLP1_alpha.npy')
    bari_lab = np.load(SCAN / 'bariatric_GLP1_labels.npy')
    delta_b = bari[bari_lab == 1].mean(axis=0) - bari[bari_lab == 0].mean(axis=0)
    db5 = float(delta_b @ AXES['PC5_IFN'])
    db9 = float(delta_b @ AXES['PC9_folate'])
    db10 = float(delta_b @ AXES['PC10_plasma'])
    n_correct = sum([db5 < 0, db9 > 0, db10 > 0])
    direction_str = f'{n_correct}/3 correct' + (' ✓' if n_correct == 3 else '')
    print(f'  {"bariatric_RYGB":<16}{"Pre→post RYGB intestine":<32}'
          f'{db5:>+13.2f}{db9:>+13.2f}{db10:>+13.2f}{direction_str:>14}',
          flush=True)
    rows.append({'cohort': 'bariatric_RYGB', 'label': 'Pre→post RYGB',
                 'dPC5': db5, 'dPC9': db9, 'dPC10': db10,
                 'directions_correct': n_correct})

    correct_cohorts = [r for r in rows if r['directions_correct'] == 3]
    print(f'\n  Cohorts with FULL 3/3 direction match: '
          f'{len(correct_cohorts)}/{len(rows)}', flush=True)
    print(f'  Cohorts ≥2/3 match: '
          f'{len([r for r in rows if r["directions_correct"] >= 2])}/{len(rows)}',
          flush=True)
    pd.DataFrame(rows).to_csv(SCAN / 'cure_trajectory_paired_cohorts.csv',
                              index=False)

    # ==========================================================
    # W6 — Anti-IL6 axis projected onto Filbin COVID Δα
    # ==========================================================
    print('\n' + '='*80, flush=True)
    print('W6 — Anti-IL6 (tocilizumab) axis ↔ Filbin COVID recovery direction',
          flush=True)
    print('='*80, flush=True)
    # Anti-IL6 direction vector (unit norm) from GSE93272
    v_il6 = np.load(SCAN / 'axis_GSE93272_antiIL6_tcz.npy').astype(np.float32)
    print(f'  Anti-IL6 axis: ||v||={np.linalg.norm(v_il6):.3f}', flush=True)

    # Filbin Δα
    filbin_t0 = np.load(CACHE / 'Filbin_alpha_t0.npy')
    filbin_t1 = np.load(CACHE / 'Filbin_alpha_t1.npy')
    n = min(filbin_t0.shape[0], filbin_t1.shape[0])
    delta_f = (filbin_t1[:n] - filbin_t0[:n]).astype(np.float32)
    # Project each patient's Δα onto anti-IL6 axis
    proj = delta_f @ v_il6
    print(f'  Filbin Δα onto anti-IL6 axis:', flush=True)
    print(f'    mean={proj.mean():+.3f}, SD={proj.std():.3f}, '
          f'n={len(proj)}', flush=True)
    pos_pct = float((proj > 0).mean() * 100)
    print(f'    % patients positive (anti-IL6 axis engaged by recovery): '
          f'{pos_pct:.1f}%', flush=True)

    # Compare to anti-TNF and MTX axes for context
    v_tnf = np.load(SCAN / 'axis_GSE93272_antiTNF_ifx.npy').astype(np.float32)
    v_mtx = np.load(SCAN / 'axis_GSE93272_MTX.npy').astype(np.float32)
    proj_tnf = delta_f @ v_tnf
    proj_mtx = delta_f @ v_mtx
    print(f'\n  Comparison — Filbin COVID Δα onto each RA-derived drug axis:',
          flush=True)
    print(f'  {"axis":<22}{"mean proj":>14}{"% positive":>14}{"interpretation":>30}',
          flush=True)
    for axis_name, projection in [('anti-IL6 (tocilizumab)', proj),
                                    ('anti-TNF (infliximab)', proj_tnf),
                                    ('MTX (methotrexate)', proj_mtx)]:
        m = projection.mean()
        pp = float((projection > 0).mean() * 100)
        interp = ('axis engaged in recovery direction' if abs(m) > 0.5 and pp > 60
                  else ('axis engaged but mixed' if abs(m) > 0.5
                         else 'axis not significantly engaged'))
        print(f'  {axis_name:<22}{m:>+14.3f}{pp:>13.1f}%{interp:>30}',
              flush=True)

    # Clinical context
    print(f'\n  Clinical context:', flush=True)
    print(f'    Anti-IL6 (tocilizumab) IS FDA-approved for severe COVID-19 '
          f'(RECOVERY trial, REMAP-CAP)', flush=True)
    print(f'    Anti-TNF is NOT used in COVID', flush=True)
    print(f'    MTX is NOT used in COVID', flush=True)
    print(f'    Framework prediction: should find anti-IL6 axis most engaged',
          flush=True)


if __name__ == '__main__':
    main()
