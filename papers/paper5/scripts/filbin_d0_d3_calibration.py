"""Filbin COVID Day-0 → Day-3 paired calibration for drug-sim operator.

Loads paired Day 0 and Day 3 Olink proteomics for the 218 patients with both
timepoints, stratifies by clinical acuity change (improved vs worsened vs
stable), and reports per-patient F shifts that ground-truth what real acute
COVID treatment biology looks like.

Outputs:
  - ‖ΔF_i‖ distribution per patient (Improved / Worsened / Stable)
  - mean β shift per group (should drop in Improved if β = disease intensity)
  - α-PC shifts per group + variability
  - Top shifted nodes (Day 0 → Day 3 in Improved patients)

Compare against drug_sim_su_covid.json's ΔF magnitudes + directions.
"""
from __future__ import annotations
import sys, json, math
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
CACHE = Path.home() / '.cache' / 'filbin_mgh_covid'

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.inference.projection import (
    build_biochem_subgraph, solve_map, decompose_beta_alpha, ModalitySetup,
)


def load_filbin_paired(day_a=0, day_b=3):
    """Load paired Day A / Day B Olink expression for COVID-positive subjects.

    Returns:
      paired_data: dict {day: {sid: {gene: npx}}}
      common_pids: list of subjects with both timepoints
      acuity_change: dict pid → (acuity_a, acuity_b)
    """
    olink_p = CACHE / 'Olink_Proteomics.xlsx'
    clin_p  = CACHE / 'Clinical_Metadata.xlsx'
    assay_p = CACHE / 'Suppl_T2_Olink_Assays_NPX.xlsx'

    # OID → gene symbol map
    t2a = pd.read_excel(assay_p, sheet_name='2A-Olink-Assay', header=1)
    oid_to_sym = dict(zip(t2a['OlinkID'].astype(str), t2a['Assay'].astype(str)))

    # Clinical metadata: COVID + Acuity at each day
    clin = pd.read_excel(clin_p, sheet_name='Subject-level metadata')
    clin['core_pid'] = clin['Public ID'].astype(str)
    pid_covid = dict(zip(clin['core_pid'], clin['COVID'].astype(int)))
    pid_acuity = {}
    for _, r in clin.iterrows():
        pid_acuity[str(r['core_pid'])] = {
            0: r.get(f'Acuity 0'),
            3: r.get(f'Acuity 3'),
            7: r.get(f'Acuity 7'),
        }

    # Olink NPX matrix
    ol = pd.read_excel(olink_p)
    ol['core_pid'] = ol['Public ID'].astype(str).str.replace(r'_D\d+$', '', regex=True)
    ol['core_pid'] = ol['core_pid'].str.replace(r'_E$', '', regex=True)
    ol['day'] = ol['Day'].astype(str)

    paired = {day_a: {}, day_b: {}}
    for _, row in ol.iterrows():
        pid = str(row['core_pid'])
        day = row['day']
        try: day_int = int(day)
        except ValueError: continue
        if day_int not in (day_a, day_b): continue
        if pid_covid.get(pid) != 1: continue  # COVID+ only

        d = {}
        for col, v in row.items():
            if not isinstance(col, str) or not col.startswith('OID'): continue
            sym = oid_to_sym.get(col)
            if not sym or sym == 'nan': continue
            try: vf = float(v)
            except (TypeError, ValueError): continue
            if pd.isna(vf): continue
            d[sym] = max(d.get(sym, -1e9), vf)
        if d:
            sid = f'{pid}_D{day_int}'
            paired[day_int][sid] = d

    # Find subjects with BOTH timepoints
    pids_a = {sid.split('_D')[0] for sid in paired[day_a]}
    pids_b = {sid.split('_D')[0] for sid in paired[day_b]}
    common_pids = sorted(pids_a & pids_b)

    # Acuity change per patient (lower acuity = milder)
    acuity_change = {}
    for pid in common_pids:
        a0 = pid_acuity.get(pid, {}).get(day_a)
        a1 = pid_acuity.get(pid, {}).get(day_b)
        if pd.notna(a0) and pd.notna(a1):
            acuity_change[pid] = (int(a0), int(a1))

    print(f'Day {day_a}/Day {day_b} paired COVID+ subjects: {len(common_pids)}')
    print(f'  With acuity scores at both days: {len(acuity_change)}')

    return paired, common_pids, acuity_change


def stratify_by_acuity(acuity_change):
    """Classify each patient as Improved / Worsened / Stable based on Acuity 0→3.

    Acuity scale (Filbin):
      1 = died
      2 = intubated/ICU
      3 = hospitalized on supp. O2
      4 = hospitalized on room air
      5 = discharged/home
    Lower number = MORE severe.
    """
    groups = {'Improved': [], 'Worsened': [], 'Stable': []}
    for pid, (a0, a1) in acuity_change.items():
        if a1 > a0:    groups['Improved'].append(pid)
        elif a1 < a0:  groups['Worsened'].append(pid)
        else:          groups['Stable'].append(pid)
    print(f'  Acuity strata: {dict((k, len(v)) for k, v in groups.items())}')
    return groups


def map_to_substrate(paired_data, sids, geom):
    all_genes = sorted({k for sid in sids for k in paired_data[sid].keys()})
    node_ids = {k: f'symbol:{k}' for k in all_genes if f'symbol:{k}' in geom.nid_idx}
    kept = [k for k in all_genes if k in node_ids]
    X = np.zeros((len(sids), len(kept)))
    for i, sid in enumerate(sids):
        for j, k in enumerate(kept):
            X[i, j] = paired_data[sid].get(k, 0.0)
    feat_node_ids = [node_ids[k] for k in kept]
    return X, feat_node_ids


def main():
    print('=' * 80)
    print('Filbin COVID Day 0 → Day 3 paired calibration')
    print('=' * 80)

    paired_d, common_pids, acuity = load_filbin_paired(day_a=0, day_b=3)
    groups = stratify_by_acuity(acuity)

    rhea_full = REPO / 'data/processed/human_full_rhea_full/graph.json'
    mg = read_json(str(rhea_full))
    geom = build_biochem_subgraph(mg, hub_cap=500)
    print(f'\nSubstrate: {len(geom.nodes)} nodes')

    out = {
        'config': {'cohort': 'Filbin_COVID_paired_D0_D3', 'hub_cap': 500},
        'n_paired': len(common_pids),
        'acuity_strata': {k: len(v) for k, v in groups.items()},
        'strata': {},
    }

    for stratum, pids in groups.items():
        if len(pids) < 10:
            print(f'\n[{stratum}] n={len(pids)} too small, skipping')
            continue
        print(f'\n--- {stratum} responders (n={len(pids)} paired) ---')

        d0_sids = [f'{p}_D0' for p in pids]
        d3_sids = [f'{p}_D3' for p in pids]
        X_d0, feat_node_ids_d0 = map_to_substrate(paired_d[0], d0_sids, geom)
        X_d3, feat_node_ids_d3 = map_to_substrate(paired_d[3], d3_sids, geom)
        print(f'  D0 mapped: {X_d0.shape[1]}, D3 mapped: {X_d3.shape[1]}')

        # Align — take intersection of mapped features
        d0_set = set(feat_node_ids_d0); d3_set = set(feat_node_ids_d3)
        common_feat = sorted(d0_set & d3_set)
        d0_keep = [feat_node_ids_d0.index(f) for f in common_feat]
        d3_keep = [feat_node_ids_d3.index(f) for f in common_feat]
        X_d0 = X_d0[:, d0_keep]; X_d3 = X_d3[:, d3_keep]
        print(f'  Aligned features: {len(common_feat)}')

        # Pooled z-score
        X_pool = np.vstack([X_d0, X_d3])
        mu = X_pool.mean(axis=0); sd = X_pool.std(axis=0) + 1e-9
        X_d0z = (X_d0 - mu) / sd; X_d3z = (X_d3 - mu) / sd

        feat_cols = [(f'feat_{k}', geom.nid_idx[common_feat[k]]) for k in range(len(common_feat))]
        all_sids = d0_sids + d3_sids
        data = {sid: {f'feat_{k}': float(X_d0z[i, k]) for k in range(len(common_feat))}
                for i, sid in enumerate(d0_sids)}
        data.update({sid: {f'feat_{k}': float(X_d3z[i, k]) for k in range(len(common_feat))}
                     for i, sid in enumerate(d3_sids)})
        main = ModalitySetup(label='main', sigma=1.0, diffusion_t=0.0,
                              feature_cols=feat_cols, data=data)
        print(f'  Solving MAP on {len(all_sids)} samples...')
        F, _ = solve_map(geom, [main], all_sids)
        F_d0 = F[:len(pids)]; F_d3 = F[len(pids):]
        dF = F_d3 - F_d0

        beta, _, alpha_pc, pca = decompose_beta_alpha(F, geom.log_pr, n_components=5)
        beta_d0 = beta[:len(pids)]; beta_d3 = beta[len(pids):]
        alpha_d0 = alpha_pc[:len(pids)]; alpha_d3 = alpha_pc[len(pids):]
        dbeta = beta_d3 - beta_d0
        dalpha = alpha_d3 - alpha_d0

        dF_norm = np.linalg.norm(dF, axis=1)
        dF_node_mean = dF.mean(axis=0)

        print(f'  ‖ΔF_i‖: mean={dF_norm.mean():.2f}, std={dF_norm.std():.2f}, '
              f'min={dF_norm.min():.2f}, max={dF_norm.max():.2f}')
        print(f'  β D0 mean = {beta_d0.mean():+.3f}, D3 mean = {beta_d3.mean():+.3f}, Δβ mean = {dbeta.mean():+.4f}, Δβ σ = {dbeta.std():.4f}')
        for k in range(5):
            v = pca.explained_variance_ratio_[k]
            print(f'  α-PC{k+1}: D0={alpha_d0[:,k].mean():+.3f} → D3={alpha_d3[:,k].mean():+.3f}  '
                  f'Δmean={dalpha[:,k].mean():+.3f}, σ={dalpha[:,k].std():.3f}  EV={v:.2%}')

        top_idx = np.argsort(-np.abs(dF_node_mean))[:15]
        top_nodes = []
        print(f'  Top shifted nodes:')
        for j in top_idx:
            nid = geom.nodes[j]
            a = mg.graph.nodes.get(nid, {})
            name = (a.get('name', '') or nid)[:50]
            print(f'    {name:<52} mean Δ={dF_node_mean[j]:+.3f}')
            top_nodes.append({'node_id': nid, 'name': a.get('name', ''),
                              'type': a.get('node_type', ''),
                              'mean_dF': float(dF_node_mean[j]),
                              'std_dF': float(dF[:, j].std())})

        out['strata'][stratum] = {
            'n': len(pids),
            'mean_dF_norm': float(dF_norm.mean()),
            'std_dF_norm': float(dF_norm.std()),
            'mean_dbeta': float(dbeta.mean()),
            'std_dbeta': float(dbeta.std()),
            'beta_D0_mean': float(beta_d0.mean()),
            'beta_D3_mean': float(beta_d3.mean()),
            'alpha_PC_D0_mean': [float(alpha_d0[:,k].mean()) for k in range(5)],
            'alpha_PC_D3_mean': [float(alpha_d3[:,k].mean()) for k in range(5)],
            'alpha_PC_dmean': [float(dalpha[:,k].mean()) for k in range(5)],
            'alpha_PC_dstd': [float(dalpha[:,k].std()) for k in range(5)],
            'top_shifted_nodes': top_nodes,
            'per_patient_dF_norm': {p: float(dF_norm[i]) for i, p in enumerate(pids)},
            'per_patient_dbeta': {p: float(dbeta[i]) for i, p in enumerate(pids)},
        }

    out_json = RESULTS / 'filbin_d0_d3_calibration.json'
    def to_native(o):
        if isinstance(o, dict): return {k: to_native(v) for k, v in o.items()}
        if isinstance(o, list): return [to_native(v) for v in o]
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        return o
    with open(out_json, 'w') as fh:
        json.dump(to_native(out), fh, indent=2)
    print(f'\nSaved: {out_json}')


if __name__ == '__main__':
    main()
