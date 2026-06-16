"""AAA — Systematic substrate-connectivity filter across cohorts.

For each cohort's top-50 z>10 PC contributors (per PC, all 15 PCs):
  - Find substrate node + 1-hop neighbors
  - Count co-loading neighbors with |z|>3 on same axis
  - Tier candidates by substrate-connectivity:
      Tier 1 (≥2 co-loaders): strong substrate-grounded
      Tier 2 (1 co-loader):   moderate
      Tier 3 (0 co-loaders):  isolated (lower-confidence)

Output: cross_disease_substrate_connectivity_atlas.csv with one row per
        (cohort × PC × contributor) annotated with z-score + co-loader count
        + tier + drug annotation lookup.
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


# Drug annotation set (compact reuse)
DRUGGED = set([
    'JAK1','JAK2','JAK3','TYK2','STAT1','STAT3','STAT4','STAT6','IFNAR1','IFNAR2',
    'IL6R','IL6','IL1B','TNF','TNFRSF1A','TNFRSF1B','IL10','IL12B','IL23R','IL17A',
    'DHFR','TYMS','MTHFR','MTR','SHMT1','SHMT2','GART','ATIC','FPGS',
    'PSMB8','PSMB9','PSMB10','MS4A1','CD19','CD20','CD38','TNFSF13B','BAFF',
    'PRDM1','XBP1','IRF4','EGFR','ERBB2','VEGFA','KDR','MET','ALK',
    'BRAF','KRAS','PARP1','PIK3CA','AKT1','MTOR','CDK4','CDK6','BCL2','TP53',
    'MDM2','HDAC1','HDAC6','EZH2','DNMT1','PDCD1','PDL1','CD274','CTLA4','LAG3',
    'TIGIT','HIF1A','EPAS1','PPARG','PPARA','HMGCR','LDLR','GLP1R','GIPR','GCG',
    'INS','INSR','LEPR','MC4R','CETP','PCSK9','AGTR1','AGT','ACE','ACE2','REN',
    'BTK','SYK','LYN','LCK','BRD4','PARP2',
])

COHORTS = [
    ('Wang_RA',         CACHE / 'Wang_RA_alpha_t0.npy'),
    ('Filbin_COVID',    CACHE / 'Filbin_alpha_t0.npy'),
    ('TB_cure',         CACHE / 'TB_alpha_t0.npy'),
    ('bariatric_GLP1',  SCAN / 'bariatric_GLP1_alpha.npy'),
    ('adipose_tumor',   SCAN / 'adipose_tumor_alpha.npy'),
    ('atherogenic',     SCAN / 'atherogenic_alpha.npy'),
    ('MDD_state',       SCAN / 'MDD_state_alpha.npy'),
    ('CPTAC_CCRCC',     SCAN / 'CPTAC_CCRCC_alpha.npy'),
    ('CPTAC_COAD',      SCAN / 'CPTAC_COAD_alpha.npy'),
    ('CPTAC_OV',        SCAN / 'CPTAC_OV_alpha.npy'),
]


def main():
    print('Loading substrate...', flush=True)
    import loocv_longitudinal_audit as audit_mod
    geom = audit_mod.geom; mg = audit_mod.mg
    G = mg.graph.to_undirected() if mg.graph.is_directed() else mg.graph
    n_dim = 86826
    mu_null = np.sqrt(2.0 / (np.pi * n_dim))
    sd_null = np.sqrt((np.pi - 2.0) / (np.pi * n_dim))

    rows = []
    for ck, apath in COHORTS:
        if not apath.exists():
            print(f'  [{ck}] missing; skip', flush=True); continue
        a = np.load(apath)
        n_c = min(15, a.shape[0] - 1)
        pca = PCA(n_components=n_c, svd_solver='randomized',
                   random_state=0).fit(a)
        del a; gc.collect()
        print(f'\n  ▼ {ck} ({n_c} PCs)', flush=True)

        for k in range(n_c):
            v = pca.components_[k]
            z = ((np.abs(v) - mu_null) / sd_null) * np.sign(v)
            # Top-50 z>10 contributors only
            cand_idx = np.where(np.abs(z) > 10)[0]
            cand_sorted = cand_idx[np.argsort(-np.abs(z[cand_idx]))][:50]
            for idx in cand_sorted:
                nid = geom.nodes[idx]
                attrs = mg.graph.nodes.get(nid, {})
                syms = attrs.get('gene_symbols', [])
                sym = syms[0] if syms else nid
                sym_clean = sym.replace('symbol:', '').strip()
                if len(sym_clean) < 2: continue
                # 1-hop substrate neighbors
                if nid not in G: continue
                neighbors = list(G.neighbors(nid))
                strong_coload = 0
                for nb in neighbors:
                    if nb not in geom.nid_idx: continue
                    nb_z = ((np.abs(v[geom.nid_idx[nb]]) - mu_null) / sd_null
                            * np.sign(v[geom.nid_idx[nb]]))
                    if abs(nb_z) > 3:
                        strong_coload += 1
                # Tier
                if strong_coload >= 2: tier = 1
                elif strong_coload == 1: tier = 2
                else: tier = 3
                is_drugged = sym_clean in DRUGGED
                rows.append({
                    'cohort': ck, 'PC': k+1,
                    'EV_pct': float(pca.explained_variance_ratio_[k] * 100),
                    'symbol': sym_clean,
                    'z_score': float(z[idx]),
                    'n_substrate_neighbors': len(neighbors),
                    'n_coloading_neighbors': strong_coload,
                    'tier': tier,
                    'in_drug_annotation': is_drugged,
                })
        gc.collect()

    df = pd.DataFrame(rows)
    df.to_csv(SCAN / 'cross_disease_substrate_connectivity_atlas.csv',
              index=False)
    print(f'\n  Atlas: {len(df)} entries; saved to '
          f'cross_disease_substrate_connectivity_atlas.csv', flush=True)

    # Summary: per-cohort tier counts of NOVEL targets
    print('\n' + '='*80, flush=True)
    print('NOVEL TARGETS per cohort by tier (z>10, not in drug annotation)',
          flush=True)
    print('='*80, flush=True)
    novel = df[~df.in_drug_annotation].drop_duplicates(
        subset=['cohort', 'symbol'])
    print(f'  {"cohort":<22}{"tier 1 (≥2)":>14}{"tier 2 (=1)":>14}'
          f'{"tier 3 (=0)":>14}', flush=True)
    for ck in novel.cohort.unique():
        sub = novel[novel.cohort == ck]
        t1 = (sub.tier == 1).sum()
        t2 = (sub.tier == 2).sum()
        t3 = (sub.tier == 3).sum()
        print(f'  {ck:<22}{t1:>14}{t2:>14}{t3:>14}', flush=True)

    # Best tier-1 novel candidates per cohort (top 5)
    print('\n  TIER-1 high-confidence novel targets per cohort (top 5):', flush=True)
    for ck in novel.cohort.unique():
        sub = novel[(novel.cohort == ck) & (novel.tier == 1)]
        sub = sub.sort_values('z_score', key=np.abs, ascending=False).head(5)
        if sub.empty: continue
        print(f'\n  {ck}:', flush=True)
        for _, r in sub.iterrows():
            print(f'    PC{int(r.PC):>2}  z={r.z_score:+6.1f}  '
                  f'{r.symbol:<14}  co-loaders={int(r.n_coloading_neighbors)}',
                  flush=True)


if __name__ == '__main__':
    main()
