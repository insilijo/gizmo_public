"""v7 Phase 2.2b: focused test of v6 §5 headline rescue under degree-preserving null.

The Phase 2.2 sweep resolved 2HG searches to the cytosolic node (R-ALL-880042),
which is NOT the node where the v6 §5 headline rescue lives. The headline
rescue is the MITOCHONDRIAL 2HG node (R-ALL-879997) at rank 26/6,406 under
z-score preprocessing.

This script tests the headline rescue against the proper null:

  TCGA_IDH_glioma mito-2HG (R-ALL-879997)
    - rank 26 under z-score (verified Phase 2.2 prior session)
    - p_random vs random metab null
    - p_degree vs degree-matched null
    - p_pagerank vs PageRank-matched null

And for completeness, the Trautwein cohort on the same node.

This is the cleanest single test of: does the v6 §5 §"GoF retraction"
narrative survive hub-bias correction?
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
import numpy as np
import networkx as nx

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.diagnostics.degree_preserving_null import (
    degree_matched_random_subsets, pagerank_matched_random_subsets,
    empirical_p_value,
)


# Explicit node IDs — bypass MetaboliteMapper which resolves to cytosolic.
ANCHORS = {
    'reactome:R-ALL-879997': 'mito-2HG',       # the v6 §5 headline rescue node
    'reactome:R-ALL-880042': 'cytosolic-2HG',  # what the mapper picked
}

COHORTS = ['IDH_glioma', 'TCGA_IDH_glioma']


def get_labels(cohort, pids):
    from per_patient_master import __dict__ as ppm
    fn = ppm[{
        'IDH_glioma': 'load_idh_glioma',
        'TCGA_IDH_glioma': 'load_tcga_idh_glioma',
    }[cohort]]
    loaded = fn()
    ylabel = loaded[2]
    return np.array([1 if ylabel.get(p) == 'active' else 0 for p in pids])


def main():
    print('=== v7 Phase 2.2b: mito-2HG vs cyto-2HG focused test ===', flush=True)
    t0 = time.time()

    print('\nLoading substrate...', flush=True)
    mg = read_json(REPO / 'data/processed/human_full/graph.json')
    g_undir = mg.graph.to_undirected() if mg.graph.is_directed() else mg.graph
    print(f'  computing PageRank...', flush=True)
    pr = nx.pagerank(g_undir)
    print(f'  PR done, {time.time() - t0:.1f}s', flush=True)

    metab_pool = [n for n, a in mg.graph.nodes(data=True)
                   if a.get('node_type') == 'metabolite']

    N_SAMPLES = 10_000

    all_results = {}
    for cohort in COHORTS:
        f_path = RESULTS / f'stage3_F_{cohort}_zscored.npz'
        if not f_path.exists():
            print(f'\nSKIP {cohort}: no F', flush=True)
            continue
        print(f'\n--- {cohort} ---', flush=True)
        npz = np.load(f_path, allow_pickle=True)
        F = npz['F']
        pids = list(npz['patient_ids'])
        nodes = list(npz['node_ids'])
        labels = get_labels(cohort, pids)
        n_act, n_ctl = int(labels.sum()), int(len(labels) - labels.sum())
        print(f'  n_active={n_act}, n_control={n_ctl}', flush=True)

        # Restrict to metabolites in F + compute Cohen's d
        metab_mask = np.array([
            mg.graph.nodes.get(n, {}).get('node_type') == 'metabolite'
            for n in nodes
        ])
        metab_in_F = [nodes[i] for i, m in enumerate(metab_mask) if m]
        F_m = F[:, metab_mask]
        mut = F_m[labels == 1]; wt = F_m[labels == 0]
        mut_m, wt_m = mut.mean(axis=0), wt.mean(axis=0)
        pooled = np.sqrt((mut.var(axis=0) + wt.var(axis=0)) / 2 + 1e-12)
        d = (mut_m - wt_m) / pooled
        order = np.argsort(-d)
        rank_of = {metab_in_F[order[r]]: r + 1 for r in range(len(metab_in_F))}

        rng = np.random.default_rng(seed=42)
        cohort_results = []
        for nid, label in ANCHORS.items():
            if nid not in rank_of:
                print(f'  {label} ({nid}): NOT in F-metab', flush=True)
                continue
            obs_rank = rank_of[nid]
            obs_pct = 100 * obs_rank / len(metab_in_F)
            obs_d = float(d[metab_in_F.index(nid)])

            # Random null
            other = [n for n in metab_in_F if n != nid]
            random_indices = rng.choice(len(other), size=N_SAMPLES, replace=True)
            random_ranks = np.array([rank_of[other[i]] for i in random_indices])

            # Degree-matched
            deg_subsets = degree_matched_random_subsets(
                g_undir, [nid], n_samples=N_SAMPLES,
                candidate_pool=metab_in_F, n_bins=10, rng=rng)
            deg_ranks = np.array([rank_of[s[0]] for s in deg_subsets])

            # PageRank-matched
            pr_subsets = pagerank_matched_random_subsets(
                g_undir, [nid], n_samples=N_SAMPLES,
                candidate_pool=metab_in_F, n_bins=10, pagerank=pr, rng=rng)
            pr_ranks = np.array([rank_of[s[0]] for s in pr_subsets])

            p_random = empirical_p_value(obs_rank, random_ranks, tail='lower')
            p_degree = empirical_p_value(obs_rank, deg_ranks, tail='lower')
            p_pagerank = empirical_p_value(obs_rank, pr_ranks, tail='lower')

            row = {
                'label': label, 'node_id': nid,
                'name': mg.graph.nodes[nid].get('name', nid),
                'observed_rank': int(obs_rank),
                'observed_pct': float(obs_pct),
                'observed_d': obs_d,
                'p_random': float(p_random),
                'p_degree': float(p_degree),
                'p_pagerank': float(p_pagerank),
                'null_random_median_rank': float(np.median(random_ranks)),
                'null_degree_median_rank': float(np.median(deg_ranks)),
                'null_pagerank_median_rank': float(np.median(pr_ranks)),
            }
            cohort_results.append(row)
            print(f'  {label:20} ({nid}):', flush=True)
            print(f'    rank {obs_rank}/{len(metab_in_F)} ({obs_pct:.2f}%)  '
                  f'd={obs_d:+.3f}', flush=True)
            print(f'    null medians: random={np.median(random_ranks):.0f}, '
                  f'degree={np.median(deg_ranks):.0f}, pagerank={np.median(pr_ranks):.0f}',
                  flush=True)
            print(f'    p_random={p_random:.4f}  p_degree={p_degree:.4f}  '
                  f'p_pagerank={p_pagerank:.4f}', flush=True)

        all_results[cohort] = {
            'n_active': n_act, 'n_control': n_ctl,
            'anchors': cohort_results,
        }

    out_path = RESULTS / 'zscored' / 'v7_mito_2hg_degree_null.json'
    out_path.write_text(json.dumps(all_results, indent=2))
    print(f'\nWrote {out_path}', flush=True)
    print(f'Total compute: {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
