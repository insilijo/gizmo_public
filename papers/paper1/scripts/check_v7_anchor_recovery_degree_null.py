"""v7 Phase 2.2: §3 anchor recovery under degree- AND PageRank-preserving nulls.

For each of the 4 horizontal-meta cohorts (Trautwein RNA-only,
TCGA_IDH RNA-only, HMP2 metab-only, GSE89408 RNA-only) under z-score
preprocessing, we test whether the cohort's literature-anchor metabolite
ranks well-vs-null in the per-patient mean |F| ranking among substrate
metabolite nodes.

Three nulls tested side by side:
  (1) Random-metabolite (v5 null — uniform over all metabolite nodes)
  (2) Degree-matched (new)
  (3) PageRank-matched (new — tighter)

For each anchor in each cohort:
  - Resolve the substrate node ID (via MetaboliteMapper)
  - Compute observed rank percentile among metabolite nodes
  - Generate 10,000 null subsets of same size, compute null rank percentiles
  - Empirical p = P(null ≤ observed)
  - Report under all three nulls

Joint p across cohorts via Stouffer's. Reports which v6 numbers survive
the degree-bias correction.
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
from gizmo.evidence.mappers import MetaboliteMapper
from gizmo.diagnostics.degree_preserving_null import (
    degree_matched_random_subsets, pagerank_matched_random_subsets,
    empirical_p_value, stouffer_combine,
)


# Cohorts with z-scored F + anchor metabolites (literature-derived).
# Searches are substring-matched against substrate node names by
# MetaboliteMapper for canonical resolution.
COHORT_ANCHORS = {
    'IDH_glioma': [
        # Trautwein RNA-only path: 2HG-mito is the canonical IDH-mut neomorphic
        # product. Validated in v6 §5 at rank 26 under z-score.
        ('2-hydroxyglutarate', '2HG'),
        ('(R)-2-hydroxyglutarate', '2HG-R'),
        ('alpha-hydroxyglutarate', '2HG-alt'),
    ],
    'TCGA_IDH_glioma': [
        ('2-hydroxyglutarate', '2HG'),
        ('(R)-2-hydroxyglutarate', '2HG-R'),
    ],
    'HMP2_IBD_CD': [
        ('propanoate', 'propionate'),
        ('propionate', 'propionate-alt'),
        ('lithocholate', 'lithocholate'),
        ('lithocholic acid', 'lithocholate-alt'),
    ],
    'GSE89408_RA': [
        ('citrulline', 'citrulline'),
        ('L-citrulline', 'L-citrulline'),
    ],
}


def get_labels(cohort, pids):
    from per_patient_master import __dict__ as ppm
    loader_map = {
        'IDH_glioma':       'load_idh_glioma',
        'TCGA_IDH_glioma':  'load_tcga_idh_glioma',
        'HMP2_IBD_CD':      'load_hmp2_ibd_cd',
        'GSE89408_RA':      'load_gse89408_ra',
    }
    fn = ppm[loader_map[cohort]]
    loaded = fn()
    ylabel = loaded[2] if len(loaded) >= 3 else {}
    labels = np.array([1 if ylabel.get(p) == 'active' else 0 for p in pids])
    return labels


def resolve_anchor_node(mapper, search_str, mg):
    """Try MetaboliteMapper first; fall back to substring scan over metabolite nodes."""
    res = mapper.map(search_str)
    nid = res[0] if isinstance(res, tuple) else res
    if nid and nid in mg.graph.nodes:
        return nid
    # Substring fallback
    search_lc = search_str.lower()
    candidates = []
    for n, attrs in mg.graph.nodes(data=True):
        if attrs.get('node_type') != 'metabolite':
            continue
        name = (attrs.get('name', '') or '').lower()
        if name == search_lc or search_lc in name:
            candidates.append((n, attrs.get('name', '')))
    return candidates[0][0] if candidates else None


def main():
    print('=== v7 Phase 2.2: §3 anchor recovery under degree- and PR-matched nulls ===',
          flush=True)
    t0 = time.time()

    print('\nLoading substrate...', flush=True)
    mg = read_json(REPO / 'data/processed/human_full/graph.json')
    mmap = MetaboliteMapper(mg)
    print(f'  substrate: {mg.graph.number_of_nodes()} nodes', flush=True)

    # Build undirected metabolite-restricted subgraph for PageRank
    metab_nodes_full = [n for n, a in mg.graph.nodes(data=True)
                         if a.get('node_type') == 'metabolite']
    print(f'  metabolite nodes in substrate: {len(metab_nodes_full)}', flush=True)

    g_undirected = mg.graph.to_undirected() if mg.graph.is_directed() else mg.graph
    print('  computing PageRank on full graph (one-time)...', flush=True)
    pr = nx.pagerank(g_undirected)
    print(f'  PageRank done, {time.time() - t0:.1f}s elapsed', flush=True)

    # Resolve all anchor node IDs up front
    print('\nResolving anchor node IDs...', flush=True)
    resolved_anchors = {}
    for cohort, anchor_list in COHORT_ANCHORS.items():
        cohort_anchors = []
        for search, label in anchor_list:
            nid = resolve_anchor_node(mmap, search, mg)
            if nid is None:
                print(f'  {cohort} {label}: NOT resolved', flush=True)
            else:
                name = mg.graph.nodes[nid].get('name', nid)
                print(f'  {cohort} {label}: {nid} ({name})', flush=True)
                cohort_anchors.append((label, nid, name))
        resolved_anchors[cohort] = cohort_anchors

    # Process each cohort
    print('\nProcessing cohorts...', flush=True)
    all_results = {}
    for cohort, anchors in resolved_anchors.items():
        if not anchors:
            print(f'\n  {cohort}: SKIP (no anchors resolved)', flush=True)
            continue
        f_path = RESULTS / f'stage3_F_{cohort}_zscored.npz'
        if not f_path.exists():
            print(f'\n  {cohort}: SKIP (z-scored F missing at {f_path})', flush=True)
            continue
        print(f'\n--- {cohort} ---', flush=True)

        npz = np.load(f_path, allow_pickle=True)
        F = npz['F']
        pids = list(npz['patient_ids'])
        nodes = list(npz['node_ids'])
        print(f'  F: {F.shape}, n_patients={len(pids)}, n_nodes={len(nodes)}',
              flush=True)

        labels = get_labels(cohort, pids)
        n_act = int(labels.sum()); n_ctl = len(labels) - n_act
        print(f'  n_active={n_act}, n_control={n_ctl}', flush=True)
        if n_act < 3 or n_ctl < 3:
            print(f'  insufficient class balance, skipping', flush=True)
            continue

        # Restrict to metabolite nodes in F
        metab_mask = np.array([
            mg.graph.nodes.get(n, {}).get('node_type') == 'metabolite'
            for n in nodes
        ])
        metab_nids_in_F = [nodes[i] for i, m in enumerate(metab_mask) if m]
        F_m = F[:, metab_mask]
        mut = F_m[labels == 1]; wt = F_m[labels == 0]
        mut_m, wt_m = mut.mean(axis=0), wt.mean(axis=0)
        pooled = np.sqrt((mut.var(axis=0) + wt.var(axis=0)) / 2 + 1e-12)
        d = (mut_m - wt_m) / pooled
        order = np.argsort(-d)  # descending
        rank_of = {metab_nids_in_F[order[r]]: r + 1 for r in range(len(metab_nids_in_F))}
        print(f'  n_metab_in_F: {len(metab_nids_in_F)}', flush=True)

        # Per-anchor empirical p-values under three nulls
        cohort_results = []
        # Restrict null sampling pool to metabolite nodes (matches the search space)
        # Use the metabolite nodes from F (not the full substrate metab pool)
        candidate_pool = metab_nids_in_F

        # Generate null subsets ONCE per cohort (shared across anchors)
        present_anchors = [(label, nid, name)
                            for label, nid, name in anchors
                            if nid in rank_of]
        if not present_anchors:
            print(f'  none of {len(anchors)} anchors present in F-metab', flush=True)
            continue

        # Treat each anchor as its own size-1 query set for proper per-anchor matching
        N_SAMPLES = 10_000
        rng = np.random.default_rng(seed=42)

        for label, nid, name in present_anchors:
            observed_rank = rank_of[nid]
            observed_d = float(d[metab_nids_in_F.index(nid)])
            observed_pct = 100 * observed_rank / len(metab_nids_in_F)

            # Random-metabolite null (uniform over candidate pool, excluding the anchor)
            other = [n for n in candidate_pool if n != nid]
            random_indices = rng.choice(len(other), size=N_SAMPLES, replace=True)
            random_ranks = np.array([rank_of[other[i]] for i in random_indices])

            # Degree-matched null
            deg_subsets = degree_matched_random_subsets(
                g_undirected, [nid], n_samples=N_SAMPLES,
                candidate_pool=candidate_pool, n_bins=10, rng=rng)
            deg_ranks = np.array([rank_of[s[0]] for s in deg_subsets])

            # PageRank-matched null
            pr_subsets = pagerank_matched_random_subsets(
                g_undirected, [nid], n_samples=N_SAMPLES,
                candidate_pool=candidate_pool, n_bins=10, pagerank=pr, rng=rng)
            pr_ranks = np.array([rank_of[s[0]] for s in pr_subsets])

            # Empirical p-values (lower tail because better rank = lower number)
            p_random = empirical_p_value(observed_rank, random_ranks, tail='lower')
            p_degree = empirical_p_value(observed_rank, deg_ranks, tail='lower')
            p_pagerank = empirical_p_value(observed_rank, pr_ranks, tail='lower')

            row = {
                'label': label, 'node_id': nid, 'name': name,
                'observed_rank': int(observed_rank),
                'observed_pct': float(observed_pct),
                'observed_d': observed_d,
                'p_random_metab_null': float(p_random),
                'p_degree_matched_null': float(p_degree),
                'p_pagerank_matched_null': float(p_pagerank),
                'null_random_median_rank': float(np.median(random_ranks)),
                'null_degree_median_rank': float(np.median(deg_ranks)),
                'null_pagerank_median_rank': float(np.median(pr_ranks)),
            }
            cohort_results.append(row)
            print(f'  {label:20} {name[:30]:30} rank {observed_rank:5d}/{len(metab_nids_in_F)} '
                  f'({observed_pct:.2f}%) d={observed_d:+.3f}', flush=True)
            print(f'    p_random={p_random:.4f}  p_degree={p_degree:.4f}  '
                  f'p_pagerank={p_pagerank:.4f}', flush=True)
            print(f'    null medians:  random={np.median(random_ranks):.0f}  '
                  f'degree={np.median(deg_ranks):.0f}  pagerank={np.median(pr_ranks):.0f}',
                  flush=True)

        all_results[cohort] = {
            'n_metab_in_F': len(metab_nids_in_F),
            'n_active': n_act, 'n_control': n_ctl,
            'anchors': cohort_results,
        }

    # Joint p via Stouffer across the best anchor per cohort under each null
    print('\n=== Joint p across cohorts (Stouffer) ===', flush=True)
    print('Using best (min-p) anchor per cohort under each null.', flush=True)
    best_p_random, best_p_degree, best_p_pagerank = [], [], []
    for cohort, r in all_results.items():
        if not r['anchors']:
            continue
        best_random = min(a['p_random_metab_null'] for a in r['anchors'])
        best_degree = min(a['p_degree_matched_null'] for a in r['anchors'])
        best_pagerank = min(a['p_pagerank_matched_null'] for a in r['anchors'])
        print(f'  {cohort:22} best p_random={best_random:.4f}  '
              f'p_degree={best_degree:.4f}  p_pagerank={best_pagerank:.4f}',
              flush=True)
        best_p_random.append(best_random)
        best_p_degree.append(best_degree)
        best_p_pagerank.append(best_pagerank)

    if best_p_random:
        joint_random = stouffer_combine(best_p_random)
        joint_degree = stouffer_combine(best_p_degree)
        joint_pagerank = stouffer_combine(best_p_pagerank)
        print(f'\n  JOINT p (Stouffer, n_cohorts={len(best_p_random)}):',
              flush=True)
        print(f'    under random-metab null:    {joint_random:.6f}', flush=True)
        print(f'    under degree-matched null:  {joint_degree:.6f}', flush=True)
        print(f'    under PageRank-matched null: {joint_pagerank:.6f}', flush=True)

    out_path = RESULTS / 'zscored' / 'v7_anchor_recovery_degree_null.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        'cohorts': all_results,
        'joint_p_random_metab': float(joint_random) if best_p_random else None,
        'joint_p_degree_matched': float(joint_degree) if best_p_random else None,
        'joint_p_pagerank_matched': float(joint_pagerank) if best_p_random else None,
        'n_samples_per_null': N_SAMPLES,
        'compute_seconds': time.time() - t0,
    }, indent=2))
    print(f'\nWrote {out_path}', flush=True)
    print(f'Total compute: {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
