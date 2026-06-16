"""v7 Phase 3b v4: rank-in-full-loading-vector interpretability test.

v2/v3 (basin top-K + reaction expansion) saturated at 3 of 81 cells with
overlap. The narrow basin connected-component metric loses information.

v4 tests a complementary, more permissive claim: across the full 38,211-D
α-PC loading vector, do source-paper key genes RANK HIGHER (by |loading|)
than random degree-matched gene sets?

Metric per cohort × PC:
  - For each gene node in substrate: compute |loading| on that PC
  - Compute AUROC of (truth-vs-non-truth membership) using |loading| as score
  - Null distribution: AUROC under degree-matched random truth-sized subset
  - Empirical p = P(null AUROC ≥ observed AUROC)

Hypotheses if results are strong:
  (a) Median truth-gene |loading| > median random-gene |loading|
  (b) AUROC > 0.5 across cohorts × PCs at p < 0.05 under degree-matched null
  (c) Strong-driver cohorts (IDH, LUAD) show AUROC > 0.7 on at least one PC

Reports per cohort × PC: AUROC, observed truth gene ranks, p_empirical.
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path
from collections import defaultdict
import numpy as np
import networkx as nx
from sklearn.metrics import roc_auc_score

REPO = Path('/home/jgardner/GIZMO')
RESULTS = REPO / 'benchmarks/results/unsupervised'
ZS_DIR = RESULTS / 'zscored'
CURATION = REPO / 'data/curation/v7_cohort_key_genes.tsv'
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'benchmarks'))

from gizmo.export.json_export import read_json
from gizmo.diagnostics.degree_preserving_null import (
    degree_matched_random_subsets,
)


CONFIDENCE_FILTER = ('HIGH', 'MEDIUM')
N_NULL = 1000


def load_curation():
    rows = []
    lines = CURATION.read_text().splitlines()
    header = lines[0].split('\t')
    for line in lines[1:]:
        cells = line.split('\t')
        if len(cells) != len(header):
            continue
        rows.append(dict(zip(header, cells)))
    return rows


def main():
    print('=== v7 Phase 3b v4: rank-in-full-loading-vector test ===', flush=True)
    t0 = time.time()

    curation = load_curation()
    filtered = [r for r in curation if r['confidence'] in CONFIDENCE_FILTER]
    by_cohort = defaultdict(list)
    for r in filtered:
        by_cohort[r['cohort']].append(r['key_gene_symbol'])

    print('\nLoading substrate...', flush=True)
    mg = read_json(REPO / 'data/processed/human_full/graph.json')
    g_undir = mg.graph.to_undirected() if mg.graph.is_directed() else mg.graph

    # Map symbol → substrate node ID (gene nodes only)
    gene_nodes = [n for n, a in mg.graph.nodes(data=True)
                   if a.get('node_type') == 'gene']
    sym_to_nid = {}
    nid_to_sym = {}
    for nid in gene_nodes:
        attrs = mg.graph.nodes.get(nid, {})
        sym = attrs.get('symbol') or attrs.get('name') or nid.replace('symbol:', '')
        if sym:
            sym_to_nid[sym] = nid
            nid_to_sym[nid] = sym
    print(f'  gene nodes: {len(gene_nodes)}, resolvable symbols: {len(sym_to_nid)}',
          flush=True)

    print('\n=== Per-cell results ===', flush=True)
    print(f'{"cohort":22} {"PC":>3} {"sign":>5} {"truth_res":>10} {"AUROC":>7} '
          f'{"med_rank_t":>11} {"med_rank_n":>11} {"p_emp":>8} {"top_genes":<40}',
          flush=True)
    print('-' * 130, flush=True)

    rng = np.random.default_rng(seed=42)
    all_cells = []
    for cohort, truth_symbols in by_cohort.items():
        # Resolve truth symbols to substrate gene-node IDs
        truth_node_ids = [sym_to_nid[s] for s in truth_symbols if s in sym_to_nid]
        if len(truth_node_ids) < 2:
            print(f'  {cohort}: only {len(truth_node_ids)} truth genes resolvable, '
                  f'skipping', flush=True)
            continue

        # Load PC loadings for this cohort
        pc_npz = ZS_DIR / cohort / 'alpha_pc_loadings.npz'
        if not pc_npz.exists():
            print(f'  {cohort}: no loadings file', flush=True)
            continue
        npz = np.load(pc_npz, allow_pickle=True)
        components = npz['components']  # (n_pcs, n_nodes)
        node_ids = list(npz['node_ids'])
        nid_to_idx = {n: i for i, n in enumerate(node_ids)}

        # Restrict to gene nodes present in this cohort's substrate
        cohort_gene_node_ids = [n for n in gene_nodes if n in nid_to_idx]
        if len(cohort_gene_node_ids) < 100:
            print(f'  {cohort}: only {len(cohort_gene_node_ids)} gene nodes in F',
                  flush=True)
            continue

        # Indices in F-space
        cohort_gene_idx = np.array([nid_to_idx[n] for n in cohort_gene_node_ids])

        # Restrict truth nodes to those present in F
        truth_in_F = [n for n in truth_node_ids if n in nid_to_idx]
        if len(truth_in_F) < 2:
            print(f'  {cohort}: only {len(truth_in_F)} truth in F-gene-nodes',
                  flush=True)
            continue
        truth_node_set = set(truth_in_F)
        truth_in_F_idx_in_cohort_genes = np.array([
            i for i, n in enumerate(cohort_gene_node_ids) if n in truth_node_set
        ])
        labels = np.zeros(len(cohort_gene_node_ids), dtype=int)
        labels[truth_in_F_idx_in_cohort_genes] = 1

        for pc in range(1, min(6, components.shape[0] + 1)):
            for sign in ('+', '-'):
                pc_idx = pc - 1
                loadings_full = components[pc_idx]
                # Restrict to gene nodes; sign filter on the loading itself
                gene_loadings = loadings_full[cohort_gene_idx]
                if sign == '+':
                    signed_scores = np.maximum(gene_loadings, 0)
                else:
                    signed_scores = np.maximum(-gene_loadings, 0)

                # Skip if degenerate
                if signed_scores.max() == signed_scores.min():
                    continue
                if labels.sum() == 0 or labels.sum() == len(labels):
                    continue

                # AUROC using signed_scores as the discriminator
                try:
                    auroc = roc_auc_score(labels, signed_scores)
                except ValueError:
                    continue

                # Truth gene ranks (1 = highest |loading| on this signed side)
                rank_order = np.argsort(-signed_scores)
                ranks_of_genes = np.empty(len(signed_scores), dtype=int)
                ranks_of_genes[rank_order] = np.arange(1, len(signed_scores) + 1)
                truth_ranks = ranks_of_genes[labels == 1]
                non_truth_ranks = ranks_of_genes[labels == 0]
                median_truth_rank = float(np.median(truth_ranks))
                median_non_truth_rank = float(np.median(non_truth_ranks))

                # Null: degree-matched random gene subsets of same size,
                # compute AUROC for each; empirical p = P(null >= observed)
                null_subsets = degree_matched_random_subsets(
                    g_undir, truth_in_F, n_samples=N_NULL,
                    candidate_pool=cohort_gene_node_ids,
                    n_bins=10, rng=rng)
                null_aurocs = []
                for sub in null_subsets:
                    sub_idx_in_cohort_genes = np.array([
                        i for i, n in enumerate(cohort_gene_node_ids)
                        if n in set(sub)
                    ])
                    if len(sub_idx_in_cohort_genes) < 2:
                        continue
                    null_labels = np.zeros(len(cohort_gene_node_ids), dtype=int)
                    null_labels[sub_idx_in_cohort_genes] = 1
                    try:
                        null_aurocs.append(roc_auc_score(null_labels, signed_scores))
                    except ValueError:
                        pass
                if not null_aurocs:
                    continue
                null_aurocs_arr = np.array(null_aurocs)
                p_emp = float((np.sum(null_aurocs_arr >= auroc) + 1)
                              / (len(null_aurocs_arr) + 1))

                # Top truth genes by signed score (for narrative)
                truth_idx = np.where(labels == 1)[0]
                truth_score = signed_scores[truth_idx]
                top_t = truth_idx[np.argsort(-truth_score)[:5]]
                top_t_syms = [nid_to_sym.get(cohort_gene_node_ids[i], '?')
                              for i in top_t]
                top_t_str = ', '.join(top_t_syms[:5])

                row = {
                    'cohort': cohort, 'pc': pc, 'sign': sign,
                    'truth_resolvable': len(truth_in_F),
                    'auroc': float(auroc),
                    'median_truth_rank': median_truth_rank,
                    'median_non_truth_rank': median_non_truth_rank,
                    'p_empirical': p_emp,
                    'top_truth_genes': top_t_syms,
                    'null_auroc_median': float(np.median(null_aurocs_arr)),
                }
                all_cells.append(row)
                print(f'{cohort:22} {pc:>3} {sign:>5} {len(truth_in_F):>10} '
                      f'{auroc:>7.3f} {median_truth_rank:>11.0f} '
                      f'{median_non_truth_rank:>11.0f} {p_emp:>8.4f} {top_t_str[:40]:<40}',
                      flush=True)

    # Aggregate
    print('\n=== Aggregate v7 §2 v4 headline ===', flush=True)
    if not all_cells:
        print('  NO cells produced', flush=True)
        return

    aurocs = np.array([c['auroc'] for c in all_cells])
    pvals = np.array([c['p_empirical'] for c in all_cells])
    n_total = len(all_cells)
    n_auroc_60 = int(np.sum(aurocs >= 0.60))
    n_auroc_70 = int(np.sum(aurocs >= 0.70))
    n_auroc_80 = int(np.sum(aurocs >= 0.80))
    n_p05 = int(np.sum(pvals < 0.05))
    n_p01 = int(np.sum(pvals < 0.01))

    print(f'  Total cells: {n_total}', flush=True)
    print(f'  Cells with AUROC ≥ 0.60: {n_auroc_60}/{n_total} '
          f'({100*n_auroc_60/n_total:.1f}%)', flush=True)
    print(f'  Cells with AUROC ≥ 0.70: {n_auroc_70}/{n_total} '
          f'({100*n_auroc_70/n_total:.1f}%)', flush=True)
    print(f'  Cells with AUROC ≥ 0.80: {n_auroc_80}/{n_total} '
          f'({100*n_auroc_80/n_total:.1f}%)', flush=True)
    print(f'  Cells with p_empirical < 0.05: {n_p05}/{n_total} '
          f'({100*n_p05/n_total:.1f}%)', flush=True)
    print(f'  Cells with p_empirical < 0.01: {n_p01}/{n_total} '
          f'({100*n_p01/n_total:.1f}%)', flush=True)
    print(f'  Median AUROC: {float(np.median(aurocs)):.3f}', flush=True)
    print(f'  Median p_empirical: {float(np.median(pvals)):.4f}', flush=True)

    # Per-cohort best
    print('\n=== Per-cohort best-cell (by p_empirical) ===', flush=True)
    by_cohort_cells = defaultdict(list)
    for c in all_cells:
        by_cohort_cells[c['cohort']].append(c)
    for cohort, cells in by_cohort_cells.items():
        best = min(cells, key=lambda c: c['p_empirical'])
        print(f'  {cohort:22}: best (PC{best["pc"]}{best["sign"]}): '
              f'AUROC={best["auroc"]:.3f}  p={best["p_empirical"]:.4f}  '
              f'top_truth={", ".join(best["top_truth_genes"][:3])}',
              flush=True)

    out = ZS_DIR / 'v7_interpretability_eval_v4.json'
    out.write_text(json.dumps({
        'cells': all_cells,
        'aggregate': {
            'n_cells': n_total,
            'n_auroc_60': n_auroc_60, 'n_auroc_70': n_auroc_70,
            'n_auroc_80': n_auroc_80,
            'n_p_05': n_p05, 'n_p_01': n_p01,
            'median_auroc': float(np.median(aurocs)),
            'median_p': float(np.median(pvals)),
        },
        'config': {
            'metric': 'AUROC of truth-vs-non-truth gene membership using |loading| as score',
            'null': 'degree-matched random gene subsets',
            'n_null': N_NULL,
        },
        'compute_seconds': time.time() - t0,
    }, indent=2))
    print(f'\nWrote {out}', flush=True)
    print(f'Compute: {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
