"""Degree- and PageRank-preserving null-distribution generators.

For substrate enrichment statistics, the appropriate null distribution
matches the degree (or PageRank centrality) of the query nodes. Naive
uniform-random null sets are systematically biased: well-studied genes
and disease-anchor metabolites are graph hubs, and any test against
uniform random nodes will look enriched by hub-membership alone.

Public API:

    degree_matched_random_subsets(graph, query_nodes, n_samples, ...)
    pagerank_matched_random_subsets(graph, query_nodes, n_samples, ...)
    empirical_p_value(observed, null_distribution, tail='lower')

References:
    The degree-preserving null is the workhorse null for
    network-enrichment tests under the multi-degree bias confound
    described in the v7 manuscript Methods §"Degree- and PageRank-
    preserving null".
"""
from __future__ import annotations

import logging
from typing import Optional, Sequence

import networkx as nx
import numpy as np

log = logging.getLogger(__name__)


def _node_degrees(graph: nx.Graph, candidate_pool: Sequence[str]) -> np.ndarray:
    """Return array of degrees for each node in candidate_pool."""
    return np.array([graph.degree(n) for n in candidate_pool], dtype=np.int64)


def _bin_assignments(values: np.ndarray, n_bins: int) -> np.ndarray:
    """Assign each value to a bin 0..n_bins-1 by quantile.

    Returns an integer array of bin indices, same shape as `values`.
    """
    # Use quantile-based bins so each bin has roughly equal candidate count.
    # This handles skewed distributions (degree, log-PR) better than fixed bins.
    if len(values) <= n_bins:
        return np.arange(len(values))
    quantile_edges = np.quantile(values, np.linspace(0, 1, n_bins + 1))
    quantile_edges[-1] = quantile_edges[-1] + 1e-9  # ensure rightmost is included
    bins = np.digitize(values, quantile_edges[1:-1], right=False)
    return bins


def degree_matched_random_subsets(
    graph: nx.Graph,
    query_nodes: Sequence[str],
    n_samples: int = 1000,
    candidate_pool: Optional[Sequence[str]] = None,
    n_bins: int = 10,
    rng: Optional[np.random.Generator] = None,
) -> list[list[str]]:
    """Generate ``n_samples`` random node subsets degree-matched to query_nodes.

    For each query node, find candidates in the same degree-quantile bin
    and sample one. This produces a random subset of the same size as
    query_nodes, with the same per-node-degree distribution.

    Parameters
    ----------
    graph : networkx.Graph
        The full substrate graph (or appropriate subgraph).
    query_nodes : sequence of node IDs
        The observed/query set whose enrichment we're testing.
    n_samples : int
        Number of random null subsets to generate.
    candidate_pool : sequence of node IDs, optional
        Restrict sampling to this pool (e.g., metabolite nodes only).
        Defaults to all nodes in the graph.
    n_bins : int
        Number of degree-quantile bins. More bins = tighter matching
        but smaller candidate pool per bin.
    rng : np.random.Generator, optional
        Random generator. Defaults to ``np.random.default_rng()``.

    Returns
    -------
    list of lists
        ``n_samples`` random subsets, each the same size as query_nodes,
        each degree-matched to query_nodes.
    """
    if rng is None:
        rng = np.random.default_rng()

    candidate_pool = list(candidate_pool) if candidate_pool is not None \
        else list(graph.nodes())
    cand_pool_set = set(candidate_pool)

    # Exclude query nodes from the candidate pool (no self-sampling)
    candidate_pool = [n for n in candidate_pool if n not in set(query_nodes)]

    cand_degrees = _node_degrees(graph, candidate_pool)
    cand_bins = _bin_assignments(cand_degrees, n_bins)

    # For each query node, determine its bin and find candidates in same bin
    query_degrees = _node_degrees(graph, query_nodes)
    # Compute bin edges from candidate pool (so query nodes can map onto them)
    quantile_edges = np.quantile(cand_degrees, np.linspace(0, 1, n_bins + 1))
    quantile_edges[-1] = quantile_edges[-1] + 1e-9
    query_bins = np.digitize(query_degrees, quantile_edges[1:-1], right=False)

    # Group candidates by bin
    cand_by_bin: dict[int, list[int]] = {}
    for i, b in enumerate(cand_bins):
        cand_by_bin.setdefault(int(b), []).append(i)

    null_subsets: list[list[str]] = []
    for s in range(n_samples):
        subset = []
        for qb in query_bins:
            qb = int(qb)
            if qb in cand_by_bin and cand_by_bin[qb]:
                idx = int(rng.choice(cand_by_bin[qb]))
                subset.append(candidate_pool[idx])
            else:
                # Fallback: nearest non-empty bin
                for offset in range(1, n_bins):
                    for candidate_b in (qb - offset, qb + offset):
                        if candidate_b in cand_by_bin and cand_by_bin[candidate_b]:
                            idx = int(rng.choice(cand_by_bin[candidate_b]))
                            subset.append(candidate_pool[idx])
                            break
                    else:
                        continue
                    break
                else:
                    # Worst case: just sample uniformly from any candidate
                    subset.append(candidate_pool[int(rng.choice(len(candidate_pool)))])
        null_subsets.append(subset)
    return null_subsets


def pagerank_matched_random_subsets(
    graph: nx.Graph,
    query_nodes: Sequence[str],
    n_samples: int = 1000,
    candidate_pool: Optional[Sequence[str]] = None,
    n_bins: int = 10,
    pagerank: Optional[dict[str, float]] = None,
    rng: Optional[np.random.Generator] = None,
) -> list[list[str]]:
    """Generate ``n_samples`` random node subsets PageRank-matched to query_nodes.

    Same as ``degree_matched_random_subsets`` but bins by log-PageRank
    instead of raw degree. Use when degree is too coarse a proxy for
    centrality (e.g., when comparing across substrate node types with
    very different degree distributions).

    Parameters
    ----------
    pagerank : dict, optional
        Precomputed PageRank dict. Defaults to ``nx.pagerank(graph)``.
        Passing in a precomputed PageRank avoids re-running it for
        every test against the same graph.
    """
    if rng is None:
        rng = np.random.default_rng()
    if pagerank is None:
        pagerank = nx.pagerank(graph)

    candidate_pool = list(candidate_pool) if candidate_pool is not None \
        else list(graph.nodes())
    candidate_pool = [n for n in candidate_pool if n not in set(query_nodes)]

    cand_logpr = np.array(
        [np.log10(pagerank.get(n, 1e-15) + 1e-15) for n in candidate_pool])
    cand_bins = _bin_assignments(cand_logpr, n_bins)

    query_logpr = np.array(
        [np.log10(pagerank.get(n, 1e-15) + 1e-15) for n in query_nodes])
    quantile_edges = np.quantile(cand_logpr, np.linspace(0, 1, n_bins + 1))
    quantile_edges[-1] = quantile_edges[-1] + 1e-9
    query_bins = np.digitize(query_logpr, quantile_edges[1:-1], right=False)

    cand_by_bin: dict[int, list[int]] = {}
    for i, b in enumerate(cand_bins):
        cand_by_bin.setdefault(int(b), []).append(i)

    null_subsets: list[list[str]] = []
    for s in range(n_samples):
        subset = []
        for qb in query_bins:
            qb = int(qb)
            if qb in cand_by_bin and cand_by_bin[qb]:
                idx = int(rng.choice(cand_by_bin[qb]))
                subset.append(candidate_pool[idx])
            else:
                for offset in range(1, n_bins):
                    for candidate_b in (qb - offset, qb + offset):
                        if candidate_b in cand_by_bin and cand_by_bin[candidate_b]:
                            idx = int(rng.choice(cand_by_bin[candidate_b]))
                            subset.append(candidate_pool[idx])
                            break
                    else:
                        continue
                    break
                else:
                    subset.append(candidate_pool[int(rng.choice(len(candidate_pool)))])
        null_subsets.append(subset)
    return null_subsets


def empirical_p_value(
    observed: float,
    null_distribution: np.ndarray,
    tail: str = 'lower',
) -> float:
    """Compute empirical p-value of ``observed`` against ``null_distribution``.

    For rank statistics (smaller = better), use tail='lower' (count
    null draws ≤ observed). For enrichment statistics (larger = better),
    use tail='upper' (count null draws ≥ observed).

    Returns the (count + 1) / (n + 1) Phipson & Smyth correction to
    avoid p = 0 when observed beats all null draws.
    """
    n = len(null_distribution)
    if tail == 'lower':
        count = int(np.sum(null_distribution <= observed))
    elif tail == 'upper':
        count = int(np.sum(null_distribution >= observed))
    else:
        raise ValueError(f"tail must be 'lower' or 'upper', got {tail!r}")
    return (count + 1) / (n + 1)


def stouffer_combine(p_values: Sequence[float], weights: Optional[Sequence[float]] = None) -> float:
    """Combine independent p-values via Stouffer's Z method.

    Returns the combined p-value. If ``weights`` is None, each test
    gets equal weight.
    """
    from scipy.stats import norm  # type: ignore
    p_arr = np.asarray(p_values, dtype=np.float64)
    # Clip to avoid -inf when p is exactly 0 or 1
    p_arr = np.clip(p_arr, 1e-300, 1.0 - 1e-15)
    z = norm.isf(p_arr)  # one-sided z-score (upper-tail inverse survival)
    if weights is None:
        w = np.ones_like(z)
    else:
        w = np.asarray(weights, dtype=np.float64)
    z_combined = (w * z).sum() / np.sqrt((w ** 2).sum())
    return float(norm.sf(z_combined))
