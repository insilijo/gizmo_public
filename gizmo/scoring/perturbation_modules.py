"""Connected-subgraph (module) discovery from Laplacian reaction scores.

Instead of returning a flat ranked list of top-perturbed reactions —
which mixes truly-coherent biology with drift — group the top-K
reactions into **connected subgraphs** (modules) and rank by aggregate
perturbation mass. Each module is evidence of a coherent phenotype
because its reactions are biochemically connected through shared
substrates, products, or catalyzing genes.

Two reactions are linked if they share at least one first-order
neighbor (a substrate, product, or catalyzing gene). Connected
components in this reaction-reaction graph are the modules.

Each module is annotated with:
  - reactions          : list of reaction node IDs
  - size               : number of reactions
  - total_perturbation : Σ|signed| across reactions
  - mean_perturbation  : avg per-reaction perturbation
  - top_pathways       : Reactome pathway IDs by count
  - genes              : catalyzing genes
  - direction          : majority sign (UP / DOWN / mixed)
  - anchor_status      : confirmed (all in anchor), novel (none),
                          bridging (mixed); only meaningful when an
                          anchor is supplied
  - n_in_anchor        : count of reactions in the curated anchor
  - anchor_fraction    : fraction in anchor

A module is the unit of phenotype evidence — much stronger than any
individual reaction. Bridging modules are the most paper-relevant
discovery output: novel reactions linked to confirmed disease biology
through specific shared-substrate/product/gene paths, with mechanistic
plausibility constraints (the novel reactions live within graph
distance of curated literature anchors).
"""
from __future__ import annotations
from collections import Counter, defaultdict
from typing import Optional

import networkx as nx


def _absval(post): return abs(post[2] - post[0])
def _signed(post): return post[2] - post[0]


def _reaction_neighbors(g, rxn_id) -> set[str]:
    """Substrates, products, and catalyzing genes of a reaction."""
    nbrs: set[str] = set()
    relevant = {"substrate", "product", "modifier"}
    for u, v, ed in g.edges(rxn_id, data=True):
        role = (ed.get("role") or ed.get("edge_type") or "").lower()
        if role in relevant or "catalys" in role:
            nbrs.add(v if u == rxn_id else u)
    for u, v, ed in g.in_edges(rxn_id, data=True):
        role = (ed.get("role") or ed.get("edge_type") or "").lower()
        if role in relevant or "catalys" in role:
            nbrs.add(u)
    return nbrs


def build_reaction_coadjacency(g, reactions: set[str]) -> nx.Graph:
    """Build a reaction-reaction graph where two reactions are linked
    iff they share at least one first-order neighbor (substrate /
    product / catalyzing gene).

    Returns an undirected ``networkx.Graph`` whose nodes are the input
    reactions. Density depends on how clustered the top-K is.
    """
    rxn_graph = nx.Graph()
    rxn_graph.add_nodes_from(reactions)

    # Index: shared_neighbor → reactions touching it
    nbr_to_rxns: dict[str, list[str]] = defaultdict(list)
    for rxn in reactions:
        for nbr in _reaction_neighbors(g, rxn):
            nbr_to_rxns[nbr].append(rxn)

    # Add an edge for each pair of reactions sharing a neighbor
    for nbr, rxn_list in nbr_to_rxns.items():
        # Skip ubiquitous neighbors (>>currency-like) that would
        # over-connect — anything touched by more than 30 of the top-K
        # reactions probably represents a hub like ATP/water and
        # shouldn't drive module structure.
        if len(rxn_list) > 30:
            continue
        for i, ra in enumerate(rxn_list):
            for rb in rxn_list[i + 1:]:
                rxn_graph.add_edge(ra, rb, shared=nbr)
    return rxn_graph


def _module_genes(g, reactions) -> set[str]:
    out: set[str] = set()
    for rxn in reactions:
        for u, v, ed in g.edges(rxn, data=True):
            role = (ed.get("role") or ed.get("edge_type") or "").lower()
            if not (role in ("modifier", "catalysis", "gene_reaction")
                    or "catalys" in role):
                continue
            other = v if u == rxn else u
            if g.nodes.get(other, {}).get("node_type") == "gene":
                gsym = (g.nodes[other].get("symbol")
                        or g.nodes[other].get("gene_symbol"))
                if gsym:
                    out.add(gsym)
        for u, v, ed in g.in_edges(rxn, data=True):
            role = (ed.get("role") or ed.get("edge_type") or "").lower()
            if not (role in ("modifier", "catalysis", "gene_reaction")
                    or "catalys" in role):
                continue
            if g.nodes.get(u, {}).get("node_type") == "gene":
                gsym = (g.nodes[u].get("symbol")
                        or g.nodes[u].get("gene_symbol"))
                if gsym:
                    out.add(gsym)
    return out


def find_perturbation_modules(
    g,
    post: dict,
    top_k_reactions: int = 200,
    min_module_size: int = 3,
    anchor: Optional[set[str]] = None,
) -> list[dict]:
    """Cluster the top-K most-perturbed reactions into connected
    subgraph modules ranked by aggregate perturbation mass.

    Parameters
    ----------
    g
        networkx Graph (mg.graph)
    post
        Reaction node_id → (p_down, p_normal, p_up) Laplacian posteriors.
    top_k_reactions
        How many top reactions (by |signed perturbation|) to consider
        as candidate module members. Larger = more inclusive but more
        dilution.
    min_module_size
        Modules with fewer reactions are dropped.
    anchor
        Optional set of curated disease-anchor reaction IDs. Modules
        get tagged with anchor_status: ``confirmed`` (all in anchor),
        ``novel`` (none), ``bridging`` (some). Bridging modules are
        the discovery candidates with mechanistic plausibility.

    Returns
    -------
    list[dict] — modules sorted by total_perturbation descending.
    """
    # Rank reactions by |signed perturbation|
    scored = []
    for nid, p in post.items():
        if g.nodes.get(nid, {}).get("node_type") != "reaction":
            continue
        scored.append((nid, _absval(p), _signed(p)))
    scored.sort(key=lambda x: -x[1])
    top = scored[:top_k_reactions]
    rxn_set = {nid for nid, _, _ in top}
    pert_map = {nid: m for nid, m, _ in top}
    sign_map = {nid: s for nid, _, s in top}

    rxn_graph = build_reaction_coadjacency(g, rxn_set)
    components = [
        c for c in nx.connected_components(rxn_graph)
        if len(c) >= min_module_size
    ]

    modules: list[dict] = []
    for comp in components:
        comp_rxns = list(comp)
        total_pert = sum(pert_map[r] for r in comp_rxns)
        mean_pert = total_pert / len(comp_rxns)

        # Direction: majority sign (UP / DOWN / mixed)
        n_up = sum(1 for r in comp_rxns if sign_map[r] > 0.01)
        n_down = sum(1 for r in comp_rxns if sign_map[r] < -0.01)
        if n_up > 0.7 * len(comp_rxns):
            direction = "UP"
        elif n_down > 0.7 * len(comp_rxns):
            direction = "DOWN"
        else:
            direction = "mixed"

        # Pathway aggregation
        path_counts: Counter = Counter()
        for r in comp_rxns:
            for p in (g.nodes[r].get("pathways") or []):
                path_counts[p] += 1

        # Catalyzing genes
        genes = _module_genes(g, comp_rxns)

        # Anchor coverage
        anchor_status = None
        n_in_anchor = 0
        anchor_fraction = 0.0
        if anchor is not None:
            n_in_anchor = sum(1 for r in comp_rxns if r in anchor)
            anchor_fraction = n_in_anchor / len(comp_rxns)
            if anchor_fraction >= 0.95:
                anchor_status = "confirmed"
            elif anchor_fraction <= 0.05:
                anchor_status = "novel"
            else:
                anchor_status = "bridging"

        modules.append({
            "reactions":          comp_rxns,
            "size":               len(comp_rxns),
            "total_perturbation": round(total_pert, 4),
            "mean_perturbation":  round(mean_pert, 4),
            "direction":          direction,
            "n_up":               n_up,
            "n_down":             n_down,
            "top_pathways":       path_counts.most_common(5),
            "genes":              sorted(genes)[:15],
            "anchor_status":      anchor_status,
            "n_in_anchor":        n_in_anchor,
            "anchor_fraction":    round(anchor_fraction, 2),
        })

    modules.sort(key=lambda m: -m["total_perturbation"])
    return modules


def module_label(g, module: dict, max_chars: int = 80) -> str:
    """Build a human-readable label for a module: top pathway name +
    representative genes."""
    if module["top_pathways"]:
        top_pid = module["top_pathways"][0][0]
        # Try to look up the pathway name from any reaction in the module
        # via its pathway_names attribute (Reactome populates this).
        # Fall back to the pathway id.
        # (For a richer label, query Reactome ContentService — too heavy.)
        pname = top_pid
    else:
        pname = "uncharacterized"
    genes = ",".join(module["genes"][:5])
    out = f"[{pname}]"
    if genes:
        out += f"  genes: {genes}"
    return out[:max_chars]
