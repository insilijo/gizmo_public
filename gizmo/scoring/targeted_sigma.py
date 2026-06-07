"""Targeted-σ evaluation: test whether REAL Laplacian scores on a curated
subset of reactions exceed a permutation null on that same subset.

Generalizes the pathway-rollup metric to use the full GrAndMA-curated
graph: not just metabolic pathway membership, but disease-pathway
anchors, druggability flags, biomarker associations, and causal-chain
neighborhoods. Each is a different "subset of reactions you'd expect
the model to recover for this cohort", and each yields its own σ.

Why this matters: the pathway-rollup top-K metric assumes signal is
concentrated in a few pathways. For systemic disease (e.g., COVID
plasma) the signal is distributed across many pathway boundaries, and
top-K fails. Disease-anchored σ replaces "find the one big pathway"
with "find the literature-curated reactions for THIS disease" — works
the same way regardless of how the disease distributes through the graph.

Anchor types
------------
- pathway       : reactions whose ``pathways`` attribute contains a
                   given Reactome pathway stID (e.g., R-HSA-9694516
                   "SARS-CoV-2 Infection")
- disease_node  : reactions reachable from a disease node within
                   ``max_hops`` via any DiseaseEdgeType
- gene_set      : reactions catalyzed by any gene in the supplied set
                   (e.g., DGIdb-druggable, OncoKB, COSMIC)
- metabolite_set: reactions whose substrate or product is in the set
                   (e.g., Orphanet biomarkers, FDA-approved markers)

All anchors yield a ``set[reaction_node_id]`` that ``targeted_sigma``
consumes uniformly.
"""
from __future__ import annotations
import statistics
from collections import deque
from typing import Optional


def _absval(post) -> float:
    return abs(post[2] - post[0])


def reactions_in_pathway(g, pathway_id: str) -> set[str]:
    """Reactions whose ``pathways`` attribute includes ``pathway_id``."""
    out: set[str] = set()
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "reaction":
            continue
        paths = attrs.get("pathways") or set()
        if pathway_id in paths:
            out.add(nid)
    return out


def reactions_in_pathways(g, pathway_ids: set[str]) -> set[str]:
    """Union of reactions in any of ``pathway_ids``."""
    pids = set(pathway_ids)
    out: set[str] = set()
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "reaction":
            continue
        paths = set(attrs.get("pathways") or set())
        if paths & pids:
            out.add(nid)
    return out


def reactions_catalyzed_by(g, gene_set: set[str]) -> set[str]:
    """Reactions catalyzed by any gene node in ``gene_set``.

    Walks ``modifier`` / ``catalysis`` / ``gene_reaction`` edges. The
    ``gene_set`` may contain either gene-node ids or gene symbols; we
    match against either.
    """
    out: set[str] = set()
    for u, v, ed in g.edges(data=True):
        role = (ed.get("role") or ed.get("edge_type") or "").lower()
        if not (role in ("modifier", "catalysis", "gene_reaction",
                          "gene_associated") or "catalys" in role):
            continue
        u_attrs = g.nodes.get(u, {})
        v_attrs = g.nodes.get(v, {})
        u_is_gene = u_attrs.get("node_type") == "gene"
        v_is_gene = v_attrs.get("node_type") == "gene"
        if u_is_gene and v_attrs.get("node_type") == "reaction":
            gene_nid, rxn_nid = u, v
        elif v_is_gene and u_attrs.get("node_type") == "reaction":
            gene_nid, rxn_nid = v, u
        else:
            continue
        # Match by gene_nid or gene symbol
        if gene_nid in gene_set:
            out.add(rxn_nid)
            continue
        gsym = g.nodes[gene_nid].get("symbol") or g.nodes[gene_nid].get("gene_symbol")
        if gsym and gsym in gene_set:
            out.add(rxn_nid)
    return out


def reactions_with_metabolite_in(g, metab_node_set: set[str]) -> set[str]:
    """Reactions whose substrate or product is in ``metab_node_set``."""
    out: set[str] = set()
    for u, v, ed in g.edges(data=True):
        role = (ed.get("role") or ed.get("edge_type") or "").lower()
        if role not in ("substrate", "product"):
            continue
        u_t = g.nodes.get(u, {}).get("node_type")
        v_t = g.nodes.get(v, {}).get("node_type")
        if u_t == "metabolite" and v_t == "reaction":
            m_nid, r_nid = u, v
        elif v_t == "metabolite" and u_t == "reaction":
            m_nid, r_nid = v, u
        else:
            continue
        if m_nid in metab_node_set:
            out.add(r_nid)
    return out


def reactions_from_disease_node(g, disease_node_id: str,
                                  max_hops: int = 2,
                                  edge_types: Optional[set[str]] = None) -> set[str]:
    """BFS from a disease node up to ``max_hops``; collect all reaction
    nodes encountered. Optionally restrict to edges in ``edge_types``
    (e.g., {'gene_associated', 'biomarker', 'pathway_associated',
    'causal'}).
    """
    if disease_node_id not in g:
        return set()
    seen = {disease_node_id}
    out: set[str] = set()
    frontier = deque([(disease_node_id, 0)])
    while frontier:
        nid, depth = frontier.popleft()
        attrs = g.nodes[nid]
        if attrs.get("node_type") == "reaction":
            out.add(nid)
        if depth >= max_hops:
            continue
        for _, nbr, ed in g.edges(nid, data=True):
            if nbr in seen:
                continue
            if edge_types is not None:
                role = (ed.get("role") or ed.get("edge_type") or "").lower()
                if role not in edge_types:
                    continue
            seen.add(nbr)
            frontier.append((nbr, depth + 1))
        # Also walk in-edges (graph is mostly directed in places)
        for nbr, _, ed in g.in_edges(nid, data=True):
            if nbr in seen:
                continue
            if edge_types is not None:
                role = (ed.get("role") or ed.get("edge_type") or "").lower()
                if role not in edge_types:
                    continue
            seen.add(nbr)
            frontier.append((nbr, depth + 1))
    return out


def targeted_sigma(post_real: dict[str, tuple],
                   perm_posts: list[dict[str, tuple]],
                   target_reactions: set[str]) -> dict:
    """Compute σ above perm null for the targeted reaction subset.

    Parameters
    ----------
    post_real
        Real-label Laplacian posteriors {reaction_id → (p_down, p_normal, p_up)}.
    perm_posts
        List of N permutation-label posterior dicts.
    target_reactions
        Subset of reaction node IDs that define the curated anchor.

    Returns
    -------
    dict with keys: n_target, real, perm_mean, perm_sd, sigma.
    """
    target = target_reactions & set(post_real.keys())
    if not target:
        return {"n_target": 0, "real": 0.0, "perm_mean": 0.0,
                "perm_sd": 0.0, "sigma": float("nan")}

    real_score = statistics.mean(_absval(post_real[r]) for r in target)
    perm_scores = []
    for pp in perm_posts:
        avail = target & set(pp.keys())
        if avail:
            perm_scores.append(statistics.mean(_absval(pp[r]) for r in avail))
    if not perm_scores:
        return {"n_target": len(target), "real": real_score,
                "perm_mean": 0.0, "perm_sd": 0.0, "sigma": float("nan")}

    mu = statistics.mean(perm_scores)
    sd = statistics.stdev(perm_scores) if len(perm_scores) > 1 else 0.0
    sigma = (real_score - mu) / sd if sd > 0 else float("nan")
    return {
        "n_target": len(target),
        "real": real_score,
        "perm_mean": mu,
        "perm_sd": sd,
        "sigma": sigma,
    }


# ---------------------------------------------------------------------------
# Convenience: cohort → anchor-builder map for the paper benchmark battery.
# ---------------------------------------------------------------------------

# Cohort anchors as multi-dimensional curated sets. Each cohort gets a
# composite anchor built from the union of:
#   - pathway IDs        (Reactome disease/biology pathways)
#   - gene_symbols       (literature-cited differential or known target genes)
#   - metabolite_names   (literature-cited differential or biomarker metabs)
#
# A reaction is in the cohort anchor if it's in any of the named pathways,
# OR catalyzed by any of the named genes, OR has any of the named
# metabolites as substrate or product. This generalizes to systemic
# disease where Reactome's pathway boundaries don't carve up the biology
# the way the curators or the published differential lists do.
COHORT_ANCHORS: dict[str, dict] = {
    "IDH": {
        "pathways": [
            "R-HSA-2978092",   # Defective IDH1 → 2-HG accumulation
            "R-HSA-71403",     # Citric acid cycle (TCA)
        ],
        "genes": ["IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "L2HGDH", "D2HGDH"],
        "metabolites": ["2-hydroxyglutarate", "isocitrate",
                          "alpha-ketoglutarate", "2-oxoglutarate"],
    },
    "RA": {
        "pathways": [
            "R-HSA-1280215",   # Cytokine signaling
            "R-HSA-446652",    # IL-1 family signaling
            "R-HSA-168249",    # Innate immune response
        ],
        "genes": ["TNF", "IL1B", "IL6", "PTGS2", "ALOX5", "NOS2",
                    "TLR2", "TLR4", "NLRP3"],
        "metabolites": ["citrullin", "histamine", "prostaglandin",
                          "leukotriene"],
    },
    "Crohn": {
        "pathways": [
            "R-HSA-168249",    # Innate immune system
            "R-HSA-1280215",   # Cytokine signaling
            "R-HSA-73847",     # Purine metabolism (thiopurine drugs)
        ],
        "genes": ["NOD2", "ATG16L1", "IL23R", "IRGM", "CARD9", "TNFAIP3"],
        "metabolites": ["nicotinic acid", "xanthine", "taurine",
                          "phenylalanine"],
    },
    "Statin": {
        "pathways": [
            "R-HSA-191273",    # Cholesterol biosynthesis (HMGCR target)
            "R-HSA-2426168",   # SREBF gene expression
            "R-HSA-1655829",   # SREBP cholesterol regulation
        ],
        "genes": ["HMGCR", "HMGCS1", "SQLE", "DHCR24", "DHCR7",
                    "LDLR", "SREBF1", "SREBF2", "ACAT2", "MVK"],
        "metabolites": ["mevalonate", "lanosterol", "cholesterol",
                          "squalene", "farnesyl-pp"],
    },
    "Su_COVID": {
        "pathways": [
            "R-HSA-9694516",   # SARS-CoV-2 Infection (Reactome disease)
            "R-HSA-913531",    # Interferon signaling
            "R-HSA-1280215",   # Cytokine signaling
            "R-HSA-71240",     # Tryptophan / kynurenine (Su's claim)
            "R-HSA-191273",    # Cholesterol biosynth (Su's steroid claim)
            "R-HSA-159418",    # Bile acid metabolism (Su's claim)
            "R-HSA-211859",    # Biological oxidations (sterol/bile)
            "R-HSA-6783783",   # IL-10 signaling
            "R-HSA-1474244",   # Extracellular matrix (clotting/fibrin)
        ],
        # From Su's S1.7 + paper text: top trans + curated COVID genes
        "genes": ["IFNG", "IL10RA", "IL6", "TNF", "CCL7", "CCL20",
                    "CXCL10", "IDO1", "TDO2", "ACE2", "TMPRSS2",
                    "C3", "C5", "F2", "F10", "VEGFA", "HGF",
                    "TNFSF11", "TNFRSF10B", "KRT19", "NADK"],
        # From our diagnostic top-30: kynurenine, steroids, bile, sterols
        "metabolites": ["kynurenine", "tryptophan", "n-acetyltryptophan",
                          "pregnanolone", "allopregnanolone",
                          "beta-sitosterol", "cholesterol",
                          "taurodeoxycholic acid", "deoxycholate",
                          "urobilin", "bilirubin", "5alpha-pregnan-diol"],
    },
    "Erawijantari": {
        "pathways": [
            "R-HSA-159418",    # Bile acid metabolism
            "R-HSA-8963743",   # Digestion & absorption
            "R-HSA-71240",     # Tryptophan metabolism (microbiome-derived)
        ],
        "genes": ["BAAT", "CYP7A1", "CYP27A1", "SLC10A2"],
        "metabolites": ["taurine", "glycine", "cholate", "deoxycholate",
                          "phenylalanine", "tyrosine"],
    },
}


def chain_back_to_anchor(g, reaction_id: str, anchor_set: set[str],
                          max_hops: int = 4) -> list[str] | None:
    """BFS from a novel reaction outward; return the shortest node-path
    to any reaction in ``anchor_set`` (within ``max_hops`` steps), or
    None if no such path exists.

    The path is a list of node IDs alternating reaction → metabolite/gene
    → reaction → ..., terminating at an anchor reaction. Tells the user
    *why* the model flagged a novel reaction as disease-relevant: through
    what curated graph chain does it connect to known biology?

    Walks substrate / product / modifier edges; does not cross
    edge-type boundaries that aren't metabolic.
    """
    if reaction_id in anchor_set:
        return [reaction_id]
    relevant_roles = {"substrate", "product", "modifier", "catalysis",
                       "gene_reaction"}
    seen = {reaction_id}
    # frontier entries: (node, path_so_far)
    frontier = [(reaction_id, [reaction_id])]
    for _ in range(max_hops):
        next_frontier = []
        for nid, path in frontier:
            for u, v, ed in g.edges(nid, data=True):
                role = (ed.get("role") or ed.get("edge_type") or "").lower()
                if role not in relevant_roles and "catalys" not in role:
                    continue
                nbr = v if u == nid else u
                if nbr in seen:
                    continue
                seen.add(nbr)
                new_path = path + [nbr]
                if (g.nodes.get(nbr, {}).get("node_type") == "reaction"
                        and nbr in anchor_set):
                    return new_path
                next_frontier.append((nbr, new_path))
            for u, v, ed in g.in_edges(nid, data=True):
                role = (ed.get("role") or ed.get("edge_type") or "").lower()
                if role not in relevant_roles and "catalys" not in role:
                    continue
                if u in seen:
                    continue
                seen.add(u)
                new_path = path + [u]
                if (g.nodes.get(u, {}).get("node_type") == "reaction"
                        and u in anchor_set):
                    return new_path
                next_frontier.append((u, new_path))
        frontier = next_frontier
        if not frontier:
            break
    return None


def chain_human_readable(g, path: list[str]) -> list[dict]:
    """Decorate a chain-back path with node names + types for display."""
    out = []
    for nid in path:
        attrs = g.nodes.get(nid, {})
        nt = attrs.get("node_type", "?")
        if nt == "reaction":
            label = attrs.get("name") or attrs.get("display_name") or nid
        elif nt == "metabolite":
            label = attrs.get("name") or attrs.get("chebi_id") or nid
        elif nt == "gene":
            label = attrs.get("symbol") or attrs.get("gene_symbol") or nid
        else:
            label = nid
        out.append({"node_id": nid, "node_type": nt, "label": label})
    return out


def composite_cohort_anchor(g, cohort_key: str) -> set[str]:
    """Build the union anchor for a cohort by traversing pathway,
    gene-catalysis, and metabolite-substrate dimensions of the graph."""
    spec = COHORT_ANCHORS.get(cohort_key)
    if not spec:
        return set()
    anchor: set[str] = set()
    if spec.get("pathways"):
        anchor |= reactions_in_pathways(g, spec["pathways"])
    if spec.get("genes"):
        anchor |= reactions_catalyzed_by(g, set(spec["genes"]))
    if spec.get("metabolites"):
        # Resolve metabolite names to graph node IDs via the mapper
        from gizmo.evidence.mappers import MetaboliteMapper
        # Build a lightweight wrapper for the bare graph
        class _Mg:
            def __init__(self, graph):
                self.graph = graph
        mm = MetaboliteMapper(_Mg(g))
        metab_node_set: set[str] = set()
        for name in spec["metabolites"]:
            for nid, _ in mm.map_all_compartments(name):
                metab_node_set.add(nid)
        anchor |= reactions_with_metabolite_in(g, metab_node_set)
    return anchor


# Backward-compat alias for the old pathway-only schema; existing scripts
# that import COHORT_PATHWAY_ANCHORS still work.
COHORT_PATHWAY_ANCHORS: dict[str, list[str]] = {
    k: v["pathways"] for k, v in COHORT_ANCHORS.items() if v.get("pathways")
}
