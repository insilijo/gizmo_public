"""
Build gene nodes from `gene_symbols` reaction-attribute lists.

Reactome's reaction nodes carry `gene_symbols: list[str]` recording the
catalysing genes, but the loader doesn't create separate gene nodes.
GeneMapper requires `node_type='gene'` nodes to resolve transcriptomic /
proteomic feature inputs to graph nodes.

This enricher:
  1. Scans every reaction node's gene_symbols
  2. Creates one GeneNode per unique symbol (canonical_id "symbol:{SYM}")
  3. Adds catalyzes-edges  gene → reaction  (role="catalysis")
  4. Re-uses the BP gene-edge couplings (positive direction)

Lightweight: no external downloads; uses the graph as-is.
"""
from __future__ import annotations

import logging
from typing import Optional

log = logging.getLogger(__name__)


def enrich_graph_genes_from_reactions(
    mg, *,
    species: str = "Homo sapiens",
    skip_existing: bool = True,
) -> dict[str, int]:
    """Add one GeneNode per unique reaction.gene_symbols entry, plus
    catalysis edges gene → reaction.

    Returns counts: {"new_genes": ..., "new_edges": ..., "skipped": ...}.
    """
    from gizmo.schema import GeneNode

    g = mg.graph

    # Pass 1 — collect all (symbol → reactions) pairs
    sym_to_rxns: dict[str, list[str]] = {}
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "reaction":
            continue
        syms = attrs.get("gene_symbols") or []
        for s in syms:
            if not s:
                continue
            s = s.strip()
            if not s:
                continue
            sym_to_rxns.setdefault(s, []).append(nid)

    log.info("enrich_graph_genes: %d unique gene symbols across %d reactions",
             len(sym_to_rxns), sum(len(v) for v in sym_to_rxns.values()))

    new_genes = 0
    new_edges = 0
    skipped = 0
    for sym, rxns in sym_to_rxns.items():
        gene_id = f"symbol:{sym}"
        if skip_existing and gene_id in g:
            # If the node exists (e.g. loaded from Orphanet/Open Targets),
            # just add catalysis edges to all its reactions if missing.
            for rxn_id in rxns:
                if not g.has_edge(gene_id, rxn_id):
                    g.add_edge(gene_id, rxn_id, role="gene_reaction",
                               edge_type="gene_reaction")
                    new_edges += 1
            skipped += 1
            continue

        # Materialise the GeneNode dataclass via add_gene to keep schema
        # validation, then add catalysis edges.
        gn = GeneNode(node_id=gene_id, symbol=sym, species=species)
        mg.add_gene(gn)
        new_genes += 1
        for rxn_id in rxns:
            if not g.has_edge(gene_id, rxn_id):
                g.add_edge(gene_id, rxn_id, role="catalysis",
                           edge_type="catalysis")
                new_edges += 1

    log.info("enrich_graph_genes: added %d gene nodes, %d catalysis edges "
             "(skipped %d existing gene nodes)",
             new_genes, new_edges, skipped)
    return {"new_genes": new_genes, "new_edges": new_edges, "skipped": skipped}
