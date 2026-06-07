"""
Core NetworkX-based metabolite graph.

Topology: directed graph with four node types:
  metabolite, reaction, disease, gene

Core bipartite edges:
  metabolite → reaction  (substrate/modifier)
  reaction   → metabolite (product)

Clinical overlay:
  disease → gene → reaction → metabolite
"""

from __future__ import annotations

from typing import Iterable

import networkx as nx

from gizmo.schema import (
    DiseaseEdge,
    DiseaseNode,
    DrugEdge,
    DrugNode,
    EdgeRole,
    GeneNode,
    MetaboliteNode,
    MicrobialEdge,
    MicrobialTaxonNode,
    PathwayEdge,
    PathwayNode,
    PhenotypeEdge,
    PhenotypeNode,
    ReactionEdge,
    ReactionNode,
    ToxEdge,
    VariantEdge,
    VariantNode,
)

_NODE_TYPES = {
    "metabolite", "reaction", "disease", "gene",
    "pathway", "phenotype", "drug", "variant",
    "microbial_taxon",
}


class GizmoGraph:
    """Directed graph integrating metabolite reactions with clinical annotations."""

    def __init__(self) -> None:
        self._g: nx.DiGraph = nx.DiGraph()

    # ------------------------------------------------------------------
    # Node management
    # ------------------------------------------------------------------

    def add_metabolite(self, node: MetaboliteNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_reaction(self, node: ReactionNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_disease(self, node: DiseaseNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_gene(self, node: GeneNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_pathway(self, node: PathwayNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_phenotype(self, node: PhenotypeNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_drug(self, node: DrugNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_variant(self, node: VariantNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_microbial_taxon(self, node: MicrobialTaxonNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_metabolites(self, nodes: Iterable[MetaboliteNode]) -> None:
        for n in nodes:
            self.add_metabolite(n)

    def add_reactions(self, nodes: Iterable[ReactionNode]) -> None:
        for n in nodes:
            self.add_reaction(n)

    def add_diseases(self, nodes: Iterable[DiseaseNode]) -> None:
        for n in nodes:
            self.add_disease(n)

    def add_genes(self, nodes: Iterable[GeneNode]) -> None:
        for n in nodes:
            self.add_gene(n)

    # ------------------------------------------------------------------
    # Edge management
    # ------------------------------------------------------------------

    def add_edge(self, edge: ReactionEdge) -> None:
        self._g.add_edge(
            edge.source,
            edge.target,
            role=edge.role.value,
            stoichiometry=edge.stoichiometry,
            compartment=edge.compartment,
        )

    def add_disease_edge(self, edge: DiseaseEdge) -> None:
        self._g.add_edge(
            edge.source,
            edge.target,
            edge_type=edge.edge_type.value,
            score=edge.score,
            evidence_count=edge.evidence_count,
            source_db=edge.source_db,
        )

    def add_edges(self, edges: Iterable[ReactionEdge]) -> None:
        for e in edges:
            self.add_edge(e)

    def add_disease_edges(self, edges: Iterable[DiseaseEdge]) -> None:
        for e in edges:
            self.add_disease_edge(e)

    def add_pathway_edge(self, edge: PathwayEdge) -> None:
        self._g.add_edge(
            edge.source, edge.target,
            edge_type=edge.edge_type,
            source_db=edge.source_db,
        )

    def add_phenotype_edge(self, edge: PhenotypeEdge) -> None:
        self._g.add_edge(
            edge.source, edge.target,
            edge_type=edge.edge_type,
            score=edge.score,
            source_db=edge.source_db,
        )

    def add_drug_edge(self, edge: DrugEdge) -> None:
        self._g.add_edge(
            edge.source, edge.target,
            edge_type=edge.edge_type,
            mechanism=edge.mechanism,
            max_phase=edge.max_phase,
            source_db=edge.source_db,
        )

    def add_variant_edge(self, edge: VariantEdge) -> None:
        self._g.add_edge(
            edge.source, edge.target,
            edge_type=edge.edge_type,
            consequence=edge.consequence,
            source_db=edge.source_db,
        )

    def add_tox_edge(self, edge: ToxEdge) -> None:
        self._g.add_edge(
            edge.source, edge.target,
            edge_type=edge.edge_type.value,
            effect_type=edge.effect_type,
            organism=edge.organism,
            p_value=edge.p_value,
            reference_count=edge.reference_count,
            direct_evidence=edge.direct_evidence,
            assay_endpoint=edge.assay_endpoint,
            assay_value=edge.assay_value,
            assay_units=edge.assay_units,
            source_db=edge.source_db,
        )

    def add_tox_edges(self, edges: Iterable[ToxEdge]) -> None:
        for e in edges:
            self.add_tox_edge(e)

    def add_microbial_edge(self, edge: MicrobialEdge) -> None:
        self._g.add_edge(
            edge.source,
            edge.target,
            edge_type=edge.edge_type,
            role=edge.role.value,
            vmh_metabolite_id=edge.vmh_metabolite_id,
            reaction_ids=edge.reaction_ids,
            abundance_weight=edge.abundance_weight,
            sample_weights=edge.sample_weights,
            source_db=edge.source_db,
        )

    # ------------------------------------------------------------------
    # Typed node accessors
    # ------------------------------------------------------------------

    @property
    def graph(self) -> nx.DiGraph:
        return self._g

    def _nodes_of_type(self, node_type: str) -> list[str]:
        return [n for n, d in self._g.nodes(data=True) if d.get("node_type") == node_type]

    def metabolite_nodes(self) -> list[str]:
        return self._nodes_of_type("metabolite")

    def reaction_nodes(self) -> list[str]:
        return self._nodes_of_type("reaction")

    def disease_nodes(self) -> list[str]:
        return self._nodes_of_type("disease")

    def gene_nodes(self) -> list[str]:
        return self._nodes_of_type("gene")

    def pathway_nodes(self) -> list[str]:
        return self._nodes_of_type("pathway")

    def phenotype_nodes(self) -> list[str]:
        return self._nodes_of_type("phenotype")

    def drug_nodes(self) -> list[str]:
        return self._nodes_of_type("drug")

    def microbial_taxon_nodes(self) -> list[str]:
        return self._nodes_of_type("microbial_taxon")

    def variant_nodes(self) -> list[str]:
        return self._nodes_of_type("variant")

    def currency_nodes(self) -> list[str]:
        """Metabolite nodes flagged as currency metabolites."""
        return [
            n
            for n, d in self._g.nodes(data=True)
            if d.get("node_type") == "metabolite" and d.get("is_currency", False)
        ]

    def reviewed_nodes(self, node_type: str | None = None) -> list[str]:
        """
        Nodes flagged as manually_reviewed=True.

        Parameters
        ----------
        node_type : optional filter — "metabolite", "reaction", "disease", "gene"
        """
        return [
            n
            for n, d in self._g.nodes(data=True)
            if d.get("manually_reviewed", False)
            and (node_type is None or d.get("node_type") == node_type)
        ]

    def reviewed_edges(self) -> list[tuple[str, str]]:
        """Edges flagged as manually_reviewed=True."""
        return [
            (u, v)
            for u, v, d in self._g.edges(data=True)
            if d.get("manually_reviewed", False)
        ]

    # ------------------------------------------------------------------
    # Traversal helpers
    # ------------------------------------------------------------------

    def neighbors_of_metabolite(self, node_id: str) -> list[str]:
        """Reaction nodes that involve this metabolite."""
        pred = [n for n in self._g.predecessors(node_id) if self._g.nodes[n].get("node_type") == "reaction"]
        succ = [n for n in self._g.successors(node_id) if self._g.nodes[n].get("node_type") == "reaction"]
        return list(set(pred + succ))

    def diseases_for_metabolite(self, node_id: str) -> list[str]:
        """Disease nodes directly linked to a metabolite (biomarker edges)."""
        return [
            n
            for n in self._g.predecessors(node_id)
            if self._g.nodes[n].get("node_type") == "disease"
        ]

    def diseases_for_reaction(self, node_id: str) -> list[str]:
        """Disease nodes linked to a reaction via pathway association."""
        return [
            n
            for n in self._g.predecessors(node_id)
            if self._g.nodes[n].get("node_type") == "disease"
        ]

    def promote_pathway_nodes(self, species: str = "Homo sapiens") -> int:
        """
        Promote the ``pathways`` stID lists stored on reaction nodes into
        first-class PathwayNode records, wiring pathway → reaction edges.

        This can be called after a Reactome build to make pathways queryable
        as graph nodes rather than plain attribute lists.

        Returns the number of pathway nodes created.
        """
        g = self._g
        pathway_names: dict[str, str] = {}

        # First pass: collect all pathway stIDs and any name hints
        for nid, attrs in g.nodes(data=True):
            if attrs.get("node_type") != "reaction":
                continue
            for stid in (attrs.get("pathways") or []):
                if stid not in pathway_names:
                    pathway_names[stid] = stid  # name unknown at this stage

        # Second pass: create nodes + edges
        created = 0
        for stid, name in pathway_names.items():
            node_id = f"reactome:{stid}" if not stid.startswith("reactome:") else stid
            if node_id not in g:
                g.add_node(node_id,
                    node_type="pathway",
                    node_id=node_id,
                    reactome_id=stid,
                    name=name,
                    species=species,
                    parent_pathways=[],
                    level=0,
                )
                created += 1

            # Wire pathway → reaction edges for all reactions in this pathway
            for nid, attrs in g.nodes(data=True):
                if attrs.get("node_type") == "reaction" and stid in (attrs.get("pathways") or []):
                    if not g.has_edge(node_id, nid):
                        g.add_edge(node_id, nid, edge_type="pathway_reaction", source_db="reactome")

        return created

    def filter_by_species(self, species: str) -> "GizmoGraph":
        """
        Return a new GizmoGraph retaining only nodes whose ``species``
        attribute matches ``species`` (case-insensitive prefix match).

        Nodes without a ``species`` attribute (metabolites, diseases,
        phenotypes, drugs, variants) are always kept.  Genes and reactions
        specific to another species are removed.

        Parameters
        ----------
        species : e.g. "Homo sapiens", "Mus musculus", "9606"

        Returns
        -------
        A new GizmoGraph filtered to the requested species.
        """
        species_lower = species.lower()
        keep: set[str] = set()
        for nid, attrs in self._g.nodes(data=True):
            node_species = (attrs.get("species") or "").lower()
            if not node_species or node_species.startswith(species_lower):
                keep.add(nid)

        sub_nx = self._g.subgraph(keep).copy()
        mg = GizmoGraph()
        mg._g = sub_nx
        return mg

    def metabolite_subgraph(self, node_ids: Iterable[str]) -> GizmoGraph:
        """
        Induced subgraph over a set of metabolite IDs, including reaction nodes
        that connect only those metabolites. Disease/gene nodes are excluded.
        """
        node_ids = set(node_ids)
        reaction_ids: set[str] = set()
        for rxn in self.reaction_nodes():
            rxn_mets = {
                n
                for n in list(self._g.predecessors(rxn)) + list(self._g.successors(rxn))
                if self._g.nodes[n].get("node_type") == "metabolite"
            }
            if rxn_mets and rxn_mets.issubset(node_ids):
                reaction_ids.add(rxn)

        sub_nx = self._g.subgraph(node_ids | reaction_ids).copy()
        mg = GizmoGraph()
        mg._g = sub_nx
        return mg

    # ------------------------------------------------------------------
    # Graph projection / collapse
    # ------------------------------------------------------------------

    def collapse_graph(
        self,
        keep_types: set[str],
        *,
        directed: bool = True,
        via_attr: str = "via",
    ) -> "GizmoGraph":
        """
        Return a new GizmoGraph retaining only nodes of ``keep_types``.

        Nodes whose ``node_type`` is **not** in ``keep_types`` are removed;
        any path that ran through them is short-circuited with a direct edge
        between the nearest kept ancestors and descendants.

        Typical projections
        -------------------
        ``{"metabolite"}``
            Metabolite co-reaction network — two metabolites are connected
            whenever they share a reaction (substrate→product).

        ``{"metabolite", "gene"}``
            Omics network — metabolites and their catalysing genes, reactions
            collapsed out.

        ``{"metabolite", "disease"}``
            Clinical biomarker network — metabolites connected to the diseases
            for which they are biomarkers (via disease→gene→reaction or
            direct disease→metabolite edges).

        ``{"metabolite", "reaction", "gene", "disease"}``
            Full graph — no collapse at all.

        Parameters
        ----------
        keep_types : node_type strings to retain in the output graph
        directed   : if True, preserve edge direction (predecessor → successor);
                     if False, add edges in both directions (co-occurrence)
        via_attr   : edge attribute name used to record which node type(s)
                     were collapsed onto this edge (list of type strings)

        Returns
        -------
        A new ``GizmoGraph`` whose nodes are exactly those in
        ``keep_types``, with edges wherever a path existed through removed
        intermediate nodes.
        """
        g = self._g
        keep_nodes: set[str] = {
            n for n, d in g.nodes(data=True) if d.get("node_type") in keep_types
        }

        # For each kept node, BFS forward (and optionally backward) through
        # *removed* nodes to find reachable kept nodes.  Collect collapsed
        # edge metadata along the way.
        #
        # collapsed_edges: (src, dst) → set of intermediate node_types
        collapsed_edges: dict[tuple[str, str], set[str]] = {}

        def _bfs_forward(start: str) -> None:
            """From `start` (kept), walk only through removed nodes to find kept successors."""
            queue = list(g.successors(start))
            visited: set[str] = set()
            while queue:
                nxt = queue.pop()
                if nxt in visited:
                    continue
                visited.add(nxt)
                ntype = g.nodes[nxt].get("node_type", "")
                if nxt in keep_nodes:
                    key = (start, nxt)
                    collapsed_edges.setdefault(key, set()).add(ntype)
                    # Don't traverse further from a kept node
                else:
                    # Removed node — carry its type and keep walking
                    collapsed_edges  # just to reference it in closure
                    for nb in g.successors(nxt):
                        if nb not in visited:
                            queue.append(nb)
                    # Also note the intermediate type for any edge we eventually create
                    # We store intermediate types per (start, *) but we won't know the
                    # final dst yet — handled below via the visited set carrying types.

        # Simpler flat BFS that accumulates intermediate types correctly:
        def _find_reachable(start: str) -> dict[str, set[str]]:
            """BFS from start through removed nodes; returns {dst_kept: {via_types}}."""
            reachable: dict[str, set[str]] = {}
            # queue entries: (node, set_of_via_types_so_far)
            queue: list[tuple[str, frozenset[str]]] = []
            for nb in g.successors(start):
                ntype = g.nodes[nb].get("node_type", "")
                if nb in keep_nodes:
                    reachable.setdefault(nb, set()).add(ntype)
                else:
                    queue.append((nb, frozenset({ntype})))

            visited: set[str] = set()
            while queue:
                cur, via = queue.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                for nb in g.successors(cur):
                    ntype = g.nodes[nb].get("node_type", "")
                    if nb in keep_nodes:
                        reachable.setdefault(nb, set()).update(via | {ntype})
                    elif nb not in visited:
                        queue.append((nb, via | frozenset({ntype})))
            return reachable

        new_g = nx.DiGraph()

        # Copy kept nodes with all their attributes
        for n in keep_nodes:
            new_g.add_node(n, **g.nodes[n])

        # Direct edges between kept nodes
        for u, v, edata in g.edges(data=True):
            if u in keep_nodes and v in keep_nodes:
                new_g.add_edge(u, v, **edata, **{via_attr: []})

        # Collapsed edges
        for start in keep_nodes:
            reachable = _find_reachable(start)
            for dst, via_types in reachable.items():
                if not new_g.has_edge(start, dst):
                    new_g.add_edge(start, dst, **{via_attr: sorted(via_types)})
                else:
                    # Merge via types into existing edge
                    existing = new_g[start][dst].get(via_attr, [])
                    new_g[start][dst][via_attr] = sorted(set(existing) | via_types)

                if not directed:
                    # Also add reverse direction for co-occurrence semantics
                    if not new_g.has_edge(dst, start):
                        new_g.add_edge(dst, start, **{via_attr: sorted(via_types)})

        mg = GizmoGraph()
        mg._g = new_g
        return mg

    # ------------------------------------------------------------------
    # In-place mutation helpers
    # ------------------------------------------------------------------

    def flag_currency(self, node_ids: Iterable[str]) -> int:
        """
        Mark nodes as currency metabolites in-place.
        Returns count of nodes flagged.
        """
        count = 0
        for nid in node_ids:
            if nid in self._g and self._g.nodes[nid].get("node_type") == "metabolite":
                self._g.nodes[nid]["is_currency"] = True
                count += 1
        return count

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------

    def summary(self) -> dict:
        return {
            "metabolites": len(self.metabolite_nodes()),
            "reactions":   len(self.reaction_nodes()),
            "pathways":    len(self.pathway_nodes()),
            "diseases":    len(self.disease_nodes()),
            "genes":       len(self.gene_nodes()),
            "phenotypes":  len(self.phenotype_nodes()),
            "drugs":       len(self.drug_nodes()),
            "variants":    len(self.variant_nodes()),
            "currency_flagged":   len(self.currency_nodes()),
            "manually_reviewed":  len(self.reviewed_nodes()),
            "reviewed_edges":     len(self.reviewed_edges()),
            "edges": self._g.number_of_edges(),
        }

    def __repr__(self) -> str:
        s = self.summary()
        parts = [
            f"{s['metabolites']} metabolites [{s['currency_flagged']} currency]",
            f"{s['reactions']} reactions",
            f"{s['diseases']} diseases",
            f"{s['genes']} genes",
        ]
        if s["phenotypes"]:
            parts.append(f"{s['phenotypes']} phenotypes")
        if s["drugs"]:
            parts.append(f"{s['drugs']} drugs")
        if s["variants"]:
            parts.append(f"{s['variants']} variants")
        parts.append(f"{s['edges']} edges")
        return f"GizmoGraph({', '.join(parts)})"
