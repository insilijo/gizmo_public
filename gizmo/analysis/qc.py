"""
Computational readiness assessment for the GIZMO metabolite graph.

Checks:
  1. Currency metabolite statistics (fraction of edge endpoints)
  2. Degree distribution (hub detection)
  3. Dead-end metabolites (appear only as substrate OR only as product)
  4. Orphan reactions (no substrates or no products in the graph)
  5. Disconnected components (weakly connected)
  6. Reaction annotation completeness (EC, gene, pathway coverage)
  7. Stoichiometric data availability (formula / charge present)
  8. Compartment coverage
  9. Clinical overlay coverage (disease / gene nodes)
 10. Metabolon compound coverage (if metabolon_name present)
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field
from typing import Optional

import networkx as nx

from gizmo.graph.network import GizmoGraph


@dataclass
class ReadinessReport:
    """Full computational readiness report for a GizmoGraph."""

    # Basic counts
    n_metabolites: int = 0
    n_reactions: int = 0
    n_diseases: int = 0
    n_genes: int = 0
    n_edges: int = 0

    # Currency
    n_currency: int = 0
    currency_edge_fraction: float = 0.0    # fraction of edges involving currency metabolites
    currency_ids: list[str] = field(default_factory=list)

    # Dead-ends
    n_dead_end_metabolites: int = 0        # only substrate or only product
    dead_end_ids: list[str] = field(default_factory=list)

    # Orphan reactions
    n_orphan_reactions: int = 0            # missing substrates or products
    orphan_ids: list[str] = field(default_factory=list)

    # Connectivity
    n_weakly_connected_components: int = 0
    largest_component_fraction: float = 0.0
    n_isolated_nodes: int = 0

    # Degree stats (metabolite reaction-degree)
    metabolite_degree_mean: float = 0.0
    metabolite_degree_max: int = 0
    metabolite_degree_p95: float = 0.0

    # Annotation completeness
    reactions_with_ec: int = 0
    reactions_with_ec_fraction: float = 0.0
    reactions_with_gene: int = 0
    reactions_with_gene_fraction: float = 0.0
    reactions_with_pathway: int = 0
    reactions_with_pathway_fraction: float = 0.0

    # Structural data
    metabolites_with_formula: int = 0
    metabolites_with_formula_fraction: float = 0.0
    metabolites_with_inchikey: int = 0
    metabolites_with_inchikey_fraction: float = 0.0
    metabolites_with_chebi: int = 0
    metabolites_with_chebi_fraction: float = 0.0

    # Compartments
    compartments: list[str] = field(default_factory=list)

    # Clinical overlay
    disease_reaction_edges: int = 0
    disease_metabolite_edges: int = 0
    disease_gene_edges: int = 0

    # Metabolon coverage
    metabolon_compounds_total: int = 0
    metabolon_with_chebi: int = 0
    metabolon_chebi_coverage: float = 0.0
    metabolon_with_hmdb: int = 0
    metabolon_hmdb_coverage: float = 0.0

    # HMDB coverage (all metabolites — needed for GeMMA/VMH interop)
    metabolites_with_hmdb: int = 0
    metabolites_with_hmdb_fraction: float = 0.0

    def print_summary(self) -> None:
        """Print a human-readable summary."""
        try:
            from rich.console import Console
            from rich.table import Table
            _rich_print(self)
        except ImportError:
            _plain_print(self)

    @property
    def is_fba_ready(self) -> bool:
        """Heuristic: graph is plausibly FBA-ready."""
        return (
            self.n_metabolites > 100
            and self.n_reactions > 50
            and self.reactions_with_ec_fraction > 0.3
            and self.metabolites_with_formula_fraction > 0.5
            and self.n_weakly_connected_components < 10
        )


def assess_readiness(mg: GizmoGraph) -> ReadinessReport:
    """Compute a full ReadinessReport for a GizmoGraph."""
    g = mg.graph
    r = ReadinessReport()

    met_ids = mg.metabolite_nodes()
    rxn_ids = mg.reaction_nodes()

    r.n_metabolites = len(met_ids)
    r.n_reactions = len(rxn_ids)
    r.n_diseases = len(mg.disease_nodes())
    r.n_genes = len(mg.gene_nodes())
    r.n_edges = g.number_of_edges()

    # --- Currency ---
    currency_set = {n for n in met_ids if g.nodes[n].get("is_currency", False)}
    r.n_currency = len(currency_set)
    r.currency_ids = list(currency_set)
    if r.n_edges > 0:
        currency_edges = sum(
            1
            for u, v in g.edges()
            if u in currency_set or v in currency_set
        )
        r.currency_edge_fraction = currency_edges / r.n_edges

    # --- Dead-end metabolites ---
    dead_ends: list[str] = []
    for nid in met_ids:
        if nid in currency_set:
            continue
        preds = [n for n in g.predecessors(nid) if g.nodes[n].get("node_type") == "reaction"]
        succs = [n for n in g.successors(nid) if g.nodes[n].get("node_type") == "reaction"]
        if (not preds) or (not succs):
            dead_ends.append(nid)
    r.n_dead_end_metabolites = len(dead_ends)
    r.dead_end_ids = dead_ends

    # --- Orphan reactions ---
    orphans: list[str] = []
    for rid in rxn_ids:
        subs = [n for n in g.predecessors(rid) if g.nodes[n].get("node_type") == "metabolite"]
        prods = [n for n in g.successors(rid) if g.nodes[n].get("node_type") == "metabolite"]
        if not subs or not prods:
            orphans.append(rid)
    r.n_orphan_reactions = len(orphans)
    r.orphan_ids = orphans

    # --- Connectivity (bipartite subgraph only) ---
    bipartite_nodes = set(met_ids) | set(rxn_ids)
    if bipartite_nodes:
        sub = g.subgraph(bipartite_nodes)
        undirected = sub.to_undirected()
        components = list(nx.connected_components(undirected))
        r.n_weakly_connected_components = len(components)
        r.largest_component_fraction = max(len(c) for c in components) / len(bipartite_nodes)
        r.n_isolated_nodes = sum(1 for c in components if len(c) == 1)

    # --- Degree distribution ---
    if met_ids:
        degrees = []
        for nid in met_ids:
            rxn_neighbors = {
                n for n in list(g.predecessors(nid)) + list(g.successors(nid))
                if g.nodes[n].get("node_type") == "reaction"
            }
            degrees.append(len(rxn_neighbors))
        degrees.sort()
        r.metabolite_degree_mean = statistics.mean(degrees) if degrees else 0.0
        r.metabolite_degree_max = max(degrees) if degrees else 0
        p95_idx = int(0.95 * len(degrees))
        r.metabolite_degree_p95 = float(degrees[p95_idx]) if degrees else 0.0

    # --- Reaction annotation completeness ---
    if rxn_ids:
        r.reactions_with_ec = sum(1 for rid in rxn_ids if g.nodes[rid].get("ec_numbers"))
        r.reactions_with_ec_fraction = r.reactions_with_ec / len(rxn_ids)
        r.reactions_with_gene = sum(1 for rid in rxn_ids if g.nodes[rid].get("gene_symbols"))
        r.reactions_with_gene_fraction = r.reactions_with_gene / len(rxn_ids)
        r.reactions_with_pathway = sum(1 for rid in rxn_ids if g.nodes[rid].get("pathways"))
        r.reactions_with_pathway_fraction = r.reactions_with_pathway / len(rxn_ids)

    # --- Metabolite structural data ---
    if met_ids:
        r.metabolites_with_formula = sum(1 for n in met_ids if g.nodes[n].get("formula"))
        r.metabolites_with_formula_fraction = r.metabolites_with_formula / len(met_ids)
        r.metabolites_with_inchikey = sum(1 for n in met_ids if g.nodes[n].get("inchikey"))
        r.metabolites_with_inchikey_fraction = r.metabolites_with_inchikey / len(met_ids)
        r.metabolites_with_chebi = sum(1 for n in met_ids if g.nodes[n].get("chebi_id"))
        r.metabolites_with_chebi_fraction = r.metabolites_with_chebi / len(met_ids)

    # --- Compartments ---
    compartments: set[str] = set()
    for nid in met_ids:
        if comp := g.nodes[nid].get("compartment"):
            compartments.add(comp)
    r.compartments = sorted(compartments)

    # --- Clinical overlay ---
    disease_ids = set(mg.disease_nodes())
    gene_ids = set(mg.gene_nodes())
    for u, v, data in g.edges(data=True):
        etype = data.get("edge_type", "")
        u_is_dis = u in disease_ids
        v_is_met = g.nodes.get(v, {}).get("node_type") == "metabolite"
        v_is_rxn = g.nodes.get(v, {}).get("node_type") == "reaction"
        v_is_gene = v in gene_ids
        if u_is_dis:
            if v_is_met:
                r.disease_metabolite_edges += 1
            elif v_is_rxn:
                r.disease_reaction_edges += 1
            elif v_is_gene:
                r.disease_gene_edges += 1

    # --- Metabolon coverage ---
    met_with_metabolon = [n for n in met_ids if g.nodes[n].get("metabolon_name")]
    r.metabolon_compounds_total = len(met_with_metabolon)
    r.metabolon_with_chebi = sum(1 for n in met_with_metabolon if g.nodes[n].get("chebi_id"))
    if r.metabolon_compounds_total > 0:
        r.metabolon_chebi_coverage = r.metabolon_with_chebi / r.metabolon_compounds_total
    r.metabolon_with_hmdb = sum(1 for n in met_with_metabolon if g.nodes[n].get("hmdb_id"))
    if r.metabolon_compounds_total > 0:
        r.metabolon_hmdb_coverage = r.metabolon_with_hmdb / r.metabolon_compounds_total

    # --- HMDB coverage (all metabolites, for GeMMA/VMH interop) ---
    if met_ids:
        r.metabolites_with_hmdb = sum(1 for n in met_ids if g.nodes[n].get("hmdb_id"))
        r.metabolites_with_hmdb_fraction = r.metabolites_with_hmdb / len(met_ids)

    return r


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

def _rich_print(r: ReadinessReport) -> None:
    from rich.console import Console
    from rich.table import Table

    console = Console()
    t = Table(title="GIZMO Computational Readiness Report", show_header=True)
    t.add_column("Check", style="bold")
    t.add_column("Value")
    t.add_column("Status")

    def _pct(f: float) -> str:
        return f"{f*100:.1f}%"

    def _ok(cond: bool) -> str:
        return "[green]OK[/green]" if cond else "[yellow]WARN[/yellow]"

    t.add_row("Metabolites", str(r.n_metabolites), _ok(r.n_metabolites > 0))
    t.add_row("Reactions", str(r.n_reactions), _ok(r.n_reactions > 0))
    t.add_row("Diseases", str(r.n_diseases), "")
    t.add_row("Genes", str(r.n_genes), "")
    t.add_row("Edges", str(r.n_edges), _ok(r.n_edges > 0))
    t.add_row("---", "", "")
    t.add_row("Currency metabolites", f"{r.n_currency}", "")
    t.add_row("Currency edge fraction", _pct(r.currency_edge_fraction),
              _ok(r.currency_edge_fraction < 0.5))
    t.add_row("Dead-end metabolites (non-currency)", str(r.n_dead_end_metabolites),
              _ok(r.n_dead_end_metabolites == 0))
    t.add_row("Orphan reactions", str(r.n_orphan_reactions),
              _ok(r.n_orphan_reactions == 0))
    t.add_row("---", "", "")
    t.add_row("Weakly connected components", str(r.n_weakly_connected_components),
              _ok(r.n_weakly_connected_components <= 5))
    t.add_row("Largest component coverage", _pct(r.largest_component_fraction),
              _ok(r.largest_component_fraction > 0.8))
    t.add_row("Isolated nodes", str(r.n_isolated_nodes), _ok(r.n_isolated_nodes == 0))
    t.add_row("---", "", "")
    t.add_row("Reactions with EC", _pct(r.reactions_with_ec_fraction),
              _ok(r.reactions_with_ec_fraction > 0.5))
    t.add_row("Reactions with gene annotation", _pct(r.reactions_with_gene_fraction),
              _ok(r.reactions_with_gene_fraction > 0.3))
    t.add_row("Reactions with pathway", _pct(r.reactions_with_pathway_fraction),
              _ok(r.reactions_with_pathway_fraction > 0.5))
    t.add_row("---", "", "")
    t.add_row("Metabolites with formula", _pct(r.metabolites_with_formula_fraction),
              _ok(r.metabolites_with_formula_fraction > 0.7))
    t.add_row("Metabolites with InChIKey", _pct(r.metabolites_with_inchikey_fraction),
              _ok(r.metabolites_with_inchikey_fraction > 0.7))
    t.add_row("Metabolites with ChEBI", _pct(r.metabolites_with_chebi_fraction),
              _ok(r.metabolites_with_chebi_fraction > 0.8))
    t.add_row("Metabolites with HMDB (VMH)", _pct(r.metabolites_with_hmdb_fraction),
              _ok(r.metabolites_with_hmdb_fraction > 0.5))
    t.add_row("Compartments", ", ".join(r.compartments) or "none", "")
    t.add_row("---", "", "")
    t.add_row("Metabolon compounds", str(r.metabolon_compounds_total), "")
    t.add_row("Metabolon → ChEBI coverage", _pct(r.metabolon_chebi_coverage),
              _ok(r.metabolon_chebi_coverage > 0.7))
    t.add_row("Metabolon → HMDB coverage", _pct(r.metabolon_hmdb_coverage),
              _ok(r.metabolon_hmdb_coverage > 0.5))
    t.add_row("---", "", "")
    t.add_row("FBA-ready heuristic", str(r.is_fba_ready),
              "[green]YES[/green]" if r.is_fba_ready else "[red]NO[/red]")

    console.print(t)


def _plain_print(r: ReadinessReport) -> None:
    lines = [
        "=== GIZMO Computational Readiness ===",
        f"Metabolites:       {r.n_metabolites}",
        f"Reactions:         {r.n_reactions}",
        f"Diseases:          {r.n_diseases}",
        f"Genes:             {r.n_genes}",
        f"Currency:          {r.n_currency} ({r.currency_edge_fraction*100:.1f}% of edges)",
        f"Dead-ends:         {r.n_dead_end_metabolites}",
        f"Orphan rxns:       {r.n_orphan_reactions}",
        f"Components:        {r.n_weakly_connected_components}",
        f"EC coverage:       {r.reactions_with_ec_fraction*100:.1f}%",
        f"Gene coverage:     {r.reactions_with_gene_fraction*100:.1f}%",
        f"ChEBI coverage:    {r.metabolites_with_chebi_fraction*100:.1f}%",
        f"HMDB coverage:     {r.metabolites_with_hmdb_fraction*100:.1f}%",
        f"Metabolon→ChEBI:   {r.metabolon_chebi_coverage*100:.1f}%",
        f"Metabolon→HMDB:    {r.metabolon_hmdb_coverage*100:.1f}%",
        f"FBA-ready:         {r.is_fba_ready}",
    ]
    print("\n".join(lines))
