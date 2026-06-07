"""
GeMMA → GIZMO integration: attach a SlimmedCommunity to a GizmoGraph.

This module bridges GeMMA's genome-scale metabolic community model with the
GIZMO knowledge graph, adding microbial taxon nodes and abundance-weighted
edges to metabolite nodes.  After calling :func:`attach_community`, the graph
gains a new layer that answers questions like:

  "Which gut bacteria produce/consume this metabolite, at what abundance?"

Node type added:  ``microbial_taxon``   (``MicrobialTaxonNode``)
Edge type added:  ``microbial_metabolite``  (``MicrobialEdge``,
                  taxon → GIZMO metabolite node)

ID translation
--------------
GeMMA uses VMH abbreviations (e.g. ``glc_D``) as metabolite identifiers.
GIZMO uses ChEBI / HMDB node IDs.  The bridge requires that the graph has been
enriched with VMH IDs first:

    >>> from gizmo.sources.vmh import enrich_graph_vmh
    >>> enrich_graph_vmh(mg)               # adds vmh_id attrs to metabolite nodes
    >>> from gizmo.integration.gemma import attach_community
    >>> attach_community(mg, community)    # wires taxon nodes + edges

Reaction role (producer / consumer / bidirectional)
----------------------------------------------------
AGORA2 SBML models encode stoichiometry, but GeMMA's lightweight streaming
parser (:func:`gemma.io.apollo.extract_model_reactions`) currently extracts
only metabolite *participation* (not the sign of the stoichiometric coefficient).
All edges are therefore marked ``bidirectional`` by default.  Pass a pre-built
``rxn_roles`` dict (``{rxn_id: {vmh_met_id: MicrobialEdgeRole}}``) to override.

Usage
-----
    from gizmo.sources.vmh import enrich_graph_vmh
    from gizmo.integration.gemma import attach_community, community_metabolites

    enrich_graph_vmh(mg)
    result = attach_community(mg, community)
    print(result.taxa_added, result.edges_added, result.unresolved_vmh_ids)

    # Query: which taxa touch glycoursodeoxycholate?
    mets = community_metabolites(mg, "CHEBI:132399")
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import pandas as pd

from gizmo.schema import MicrobialEdge, MicrobialEdgeRole, MicrobialTaxonNode

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result object
# ---------------------------------------------------------------------------

@dataclass
class AttachResult:
    """Summary of an :func:`attach_community` call."""
    taxa_added:        int = 0
    taxa_updated:      int = 0
    edges_added:       int = 0
    edges_updated:     int = 0
    unresolved_vmh_ids: list[str] = field(default_factory=list)
    """VMH metabolite IDs in the community that could not be matched to a GIZMO node."""


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def attach_community(
    mg,
    community,
    *,
    sample_id: str | None = None,
    rxn_roles: dict[str, dict[str, MicrobialEdgeRole]] | None = None,
    min_abundance: float = 0.0,
) -> AttachResult:
    """Attach a GeMMA :class:`~gemma.models.slimmer.SlimmedCommunity` to *mg*.

    Parameters
    ----------
    mg:
        A :class:`~gizmo.graph.network.GizmoGraph` enriched with VMH IDs
        (call :func:`~gizmo.sources.vmh.enrich_graph_vmh` first).
    community:
        A :class:`~gemma.models.slimmer.SlimmedCommunity` returned by
        :func:`~gemma.models.slimmer.slim`.
    sample_id:
        If given, per-sample abundance weights from
        ``community.reaction_weights[sample_id]`` are stored on edges.
    rxn_roles:
        Optional ``{rxn_id: {vmh_met_id: MicrobialEdgeRole}}`` dict to
        override the default ``bidirectional`` role assignment.
    min_abundance:
        Skip taxa whose mean abundance (across all samples) is below this
        threshold.  Default 0.0 keeps everything already selected by
        :func:`~gemma.models.slimmer.slim`.

    Returns
    -------
    :class:`AttachResult` with counts of added/updated nodes and edges.
    """
    result = AttachResult()

    # Build VMH → GIZMO node_id index from graph attributes
    vmh_to_node = _build_vmh_index(mg)
    log.debug("VMH index: %d entries", len(vmh_to_node))

    # Build taxon → mean abundance map from reaction_weights
    taxon_mean_abund = _taxon_mean_abundances(community)

    # Per-taxon per-metabolite: collect exchangeable metabolites + weights.
    #
    # Preferred source: community.taxon_secretome (AGORA2 R_EX_* boundary
    # exchange reactions only) — these are metabolites the organism actually
    # secretes into or imports from the gut lumen.
    #
    # Fallback (old communities without secretome): iterate all reaction
    # metabolites (original behaviour, kept for backward compatibility).
    #
    # Structure: taxon → vmh_met_id → {"rxn_ids": [...], "weight": float,
    #                                   "sample_weights": {...}}
    taxon_met_data: dict[str, dict[str, dict]] = {}

    use_secretome = bool(getattr(community, "taxon_secretome", None))

    if use_secretome:
        log.debug("attach_community: using taxon_secretome (exchange reactions only)")
        for taxon, vmh_mets in community.taxon_secretome.items():
            if taxon_mean_abund.get(taxon, 0.0) < min_abundance:
                continue
            taxon_abund_weight = taxon_mean_abund.get(taxon, 0.0)
            taxon_sample_w = _taxon_sample_abundances(community, taxon)
            taxon_met_data.setdefault(taxon, {})
            for vmh_met in vmh_mets:
                taxon_met_data[taxon].setdefault(vmh_met, {
                    "rxn_ids": [],
                    "weight": taxon_abund_weight,
                    "sample_weights": taxon_sample_w,
                })
    else:
        log.debug(
            "attach_community: taxon_secretome absent — falling back to all "
            "reaction metabolites (re-run slim() to populate secretome)"
        )
        for rxn_id, taxa in community.reaction_taxa.items():
            mets = community.reaction_metabolites.get(rxn_id, set())
            rxn_weight = _rxn_mean_weight(community, rxn_id)
            sample_weights = _rxn_sample_weights(community, rxn_id, sample_id)

            for taxon in taxa:
                if taxon_mean_abund.get(taxon, 0.0) < min_abundance:
                    continue
                taxon_met_data.setdefault(taxon, {})
                for vmh_met in mets:
                    vmh_clean = _strip_compartment(vmh_met)
                    entry = taxon_met_data[taxon].setdefault(vmh_clean, {
                        "rxn_ids": [], "weight": 0.0, "sample_weights": {},
                    })
                    if rxn_id not in entry["rxn_ids"]:
                        entry["rxn_ids"].append(rxn_id)
                    entry["weight"] = max(entry["weight"], rxn_weight)
                    for sid, w in sample_weights.items():
                        entry["sample_weights"][sid] = max(
                            entry["sample_weights"].get(sid, 0.0), w
                        )

        # Also handle community.metabolite_taxa directly
        for vmh_met, taxa in community.metabolite_taxa.items():
            vmh_clean = _strip_compartment(vmh_met)
            for taxon in taxa:
                if taxon_mean_abund.get(taxon, 0.0) < min_abundance:
                    continue
                taxon_met_data.setdefault(taxon, {})
                taxon_met_data[taxon].setdefault(vmh_clean, {
                    "rxn_ids": [], "weight": 0.0, "sample_weights": {},
                })

    # ── Add / update MicrobialTaxonNode for each taxon ────────────────────
    for taxon in taxon_met_data:
        node_id = _taxon_node_id(taxon, community)
        model_ids = community.taxon_models.get(taxon, [])
        mean_abund = taxon_mean_abund.get(taxon, 0.0)

        # Per-sample abundances from reaction_weights columns
        sample_abunds = _taxon_sample_abundances(community, taxon)

        if node_id in mg.graph.nodes:
            # Update existing node attributes in-place (graph nodes are mutable dicts)
            mg.graph.nodes[node_id]["agora2_model_ids"] = model_ids
            mg.graph.nodes[node_id]["mean_abundance"] = mean_abund
            if sample_abunds:
                mg.graph.nodes[node_id]["sample_abundances"] = sample_abunds
            result.taxa_updated += 1
        else:
            taxon_node = MicrobialTaxonNode(
                node_id=node_id,
                name=taxon,
                rank=_infer_rank(taxon, community),
                agora2_model_ids=model_ids,
                mean_abundance=mean_abund,
                sample_abundances=sample_abunds,
            )
            mg.add_microbial_taxon(taxon_node)
            result.taxa_added += 1

    # ── Add / update MicrobialEdge for each (taxon, metabolite) pair ──────
    unresolved: set[str] = set()

    for taxon, met_map in taxon_met_data.items():
        taxon_nid = _taxon_node_id(taxon, community)

        for vmh_clean, entry in met_map.items():
            gizmo_nid = vmh_to_node.get(vmh_clean)
            if gizmo_nid is None:
                unresolved.add(vmh_clean)
                continue

            role = MicrobialEdgeRole.BIDIRECTIONAL
            if rxn_roles:
                roles_for_met = {
                    rxn_roles[r][vmh_clean]
                    for r in entry["rxn_ids"]
                    if r in rxn_roles and vmh_clean in rxn_roles[r]
                }
                if roles_for_met == {MicrobialEdgeRole.PRODUCER}:
                    role = MicrobialEdgeRole.PRODUCER
                elif roles_for_met == {MicrobialEdgeRole.CONSUMER}:
                    role = MicrobialEdgeRole.CONSUMER

            edge_key = (taxon_nid, gizmo_nid)
            if mg.graph.has_edge(*edge_key):
                # Update weight if higher
                existing = mg.graph.edges[edge_key]
                if entry["weight"] > existing.get("abundance_weight", 0.0):
                    existing["abundance_weight"] = entry["weight"]
                    existing["reaction_ids"] = entry["rxn_ids"]
                result.edges_updated += 1
            else:
                edge = MicrobialEdge(
                    source=taxon_nid,
                    target=gizmo_nid,
                    role=role,
                    vmh_metabolite_id=vmh_clean,
                    reaction_ids=entry["rxn_ids"],
                    abundance_weight=entry["weight"],
                    sample_weights=entry["sample_weights"],
                )
                mg.add_microbial_edge(edge)
                result.edges_added += 1

    result.unresolved_vmh_ids = sorted(unresolved)
    log.info(
        "attach_community: +%d taxa (%d updated), +%d edges (%d updated), "
        "%d VMH IDs unresolved",
        result.taxa_added, result.taxa_updated,
        result.edges_added, result.edges_updated,
        len(result.unresolved_vmh_ids),
    )
    return result


# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------

def community_metabolites(mg, metabolite_node_id: str) -> list[dict]:
    """Return taxa connected to *metabolite_node_id* as a list of dicts.

    Each dict has: taxon_name, role, abundance_weight, reaction_ids, vmh_id.
    Results are sorted by abundance_weight descending.
    """
    results = []
    for src, tgt, attrs in mg.graph.in_edges(metabolite_node_id, data=True):
        if attrs.get("edge_type") != "microbial_metabolite":
            continue
        src_attrs = mg.graph.nodes.get(src, {})
        results.append({
            "taxon_name":      src_attrs.get("name", src),
            "role":            attrs.get("role", "bidirectional"),
            "abundance_weight": attrs.get("abundance_weight", 0.0),
            "reaction_ids":    attrs.get("reaction_ids", []),
            "vmh_id":          attrs.get("vmh_metabolite_id"),
            "mean_abundance":  src_attrs.get("mean_abundance", 0.0),
        })
    return sorted(results, key=lambda x: x["abundance_weight"], reverse=True)


def taxon_metabolites(mg, taxon_name: str) -> list[dict]:
    """Return metabolites a taxon is connected to.

    Each dict has: metabolite_node_id, name, vmh_id, role, abundance_weight.
    """
    taxon_nid = next(
        (n for n, d in mg.graph.nodes(data=True)
         if d.get("node_type") == "microbial_taxon" and d.get("name") == taxon_name),
        None,
    )
    if taxon_nid is None:
        return []
    results = []
    for src, tgt, attrs in mg.graph.out_edges(taxon_nid, data=True):
        if attrs.get("edge_type") != "microbial_metabolite":
            continue
        tgt_attrs = mg.graph.nodes.get(tgt, {})
        results.append({
            "metabolite_node_id": tgt,
            "name":               tgt_attrs.get("name", tgt),
            "vmh_id":             attrs.get("vmh_metabolite_id"),
            "role":               attrs.get("role", "bidirectional"),
            "abundance_weight":   attrs.get("abundance_weight", 0.0),
        })
    return sorted(results, key=lambda x: x["abundance_weight"], reverse=True)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _build_vmh_index(mg) -> dict[str, str]:
    """Return {vmh_abbreviation_lower: node_id} for all metabolite nodes."""
    index: dict[str, str] = {}
    for nid, attrs in mg.graph.nodes(data=True):
        if attrs.get("node_type") != "metabolite":
            continue
        vmh = attrs.get("vmh_id")
        if vmh:
            index[str(vmh).strip().lower()] = nid
    return index


def _strip_compartment(vmh_met: str) -> str:
    """Strip SBML compartment suffix from VMH metabolite ID.

    ``M_glc_D[c]`` → ``glc_D``
    ``M_glc_D__91__c__93__`` → ``glc_D``  (URL-encoded brackets)
    """
    s = vmh_met.strip()
    # Strip leading M_ prefix
    if s.startswith("M_"):
        s = s[2:]
    # URL-encoded brackets: __91__ = '[', __93__ = ']'
    s = s.replace("__91__", "[").replace("__93__", "]")
    # Strip compartment [c], [e], [p] etc.
    if s.endswith("]") and "[" in s:
        s = s[:s.rfind("[")]
    return s.lower()


def _taxon_node_id(taxon: str, community) -> str:
    rank = _infer_rank(taxon, community)
    return f"taxon:{rank}__{taxon}"


def _infer_rank(taxon: str, community) -> str:
    """Infer taxonomic rank from clade string or community metadata."""
    # MetaPhlAn clade format: "k__...|g__Prevotella"
    if "g__" in taxon:
        return "genus"
    if "s__" in taxon:
        return "species"
    # Default: assume genus-level (most common for 16S)
    return "genus"


def _taxon_mean_abundances(community) -> dict[str, float]:
    """Return {taxon: mean_relative_abundance} from reaction_weights."""
    rw = community.reaction_weights
    if rw is None or rw.empty:
        return {}
    # reaction_weights is (reactions × samples); mean across samples per reaction
    # We want taxon-level mean abundance — use taxon_models keys and reaction_taxa
    result: dict[str, float] = {}
    for taxon, model_ids in community.taxon_models.items():
        # Find reactions attributed to this taxon
        taxon_rxns = [r for r, taxa in community.reaction_taxa.items() if taxon in taxa]
        if not taxon_rxns:
            result[taxon] = 0.0
            continue
        available = [r for r in taxon_rxns if r in rw.index]
        if available:
            result[taxon] = float(rw.loc[available].mean(axis=None))
        else:
            result[taxon] = 0.0
    return result


def _taxon_sample_abundances(community, taxon: str) -> dict[str, float]:
    """Return {sample_id: mean_reaction_weight} for a taxon across samples."""
    rw = community.reaction_weights
    if rw is None or rw.empty:
        return {}
    taxon_rxns = [r for r, taxa in community.reaction_taxa.items() if taxon in taxa]
    available = [r for r in taxon_rxns if r in rw.index]
    if not available:
        return {}
    series = rw.loc[available].mean(axis=0)
    return {str(col): float(val) for col, val in series.items()}


def _rxn_mean_weight(community, rxn_id: str) -> float:
    rw = community.reaction_weights
    if rw is None or rw.empty or rxn_id not in rw.index:
        return 0.0
    return float(rw.loc[rxn_id].mean())


def _rxn_sample_weights(community, rxn_id: str, sample_id: str | None) -> dict[str, float]:
    rw = community.reaction_weights
    if rw is None or rw.empty or rxn_id not in rw.index:
        return {}
    if sample_id and sample_id in rw.columns:
        return {sample_id: float(rw.loc[rxn_id, sample_id])}
    return {str(col): float(val) for col, val in rw.loc[rxn_id].items()}


# ---------------------------------------------------------------------------
# GIZMO pathway sets — Reactome graph → {pathway_name: set[hmdb_ids]}
# ---------------------------------------------------------------------------

def build_pathway_sets_from_gizmo(
    mg,
    *,
    min_size: int = 3,
    max_size: int = 500,
    include_microbial_only: bool = False,
) -> dict[str, set[str]]:
    """Extract Reactome pathway → HMDB ID sets from a GizmoGraph.

    Traverses ReactionNode.pathways → neighbouring MetaboliteNodes →
    ``hmdb_id`` attribute.  Returns only pathways that have at least
    *min_size* members with HMDB IDs.

    This replaces the AGORA2-subsystem-based :func:`~gemma.enrichment.hypergeometric.build_pathway_sets`
    with the full Reactome human metabolism graph, giving coverage of
    sphingolipid synthesis, steroid metabolism, and other host pathways
    that AGORA2 models do not include.

    Parameters
    ----------
    mg:
        GizmoGraph — must have Reactome reactions loaded (via ReactomeLoader).
    min_size:
        Minimum number of HMDB-mapped metabolites per pathway.
    max_size:
        Maximum pathway size (avoids generic "Metabolism" catch-alls).
    include_microbial_only:
        If True, only return pathways that contain at least one metabolite
        connected to a microbial taxon node (requires prior
        :func:`attach_community` call).

    Returns
    -------
    ``{pathway_name: set[hmdb_ids]}``
    """
    g = mg.graph

    # Build pathway name lookup from PathwayNode entries
    pathway_names: dict[str, str] = {}
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") == "pathway":
            stid = attrs.get("reactome_id") or nid
            name = attrs.get("name") or stid
            pathway_names[stid] = name

    # reaction_id → list[pathway_stid]
    rxn_to_pathways: dict[str, list[str]] = {}
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "reaction":
            continue
        pws = attrs.get("pathways") or []
        if pws:
            rxn_to_pathways[nid] = list(pws)

    # pathway_stid → set[metabolite_node_ids]
    pw_met_nodes: dict[str, set[str]] = {}
    for rxn_id, pws in rxn_to_pathways.items():
        # Collect metabolite neighbours (substrates and products)
        nbrs = list(g.successors(rxn_id)) + list(g.predecessors(rxn_id))
        for nbr in nbrs:
            if g.nodes[nbr].get("node_type") == "metabolite":
                for pw in pws:
                    pw_met_nodes.setdefault(pw, set()).add(nbr)

    # Convert to HMDB ID sets
    def _node_hmdb(nid: str) -> str | None:
        attrs = g.nodes.get(nid, {})
        hmdb = attrs.get("hmdb_id") or ""
        if hmdb and str(hmdb).upper().startswith("HMDB"):
            return _normalise_hmdb(str(hmdb))
        return None

    # Microbial node set for filtering
    microbial_met_nodes: set[str] = set()
    if include_microbial_only:
        for src, tgt, attrs in g.edges(data=True):
            if attrs.get("edge_type") == "microbial_metabolite":
                microbial_met_nodes.add(tgt)

    pathway_sets: dict[str, set[str]] = {}
    for pw_stid, met_nodes in pw_met_nodes.items():
        if include_microbial_only:
            met_nodes = met_nodes & microbial_met_nodes
            if not met_nodes:
                continue
        hmdb_ids = {h for n in met_nodes if (h := _node_hmdb(n)) is not None}
        if min_size <= len(hmdb_ids) <= max_size:
            name = pathway_names.get(pw_stid, pw_stid)
            pathway_sets[name] = hmdb_ids

    log.info(
        "build_pathway_sets_from_gizmo: %d pathways (HMDB-mapped, size %d–%d)",
        len(pathway_sets), min_size, max_size,
    )
    return pathway_sets


def build_microbial_pathway_weights(
    mg,
    pathway_sets: dict[str, set[str]],
) -> dict[str, float]:
    """Compute community abundance weight for each pathway.

    For each pathway, averages the ``abundance_weight`` of incoming
    ``microbial_metabolite`` edges for the pathway's member metabolites.
    Pathways whose members have no microbial connections receive weight 0.

    Requires a prior :func:`attach_community` call.

    Parameters
    ----------
    mg:
        GizmoGraph with microbial edges attached.
    pathway_sets:
        ``{pathway_name: set[hmdb_ids]}`` — as returned by
        :func:`build_pathway_sets_from_gizmo`.

    Returns
    -------
    ``{pathway_name: float}`` — higher = more community metabolic capacity
    in this pathway.
    """
    g = mg.graph

    # HMDB → max microbial abundance weight across all incoming microbial edges
    hmdb_weight: dict[str, float] = {}
    for src, tgt, attrs in g.edges(data=True):
        if attrs.get("edge_type") != "microbial_metabolite":
            continue
        tgt_attrs = g.nodes.get(tgt, {})
        hmdb = tgt_attrs.get("hmdb_id")
        if hmdb:
            h = _normalise_hmdb(str(hmdb))
            w = attrs.get("abundance_weight", 0.0)
            hmdb_weight[h] = max(hmdb_weight.get(h, 0.0), w)

    weights: dict[str, float] = {}
    for pw_name, hmdb_ids in pathway_sets.items():
        member_weights = [hmdb_weight.get(h, 0.0) for h in hmdb_ids]
        weights[pw_name] = float(sum(member_weights) / len(member_weights)) if member_weights else 0.0

    return weights


# ---------------------------------------------------------------------------
# End-to-end GIZMO-aware ORA
# ---------------------------------------------------------------------------
# InChIKey-based HMDB augmentation helpers
# ---------------------------------------------------------------------------

def build_inchikey_to_hmdb(mg) -> dict[str, str]:
    """Return ``{inchikey_14char: hmdb_id}`` from graph nodes that have both.

    Uses only the first 14 characters of the InChIKey (connectivity layer),
    which is robust to stereochemistry differences between databases.
    The first HMDB hit per InChIKey prefix is kept.
    """
    result: dict[str, str] = {}
    for _nid, attrs in mg.graph.nodes(data=True):
        if attrs.get("node_type") != "metabolite":
            continue
        ik = attrs.get("inchikey")
        hmdb = attrs.get("hmdb_id")
        if ik and hmdb:
            ik14 = str(ik)[:14]
            if ik14 not in result:
                result[ik14] = _normalise_hmdb(str(hmdb))
    return result


def build_inchikey_to_vmh(mg) -> dict[str, str]:
    """Return ``{inchikey_14char: vmh_id}`` from graph nodes that have both.

    Useful for extending the ``hmdb_to_vmh`` dict in punch-above-abundance
    analysis when HMDB IDs are missing but InChIKeys are available.
    """
    result: dict[str, str] = {}
    for _nid, attrs in mg.graph.nodes(data=True):
        if attrs.get("node_type") != "metabolite":
            continue
        ik = attrs.get("inchikey")
        vmh = attrs.get("vmh_id")
        if ik and vmh:
            ik14 = str(ik)[:14]
            if ik14 not in result:
                result[ik14] = str(vmh)
    return result


def _augment_metabolomics_hmdb(metabolomics, mg):
    """Return a copy of *metabolomics* with HMDB IDs filled in via InChIKey.

    For features that have an ``INCHIKEY`` (or ``InChIKey``) column in their
    metadata but lack an ``HMDB`` value, this function looks up the HMDB ID
    from GIZMO graph nodes that share the same InChIKey connectivity layer
    (first 14 characters).

    If no ``INCHIKEY`` column is present or no new matches are found, the
    original object is returned unchanged.
    """
    from gemma.io.metabolomics import MetabolomicsTable

    meta = metabolomics.metadata.copy()

    # Normalise column names
    ik_col = next((c for c in meta.columns if c.upper() == "INCHIKEY"), None)
    hmdb_col = next((c for c in meta.columns if c.upper() == "HMDB"), None)

    if ik_col is None:
        log.debug("_augment_metabolomics_hmdb: no INCHIKEY column found in metadata")
        return metabolomics

    if hmdb_col is None:
        meta["HMDB"] = None
        hmdb_col = "HMDB"

    ik14_to_hmdb = build_inchikey_to_hmdb(mg)
    if not ik14_to_hmdb:
        log.debug("_augment_metabolomics_hmdb: graph has no nodes with both inchikey+hmdb_id")
        return metabolomics

    filled = 0
    for idx in meta.index:
        if pd.notna(meta.at[idx, hmdb_col]) and str(meta.at[idx, hmdb_col]).strip():
            continue  # already has HMDB
        ik = meta.at[idx, ik_col]
        if pd.isna(ik):
            continue
        hmdb = ik14_to_hmdb.get(str(ik)[:14])
        if hmdb:
            meta.at[idx, hmdb_col] = hmdb
            filled += 1

    if filled == 0:
        return metabolomics

    log.info(
        "_augment_metabolomics_hmdb: assigned HMDB IDs to %d features via InChIKey matching",
        filled,
    )
    return MetabolomicsTable(data=metabolomics.data, metadata=meta)


# ---------------------------------------------------------------------------

def run_ora_gizmo(
    mg,
    community,
    metabolomics,
    *,
    min_size: int = 3,
    max_size: int = 500,
    min_prevalence: float = 0.0,
    microbial_weight_boost: float = 1.5,
    alpha: float = 0.05,
) -> "ORAGizmoResult":
    """Run pathway ORA using GIZMO's Reactome graph + community abundance weights.

    This replaces the AGORA2-subsystem ORA with a host+microbiome integrated
    analysis:

    * **Background**: all HMDB IDs present in the GIZMO graph (thousands vs
      the ~1600 in AGORA2 alone)
    * **Pathway sets**: Reactome human pathways (sphingolipids, steroids,
      amino acids, …) extracted from GIZMO
    * **Microbial weighting**: pathways with high community metabolic
      activity receive a fold-enrichment boost of up to *microbial_weight_boost*×

    Parameters
    ----------
    mg:
        GizmoGraph enriched with VMH IDs and community attached
        (call :func:`~gizmo.sources.vmh.enrich_graph_vmh` and
        :func:`attach_community` first).
    community:
        :class:`~gemma.models.slimmer.SlimmedCommunity` — used to compute
        microbial pathway weights.
    metabolomics:
        :class:`~gemma.io.metabolomics.MetabolomicsTable`.
    min_size / max_size:
        Pathway size filter.
    min_prevalence:
        Minimum sample prevalence for a feature to enter the hit set.
    microbial_weight_boost:
        Maximum fold-enrichment multiplier applied to pathways with
        maximum microbial activity (linear scale between 0 and this value).
    alpha:
        FDR threshold for reporting.

    Returns
    -------
    :class:`ORAGizmoResult`
    """
    from gemma.enrichment.hypergeometric import run_ora_per_sample, EnrichmentResult

    # Build pathway sets from Reactome graph
    pathway_sets = build_pathway_sets_from_gizmo(mg, min_size=min_size, max_size=max_size)
    if not pathway_sets:
        raise RuntimeError(
            "build_pathway_sets_from_gizmo returned 0 pathways. "
            "Ensure the GizmoGraph has Reactome reactions loaded."
        )

    # Background: all HMDB IDs in the graph
    background: set[str] = set()
    for _nid, attrs in mg.graph.nodes(data=True):
        if attrs.get("node_type") == "metabolite":
            hmdb = attrs.get("hmdb_id")
            if hmdb:
                background.add(_normalise_hmdb(str(hmdb)))
    log.info("GIZMO background universe: %d HMDB IDs", len(background))

    # InChIKey-augmented metabolomics: assign HMDB IDs to features that lack
    # them but have an InChIKey matching a graph node.
    metabolomics = _augment_metabolomics_hmdb(metabolomics, mg)

    # Per-sample ORA using HMDB IDs
    from gemma.io.metabolomics import _normalise_hmdb as _gemma_normalise_hmdb
    ora_results = run_ora_per_sample(
        metabolomics,
        pathway_sets,
        background,
        min_prevalence=min_prevalence,
        min_pathway_size=min_size,
        max_pathway_size=max_size,
    )

    # Compute microbial pathway weights
    microbial_weights = build_microbial_pathway_weights(mg, pathway_sets)
    max_w = max(microbial_weights.values()) if microbial_weights else 1.0

    # Apply microbial weight boost to fold_enrichment and re-rank
    boosted_results: dict[str, "pd.DataFrame"] = {}
    for sample_id, er in ora_results.items():
        df = er.table.copy()
        if not df.empty and max_w > 0:
            boost = df["pathway"].map(
                lambda pw: 1.0 + (microbial_weight_boost - 1.0) * microbial_weights.get(pw, 0.0) / max_w
            )
            df["microbial_weight"] = df["pathway"].map(
                lambda pw: microbial_weights.get(pw, 0.0)
            )
            df["fold_enrichment_boosted"] = df["fold_enrichment"] * boost
        else:
            df["microbial_weight"] = 0.0
            df["fold_enrichment_boosted"] = df.get("fold_enrichment", 0.0)
        boosted_results[sample_id] = df

    return ORAGizmoResult(
        pathway_sets=pathway_sets,
        background=background,
        microbial_weights=microbial_weights,
        sample_results=boosted_results,
        n_pathways=len(pathway_sets),
        n_background=len(background),
    )


import pandas as pd
from dataclasses import dataclass as _dataclass, field as _field


@_dataclass
class ORAGizmoResult:
    """Result of :func:`run_ora_gizmo`.

    Attributes
    ----------
    pathway_sets:
        ``{pathway_name: set[hmdb_ids]}`` used for testing.
    background:
        HMDB ID universe drawn from the GIZMO graph.
    microbial_weights:
        ``{pathway_name: float}`` — community abundance weight per pathway.
    sample_results:
        ``{sample_id: pd.DataFrame}`` — per-sample ORA tables with added
        ``microbial_weight`` and ``fold_enrichment_boosted`` columns.
    n_pathways:
        Number of Reactome pathways tested.
    n_background:
        Size of HMDB background universe.
    """

    pathway_sets:      dict
    background:        set
    microbial_weights: dict
    sample_results:    dict
    n_pathways:        int = 0
    n_background:      int = 0

    def significant(
        self,
        alpha: float = 0.05,
        min_microbial_weight: float = 0.0,
    ) -> "pd.DataFrame":
        """Return significant pathways pooled across samples.

        Parameters
        ----------
        alpha:
            FDR threshold.
        min_microbial_weight:
            If > 0, only return pathways with microbial weight above this value
            (filters to microbiome-attributable signals).
        """
        frames = []
        for sid, df in self.sample_results.items():
            if df.empty:
                continue
            sig = df[df["p_adj"] < alpha].copy()
            if min_microbial_weight > 0:
                sig = sig[sig["microbial_weight"] >= min_microbial_weight]
            sig.insert(0, "sample_id", sid)
            frames.append(sig)
        if not frames:
            return pd.DataFrame()
        out = pd.concat(frames, ignore_index=True)
        return out.sort_values(["p_adj", "fold_enrichment_boosted"], ascending=[True, False])

    def summary(self) -> str:
        n_samples = len(self.sample_results)
        n_sig = sum(
            (df["p_adj"] < 0.05).sum()
            for df in self.sample_results.values()
            if not df.empty
        )
        mw_nonzero = sum(1 for w in self.microbial_weights.values() if w > 0)
        return (
            f"ORAGizmoResult: {self.n_pathways} Reactome pathways × "
            f"{n_samples} samples\n"
            f"  Background: {self.n_background} HMDB IDs\n"
            f"  Microbial-active pathways: {mw_nonzero}/{self.n_pathways}\n"
            f"  Significant pathway-sample pairs (FDR<0.05): {n_sig}"
        )


# ---------------------------------------------------------------------------
# GIZMO-aware concordance
# ---------------------------------------------------------------------------

def build_gizmo_capacity(
    mg,
    pathway_sets: dict[str, set[str]],
) -> "pd.DataFrame":
    """Build a (pathways × samples) community capacity matrix from microbial edges.

    For each Reactome pathway and each sample, capacity is the mean
    ``abundance_weight`` of incoming microbial edges on that pathway's
    member metabolites.  Zero means the community has no modelled capacity
    for those metabolites in that sample.

    Requires a prior :func:`attach_community` call so that microbial edges
    with ``sample_weights`` are present in the graph.

    Parameters
    ----------
    mg:
        GizmoGraph with microbial edges (``edge_type == "microbial_metabolite"``).
    pathway_sets:
        ``{pathway_name: set[hmdb_ids]}`` — as returned by
        :func:`build_pathway_sets_from_gizmo`.

    Returns
    -------
    DataFrame with pathways as rows and sample IDs as columns, values in [0, ∞).
    """
    g = mg.graph

    # Build HMDB → max sample_weights across all incoming microbial edges
    # {hmdb_id: {sample_id: max_weight}}
    hmdb_sample_weights: dict[str, dict[str, float]] = {}
    for src, tgt, attrs in g.edges(data=True):
        if attrs.get("edge_type") != "microbial_metabolite":
            continue
        tgt_attrs = g.nodes.get(tgt, {})
        hmdb = tgt_attrs.get("hmdb_id")
        if not hmdb:
            continue
        h = _normalise_hmdb(str(hmdb))
        sample_w = attrs.get("sample_weights") or {}
        entry = hmdb_sample_weights.setdefault(h, {})
        for sid, w in sample_w.items():
            entry[sid] = max(entry.get(sid, 0.0), float(w))

    # Collect all sample IDs
    all_samples: set[str] = set()
    for sw in hmdb_sample_weights.values():
        all_samples.update(sw.keys())
    samples = sorted(all_samples)

    if not samples:
        log.warning("build_gizmo_capacity: no sample_weights found on microbial edges. "
                    "Re-run attach_community with a sample_id or re-run slim() with per-sample data.")
        return pd.DataFrame()

    # Build capacity matrix
    rows: dict[str, list[float]] = {}
    for pw_name, hmdb_ids in pathway_sets.items():
        weights_per_sample: dict[str, list[float]] = {s: [] for s in samples}
        for h in hmdb_ids:
            if h in hmdb_sample_weights:
                for sid in samples:
                    w = hmdb_sample_weights[h].get(sid, 0.0)
                    weights_per_sample[sid].append(w)
        row = [
            float(sum(weights_per_sample[s]) / len(weights_per_sample[s]))
            if weights_per_sample[s] else 0.0
            for s in samples
        ]
        rows[pw_name] = row

    cap = pd.DataFrame(rows, index=samples).T   # pathways × samples
    log.info("build_gizmo_capacity: %d pathways × %d samples", cap.shape[0], cap.shape[1])
    return cap


def compute_concordance_gizmo(
    mg,
    community,
    ora_result: "ORAGizmoResult",
    *,
    capacity_threshold: float = 0.3,
    activity_threshold: float = 0.1,
) -> "ConcordanceResult":
    """Run GeMMA concordance using GIZMO Reactome pathways.

    Compares:

    * **Capacity (C)** — community metabolic capacity per Reactome pathway per
      sample, derived from microbial edge ``sample_weights`` in the GIZMO graph.
    * **Activity (A)** — metabolite detection fraction per pathway per sample,
      from the GIZMO ORA results (``n_overlap / n_pathway``).

    Labels each (pathway, sample) cell as one of:

    * ``active``    — community has capacity AND metabolites are detected
    * ``silent``    — capacity present but metabolites not detected
                      (regulation, export, host suppression?)
    * ``exogenous`` — metabolites detected but no community capacity
                      (host synthesis, diet, or AGORA2 model gap)
    * ``absent``    — neither signal

    Parameters
    ----------
    mg:
        GizmoGraph enriched with VMH IDs and microbial edges attached.
    community:
        :class:`~gemma.models.slimmer.SlimmedCommunity` used in the original
        :func:`attach_community` call.
    ora_result:
        :class:`ORAGizmoResult` from :func:`run_ora_gizmo`.
    capacity_threshold / activity_threshold:
        Normalised thresholds for labelling (both 0–1 after min-max scaling).

    Returns
    -------
    :class:`~gemma.integration.concordance.ConcordanceResult`
    """
    from gemma.integration.concordance import (
        ConcordanceResult, _minmax_norm, _assign_labels,
    )

    pathway_sets = ora_result.pathway_sets

    # Capacity matrix from microbial edges
    cap_raw = build_gizmo_capacity(mg, pathway_sets)
    if cap_raw.empty:
        log.warning("compute_concordance_gizmo: empty capacity matrix — no per-sample weights.")
        return ConcordanceResult(
            capacity=pd.DataFrame(), activity=pd.DataFrame(),
            labels=pd.DataFrame(), rank="gizmo",
        )
    cap = _minmax_norm(cap_raw)   # pathways × samples, 0–1

    # Activity matrix from ORA detection fractions
    samples = list(ora_result.sample_results.keys())
    all_pathways = list(pathway_sets.keys())

    act_data: dict[str, list[float]] = {pw: [] for pw in all_pathways}
    for sample in samples:
        df = ora_result.sample_results.get(sample, pd.DataFrame())
        if df.empty:
            score_map: dict[str, float] = {}
        else:
            score_map = {
                row["pathway"]: (row["n_overlap"] / row["n_pathway"]
                                 if row.get("n_pathway", 0) > 0 else 0.0)
                for _, row in df.iterrows()
            }
        for pw in all_pathways:
            act_data[pw].append(score_map.get(pw, 0.0))

    act_raw = pd.DataFrame(act_data, index=samples).T   # pathways × samples
    act = _minmax_norm(act_raw)

    # Align
    cap = cap.reindex(index=all_pathways, columns=samples, fill_value=0.0)
    act = act.reindex(index=all_pathways, columns=samples, fill_value=0.0)

    labels = _assign_labels(cap, act, capacity_threshold, activity_threshold)

    log.info(
        "compute_concordance_gizmo: %d pathways × %d samples — "
        "active=%d, silent=%d, exogenous=%d, absent=%d",
        len(all_pathways), len(samples),
        (labels == "active").sum().sum(),
        (labels == "silent").sum().sum(),
        (labels == "exogenous").sum().sum(),
        (labels == "absent").sum().sum(),
    )
    return ConcordanceResult(
        capacity=cap, activity=act, labels=labels, rank="gizmo",
    )


# ---------------------------------------------------------------------------
# Private normalisation helper (mirrors gemma.io.metabolomics._normalise_hmdb)
# ---------------------------------------------------------------------------

def _normalise_hmdb(s: str) -> str:
    s = str(s).strip().upper()
    if s.startswith("HMDB"):
        return f"HMDB{s[4:].zfill(7)}"
    return s
