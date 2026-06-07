"""
BuildPipeline — reproducible, versioned GIZMO graph construction.

Orchestrates all source adapters, runs QC, and writes a graph bundle:

    data/processed/<graph_name>/
        graph.json
        graph.graphml
        qc_report.json
        manifest.json

Usage::

    from gizmo.build.pipeline import BuildPipeline

    pipe = BuildPipeline("human_iem", cache_dir="data/raw")
    pipe.add_reactome(species="Homo sapiens")
    pipe.add_mondo_iem()
    pipe.add_orphanet()
    pipe.add_open_targets(disease_ids=["MONDO:0019052"])
    pipe.add_metabolon("data/resources/gizmo/sources/metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv",
                       metanetx_prop="data/resources/gizmo/metanetx/chem_prop.tsv",
                       metanetx_xref="data/resources/gizmo/metanetx/chem_xref.tsv")
    pipe.add_stringdb(min_score=0.4)
    mg, manifest = pipe.run()
    pipe.save_bundle(mg, manifest, output_dir="data/processed")
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from gizmo.build.manifest import GraphManifest, SourceRecord
from gizmo.graph.network import GizmoGraph

log = logging.getLogger(__name__)


@dataclass
class _Step:
    name: str
    fn: Any                    # callable(mg) -> SourceRecord
    enabled: bool = True


class BuildPipeline:
    """
    Declarative, ordered graph build pipeline.

    Each ``add_*`` method registers a build step.  Calling ``.run()``
    executes them in order, then flags currency metabolites and returns
    (GizmoGraph, GraphManifest).
    """

    def __init__(
        self,
        graph_name: str,
        cache_dir: str | Path = "data/raw",
        notes: str = "",
    ) -> None:
        self.graph_name = graph_name
        self.cache_dir  = Path(cache_dir)
        self.notes      = notes
        self._steps: list[_Step] = []
        self._params: dict[str, Any] = {}

    # ------------------------------------------------------------------
    # Build step registration
    # ------------------------------------------------------------------

    def add_reactome(
        self,
        species: str = "Homo sapiens",
        pathway_stids: list[str] | None = None,
        max_workers: int = 20,
    ) -> "BuildPipeline":
        """
        Load Reactome reactions.

        Parameters
        ----------
        species      : load full species (ignored if pathway_stids given)
        pathway_stids: load only these top-level pathway stIds
        max_workers  : parallel fetch threads
        """
        self._params["reactome"] = {
            "species": species,
            "pathway_stids": pathway_stids,
        }

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.reactome import ReactomeLoader
            loader = ReactomeLoader(
                cache_dir=str(self.cache_dir / "reactome"),
                max_workers=max_workers,
            )
            if pathway_stids:
                log.info("Reactome: loading %d pathway(s)…", len(pathway_stids))
                sub = loader.load_pathways(pathway_stids)
            else:
                log.info("Reactome: loading full species '%s'…", species)
                sub = loader.load_species(species)

            # Merge sub into mg
            for nid, attrs in sub.graph.nodes(data=True):
                mg.graph.add_node(nid, **attrs)
            for u, v, attrs in sub.graph.edges(data=True):
                mg.graph.add_edge(u, v, **attrs)

            n_rxns = len(sub.reaction_nodes())
            log.info("Reactome: merged %d reactions.", n_rxns)
            return SourceRecord(
                name="reactome",
                version="current",
                url="https://reactome.org",
                license="CC BY 4.0",
                n_records=n_rxns,
            )

        self._steps.append(_Step("reactome", _step))
        return self

    def add_mondo_iem(self) -> "BuildPipeline":
        """Add MONDO inborn errors of metabolism disease nodes."""
        self._params["mondo_iem"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.mondo import MondoClient
            client = MondoClient(cache_dir=str(self.cache_dir / "mondo"))
            client.download()
            diseases = client.load_iem_subset()
            mg.add_diseases(diseases)
            log.info("MONDO IEM: added %d disease nodes.", len(diseases))
            return SourceRecord(
                name="mondo",
                version="current",
                url="https://mondo.monarchinitiative.org",
                license="CC BY 4.0",
                n_records=len(diseases),
            )

        self._steps.append(_Step("mondo_iem", _step))
        return self

    def add_mondo_all(self) -> "BuildPipeline":
        """Add all MONDO disease nodes (not just IEM)."""
        self._params["mondo_all"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.mondo import MondoClient
            client = MondoClient(cache_dir=str(self.cache_dir / "mondo"))
            client.download()
            diseases = client.load_all()
            mg.add_diseases(diseases)
            log.info("MONDO all: added %d disease nodes.", len(diseases))
            return SourceRecord(
                name="mondo",
                version="current",
                url="https://mondo.monarchinitiative.org",
                license="CC BY 4.0",
                n_records=len(diseases),
            )

        self._steps.append(_Step("mondo_all", _step))
        return self

    def add_orphanet(self, iem_only: bool = True) -> "BuildPipeline":
        """
        Add Orphanet disease nodes and gene–disease edges.

        Parameters
        ----------
        iem_only : if True, restrict to inborn errors of metabolism
        """
        self._params["orphanet"] = {"iem_only": iem_only}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.orphanet import OrphanetClient
            from gizmo.schema import DiseaseEdge, DiseaseEdgeType

            client = OrphanetClient(cache_dir=str(self.cache_dir / "orphanet"))
            client.download()
            diseases = (
                client.load_iem_diseases() if iem_only
                else client.load_diseases()
            )
            mg.add_diseases(diseases)

            # Gene–disease edges (returns tuple: gene_nodes, edges)
            gene_nodes, edges = client.load_gene_associations()
            for gene in gene_nodes:
                mg.add_gene(gene)
            n_orpha_edges = 0
            for edge in edges:
                try:
                    mg.add_disease_edge(edge)
                    n_orpha_edges += 1
                except Exception:
                    pass

            # Cross-reference: propagate Orphanet gene edges to matching MONDO nodes.
            # MONDO nodes store xref_orphanet in attrs; use that to add duplicate
            # disease→gene edges from the MONDO node so MONDO diseases are connected.
            orpha_gene_edges: dict[str, list[DiseaseEdge]] = {}
            for edge in edges:
                orpha_gene_edges.setdefault(edge.source, []).append(edge)

            n_mondo_edges = 0
            for nid, attrs in list(mg.graph.nodes(data=True)):
                if attrs.get("node_type") != "disease":
                    continue
                if not nid.startswith("MONDO:"):
                    continue
                xref_orphanet = attrs.get("xref_orphanet") or []
                if isinstance(xref_orphanet, str):
                    xref_orphanet = [xref_orphanet]
                for orpha_ref in xref_orphanet:
                    # normalise to "Orphanet:XXXX" key
                    orpha_key = orpha_ref if orpha_ref.startswith("Orphanet:") else f"Orphanet:{orpha_ref}"
                    for orig_edge in orpha_gene_edges.get(orpha_key, []):
                        rewired = DiseaseEdge(
                            source=nid,
                            target=orig_edge.target,
                            edge_type=DiseaseEdgeType.GENE_ASSOCIATED,
                            score=orig_edge.score,
                            source_db="orphanet",
                        )
                        try:
                            mg.add_disease_edge(rewired)
                            n_mondo_edges += 1
                        except Exception:
                            pass

            log.info(
                "Orphanet: added %d disease nodes, %d genes, %d Orphanet→gene edges, "
                "%d MONDO→gene edges (via xref).",
                len(diseases), len(gene_nodes), n_orpha_edges, n_mondo_edges,
            )
            return SourceRecord(
                name="orphanet",
                version="current",
                url="https://www.orphadata.com",
                license="CC BY 4.0",
                n_records=len(diseases),
            )

        self._steps.append(_Step("orphanet", _step))
        return self

    def add_open_targets(
        self,
        disease_ids: list[str] | None = None,
        min_score: float = 0.1,
    ) -> "BuildPipeline":
        """
        Add gene–disease associations from Open Targets.

        Parameters
        ----------
        disease_ids : MONDO IDs to query (queries disease nodes already in graph if None)
        min_score   : minimum association score [0, 1]
        """
        self._params["open_targets"] = {
            "disease_ids": disease_ids,
            "min_score": min_score,
        }

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.open_targets import OpenTargetsClient
            from gizmo.schema import DiseaseEdge, DiseaseEdgeType, GeneNode

            client = OpenTargetsClient()

            # Determine disease IDs to query
            ids_to_query: list[str] = disease_ids or [
                nid for nid in mg.disease_nodes()
                if nid.startswith("MONDO:")
            ]
            if not ids_to_query:
                log.warning("Open Targets: no MONDO disease IDs to query.")
                return SourceRecord(
                    name="open_targets", license="CC BY 4.0",
                    url="https://platform.opentargets.org",
                )

            n_edges = 0
            seen_genes: set[str] = set()
            for did in ids_to_query[:50]:   # practical API limit
                try:
                    gene_nodes, disease_edges = client.gene_associations_for_disease(
                        did, min_score=min_score,
                    )
                except Exception as exc:
                    log.debug("Open Targets skip %s: %s", did, exc)
                    continue
                for gene in gene_nodes:
                    if gene.node_id not in seen_genes:
                        mg.add_gene(gene)
                        seen_genes.add(gene.node_id)
                for edge in disease_edges:
                    try:
                        mg.add_disease_edge(edge)
                        n_edges += 1
                    except Exception:
                        pass

            log.info("Open Targets: added %d gene-disease edges.", n_edges)
            return SourceRecord(
                name="open_targets",
                version="current",
                url="https://platform.opentargets.org",
                license="CC BY 4.0",
                n_records=n_edges,
            )

        self._steps.append(_Step("open_targets", _step))
        return self

    def add_metabolon(
        self,
        csv_path: str | Path,
        metanetx_prop: str | Path | None = None,
        metanetx_xref: str | Path | None = None,
        overrides_path: str | Path | None = None,
    ) -> "BuildPipeline":
        """
        Map Metabolon compounds and merge into the graph.

        Parameters
        ----------
        csv_path       : Metabolon data dictionary CSV
        metanetx_prop  : chem_prop.tsv path (for fast local InChIKey lookup)
        metanetx_xref  : chem_xref.tsv path
        overrides_path : curated overrides JSON from MetabolonCurator
        """
        self._params["metabolon"] = {"csv_path": str(csv_path)}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.metabolon import MetabolonLoader

            loader = MetabolonLoader(str(csv_path))

            prop = Path(metanetx_prop) if metanetx_prop else None
            xref = Path(metanetx_xref) if metanetx_xref else None
            if prop and prop.exists() and xref and xref.exists():
                log.info("Metabolon: building MetaNetX index…")
                loader.load_metanetx_index(str(prop), str(xref))

            # Apply curator overrides if available
            if overrides_path and Path(overrides_path).exists():
                from gizmo.curation.metabolon_curator import MetabolonCurator
                curator = MetabolonCurator(
                    loader, graph=mg, overrides_path=overrides_path,
                )
                n_applied = curator.apply()
                log.info("Metabolon: applied %d curation overrides.", n_applied)

            nodes, report = loader.to_metabolite_nodes()
            mg.add_metabolites(nodes)
            log.info("Metabolon: added %d compound nodes.  %s", len(nodes), report)
            return SourceRecord(
                name="metabolon",
                url="https://metabolon.com",
                license="open subset",
                n_records=len(nodes),
            )

        self._steps.append(_Step("metabolon", _step))
        return self

    def add_stringdb(self, min_score: float = 0.4) -> "BuildPipeline":
        """
        Add STRING protein–protein interaction edges between existing gene nodes.

        Parameters
        ----------
        min_score : combined confidence threshold [0, 1]
        """
        self._params["stringdb"] = {"min_score": min_score}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.stringdb import StringDBLoader
            loader = StringDBLoader(min_score=min_score)
            n = loader.enrich(mg)
            return SourceRecord(
                name="stringdb",
                version="v12",
                url="https://string-db.org",
                license="CC BY 4.0",
                n_records=n,
            )

        self._steps.append(_Step("stringdb", _step))
        return self

    def add_chemical_enrichment(
        self,
        chem_prop_path: str | Path | None = None,
        chem_xref_path: str | Path | None = None,
    ) -> "BuildPipeline":
        """
        Enrich all metabolite nodes with SMILES, InChI, formula, and MetaNetX IDs
        from the MetaNetX flat files.  Fills in only missing fields — never overwrites.

        Metabolon-sourced nodes already have these from load_metanetx_index, but
        Reactome metabolites (identified by ChEBI ID or InChIKey) benefit most.

        Parameters
        ----------
        chem_prop_path : path to chem_prop.tsv (defaults to self.cache_dir / 'metanetx/chem_prop.tsv')
        chem_xref_path : path to chem_xref.tsv (defaults to self.cache_dir / 'metanetx/chem_xref.tsv')
        """
        self._params["chemical_enrichment"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.metanetx import _mnx_header_info
            import re, math

            prop_path = Path(chem_prop_path) if chem_prop_path else self.cache_dir / "metanetx" / "chem_prop.tsv"
            xref_path = Path(chem_xref_path) if chem_xref_path else self.cache_dir / "metanetx" / "chem_xref.tsv"

            if not prop_path.exists() or not xref_path.exists():
                log.warning(
                    "Chemical enrichment skipped: MetaNetX files not found (%s, %s).",
                    prop_path, xref_path,
                )
                return SourceRecord(name="chemical_enrichment", license="CC BY 4.0",
                                    url="https://metanetx.org")

            # --- Step 0: collect target InChIKeys and ChEBI IDs from the graph ---
            # Only index compounds that are actually in the graph (<<1M).
            target_inchikeys: set[str] = set()
            target_chebi: set[str] = set()
            for nid, data in mg.graph.nodes(data=True):
                if data.get("node_type") != "metabolite":
                    continue
                ik = (data.get("inchikey") or "").strip()
                ch = (data.get("chebi_id") or "").strip()
                if ik:
                    target_inchikeys.add(ik)
                if ch:
                    target_chebi.add(ch)

            if not target_inchikeys and not target_chebi:
                log.info("Chemical enrichment: no metabolite nodes to enrich.")
                return SourceRecord(name="chemical_enrichment", license="CC BY 4.0",
                                    url="https://metanetx.org")

            # --- Pass 1: chem_xref → {ChEBI_ID: MNX_ID} for our target ChEBI set ---
            xref_cols, xref_start = _mnx_header_info(xref_path)
            mnx_idx = xref_cols.index("ID") if "ID" in xref_cols else 1
            chebi_to_mnx: dict[str, str] = {}
            mnx_target: set[str] = set()  # MNX IDs for our target ChEBI nodes
            with open(xref_path, encoding="utf-8") as fh:
                for lineno, line in enumerate(fh):
                    if lineno < xref_start:
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) <= mnx_idx:
                        continue
                    source = parts[0]
                    mnx_id = parts[mnx_idx]
                    if source.startswith("chebi:"):
                        chebi_id = "CHEBI:" + re.sub(r"^chebi:(?:CHEBI:)?", "", source)
                        if chebi_id in target_chebi:
                            chebi_to_mnx.setdefault(chebi_id, mnx_id)
                            mnx_target.add(mnx_id)

            # --- Pass 2: chem_prop → only index rows matching our InChIKeys or MNX IDs ---
            prop_cols, prop_start = _mnx_header_info(prop_path)

            def _ci(name: str) -> int:
                try:
                    return next(i for i, c in enumerate(prop_cols) if c.lower() == name.lower())
                except StopIteration:
                    return -1

            id_idx = _ci("ID") if _ci("ID") >= 0 else 0
            ik_idx = _ci("InChIKey")
            sm_idx = _ci("SMILES")
            in_idx = _ci("InChI")
            fo_idx = _ci("formula")
            ma_idx = _ci("mass")

            # {inchikey: props}  and  {mnx_id: props} — only for our target compounds
            ik_props:  dict[str, dict] = {}
            mnx_props: dict[str, dict] = {}

            with open(prop_path, encoding="utf-8") as fh:
                for lineno, line in enumerate(fh):
                    if lineno < prop_start:
                        continue
                    parts = line.rstrip("\n").split("\t")
                    mnx_id   = parts[id_idx] if id_idx < len(parts) else ""
                    inchikey  = (parts[ik_idx].strip() if ik_idx >= 0 and ik_idx < len(parts) else "")

                    # Only process rows relevant to our graph
                    if inchikey not in target_inchikeys and mnx_id not in mnx_target:
                        continue

                    def _get(idx: int) -> str:
                        return parts[idx].strip() if idx >= 0 and idx < len(parts) else ""

                    mass_s = _get(ma_idx)
                    try:
                        mass_v: float | None = float(mass_s) if mass_s else None
                    except ValueError:
                        mass_v = None
                    if mass_v is not None and math.isnan(mass_v):
                        mass_v = None

                    props = {
                        "metanetx_id": mnx_id,
                        "smiles":  _get(sm_idx),
                        "inchi":   _get(in_idx),
                        "formula": _get(fo_idx),
                        "mass":    mass_v,
                    }
                    if inchikey:
                        ik_props.setdefault(inchikey, props)
                    if mnx_id:
                        mnx_props.setdefault(mnx_id, props)

            # --- Enrich metabolite nodes ---
            enriched = 0
            total_met = 0
            for nid, data in mg.graph.nodes(data=True):
                if data.get("node_type") != "metabolite":
                    continue
                total_met += 1
                inchikey = (data.get("inchikey") or "").strip()
                chebi_id = (data.get("chebi_id") or "").strip()

                props = (
                    ik_props.get(inchikey)
                    or (mnx_props.get(chebi_to_mnx[chebi_id]) if chebi_id in chebi_to_mnx else None)
                )
                if not props:
                    continue

                changed = False
                for field in ("smiles", "inchi", "formula", "metanetx_id"):
                    if props.get(field) and not data.get(field):
                        mg.graph.nodes[nid][field] = props[field]
                        changed = True
                if props.get("mass") is not None and not data.get("mass"):
                    mg.graph.nodes[nid]["mass"] = props["mass"]
                    changed = True
                if not data.get("metanetx_id") and chebi_id in chebi_to_mnx:
                    mg.graph.nodes[nid]["metanetx_id"] = chebi_to_mnx[chebi_id]
                    changed = True
                if changed:
                    enriched += 1

            log.info("Chemical enrichment: enriched %d / %d metabolite nodes with MetaNetX data.",
                     enriched, total_met)
            return SourceRecord(
                name="chemical_enrichment",
                version="4.5",
                url="https://metanetx.org",
                license="CC BY 4.0",
                n_records=enriched,
            )

        self._steps.append(_Step("chemical_enrichment", _step))
        return self

    def collapse_orphan_metab_twins(self) -> "BuildPipeline":
        """
        Merge orphan metabolite nodes into their connected Reactome twins.

        ROOT CAUSE this fixes: ChEBI / Metabolon / PubChem / HMDB adapters
        each add their own metabolite nodes. Reactome metabolites are
        connected to reactions via substrate/product edges, but the other
        sources' nodes have ZERO edges (orphans). The MetaboliteMapper
        prefers exact-ID hits and lands on orphans, so the downstream
        graph signal vanishes silently.

        FIX: For each orphan metabolite, find its connected Reactome twin
        by InChIKey14, ChEBI ID, or normalized name. Merge the orphan's
        IDs/metadata into the twin, rewrite any edges touching the orphan
        to point to the twin, delete the orphan.

        Must run AFTER add_chemical_enrichment() so Reactome twins have
        their chebi_id / inchikey populated for matching.

        Orphans that don't match any Reactome twin are kept as-is
        (genuinely not in Reactome — many specialized lipids, xenobiotics,
        Tier-2/3 isomer mixes).
        """
        self._params["collapse_orphan_metab_twins"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            import re
            from collections import defaultdict, Counter

            _compartment_re = re.compile(r"\s*\[[^\]]+\]\s*$")
            _paren_re = re.compile(r"\s*\([^)]*\)\s*$")
            def _normalize_name(s):
                if not s: return ""
                return _compartment_re.sub("", s.lower().strip()).strip()
            def _loose_normalize(s):
                t = _normalize_name(s)
                t = _paren_re.sub("", t).strip()
                return re.sub(r"[\s\-_,;]+", "", t)

            g = mg.graph
            # Index Reactome metab nodes
            reactome_by_name = defaultdict(list)
            reactome_by_loose = defaultdict(list)
            reactome_by_ik14 = defaultdict(list)
            reactome_by_chebi = defaultdict(list)
            for nid, attrs in g.nodes(data=True):
                if attrs.get("node_type") != "metabolite": continue
                if not nid.startswith("reactome:"): continue
                nm = _normalize_name(attrs.get("name") or "")
                if nm: reactome_by_name[nm].append(nid)
                ln = _loose_normalize(attrs.get("name") or "")
                if ln: reactome_by_loose[ln].append(nid)
                ik = attrs.get("inchikey")
                if ik and len(ik) >= 14: reactome_by_ik14[ik[:14]].append(nid)
                ch = attrs.get("chebi_id")
                if ch:
                    cb = str(ch).replace("CHEBI:", "").replace("chebi:", "")
                    reactome_by_chebi[cb].append(nid)
                    reactome_by_chebi[f"CHEBI:{cb}"].append(nid)
                for alt_ch in (attrs.get("all_chebi_ids") or []):
                    cb = str(alt_ch).replace("CHEBI:", "").replace("chebi:", "")
                    if nid not in reactome_by_chebi[cb]:
                        reactome_by_chebi[cb].append(nid)
                    if nid not in reactome_by_chebi[f"CHEBI:{cb}"]:
                        reactome_by_chebi[f"CHEBI:{cb}"].append(nid)

            # Build orphan → twin
            orphan_to_twin = {}
            strat = Counter()
            for nid, attrs in g.nodes(data=True):
                if attrs.get("node_type") != "metabolite": continue
                if nid.startswith("reactome:"): continue
                twin = None; strategy = None
                ik = attrs.get("inchikey")
                if not twin and ik and len(ik) >= 14:
                    cands = reactome_by_ik14.get(ik[:14], [])
                    if cands: twin = cands[0]; strategy = "inchikey14"
                if not twin:
                    ch = attrs.get("chebi_id")
                    if not ch and nid.startswith("CHEBI:"):
                        ch = nid
                    if ch:
                        cb = str(ch).replace("CHEBI:", "").replace("chebi:", "")
                        for key in (f"CHEBI:{cb}", cb):
                            cands = reactome_by_chebi.get(key, [])
                            if cands: twin = cands[0]; strategy = "chebi_id"; break
                if not twin:
                    nm = _normalize_name(attrs.get("name") or "")
                    cands = reactome_by_name.get(nm, []) if nm else []
                    if cands: twin = cands[0]; strategy = "name_exact"
                if not twin:
                    ln = _loose_normalize(attrs.get("name") or "")
                    cands = reactome_by_loose.get(ln, []) if ln else []
                    if cands: twin = cands[0]; strategy = "name_loose"
                if twin:
                    orphan_to_twin[nid] = twin
                    strat[strategy] += 1

            # Merge metadata + rewrite edges + drop orphans
            n_edges_rewritten = 0
            for oid, tid in orphan_to_twin.items():
                # Copy IDs/metadata from orphan onto twin if twin has None
                o_attrs = g.nodes[oid]
                t_attrs = g.nodes[tid]
                for k in ("chebi_id", "hmdb_id", "pubchem_cid", "metabolon_name",
                          "vmh_id", "cas_id", "inchi", "inchikey", "formula",
                          "smiles", "mass", "charge"):
                    if not t_attrs.get(k) and o_attrs.get(k):
                        t_attrs[k] = o_attrs[k]
                t_attrs.setdefault("alt_node_ids", [])
                if oid not in t_attrs["alt_node_ids"]:
                    t_attrs["alt_node_ids"].append(oid)
                # Merge synonyms
                syns_t = t_attrs.get("synonyms") or []
                syns_o = o_attrs.get("synonyms") or []
                if syns_t or syns_o:
                    seen = set(s.lower() if isinstance(s, str) else s for s in syns_t)
                    for s in syns_o:
                        if isinstance(s, str) and s.lower() not in seen:
                            syns_t.append(s); seen.add(s.lower())
                    t_attrs["synonyms"] = syns_t
                # Rewrite edges
                for pred in list(g.predecessors(oid)):
                    if pred == tid: continue
                    if not g.has_edge(pred, tid):
                        edata = g.get_edge_data(pred, oid) or {}
                        g.add_edge(pred, tid, **edata)
                        n_edges_rewritten += 1
                for succ in list(g.successors(oid)):
                    if succ == tid: continue
                    if not g.has_edge(tid, succ):
                        edata = g.get_edge_data(oid, succ) or {}
                        g.add_edge(tid, succ, **edata)
                        n_edges_rewritten += 1
                g.remove_node(oid)

            log.info("collapse_orphan_metab_twins: collapsed %d orphans "
                     "(by_inchikey14=%d, by_chebi=%d, by_name_exact=%d, "
                     "by_name_loose=%d), rewrote %d edges.",
                     len(orphan_to_twin),
                     strat["inchikey14"], strat["chebi_id"],
                     strat["name_exact"], strat["name_loose"],
                     n_edges_rewritten)
            return SourceRecord(
                name="collapse_orphan_metab_twins",
                license="internal",
                n_records=len(orphan_to_twin),
            )

        self._steps.append(_Step("collapse_orphan_metab_twins", _step))
        return self

    def mark_unannotatable_metabolites(
        self,
        metabolon_csv: str | Path | None = None,
    ) -> "BuildPipeline":
        """
        Mark Metabolon-derived metabolite nodes with no canonical structure
        as ``unannotatable=True``.

        ROOT CAUSE: Metabolon Tier-2/3 entries (e.g. "(11 or 12)-methyl-
        tridecanoate (a14:0 or i14:0)") report an ambiguous isomer MIX —
        no single canonical InChIKey exists by Metabolon's own annotation.
        Downstream analysis that tries to propagate signal from these
        nodes silently fails because they have no chemical identity.

        FIX: For each metabolite node matching a Metabolon biochemical
        whose dictionary entry has NEITHER InChIKey NOR PubChem CID,
        flag it ``unannotatable=True`` with reason
        ``metabolon_no_canonical_structure``. Downstream MAP reconstruction
        and graph signal computation should honor this flag and skip these
        nodes' contribution rather than attempting to propagate noise.

        Must run AFTER add_metabolon().
        """
        self._params["mark_unannotatable_metabolites"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            import csv, re
            csv_path = metabolon_csv
            if csv_path is None:
                from gizmo.resources import _project_root
                csv_path = (_project_root() / "data/resources/gizmo/sources"
                              / "metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv")
            csv_path = Path(csv_path)
            if not csv_path.exists():
                log.warning("mark_unannotatable_metabolites: dict CSV not found at %s; skipping",
                             csv_path)
                return SourceRecord(name="mark_unannotatable_metabolites",
                                     license="internal")

            def _norm(name): return re.sub(r"\s+", "_", name.strip().lower())

            by_norm = {}
            by_lower = {}
            with open(csv_path, encoding="utf-8-sig") as f:
                for r in csv.DictReader(f):
                    name = (r.get("BIOCHEMICAL") or "").strip()
                    if not name: continue
                    entry = {
                        "inchikey": (r.get("INCHIKEY") or "").strip() or None,
                        "pubchem_cid": (r.get("PUBCHEM") or "").strip() or None,
                    }
                    by_norm[_norm(name)] = entry
                    by_lower[name.lower()] = entry

            n_flagged = 0
            g = mg.graph
            for nid, attrs in g.nodes(data=True):
                if attrs.get("node_type") != "metabolite": continue
                if attrs.get("inchikey") or attrs.get("pubchem_cid"):
                    continue  # has structure, not unannotatable
                # Try to find a Metabolon dict entry for this node
                entry = None
                if nid.startswith("metabolon:"):
                    key = nid[len("metabolon:"):]
                    entry = by_norm.get(key)
                nm = (attrs.get("name") or "").strip()
                if not entry and nm:
                    entry = by_lower.get(nm.lower()) or by_norm.get(_norm(nm))
                # Flag if dict entry exists but has no canonical structure
                if entry and not entry["inchikey"] and not entry["pubchem_cid"]:
                    attrs["unannotatable"] = True
                    attrs["unannotatable_reason"] = "metabolon_no_canonical_structure"
                    n_flagged += 1

            log.info("mark_unannotatable_metabolites: flagged %d nodes "
                     "as Metabolon Tier-2/3 (no canonical structure)",
                     n_flagged)
            return SourceRecord(
                name="mark_unannotatable_metabolites",
                license="internal",
                n_records=n_flagged,
            )

        self._steps.append(_Step("mark_unannotatable_metabolites", _step))
        return self

    def enrich_metab_synonyms_from_parquet(
        self,
        parquet_path: str | Path,
        max_synonyms_per_node: int = 30,
    ) -> "BuildPipeline":
        """
        Populate ``attrs['synonyms']`` on metabolite nodes from a local
        PubChem CID→synonym parquet (long format, columns
        ``pubchem_cid``, ``synonym``).

        Used in lieu of the PUG-REST live fetch (see
        ``enrich_pubchem_synonyms`` in gizmo/sources/pubchem_synonyms.py)
        when a local cached parquet is available. Reads by row group to
        avoid loading 1.2GB in one shot.

        Synonyms become additional lookup keys in MetaboliteMapper's
        confidence ladder (exact-synonym 0.8 > substring 0.7 > fuzzy).
        This fixes name-resolution failures like "Adipic acid" not
        matching the Reactome node "adipate", and reduces wrong fuzzy
        matches like "2-aminoisobutyric acid" → "L-Leu" by giving the
        canonical-synonym path a higher-confidence match before fuzzy
        fallback fires.
        """
        self._params["enrich_metab_synonyms_from_parquet"] = {
            "parquet_path": str(parquet_path),
            "max_synonyms_per_node": int(max_synonyms_per_node),
        }

        def _step(mg: GizmoGraph) -> SourceRecord:
            p = Path(parquet_path)
            if not p.exists():
                log.warning("enrich_metab_synonyms_from_parquet: parquet not found at %s; skipping",
                             p)
                return SourceRecord(name="enrich_metab_synonyms_from_parquet",
                                     license="PubChem (NLM)")
            try:
                import pyarrow.parquet as pq
            except ImportError:
                log.warning("enrich_metab_synonyms_from_parquet: pyarrow not installed; skipping")
                return SourceRecord(name="enrich_metab_synonyms_from_parquet",
                                     license="PubChem (NLM)")

            from collections import defaultdict
            g = mg.graph
            cid_to_nodes = defaultdict(list)
            for nid, attrs in g.nodes(data=True):
                if attrs.get("node_type") != "metabolite": continue
                cid = attrs.get("pubchem_cid")
                if not cid: continue
                cid_s = str(cid).replace("CID:", "").replace("cid:", "").strip()
                if cid_s.endswith(".0"): cid_s = cid_s[:-2]
                if cid_s.isdigit():
                    cid_to_nodes[cid_s].append(nid)
            if not cid_to_nodes:
                return SourceRecord(name="enrich_metab_synonyms_from_parquet",
                                     license="PubChem (NLM)")

            target = set(cid_to_nodes.keys())
            cid_synonyms = defaultdict(list)
            f = pq.ParquetFile(p)
            for rg_i in range(f.num_row_groups):
                tbl = f.read_row_group(rg_i, columns=["pubchem_cid", "synonym"])
                cid_col = tbl["pubchem_cid"].to_pylist()
                syn_col = tbl["synonym"].to_pylist()
                for c, s in zip(cid_col, syn_col):
                    if c in target and s and len(cid_synonyms[c]) < max_synonyms_per_node:
                        cid_synonyms[c].append(s)

            n_nodes_enriched = 0
            n_total = 0
            for cid_s, nids in cid_to_nodes.items():
                syns = cid_synonyms.get(cid_s)
                if not syns: continue
                # Dedup case-insensitively
                seen = set()
                uniq = []
                for s in syns:
                    sl = s.casefold().strip()
                    if sl and sl not in seen:
                        seen.add(sl); uniq.append(s.strip())
                for nid in nids:
                    attrs = g.nodes[nid]
                    existing = attrs.get("synonyms") or []
                    existing_set = set(s.casefold() for s in existing if isinstance(s, str))
                    additions = [s for s in uniq if s.casefold() not in existing_set]
                    if additions:
                        attrs["synonyms"] = existing + additions
                        n_nodes_enriched += 1
                        n_total += len(additions)

            log.info("enrich_metab_synonyms_from_parquet: enriched %d nodes "
                     "with %d new synonyms (%.1f per node, from %d unique CIDs in source)",
                     n_nodes_enriched, n_total,
                     n_total / max(n_nodes_enriched, 1), len(cid_synonyms))
            return SourceRecord(
                name="enrich_metab_synonyms_from_parquet",
                url="https://pubchem.ncbi.nlm.nih.gov/",
                license="PubChem (NLM)",
                n_records=n_nodes_enriched,
            )

        self._steps.append(_Step("enrich_metab_synonyms_from_parquet", _step))
        return self

    def add_pathway_nodes(
        self,
        species: str = "Homo sapiens",
    ) -> "BuildPipeline":
        """
        Promote Reactome pathway stIDs (stored as attributes on reaction nodes)
        into first-class PathwayNode records with pathway→reaction edges.

        Call this after add_reactome() to make pathways queryable as nodes.
        """
        self._params["pathway_nodes"] = {"species": species}

        def _step(mg: GizmoGraph) -> SourceRecord:
            n = mg.promote_pathway_nodes(species=species)
            log.info("Pathway nodes: promoted %d pathway nodes.", n)
            return SourceRecord(
                name="pathway_nodes",
                version="current",
                url="https://reactome.org",
                license="CC BY 4.0",
                n_records=n,
            )

        self._steps.append(_Step("pathway_nodes", _step))
        return self

    def add_hpo(
        self,
        metabolic_only: bool = False,
    ) -> "BuildPipeline":
        """
        Add HPO phenotype nodes and phenotype→disease / phenotype→gene edges.

        Parameters
        ----------
        metabolic_only : if True, only load phenotype terms under
                         HP:0001939 "Abnormality of metabolism or catabolism"
        """
        self._params["hpo"] = {"metabolic_only": metabolic_only}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.hpo import HPOClient, build_mondo_id_map

            client = HPOClient(cache_dir=str(self.cache_dir / "hpo"))
            client.download()

            phenotypes = client.load_phenotypes()
            if metabolic_only:
                phenotypes = [p for p in phenotypes if p.is_metabolic]
            for p in phenotypes:
                mg.add_phenotype(p)

            gene_edges = client.load_gene_edges()
            for e in gene_edges:
                # Only add edge if target gene node exists in graph
                if e.target in mg.graph or e.target.startswith("symbol:"):
                    try:
                        mg.add_phenotype_edge(e)
                    except Exception:
                        pass

            mondo_map = build_mondo_id_map(mg)
            disease_edges = client.load_disease_edges(mondo_map)
            for e in disease_edges:
                if e.target in mg.graph or e.target.startswith(("OMIM:", "ORPHA:", "MONDO:")):
                    try:
                        mg.add_phenotype_edge(e)
                    except Exception:
                        pass

            log.info(
                "HPO: added %d phenotype nodes, %d gene edges, %d disease edges.",
                len(phenotypes), len(gene_edges), len(disease_edges),
            )
            return SourceRecord(
                name="hpo",
                version="current",
                url="https://hpo.jax.org",
                license="CC BY 4.0",
                n_records=len(phenotypes),
            )

        self._steps.append(_Step("hpo", _step))
        return self

    def add_gtex(
        self,
        min_tpm: float = 1.0,
    ) -> "BuildPipeline":
        """
        Enrich gene nodes with GTEx tissue expression data.

        Parameters
        ----------
        min_tpm : minimum median TPM to store for a tissue (0 = all non-zero)
        """
        self._params["gtex"] = {"min_tpm": min_tpm}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.gtex import GTExClient

            client = GTExClient(cache_dir=str(self.cache_dir / "gtex"))
            client.download()
            n = client.enrich_graph(mg, min_tpm=min_tpm)
            return SourceRecord(
                name="gtex",
                version="v8",
                url="https://gtexportal.org",
                license="dbGaP phs000424",
                n_records=n,
            )

        self._steps.append(_Step("gtex", _step))
        return self

    def add_clinvar(
        self,
        min_stars: int = 1,
        assembly: str = "GRCh38",
    ) -> "BuildPipeline":
        """
        Add ClinVar pathogenic/likely-pathogenic variant nodes.

        Parameters
        ----------
        min_stars : minimum ClinVar review-star rating (0–4)
        assembly  : genome assembly to filter ("GRCh38" | "GRCh37" | None)
        """
        self._params["clinvar"] = {"min_stars": min_stars, "assembly": assembly}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.clinvar import ClinVarClient

            client = ClinVarClient(cache_dir=str(self.cache_dir / "clinvar"))
            client.download()
            nodes, edges = client.load_pathogenic(min_stars=min_stars, assembly=assembly)
            for node in nodes:
                mg.add_variant(node)
            for edge in edges:
                try:
                    mg.add_variant_edge(edge)
                except Exception:
                    pass

            log.info("ClinVar: added %d variant nodes, %d edges.", len(nodes), len(edges))
            return SourceRecord(
                name="clinvar",
                version="current",
                url="https://www.ncbi.nlm.nih.gov/clinvar/",
                license="public domain",
                n_records=len(nodes),
            )

        self._steps.append(_Step("clinvar", _step))
        return self

    def add_drugs(
        self,
        min_phase: int = 2,
    ) -> "BuildPipeline":
        """
        Add first-class drug nodes and drug→gene edges from ChEMBL.

        Parameters
        ----------
        min_phase : minimum clinical phase (2 = Phase 2+ / approved)
        """
        self._params["drugs"] = {"min_phase": min_phase}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.chembl_drugs import ChEMBLDrugClient

            client = ChEMBLDrugClient(cache_dir=str(self.cache_dir / "chembl_drugs"))
            nodes, edges = client.load_for_graph(mg, min_phase=min_phase)
            for node in nodes:
                mg.add_drug(node)
            for edge in edges:
                try:
                    mg.add_drug_edge(edge)
                except Exception:
                    pass

            log.info("Drugs: added %d drug nodes, %d drug-gene edges.", len(nodes), len(edges))
            return SourceRecord(
                name="chembl_drugs",
                version="current",
                url="https://www.ebi.ac.uk/chembl/",
                license="CC BY-SA 3.0",
                n_records=len(nodes),
            )

        self._steps.append(_Step("drugs", _step))
        return self

    # ------------------------------------------------------------------
    # Toxicology sources
    # ------------------------------------------------------------------

    def add_ctd(
        self,
        direct_evidence_only: bool = True,
        organisms: Optional[set[str]] = None,
        chem_xref_path: Optional[str | Path] = None,
    ) -> "BuildPipeline":
        """
        Add CTD chemical–gene and chemical–disease toxicology edges.

        Parameters
        ----------
        direct_evidence_only : only include curated CTD evidence (not inferred)
        organisms            : organism filter (default: {"Homo sapiens"})
        chem_xref_path       : MetaNetX chem_xref.tsv for CAS→ChEBI mapping
                               (defaults to cache_dir/metanetx/chem_xref.tsv)
        """
        self._params["ctd"] = {
            "direct_evidence_only": direct_evidence_only,
            "organisms": list(organisms) if organisms else ["Homo sapiens"],
        }

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.ctd import CTDClient

            xref = (
                Path(chem_xref_path)
                if chem_xref_path
                else self.cache_dir / "metanetx" / "chem_xref.tsv"
            )

            client = CTDClient(
                cache_dir=str(self.cache_dir / "ctd"),
                organisms=organisms,
            )
            client.download()
            client.build_cas_index(mg, chem_xref_path=xref if xref.exists() else None)

            gene_edges    = client.load_chem_gene_interactions(mg, direct_evidence_only)
            disease_edges = client.load_chem_disease_associations(mg, direct_evidence_only)

            mg.add_tox_edges(gene_edges)
            mg.add_tox_edges(disease_edges)

            n = len(gene_edges) + len(disease_edges)
            log.info("CTD: added %d tox edges (%d gene, %d disease).",
                     n, len(gene_edges), len(disease_edges))
            return SourceRecord(
                name="ctd",
                version="current",
                url="https://ctdbase.org",
                license="non-commercial",
                n_records=n,
            )

        self._steps.append(_Step("ctd", _step))
        return self

    def add_comptox(
        self,
        force: bool = False,
    ) -> "BuildPipeline":
        """
        Enrich metabolite nodes with EPA CompTox/DSSTox data.

        Adds ``cas_id``, ``dtxsid``, and ``epa_hazard_flags`` attributes to
        metabolite nodes that have an InChIKey.  Uses the CCTE REST API.

        Parameters
        ----------
        force : re-query even if cas_id is already set
        """
        self._params["comptox"] = {"force": force}

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.comptox import CompToxClient

            client = CompToxClient(cache_dir=str(self.cache_dir / "comptox"))
            n = client.enrich_graph(mg, force=force)
            return SourceRecord(
                name="comptox",
                version="current",
                url="https://comptox.epa.gov",
                license="CC BY 4.0",
                n_records=n,
            )

        self._steps.append(_Step("comptox", _step))
        return self

    def add_t3db(self) -> "BuildPipeline":
        """
        Add T3DB (Toxic Exposome Database) toxin–target edges.

        Toxins are matched to existing metabolite nodes via HMDB, ChEBI,
        CAS, and PubChem IDs.  Running add_comptox() first improves
        CAS coverage and therefore T3DB match rate.
        """
        self._params["t3db"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.t3db import T3DBClient

            client = T3DBClient(cache_dir=str(self.cache_dir / "t3db"))
            client.download()
            edges = client.load_tox_edges(mg)
            mg.add_tox_edges(edges)

            log.info("T3DB: added %d tox-gene edges.", len(edges))
            return SourceRecord(
                name="t3db",
                version="current",
                url="https://t3db.ca",
                license="free academic",
                n_records=len(edges),
            )

        self._steps.append(_Step("t3db", _step))
        return self

    def add_chembl_tox(self) -> "BuildPipeline":
        """
        Add ChEMBL ADMET / toxicology assay annotations.

        For metabolite nodes that have a ``chembl_id`` attribute (populated by
        add_drugs()), queries the ChEMBL ADMET assay endpoint and:
          - Adds hERG inhibition as a ToxEdge → KCNH2/hERG gene node
          - Annotates nodes with other ADMET values (LD50, AMES, etc.)

        Should be called after add_drugs().
        """
        self._params["chembl_tox"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.chembl import ChEMBLClient

            client = ChEMBLClient(cache_dir=str(self.cache_dir / "chembl"))
            tox_edges = client.tox_edges_for_graph(mg)
            mg.add_tox_edges([e for e in tox_edges if hasattr(e, "edge_type")])

            log.info("ChEMBL tox: added %d tox edges.", len(tox_edges))
            return SourceRecord(
                name="chembl_tox",
                version="current",
                url="https://www.ebi.ac.uk/chembl/",
                license="CC BY-SA 3.0",
                n_records=len(tox_edges),
            )

        self._steps.append(_Step("chembl_tox", _step))
        return self

    def add_hmdb_enrichment(
        self,
        metanetx_xref: str | Path | None = None,
    ) -> "BuildPipeline":
        """
        Back-fill VMH-compatible HMDB IDs for metabolite nodes that have a ChEBI
        ID but no hmdb_id (typically all Reactome-sourced nodes).

        Should be called after add_reactome() and optionally add_metabolon().
        Uses MetaNetX chem_xref ChEBI→MNX→HMDB mapping.

        Parameters
        ----------
        metanetx_xref : path to chem_xref.tsv; defaults to <cache_dir>/metanetx/chem_xref.tsv
        """
        self._params["hmdb_enrichment"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.sources.metanetx import MetaNetXClient
            xref_path = (
                Path(metanetx_xref) if metanetx_xref
                else self.cache_dir / "metanetx" / "chem_xref.tsv"
            )
            if not xref_path.exists():
                log.warning(
                    "HMDB enrichment skipped — chem_xref.tsv not found at %s", xref_path
                )
                return SourceRecord(name="hmdb_enrichment", license="CC BY 4.0 (MetaNetX)")
            mnx = MetaNetXClient(cache_dir=xref_path.parent)
            n = mnx.enrich_graph_hmdb(mg)
            log.info("HMDB enrichment: back-filled %d metabolite nodes.", n)
            return SourceRecord(
                name="hmdb_enrichment",
                url="https://www.metanetx.org",
                license="CC BY 4.0",
                n_records=n,
            )

        self._steps.append(_Step("hmdb_enrichment", _step))
        return self

    def add_currency_flags(self) -> "BuildPipeline":
        """Flag currency metabolites (ATP, NAD+, water, etc.) in-place."""
        self._params["currency_flags"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.analysis.currency import flag_currency_metabolites
            n = flag_currency_metabolites(mg)
            log.info("Currency: flagged %d metabolites.", n)
            return SourceRecord(name="currency_flags", license="internal")

        self._steps.append(_Step("currency_flags", _step))
        return self

    def link_genes_to_reactions(self) -> "BuildPipeline":
        """
        Wire gene nodes to reaction nodes via GENE_REACTION edges.

        Reactome stores gene symbols as an attribute on reaction nodes
        (``gene_symbols``), while Open Targets / Orphanet create separate
        gene nodes identified by Ensembl ID with a ``symbol`` attribute.
        This step cross-references the two and adds directed
        gene → reaction edges so the full path
        disease → gene → reaction → metabolite is navigable in the graph.

        Call after both ``add_reactome()`` and any gene-adding step
        (``add_open_targets``, ``add_orphanet``).
        """
        self._params["link_genes_to_reactions"] = True

        def _step(mg: GizmoGraph) -> SourceRecord:
            from gizmo.schema import DiseaseEdgeType

            # Build symbol → node_id index from existing gene nodes
            symbol_to_nid: dict[str, str] = {}
            for nid, attrs in mg.graph.nodes(data=True):
                if attrs.get("node_type") == "gene":
                    sym = attrs.get("symbol", "")
                    if sym:
                        symbol_to_nid[sym.upper()] = nid

            if not symbol_to_nid:
                log.warning("link_genes_to_reactions: no gene nodes found; skipping.")
                return SourceRecord(name="link_genes_to_reactions", license="internal")

            n_added = 0
            for rxn_id, attrs in list(mg.graph.nodes(data=True)):
                if attrs.get("node_type") != "reaction":
                    continue
                for sym in attrs.get("gene_symbols") or []:
                    gene_nid = symbol_to_nid.get(sym.upper())
                    if gene_nid and not mg.graph.has_edge(gene_nid, rxn_id):
                        mg.graph.add_edge(
                            gene_nid, rxn_id,
                            edge_type=DiseaseEdgeType.GENE_REACTION.value,
                            gene_symbol=sym,
                        )
                        n_added += 1

            log.info("link_genes_to_reactions: added %d gene→reaction edges.", n_added)
            return SourceRecord(name="link_genes_to_reactions", license="internal")

        self._steps.append(_Step("link_genes_to_reactions", _step))
        return self

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self) -> tuple[GizmoGraph, GraphManifest]:
        """
        Execute all registered steps and return (GizmoGraph, GraphManifest).

        Steps run in registration order.  Errors in individual steps are
        logged and skipped so that a partial build still completes.
        """
        if not self._steps:
            raise RuntimeError("No build steps registered.  Call add_reactome() etc. first.")

        mg = GizmoGraph()
        sources: list[SourceRecord] = []

        import time as _time
        total_steps = sum(1 for s in self._steps if s.enabled)
        for si, step in enumerate(s for s in self._steps if s.enabled):
            log.info("[%d/%d] Build step: %s …", si + 1, total_steps, step.name)
            t0 = _time.time()
            try:
                rec = step.fn(mg)
                if rec is not None:
                    sources.append(rec)
                log.info("[%d/%d] %s done (%.1fs) — graph now %d nodes, %d edges",
                         si + 1, total_steps, step.name, _time.time() - t0,
                         mg.graph.number_of_nodes(), mg.graph.number_of_edges())
            except Exception as exc:
                log.error("[%d/%d] Build step '%s' failed (%.1fs): %s",
                          si + 1, total_steps, step.name, _time.time() - t0,
                          exc, exc_info=True)

        manifest = GraphManifest.from_graph(
            mg,
            graph_name=self.graph_name,
            build_params=self._params,
            sources=sources,
            notes=self.notes,
        )

        log.info(
            "Build complete: %d nodes, %d edges.",
            mg.graph.number_of_nodes(),
            mg.graph.number_of_edges(),
        )
        return mg, manifest

    # ------------------------------------------------------------------
    # Bundle output
    # ------------------------------------------------------------------

    def save_bundle(
        self,
        mg: GizmoGraph,
        manifest: GraphManifest,
        output_dir: str | Path = "data/processed",
    ) -> Path:
        """
        Write the graph bundle to ``output_dir/<graph_name>/``.

        Bundle contents:
            graph.json       — node-link JSON
            graph.graphml    — GATOM/igraph compatible GraphML
            qc_report.json   — full QC report
            manifest.json    — provenance manifest

        Returns the bundle directory path.
        """
        bundle_dir = Path(output_dir) / self.graph_name
        bundle_dir.mkdir(parents=True, exist_ok=True)

        # graph.json
        from gizmo.export.json_export import write_json
        write_json(mg, bundle_dir / "graph.json")
        log.info("Saved graph.json")

        # graph.graphml
        try:
            from gizmo.export.graphml import write_graphml
            write_graphml(mg, bundle_dir / "graph.graphml")
            log.info("Saved graph.graphml")
        except Exception as exc:
            log.warning("GraphML export failed: %s", exc)

        # qc_report.json
        try:
            from gizmo.analysis.qc import assess_readiness
            r = assess_readiness(mg)
            import dataclasses
            qc_path = bundle_dir / "qc_report.json"
            qc_path.write_text(json.dumps(dataclasses.asdict(r), indent=2, default=str))
            log.info("Saved qc_report.json")
        except Exception as exc:
            log.warning("QC report failed: %s", exc, exc_info=True)

        # manifest.json
        manifest.save(bundle_dir / "manifest.json")
        log.info("Saved manifest.json")

        log.info("Bundle written to %s", bundle_dir)
        return bundle_dir
