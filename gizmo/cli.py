"""
GIZMO command-line interface.

Usage:
  gizmo build  --name human_iem [--pathways R-HSA-XXX] [--metabolon-csv …]
               [--hpo] [--gtex] [--clinvar] [--drugs]
               [--ctd] [--comptox] [--t3db] [--chembl-tox]
               [--chemical-enrichment] [--pathway-nodes] [--link-genes]
  gizmo qc     --graph graph.json
  gizmo score  --graph graph.json --metabolomics met.csv [--transcriptomics rna.csv]
  gizmo app    --graph graph.json [--port 8050]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gizmo.resources import (
    metabolon_csv_default,
    metanetx_prop_default,
    metanetx_xref_default,
    overrides_default,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gizmo",
        description="Graph-Integrated Zone of Metabolite Operations",
    )
    sub = parser.add_subparsers(dest="command")

    build_p = sub.add_parser(
        "build", help="Build and save a reproducible graph bundle"
    )
    build_p.add_argument("--name", required=True, help="Graph bundle name (used as output directory)")
    build_p.add_argument("--output-dir", default="data/processed",
                         help="Parent directory for the graph bundle (default: data/processed)")
    build_p.add_argument("--cache-dir",  default="data/raw",
                         help="Raw data cache root (default: data/raw)")
    build_p.add_argument("--species",  default="Homo sapiens",
                         help="Species for full Reactome load")
    build_p.add_argument("--pathways", default=None,
                         help="Comma-separated Reactome stIds (overrides --species)")
    build_p.add_argument("--no-reactome", action="store_true")
    build_p.add_argument("--no-mondo",        action="store_true")
    build_p.add_argument("--no-orphanet",     action="store_true")
    build_p.add_argument("--no-open-targets", action="store_true")
    build_p.add_argument("--metabolon-csv",  default=None)
    build_p.add_argument("--metanetx-prop",  default=metanetx_prop_default())
    build_p.add_argument("--metanetx-xref",  default=metanetx_xref_default())
    build_p.add_argument("--overrides",      default=overrides_default(),
                         help="Metabolon curation overrides JSON")
    build_p.add_argument("--no-stringdb",    action="store_true")
    build_p.add_argument("--string-min-score", type=float, default=0.4)
    # Graph enrichment
    build_p.add_argument("--chemical-enrichment", action="store_true",
                         help="Enrich metabolite nodes with MetaNetX structural data (SMILES, formula, mass)")
    build_p.add_argument("--pathway-nodes", action="store_true",
                         help="Promote Reactome pathway stIDs to first-class graph nodes")
    build_p.add_argument("--link-genes",    action="store_true",
                         help="Wire gene nodes to reactions via GENE_REACTION edges")
    build_p.add_argument("--hpo",           action="store_true",
                         help="Add HPO phenotype nodes and phenotype–disease/gene edges")
    build_p.add_argument("--hpo-metabolic-only", action="store_true",
                         help="With --hpo: restrict to metabolic phenotype terms (HP:0001939)")
    build_p.add_argument("--gtex",          action="store_true",
                         help="Enrich gene nodes with GTEx tissue expression (min TPM 1.0)")
    build_p.add_argument("--gtex-min-tpm",  type=float, default=1.0)
    build_p.add_argument("--clinvar",       action="store_true",
                         help="Add ClinVar pathogenic/likely-pathogenic variant nodes")
    build_p.add_argument("--clinvar-min-stars", type=int, default=1)
    build_p.add_argument("--drugs",         action="store_true",
                         help="Add ChEMBL drug nodes and drug→gene edges (Phase 2+)")
    build_p.add_argument("--drugs-min-phase", type=int, default=2)
    # Toxicology sources
    build_p.add_argument("--ctd",           action="store_true",
                         help="Add CTD chemical–gene/disease tox edges")
    build_p.add_argument("--ctd-all-evidence", action="store_true",
                         help="Include inferred (not just curated) CTD evidence")
    build_p.add_argument("--comptox",       action="store_true",
                         help="Enrich metabolites with EPA CompTox CAS/hazard data")
    build_p.add_argument("--t3db",          action="store_true",
                         help="Add T3DB toxin–target edges")
    build_p.add_argument("--chembl-tox",    action="store_true",
                         help="Add ChEMBL ADMET/hERG tox assay annotations")
    build_p.add_argument("--notes", default="")

    qc_p = sub.add_parser("qc", help="Run computational readiness report on a saved graph")
    qc_p.add_argument("--graph", required=True, help="Path to graph JSON file")

    score_p = sub.add_parser("score", help="Score reactions from multiomic evidence")
    score_p.add_argument("--graph",          required=True, help="Path to graph JSON")
    score_p.add_argument("--metabolomics",   default=None,
                         help="Metabolomics CSV (FEATURE_ID, EFFECT_SIZE, P_VALUE, FDR)")
    score_p.add_argument("--transcriptomics",default=None,
                         help="Transcriptomics CSV (FEATURE_ID, EFFECT_SIZE, P_VALUE, FDR)")
    score_p.add_argument("--sample-id",      default=None)
    score_p.add_argument("--top-n",          type=int, default=20,
                         help="Show top N ranked reactions (default 20)")
    score_p.add_argument("--output",         default=None,
                         help="Save ranked reactions to this JSON file")
    score_p.add_argument("--min-evidence",   type=int, default=1,
                         help="Min mapped evidence features per reaction (default 1)")
    score_p.add_argument("--pathways",       action="store_true",
                         help="Also print a pathway-level summary")
    score_p.add_argument("--chains",         action="store_true",
                         help="Also print top causal chain hypotheses")
    score_p.add_argument("--propagation",    action="store_true",
                         help="Include metabolite→reaction→metabolite chains (with --chains)")
    score_p.add_argument("--actionability",  action="store_true",
                         help="Score druggability (ChEMBL) + perturbability for each reaction")
    score_p.add_argument("--chembl-cache",   default=None,
                         help="Directory for ChEMBL response cache (default: data/raw/chembl)")

    app_p = sub.add_parser("app", help="Launch the GIZMO Dash explorer")
    app_p.add_argument("--graph", required=True, help="Path to graph JSON file")
    app_p.add_argument("--port", type=int, default=8050)
    app_p.add_argument("--host", default="127.0.0.1")
    app_p.add_argument("--reactome-cache", default="data/raw/reactome")
    app_p.add_argument("--metabolon-csv", default=metabolon_csv_default())
    app_p.add_argument("--metanetx-prop", default=metanetx_prop_default())
    app_p.add_argument("--metanetx-xref", default=metanetx_xref_default())
    app_p.add_argument("--overrides",     default=overrides_default())
    app_p.add_argument("--max-nodes",     type=int, default=500)
    app_p.add_argument("--debug",         action="store_true")

    met_p = sub.add_parser("metabolon", help="Report ChEBI coverage for a Metabolon CSV")
    met_p.add_argument("--csv",          required=True)
    met_p.add_argument("--metanetx-prop", default=metanetx_prop_default())
    met_p.add_argument("--metanetx-xref", default=metanetx_xref_default())

    args = parser.parse_args(argv)

    if args.command == "build":
        import logging
        from pathlib import Path
        logging.basicConfig(level=logging.INFO,
                            format="%(asctime)s  %(levelname)-8s  %(message)s")

        from gizmo.build.pipeline import BuildPipeline

        pipe = BuildPipeline(
            graph_name=args.name,
            cache_dir=args.cache_dir,
            notes=args.notes,
        )

        if not args.no_reactome:
            stids = [s.strip() for s in args.pathways.split(",")] if args.pathways else None
            pipe.add_reactome(
                species=args.species,
                pathway_stids=stids,
            )

        if not args.no_mondo:
            pipe.add_mondo_iem()

        if not args.no_orphanet:
            pipe.add_orphanet()

        if not args.no_open_targets:
            pipe.add_open_targets()

        metabolon_csv = args.metabolon_csv
        if not metabolon_csv:
            discovered_csv = Path(metabolon_csv_default())
            if discovered_csv.exists():
                metabolon_csv = str(discovered_csv)

        if metabolon_csv:
            pipe.add_metabolon(
                csv_path=metabolon_csv,
                metanetx_prop=args.metanetx_prop,
                metanetx_xref=args.metanetx_xref,
                overrides_path=args.overrides,
            )

        if not args.no_stringdb:
            pipe.add_stringdb(min_score=args.string_min_score)

        if args.link_genes:
            pipe.link_genes_to_reactions()

        # Back-fill VMH-compatible HMDB IDs for all ChEBI-mapped metabolite nodes
        pipe.add_hmdb_enrichment(metanetx_xref=args.metanetx_xref)

        if args.chemical_enrichment:
            pipe.add_chemical_enrichment(
                metanetx_prop=args.metanetx_prop,
                metanetx_xref=args.metanetx_xref,
            )

        if args.pathway_nodes:
            pipe.add_pathway_nodes(species=args.species)

        if args.hpo:
            pipe.add_hpo(metabolic_only=args.hpo_metabolic_only)

        if args.gtex:
            pipe.add_gtex(min_tpm=args.gtex_min_tpm)

        if args.clinvar:
            pipe.add_clinvar(min_stars=args.clinvar_min_stars)

        if args.drugs:
            pipe.add_drugs(min_phase=args.drugs_min_phase)

        # Toxicology — recommended order: comptox (CAS enrichment) → ctd → t3db → chembl_tox
        if args.comptox:
            pipe.add_comptox()

        if args.ctd:
            pipe.add_ctd(
                direct_evidence_only=not args.ctd_all_evidence,
                chem_xref_path=args.metanetx_xref,
            )

        if args.t3db:
            pipe.add_t3db()

        if args.chembl_tox:
            pipe.add_chembl_tox()

        pipe.add_currency_flags()

        mg, manifest = pipe.run()
        bundle_dir = pipe.save_bundle(mg, manifest, output_dir=args.output_dir)
        print(f"\nBundle saved to: {bundle_dir}")
        manifest.print_summary()

    elif args.command == "qc":
        from gizmo.export.json_export import read_json
        from gizmo.analysis.qc import assess_readiness

        mg = read_json(args.graph)
        report = assess_readiness(mg)
        report.print_summary()

    elif args.command == "score":
        import json as _json
        from gizmo.export.json_export import read_json
        from gizmo.evidence.ingest import load_metabolomics_csv, load_transcriptomics_csv
        from gizmo.evidence.model import SampleContext
        from gizmo.scoring.reaction_scorer import score_reactions

        print(f"Loading graph from {args.graph}…")
        mg = read_json(args.graph)

        ctx = SampleContext(sample_id=args.sample_id or "sample")

        if args.metabolomics:
            print(f"Loading metabolomics from {args.metabolomics}…")
            load_metabolomics_csv(args.metabolomics, mg, ctx)
        if args.transcriptomics:
            print(f"Loading transcriptomics from {args.transcriptomics}…")
            load_transcriptomics_csv(args.transcriptomics, mg, ctx)

        s = ctx.summary()
        print(f"Evidence: {s['n_metabolomics']} metabolomics  "
              f"{s['n_transcriptomics']} transcriptomics  "
              f"({s['n_mapped_nodes']} nodes mapped)")

        scores = score_reactions(mg, ctx, min_evidence=args.min_evidence)
        scores.sort(key=lambda r: abs(r.score), reverse=True)

        print(f"\nTop {args.top_n} reactions by |score|:\n")
        print(f"{'Rank':<5} {'Score':>8}  {'Dir':>5}  {'Evid':>5}  {'Reaction'}")
        print("-" * 72)
        for i, r in enumerate(scores[:args.top_n], 1):
            direction = "▲" if r.direction > 0.1 else ("▼" if r.direction < -0.1 else "~")
            print(f"{i:<5} {r.score:>8.3f}  {direction:>5}  {r.evidence_count:>5}  {r.reaction_id}")

        if args.pathways:
            from gizmo.scoring.pathway_scorer import summarise_pathways, print_pathway_report
            summaries = summarise_pathways(scores)
            print_pathway_report(summaries, top_n=args.top_n)

        if args.chains:
            from gizmo.scoring.chain_ranker import rank_chains, print_chain_report
            chains = rank_chains(
                mg, scores, ctx=ctx,
                include_metabolite_propagation=args.propagation,
            )
            print(f"\nTop {args.top_n} causal chain hypotheses ({len(chains)} total):")
            print_chain_report(chains, top_n=args.top_n)

        if args.actionability:
            from gizmo.actionability import (
                score_druggability, score_perturbability,
                combine_actionability, print_actionability_report,
            )
            from gizmo.sources.chembl import ChEMBLClient

            cache_dir = args.chembl_cache or "data/raw/chembl"
            chembl = ChEMBLClient(cache_dir=cache_dir)
            print("\nScoring druggability via ChEMBL (network calls for uncached genes)…")
            drug_scores = score_druggability(mg, chembl)
            pert_scores = score_perturbability(mg, scores)
            combined    = combine_actionability(drug_scores, pert_scores, scores)
            print_actionability_report(combined, top_n=args.top_n)

        if args.output:
            import dataclasses
            out = [dataclasses.asdict(r) for r in scores]
            _json.dump(out, open(args.output, "w"), indent=2, default=str)
            print(f"\nSaved {len(scores)} scored reactions → {args.output}")

    elif args.command == "app":
        from pathlib import Path
        from gizmo.export.json_export import read_json
        from gizmo.app.dash_app import create_app

        print(f"Loading graph from {args.graph}…")
        mg = read_json(args.graph)

        met_loader = None
        if args.metabolon_csv and Path(args.metabolon_csv).exists():
            from gizmo.sources.metabolon import MetabolonLoader
            from gizmo.curation.metabolon_curator import MetabolonCurator
            print(f"Loading Metabolon CSV from {args.metabolon_csv}…")
            met_loader = MetabolonLoader(args.metabolon_csv)
            if Path(args.metanetx_prop).exists() and Path(args.metanetx_xref).exists():
                print("Building MetaNetX lookup index (streaming, may take ~30-60 s)…")
                met_loader.load_metanetx_index(args.metanetx_prop, args.metanetx_xref)
            curator = MetabolonCurator(met_loader, graph=mg, overrides_path=args.overrides)
            curator.apply()
        else:
            curator = None

        app = create_app(
            mg,
            met_loader=met_loader,
            curator=curator,
            reactome_cache_dir=args.reactome_cache,
            max_nodes=args.max_nodes,
        )
        app.run(debug=args.debug, host=args.host, port=args.port)

    elif args.command == "metabolon":
        from gizmo.sources.metabolon import MetabolonLoader
        print(f"Loading Metabolon CSV: {args.csv}")
        loader = MetabolonLoader(args.csv)
        print("Building MetaNetX index…")
        loader.load_metanetx_index(args.metanetx_prop, args.metanetx_xref)
        rep = loader.coverage_report()
        print(f"Total compounds: {rep['total']}")
        print(f"Mapped to ChEBI: {rep['mapped']} ({rep['coverage']:.1f}%)")
        print("\nTop unmatched compounds:")
        for name, n in rep["top_unmatched"][:20]:
            print(f"  {name} ({n})")

    else:
        parser.print_help(sys.stderr)
        return 1

    return 0
