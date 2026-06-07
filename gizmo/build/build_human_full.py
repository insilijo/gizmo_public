"""Build the canonical multi-source human_full GIZMO graph bundle.

This is the production build entry point. Promoted from
benchmarks/build_full_human_graph.py.

Sources integrated:
  - Reactome reactions + pathway hierarchy (cached)
  - Mondo all-disease ontology (cached)
  - Orphanet rare-disease catalog (cached)
  - Open Targets gene-disease associations (queries existing MONDO IDs)
  - StringDB PPI at min_score=0.7 (downloads ~50MB on first run)
  - Metabolon HD4 compounds (4,851 compounds via local CSV)
  - PubChem synonym enrichment (cached)
  - MetaNetX chemical enrichment
  - Reactome gene-symbol enrichment (~12k catalysis edges)
  - Conditional currency edges (EC-class + pathway overrides)

Output: data/processed/human_full/{graph.json, graph.graphml,
                                    manifest.json, qc_report.json}

Usage::

    from gizmo.build.build_human_full import build_human_full
    mg, manifest = build_human_full()

CLI::

    python -m gizmo.build.build_human_full --output data/processed
"""
from __future__ import annotations
import sys, time, argparse, logging, hashlib, json, os
from pathlib import Path
from datetime import datetime, timezone

log = logging.getLogger(__name__)


def _read_bundle_version() -> str:
    """Bundle semver, sourced from gizmo/build/_VERSION (so it can be bumped
    without editing code)."""
    f = Path(__file__).parent / "_VERSION"
    if f.exists():
        return f.read_text().strip()
    return "0.0.0-dev"


BUNDLE_VERSION = _read_bundle_version()


def bump_version(part: str = "patch") -> str:
    """Bump the version file by part ∈ {major, minor, patch}.
    Writes to _VERSION and returns the new version string."""
    parts = BUNDLE_VERSION.split("-")[0].split(".")
    if len(parts) != 3:
        raise ValueError(f"Invalid version: {BUNDLE_VERSION}")
    major, minor, patch = (int(x) for x in parts)
    if part == "major": major, minor, patch = major + 1, 0, 0
    elif part == "minor": minor, patch = minor + 1, 0
    elif part == "patch": patch += 1
    else: raise ValueError(f"Invalid part: {part}")
    new_v = f"{major}.{minor}.{patch}"
    (Path(__file__).parent / "_VERSION").write_text(new_v + "\n")
    return new_v


def _sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            buf = f.read(chunk)
            if not buf: break
            h.update(buf)
    return h.hexdigest()


def _stamp_version_metadata(bundle_dir: Path, version: str):
    """Append version + content hashes to manifest.json after save_bundle."""
    manifest_path = bundle_dir / "manifest.json"
    if not manifest_path.exists(): return
    manifest = json.loads(manifest_path.read_text())
    manifest["bundle_version"] = version
    manifest["versioned_at"] = datetime.now(timezone.utc).isoformat()
    # Content hashes for graph artifacts
    hashes = {}
    for fname in ("graph.json", "graph.graphml"):
        fpath = bundle_dir / fname
        if fpath.exists():
            hashes[fname] = {
                "sha256": _sha256_file(fpath),
                "size_bytes": fpath.stat().st_size,
            }
    manifest["artifact_hashes"] = hashes
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str))


def build_human_full(
    *,
    cache_dir: str | Path = "data/raw",
    output_dir: str | Path = "data/processed",
    metabolon_csv: str | Path | None = None,
    stringdb_min_score: float = 0.7,
    open_targets_min_score: float = 0.1,
    bundle_name: str = "human_full",
):
    """Build human_full and save the bundle.

    Returns (GizmoGraph, GraphManifest).
    """
    from gizmo.build.pipeline import BuildPipeline

    cache_dir = Path(cache_dir)
    pipe = BuildPipeline(bundle_name, cache_dir=str(cache_dir))

    log.info("[1/12] Reactome reactions")
    pipe.add_reactome(species="Homo sapiens")

    log.info("[2/12] Reactome pathway hierarchy")
    pipe.add_pathway_nodes(species="Homo sapiens")

    log.info("[3/12] Metabolon HD4 compounds")
    if metabolon_csv is None:
        # Default to repo's bundled CSV
        from gizmo.resources import _project_root
        metabolon_csv = (_project_root() / "data/resources/gizmo/sources"
                          / "metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv")
    if Path(metabolon_csv).exists():
        chem_prop = cache_dir / "metanetx" / "chem_prop.tsv"
        chem_xref = cache_dir / "metanetx" / "chem_xref.tsv"
        kw = {"csv_path": str(metabolon_csv)}
        if chem_prop.exists() and chem_xref.exists():
            kw["metanetx_prop"] = str(chem_prop)
            kw["metanetx_xref"] = str(chem_xref)
        pipe.add_metabolon(**kw)
    else:
        log.warning(f"  Metabolon CSV not found at {metabolon_csv}, skipping")

    log.info("[4/12] Mondo all-disease ontology")
    pipe.add_mondo_all()

    log.info("[5/12] Orphanet rare-disease catalog")
    pipe.add_orphanet(iem_only=False)

    log.info("[6/12] Open Targets gene-disease")
    pipe.add_open_targets(disease_ids=None, min_score=open_targets_min_score)

    log.info("[7/12] StringDB PPI")
    pipe.add_stringdb(min_score=stringdb_min_score)

    log.info("[8/12] Chemical enrichment")
    chem_prop = cache_dir / "metanetx" / "chem_prop.tsv"
    chem_xref = cache_dir / "metanetx" / "chem_xref.tsv"
    if chem_prop.exists() and chem_xref.exists():
        pipe.add_chemical_enrichment(chem_prop_path=str(chem_prop),
                                      chem_xref_path=str(chem_xref))

    log.info("[8a/12] Mark Metabolon Tier-2/3 isomer mixes as unannotatable")
    pipe.mark_unannotatable_metabolites(metabolon_csv=metabolon_csv)

    log.info("[8b/12] Enrich metabolite synonyms from PubChem parquet (if available)")
    pubchem_syn_parquet = (cache_dir / "pubchem" / "pubchem_cid_synonym.parquet")
    if not pubchem_syn_parquet.exists():
        # Fall back to SQuID-INC's cached parquet if available
        squid_inc_parquet = Path("/home/jgardner/SQuID-INC/data/processed/parquet/raw_pubchem/pubchem_cid_synonym.parquet")
        if squid_inc_parquet.exists():
            pubchem_syn_parquet = squid_inc_parquet
    if pubchem_syn_parquet.exists():
        pipe.enrich_metab_synonyms_from_parquet(str(pubchem_syn_parquet))
    else:
        log.warning(f"  PubChem synonym parquet not found; skipping "
                    f"(checked {cache_dir / 'pubchem' / 'pubchem_cid_synonym.parquet'} "
                    f"and SQuID-INC cache)")

    log.info("[8c/12] Collapse orphan metabolite twins into Reactome connected equivalents")
    pipe.collapse_orphan_metab_twins()

    log.info("[9/12] link_genes_to_reactions (OT/Orphanet → reactions)")
    pipe.link_genes_to_reactions()

    log.info("[10/12] Currency flags")
    pipe.add_currency_flags()

    t0 = time.time()
    log.info("Running BuildPipeline...")
    mg, manifest = pipe.run()
    log.info(f"Pipeline run: {time.time()-t0:.1f}s")

    # Post-pipeline steps
    log.info("[11/12] Reactome gene-symbol enrichment")
    from gizmo.sources.gene_enrichment import enrich_graph_genes_from_reactions
    counts = enrich_graph_genes_from_reactions(mg)
    log.info(f"  added {counts['new_genes']} gene nodes, "
             f"{counts['new_edges']} catalysis edges")

    log.info("[11b/12] PubChem synonym enrichment")
    from gizmo.sources.pubchem_synonyms import enrich_pubchem_synonyms
    syn_cache = cache_dir / "pubchem" / "synonyms.json"
    if syn_cache.exists():
        n_syn = enrich_pubchem_synonyms(mg, cache_path=syn_cache)
        log.info(f"  enriched {n_syn} metabolite nodes with synonyms")

    log.info("[12/12] Conditional currency edges")
    from gizmo.analysis.currency import compute_conditional_currency_edges
    compute_conditional_currency_edges(mg)

    # Save bundle
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bundle_dir = output_dir / bundle_name
    log.info(f"Saving bundle to {bundle_dir}")
    pipe.save_bundle(mg, manifest, output_dir=str(output_dir))

    # Stamp bundle_version + artifact hashes into manifest.json
    _stamp_version_metadata(bundle_dir, BUNDLE_VERSION)
    log.info(f"Bundle version stamped: {BUNDLE_VERSION}")

    return mg, manifest


def verify_bundle(bundle_path: str | Path) -> dict:
    """Verify a saved bundle has expected structure + sources.

    Returns dict with status: ok | stale | missing_source | malformed.
    """
    bundle = Path(bundle_path)
    if not bundle.exists() or not bundle.is_dir():
        return {"status": "missing_path", "path": str(bundle)}
    needed = ["graph.json", "manifest.json"]
    missing = [f for f in needed if not (bundle / f).exists()]
    if missing:
        return {"status": "malformed", "missing_files": missing}

    import json
    manifest = json.loads((bundle / "manifest.json").read_text())
    sources = manifest.get("sources", [])
    src_names = {s.get("name") if isinstance(s, dict) else s for s in sources}

    expected_sources = {
        "reactome", "metabolon", "mondo", "orphanet", "open_targets",
        "stringdb", "currency_flags",
    }
    missing_sources = expected_sources - src_names
    if missing_sources:
        return {
            "status": "stale",
            "found_sources": sorted(src_names),
            "missing_sources": sorted(missing_sources),
        }

    qc = {}
    qc_path = bundle / "qc_report.json"
    if qc_path.exists():
        qc = json.loads(qc_path.read_text())

    # Re-verify content hashes if present
    hash_status = "not_versioned"
    saved_hashes = manifest.get("artifact_hashes") or {}
    if saved_hashes:
        all_match = True
        for fname, expected in saved_hashes.items():
            fpath = bundle / fname
            if not fpath.exists():
                all_match = False
                continue
            actual = _sha256_file(fpath)
            if actual != expected.get("sha256"):
                all_match = False
        hash_status = "match" if all_match else "drift"

    return {
        "status": "ok",
        "bundle_version": manifest.get("bundle_version", "unversioned"),
        "versioned_at": manifest.get("versioned_at"),
        "artifact_hash_status": hash_status,
        "found_sources": sorted(src_names),
        "node_count": qc.get("n_metabolites", 0) + qc.get("n_reactions", 0)
                       + qc.get("n_genes", 0) + qc.get("n_diseases", 0),
        "edge_count": qc.get("n_edges", 0),
        "built_at": manifest.get("built_at"),
        "graph_name": manifest.get("graph_name"),
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Build/verify GIZMO human_full bundle")
    sub = parser.add_subparsers(dest="cmd", required=False)

    p_build = sub.add_parser("build", help="Build the bundle (default)")
    p_build.add_argument("--cache-dir", default="data/raw")
    p_build.add_argument("--output", default="data/processed")
    p_build.add_argument("--metabolon-csv", default=None)
    p_build.add_argument("--stringdb-min-score", type=float, default=0.7)
    p_build.add_argument("--bundle-name", default="human_full")

    p_verify = sub.add_parser("verify", help="Verify a saved bundle")
    p_verify.add_argument("bundle_path")

    p_version = sub.add_parser("version", help="Print current bundle version")

    p_bump = sub.add_parser("bump", help="Bump the bundle version")
    p_bump.add_argument("part", choices=("major", "minor", "patch"),
                         help="Version part to bump")

    p_stamp = sub.add_parser("stamp", help="Re-stamp version into existing bundle")
    p_stamp.add_argument("bundle_path")

    args = parser.parse_args()
    cmd = args.cmd or "build"

    logging.basicConfig(level=logging.INFO,
                         format="%(asctime)s  %(levelname)-8s  %(message)s")

    if cmd == "version":
        print(BUNDLE_VERSION)
        return

    if cmd == "bump":
        new_v = bump_version(args.part)
        print(f"Bumped bundle version → {new_v}")
        return

    if cmd == "stamp":
        bundle = Path(args.bundle_path)
        _stamp_version_metadata(bundle, BUNDLE_VERSION)
        print(f"Stamped {bundle} with version {BUNDLE_VERSION}")
        return

    if cmd == "verify":
        result = verify_bundle(args.bundle_path)
        import json
        print(json.dumps(result, indent=2, default=str))
        sys.exit(0 if result.get("status") == "ok" else 1)
    else:
        build_human_full(
            cache_dir=args.cache_dir,
            output_dir=args.output,
            metabolon_csv=args.metabolon_csv,
            stringdb_min_score=args.stringdb_min_score,
            bundle_name=args.bundle_name,
        )


if __name__ == "__main__":
    main()
