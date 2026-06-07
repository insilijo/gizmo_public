"""make_zenodo_bundle.py — assemble per-cohort Zenodo deposit bundles.

For each cohort GIZMO has been applied to, produces a self-contained
directory with the full analytical bundle:

  cohorts/<cohort>/
    README.md                   case-study writeup (auto-populated template,
                                 hand-edit before final submission)
    F.npz                       per-patient state vector in substrate coords
    beta_alpha_per_patient.tsv  β, α-PC1..5, ‖α‖₂ per patient
    signed_basins/              per-α-PC signed basin decomposition
      alpha_pc<k>_basins.tsv    basin node IDs, loadings, mass per PC
      summary.tsv               cross-PC summary (basin sizes, masses, AUC)
    mofa_weights.json           cached MOFA+ factor weights (if available)
    mofa_weights_sm.json        substrate-matched MOFA+ (if available)
    metadata.tsv                per-patient clinical / phenotype metadata
    cohort_provenance.json      acquisition source, n, modality, license
    checksums.txt               SHA256 of every file above

Plus a top-level deposit/:
  README.md            paper series intro + deposit manifest
  substrate/           graph.json + license + hub annotation
  manifest.json        programmatic index of every cohort bundle
  zenodo_metadata.json template metadata for Zenodo API

Invocation:
  python3 benchmarks/make_zenodo_bundle.py [--out PATH] [--cohorts c1,c2]
  python3 benchmarks/make_zenodo_bundle.py --dry-run    # report what would
                                                          be written

Does NOT upload to Zenodo. The deposit script + bundle is ready; manual
upload via web UI or `zenodo_get` / `zenodo-python` after curation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))

RESULTS = REPO / "benchmarks" / "results"
UR = RESULTS / "unsupervised"
MOFA_DIR = UR / "mofa_weights"
SUBSTRATE_FILE = REPO / "data" / "processed" / "human_full" / "graph.json"


# ---------------------------------------------------------------------------
# Cohort registry — provenance metadata for each cohort GIZMO has touched.
# Edit this when adding new cohorts; the registry drives bundle assembly.
# ---------------------------------------------------------------------------

COHORT_REGISTRY = {
    # 17 panel cohorts (Paper 1 main)
    "IDH_glioma": {
        "disease": "Diffuse glioma (IDH-mut vs WT)",
        "n": 88, "modality": "paired NMR + RNA-seq",
        "acquisition": "Trautwein 2022 JCI Insight; GSE190504 (RNA) + MTBLS3873 (NMR)",
        "license": "GEO + MetaboLights public",
        "panel": True,
    },
    "TCGA_IDH_glioma": {
        "disease": "TCGA diffuse glioma (IDH-mut vs WT)",
        "n": 458, "modality": "RNA-seq",
        "acquisition": "TCGA via Broad Firehose RSEM",
        "license": "TCGA open-access",
        "panel": True,
    },
    "TCGA_LUAD": {
        "disease": "Lung adenocarcinoma (KRAS / EGFR / etc. subtypes)",
        "n": 508, "modality": "RNA-seq",
        "acquisition": "TCGA Firehose + cBioPortal mutation calls",
        "license": "TCGA open-access + cBioPortal CC-BY",
        "panel": True,
    },
    "KMPLOT_BRCA": {
        "disease": "Breast cancer (grade / subtype)",
        "n": 645, "modality": "mRNA",
        "acquisition": "KMPlot.com aggregate",
        "license": "Academic use, with attribution",
        "panel": True,
    },
    "CPTAC_CCRCC": {
        "disease": "Clear-cell renal cell carcinoma (tumor vs normal)",
        "n": 185, "modality": "RNA + proteome + phosphoproteome",
        "acquisition": "CPTAC LinkedOmics",
        "license": "CPTAC open-access",
        "panel": True,
    },
    "CPTAC_COAD": {
        "disease": "Colon adenocarcinoma (tumor vs normal)",
        "n": 207, "modality": "RNA + proteome",
        "acquisition": "CPTAC LinkedOmics",
        "license": "CPTAC open-access",
        "panel": True,
    },
    "CPTAC_OV": {
        "disease": "Ovarian serous cystadenocarcinoma (tumor vs normal)",
        "n": 103, "modality": "RNA + proteome",
        "acquisition": "CPTAC LinkedOmics",
        "license": "CPTAC open-access",
        "panel": True,
    },
    "Su_COVID": {
        "disease": "COVID-19 severity",
        "n": 270, "modality": "Olink plasma proteomics + plasma metab",
        "acquisition": "Su et al. 2020 Cell",
        "license": "Public per publication",
        "panel": True,
    },
    "Filbin_COVID": {
        "disease": "COVID-19 severity (D0 acute)",
        "n": 383, "modality": "Olink plasma proteomics",
        "acquisition": "Filbin et al. 2021 Cell Reports Medicine",
        "license": "Public per publication",
        "panel": True,
    },
    "GSE65391_SLE": {
        "disease": "Pediatric systemic lupus erythematosus (longitudinal)",
        "n": 996, "modality": "Illumina HT-12 microarray RNA",
        "acquisition": "GEO GSE65391 (Banchereau et al. 2016 Cell)",
        "license": "GEO public",
        "panel": True,
    },
    "GSE65682_sepsis": {
        "disease": "Sepsis (community-acquired pneumonia + healthy)",
        "n": 802, "modality": "Affymetrix HG-U219 microarray RNA",
        "acquisition": "GEO GSE65682 (Scicluna et al. 2017 Lancet Resp Med)",
        "license": "GEO public",
        "panel": True,
    },
    "Gao_RA": {
        "disease": "Rheumatoid arthritis (active vs healthy)",
        "n": 52, "modality": "Olink proteomics + plasma metabolomics",
        "acquisition": "Gao et al. 2024 (paired multi-omic, RA cohort)",
        "license": "Public per publication",
        "panel": True,
    },
    "GSE89408_RA": {
        "disease": "Rheumatoid arthritis vs osteoarthritis (synovial biopsy)",
        "n": 174, "modality": "RNA-seq",
        "acquisition": "GEO GSE89408 (Guo et al. 2017 Arthritis Rheumatol)",
        "license": "GEO public",
        "panel": True,
    },
    "Crohn": {
        "disease": "Crohn disease (thiopurine response)",
        "n": 33, "modality": "Olink proteomics + plasma metabolomics",
        "acquisition": "Koopman et al. 2025 (paired multi-omic, Crohn cohort)",
        "license": "Public per publication",
        "panel": True,
    },
    "HMP2_IBD_CD": {
        "disease": "IBD (CD subset)",
        "n": 399, "modality": "Stool metabolomics + 16S microbiome",
        "acquisition": "HMP2 / iHMP IBDMDB",
        "license": "iHMP public",
        "panel": True,
    },
    "Erawijantari": {
        "disease": "Gastric microbiome × metabolome",
        "n": 96, "modality": "Fecal metabolomics",
        "acquisition": "Erawijantari et al. 2020 Gut",
        "license": "Public per publication",
        "panel": True,
    },

    # LOOCV-validation cohorts (Paper 1 supp; documented in Methods).
    # `f_alias` overrides the default F filename pattern when the cohort
    # was saved under a different name.
    "NEPTUNE_kidney": {
        "disease": "Nephrotic syndrome / CKD discovery cohort",
        "n": 276, "modality": "Glomerular RNA-seq",
        "acquisition": "NEPTUNE-portal (NIDDK Nephrotic Syndrome Study Network)",
        "license": "Data-use agreement required",
        "panel": False, "loocv_validation": True,
        "f_alias": "NEPTUNE",
    },
    "Wang_RA": {
        "disease": "RA baseline + post-treatment",
        "n": 364, "modality": "Peripheral blood RNA-seq",
        "acquisition": "GEO GSE176051 (Wang et al.)",
        "license": "GEO public",
        "panel": False, "loocv_validation": True,
        "f_status": "computed-inline-not-saved",
    },
    "TB_DX": {
        "disease": "Tuberculosis diagnosis vs healthy donor",
        "n": 106, "modality": "Whole-blood RNA-seq",
        "acquisition": "GEO GSE89403 baseline subset",
        "license": "GEO public",
        "panel": False, "loocv_validation": True,
        "f_status": "computed-inline-not-saved",
    },
}


# ---------------------------------------------------------------------------
# Bundle assembly
# ---------------------------------------------------------------------------

def find_F_path(cohort):
    """F matrices may be saved with a propagation-variant suffix and/or under
    an alias name (see COHORT_REGISTRY['f_alias'])."""
    meta = COHORT_REGISTRY.get(cohort, {})
    candidates = [meta.get("f_alias"), cohort]
    for base in candidates:
        if not base: continue
        for suffix in ("", "_edge_informed", "_combined", "_node_informed"):
            p = UR / f"stage3_F_{base}{suffix}.npz"
            if p.exists():
                return p
    return None


def load_substrate_for_decomposition():
    """Lazy substrate load — only needed for β/α + basin computation."""
    from gizmo.export.json_export import read_json
    from per_patient_wlsp_v2 import biochem_subgraph
    import networkx as nx
    mg = read_json(SUBSTRATE_FILE)
    sub_dir, nodes, _ = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)
    return mg, sub, nodes, log_pr


def compute_beta_alpha_pcs(F, log_pr, n_components=5):
    """Returns β per patient, α-PC1..n loadings (n_components × n_nodes),
    α-PC scores per patient, and ‖α‖₂ per patient."""
    from sklearn.decomposition import PCA
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    alpha_norm = np.linalg.norm(alpha, axis=1)
    n = min(n_components, alpha.shape[0])
    pca = PCA(n_components=n, random_state=0)
    scores = pca.fit_transform(alpha)
    loadings = pca.components_   # n × n_nodes
    return beta, loadings, scores, alpha_norm


def compute_signed_basins(loadings, sub, nodes, mg, max_top=200):
    """For each α-PC, partition nodes by loading sign and identify the
    largest connected component in each sign-class. Returns a dict
    {pc_idx (1-based): {'pos_nodes': [...], 'neg_nodes': [...],
                          'pos_mass': float, 'neg_mass': float}}.
    Node entries are (node_id, loading, abs_loading) sorted by |loading|."""
    import networkx as nx
    out = {}
    for k in range(loadings.shape[0]):
        pc = loadings[k]
        pos_nids = [nodes[i] for i in range(len(nodes)) if pc[i] > 0]
        neg_nids = [nodes[i] for i in range(len(nodes)) if pc[i] < 0]
        pos_sub = sub.subgraph(pos_nids)
        neg_sub = sub.subgraph(neg_nids)
        pos_ccs = sorted(nx.connected_components(pos_sub), key=len, reverse=True)
        neg_ccs = sorted(nx.connected_components(neg_sub), key=len, reverse=True)
        pos_cc = pos_ccs[0] if pos_ccs else set()
        neg_cc = neg_ccs[0] if neg_ccs else set()
        total_mass = float((pc ** 2).sum())
        pos_mass = float(sum(pc[nodes.index(n)] ** 2 for n in pos_cc)) / max(total_mass, 1e-12)
        neg_mass = float(sum(pc[nodes.index(n)] ** 2 for n in neg_cc)) / max(total_mass, 1e-12)

        def _top_records(cc, sign):
            rows = []
            for n in cc:
                i = nodes.index(n)
                load = float(pc[i])
                if sign == "+" and load <= 0: continue
                if sign == "-" and load >= 0: continue
                attrs = mg.graph.nodes.get(n, {})
                sym = attrs.get("symbol") or attrs.get("name") or attrs.get("display_name") or n
                rows.append({
                    "node_id": n, "symbol": str(sym),
                    "node_type": attrs.get("node_type", "?"),
                    "loading": load, "abs_loading": abs(load),
                })
            rows.sort(key=lambda r: -r["abs_loading"])
            return rows[:max_top]

        out[k + 1] = {
            "pos_basin": _top_records(pos_cc, "+"),
            "neg_basin": _top_records(neg_cc, "-"),
            "pos_basin_size": len(pos_cc),
            "neg_basin_size": len(neg_cc),
            "pos_mass": pos_mass,
            "neg_mass": neg_mass,
        }
    return out


def _load_cohort_metadata(cohort):
    """Dispatch to the right per-cohort loader in
    benchmarks.diagnostics.axis_metadata_extended. Returns DataFrame or None."""
    try:
        from benchmarks.diagnostics.axis_metadata_extended import (
            load_kmplot_metadata, load_cptac_metadata,
            load_tcga_idh_glioma_metadata, load_idh_glioma_trautwein_metadata,
            load_gao_ra_metadata, load_crohn_metadata, load_filbin_metadata,
            load_erawijantari_metadata, load_hmp2_metadata,
            load_tcga_luad_metadata, load_su_covid_metadata,
            load_gse_series_metadata,
        )
    except Exception as e:
        print(f"      could not import metadata loaders: {e}", flush=True)
        return None
    handlers = {
        "KMPLOT_BRCA":     load_kmplot_metadata,
        "TCGA_IDH_glioma": load_tcga_idh_glioma_metadata,
        "IDH_glioma":      load_idh_glioma_trautwein_metadata,
        "Gao_RA":          load_gao_ra_metadata,
        "Crohn":           load_crohn_metadata,
        "Filbin_COVID":    load_filbin_metadata,
        "Erawijantari":    load_erawijantari_metadata,
        "HMP2_IBD_CD":     load_hmp2_metadata,
        "TCGA_LUAD":       load_tcga_luad_metadata,
        "Su_COVID":        load_su_covid_metadata,
        "CPTAC_CCRCC":     lambda: load_cptac_metadata("CPTAC_CCRCC"),
        "CPTAC_COAD":      lambda: load_cptac_metadata("CPTAC_COAD"),
        "CPTAC_OV":        lambda: load_cptac_metadata("CPTAC_OV"),
        "GSE65391_SLE":    lambda: load_gse_series_metadata(
            REPO / "data/cohorts/GSE65391_SLE/GSE65391_series_matrix.txt.gz"),
        "GSE65682_sepsis": lambda: load_gse_series_metadata(
            REPO / "data/cohorts/GSE65682_sepsis/GSE65682_series_matrix.txt.gz"),
    }
    if cohort not in handlers:
        return None
    try:
        return handlers[cohort]()
    except Exception as e:
        print(f"      metadata loader failed for {cohort}: {e}", flush=True)
        return None


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def write_cohort_readme(cohort, meta, F_shape, basin_summary, has_mofa, has_mofa_sm,
                          out_path):
    """Auto-populated case-study writeup template. Edit before final submission."""
    lines = [
        f"# {cohort}",
        "",
        f"**Disease:** {meta.get('disease', 'TBD')}",
        f"**n samples:** {F_shape[0]}",
        f"**Modality:** {meta.get('modality', 'TBD')}",
        f"**Acquisition:** {meta.get('acquisition', 'TBD')}",
        f"**License:** {meta.get('license', 'TBD')}",
        "",
        "## Bundle contents",
        "",
        f"- `F.npz` — per-patient state vector ({F_shape[0]} patients × {F_shape[1]} substrate nodes)",
        "- `beta_alpha_per_patient.tsv` — β, α-PC1..5 scores, ‖α‖₂ per patient",
        "- `signed_basins/alpha_pc<k>_basins.tsv` — per-PC signed-basin decomposition (basin node IDs, loadings, mass)",
        "- `signed_basins/summary.tsv` — cross-PC summary (basin sizes, masses)",
    ]
    if has_mofa:
        lines.append("- `mofa_weights.json` — cached MOFA+ factor weights (full input universe)")
    if has_mofa_sm:
        lines.append("- `mofa_weights_sm.json` — substrate-matched MOFA+ factor weights (input restricted to GIZMO's substrate-mappable universe)")
    lines += [
        "- `metadata.tsv` — per-patient clinical / phenotype metadata",
        "- `cohort_provenance.json` — acquisition source, license, modality, panel/validation status",
        "- `checksums.txt` — SHA256 of every file above",
        "",
        "## Per-PC signed-basin summary",
        "",
        "| α-PC | + basin nodes | + mass | − basin nodes | − mass | Top + node | Top − node |",
        "|---|---|---|---|---|---|---|",
    ]
    for k in sorted(basin_summary.keys()):
        b = basin_summary[k]
        top_pos = b["pos_basin"][0]["symbol"] if b["pos_basin"] else "—"
        top_neg = b["neg_basin"][0]["symbol"] if b["neg_basin"] else "—"
        lines.append(
            f"| α-PC{k} | {b['pos_basin_size']:,} | {b['pos_mass']*100:.1f}% "
            f"| {b['neg_basin_size']:,} | {b['neg_mass']*100:.1f}% | {top_pos} | {top_neg} |"
        )
    lines += [
        "",
        "## Case-study notes",
        "",
        "*Hand-edit before final submission: cross-reference top-loaded basin members against the canonical disease mechanism (CTD / Reactome / KEGG / cohort publication); note which α-PCs map to which clinical metadata field; flag any cohort-specific scope conditions (small n, modality budget, missing controls, etc.).*",
        "",
        "## Citation",
        "",
        "If you use this cohort bundle, please cite:",
        "- The cohort source publication (see Acquisition above)",
        "- GIZMO Paper 1: *A biochemistry substrate as a fixed coordinate system for multi-omic per-patient projection* (Zenodo DOI on acceptance)",
        f"- This bundle: `{cohort}` per-cohort case study at Zenodo DOI on acceptance",
        "",
        "",
    ]
    out_path.write_text("\n".join(lines))


def assemble_cohort_bundle(cohort, meta, state, deposit_root, dry_run=False):
    """Build cohorts/<cohort>/ deposit directory.
    state = (mg, sub, nodes, log_pr) — pre-loaded substrate context."""
    bundle = deposit_root / "cohorts" / cohort
    if dry_run:
        print(f"  [dry] would write {bundle}/", flush=True)
        return None

    bundle.mkdir(parents=True, exist_ok=True)

    F_path = find_F_path(cohort)
    if F_path is None:
        # Some LOOCV-validation cohorts were computed inline and never saved
        # to disk. Write a stub bundle that documents this so the cohort is
        # still represented in the manifest, with instructions for
        # downstream consumers to recompute.
        if meta.get("f_status") == "computed-inline-not-saved":
            (bundle / "F_NOT_SAVED.md").write_text(
                f"# {cohort} — F matrix not deposited\n\n"
                "This cohort's per-patient F matrix was computed inline during "
                "the LOOCV analysis pass and not saved to disk. To regenerate, "
                "load the cohort via the appropriate loader in "
                "`benchmarks/per_patient_master.py` or `benchmarks/drug_sim_multi_cohort.py`, "
                "then run `per_patient_master.fit_per_patient_F(...)` against the "
                "v1 substrate. See the cohort's source in "
                "`cohort_provenance.json`.\n"
            )
            (bundle / "cohort_provenance.json").write_text(json.dumps(meta, indent=2))
            print(f"  ⚠ {cohort}: F not saved — wrote stub bundle", flush=True)
            return {
                "cohort": cohort,
                "bundle_dir": str(bundle.relative_to(deposit_root)),
                "n_patients": meta.get("n", None),
                "n_features": None, "n_pcs": 0,
                "has_mofa_full": False, "has_mofa_sm": False,
                "panel": meta.get("panel", False),
                "loocv_validation": meta.get("loocv_validation", False),
                "disease": meta.get("disease", "TBD"),
                "n_files": 2,
                "stub_only": True,
            }
        print(f"  ⚠ {cohort}: no F matrix found — skipping bundle", flush=True)
        return None

    # F matrix
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"]
    patient_ids = list(fd["patient_ids"])
    np.savez_compressed(bundle / "F.npz", F=F,
                          patient_ids=np.array(patient_ids, dtype=object))
    print(f"    F.npz {F.shape}", flush=True)

    # β/α per patient + α-PC loadings
    mg, sub, nodes, log_pr = state
    beta, loadings, scores, alpha_norm = compute_beta_alpha_pcs(
        F.astype(np.float64), log_pr, n_components=5)
    per_pt = pd.DataFrame({
        "patient_id": patient_ids,
        "beta": beta,
        "alpha_norm": alpha_norm,
        **{f"alpha_pc{k+1}_score": scores[:, k] for k in range(scores.shape[1])},
    })
    per_pt.to_csv(bundle / "beta_alpha_per_patient.tsv", sep="\t", index=False)
    print(f"    beta_alpha_per_patient.tsv {per_pt.shape}", flush=True)

    # Signed-basin decomposition per α-PC
    basin_dir = bundle / "signed_basins"
    basin_dir.mkdir(exist_ok=True)
    basins = compute_signed_basins(loadings, sub, nodes, mg, max_top=200)
    summary_rows = []
    for k, b in basins.items():
        rows = []
        for r in b["pos_basin"]:
            rows.append({**r, "basin": "+"})
        for r in b["neg_basin"]:
            rows.append({**r, "basin": "-"})
        pd.DataFrame(rows).to_csv(basin_dir / f"alpha_pc{k}_basins.tsv",
                                    sep="\t", index=False)
        summary_rows.append({
            "alpha_pc": k,
            "pos_basin_size": b["pos_basin_size"],
            "pos_basin_mass_pct": round(100 * b["pos_mass"], 2),
            "neg_basin_size": b["neg_basin_size"],
            "neg_basin_mass_pct": round(100 * b["neg_mass"], 2),
            "top_pos_node": b["pos_basin"][0]["symbol"] if b["pos_basin"] else "",
            "top_neg_node": b["neg_basin"][0]["symbol"] if b["neg_basin"] else "",
        })
    pd.DataFrame(summary_rows).to_csv(basin_dir / "summary.tsv",
                                         sep="\t", index=False)
    print(f"    signed_basins/ {len(basins)} PCs", flush=True)

    # MOFA+ weights (if cached)
    has_mofa = False; has_mofa_sm = False
    mofa_src = MOFA_DIR / f"mofa_weights_{cohort}.json"
    if mofa_src.exists():
        shutil.copy2(mofa_src, bundle / "mofa_weights.json")
        has_mofa = True
        print(f"    mofa_weights.json copied", flush=True)
    mofa_sm_src = MOFA_DIR / f"mofa_weights_{cohort}_sm.json"
    if mofa_sm_src.exists():
        shutil.copy2(mofa_sm_src, bundle / "mofa_weights_sm.json")
        has_mofa_sm = True
        print(f"    mofa_weights_sm.json copied", flush=True)

    # Metadata
    md = _load_cohort_metadata(cohort)
    if md is not None and not md.empty:
        md.to_csv(bundle / "metadata.tsv", sep="\t", index=False)
        print(f"    metadata.tsv {md.shape}", flush=True)
    else:
        print(f"    metadata unavailable for {cohort}", flush=True)

    # Cohort provenance JSON
    (bundle / "cohort_provenance.json").write_text(json.dumps(meta, indent=2))

    # Auto-populated README case-study
    write_cohort_readme(cohort, meta, F.shape, basins, has_mofa, has_mofa_sm,
                          bundle / "README.md")

    # SHA256 checksums
    checksums = []
    for f in sorted(bundle.rglob("*")):
        if f.is_file() and f.name != "checksums.txt":
            rel = f.relative_to(bundle)
            checksums.append(f"{file_sha256(f)}  {rel}")
    (bundle / "checksums.txt").write_text("\n".join(checksums) + "\n")
    print(f"    checksums.txt ({len(checksums)} files)", flush=True)

    return {
        "cohort": cohort,
        "bundle_dir": str(bundle.relative_to(deposit_root)),
        "n_patients": F.shape[0],
        "n_features": F.shape[1],
        "n_pcs": loadings.shape[0],
        "has_mofa_full": has_mofa,
        "has_mofa_sm": has_mofa_sm,
        "panel": meta.get("panel", False),
        "loocv_validation": meta.get("loocv_validation", False),
        "disease": meta.get("disease", "TBD"),
        "n_files": len(checksums) + 1,   # +1 for checksums.txt itself
    }


def assemble_substrate_dir(deposit_root, dry_run=False):
    """Copy substrate file + write top-level substrate README."""
    sub_dir = deposit_root / "substrate"
    if dry_run:
        print(f"  [dry] would write {sub_dir}/", flush=True)
        return
    sub_dir.mkdir(parents=True, exist_ok=True)
    if SUBSTRATE_FILE.exists():
        shutil.copy2(SUBSTRATE_FILE, sub_dir / "graph.json")
    (sub_dir / "README.md").write_text(
        "# GIZMO biochemistry substrate (v1)\n\n"
        "38,148-node merged graph (16,343 genes + 6,406 metabolites + "
        "15,399 reactions + ancillary). Sources: Reactome, StringDB, HMDB, "
        "KEGG. License: CC-BY 4.0.\n\n"
        "Format: nx node-link JSON. Use `gizmo.export.json_export.read_json` to load.\n\n"
        "## Citation\n\n"
        "If you use the substrate, please cite GIZMO Paper 1 "
        "(Zenodo DOI on acceptance) and the four underlying source databases "
        "(Reactome, StringDB, HMDB, KEGG).\n"
    )


def write_top_readme(deposit_root, bundle_records):
    """Top-level paper-series intro + deposit manifest."""
    lines = [
        "# GIZMO substrate + per-cohort case studies — Zenodo deposit",
        "",
        "Companion deposit to *A biochemistry substrate as a fixed coordinate "
        "system for multi-omic per-patient projection* (Paper 1 of the GIZMO series).",
        "",
        "## Contents",
        "",
        "- `substrate/` — the 38,148-node GIZMO biochemistry substrate "
        "(CC-BY 4.0). The primary deposit of this paper.",
        f"- `cohorts/<cohort>/` — {len(bundle_records)} per-cohort case-study "
        "bundles. Each contains the F matrix, β/α decomposition, signed-basin "
        "biology, MOFA+ comparison weights, metadata, and a case-study writeup.",
        "- `manifest.json` — programmatic index of every bundle.",
        "- `zenodo_metadata.json` — template metadata for the Zenodo deposit record.",
        "",
        "## Cohort manifest",
        "",
        "| Cohort | Disease | n | Modality | Role | MOFA+ | MOFA+_sm |",
        "|---|---|---|---|---|---|---|",
    ]
    for b in bundle_records:
        role = "panel" if b["panel"] else ("LOOCV-val" if b["loocv_validation"] else "supplementary")
        mof = "✓" if b["has_mofa_full"] else "—"
        mof_sm = "✓" if b["has_mofa_sm"] else "—"
        lines.append(f"| {b['cohort']} | {b['disease']} | {b['n_patients']} | "
                     f"see bundle | {role} | {mof} | {mof_sm} |")
    lines += [
        "",
        "## Reproducing per-cohort projections",
        "",
        "Each cohort bundle is self-contained. Reload the F matrix and recompute "
        "β / α / signed basins via the published GIZMO codebase (commit hash + tag "
        "linked from the manuscript). Cross-cohort cosine comparison is direct on "
        "the substrate-coordinate α-PC loadings shipped per bundle.",
        "",
        "## License",
        "",
        "- Substrate (`substrate/graph.json`): CC-BY 4.0.",
        "- Per-cohort F matrices, β/α, basin output, MOFA+ weights: CC-BY 4.0.",
        "- Per-cohort raw input data: governed by each cohort's original source "
        "license (see `cohorts/<cohort>/cohort_provenance.json`).",
        "",
    ]
    (deposit_root / "README.md").write_text("\n".join(lines))


def write_manifest_and_zenodo_metadata(deposit_root, bundle_records):
    manifest = {
        "deposit_version": "v1",
        "substrate_nodes": 38148,
        "n_cohorts": len(bundle_records),
        "cohorts": bundle_records,
    }
    (deposit_root / "manifest.json").write_text(json.dumps(manifest, indent=2))

    zen_md = {
        "metadata": {
            "title": "GIZMO biochemistry substrate + per-cohort case-study bundles (Paper 1)",
            "upload_type": "dataset",
            "description": (
                "Companion deposit to GIZMO Paper 1. Contains the 38,148-node "
                "biochemistry substrate (Reactome + StringDB + HMDB + KEGG; "
                "CC-BY 4.0) plus per-cohort case-study bundles for every "
                "cohort GIZMO has been applied to. Each cohort bundle ships "
                "the per-patient F matrix in substrate coordinates, β/α "
                "decomposition, signed-basin biology decomposition per α-PC, "
                "MOFA+ comparison weights (both full-input and substrate-matched), "
                "per-patient metadata, and a per-cohort case-study writeup."
            ),
            "creators": [
                {"name": "Gardner, Joseph J.",
                 "affiliation": "Insilijo (independent)",
                 "orcid": "TBD"}
            ],
            "keywords": [
                "multi-omic integration", "biochemistry substrate",
                "MAP reconstruction", "graph Laplacian", "α-PC decomposition",
                "signed-basin biology", "cross-cohort meta-analysis",
                "GIZMO", "Reactome", "HMDB", "KEGG", "StringDB",
            ],
            "license": "cc-by-4.0",
            "communities": [{"identifier": "biomedical-data"}],
            "related_identifiers": [
                {"identifier": "https://github.com/insilijo/gizmo",
                 "relation": "isSupplementTo", "scheme": "url"},
            ],
        }
    }
    (deposit_root / "zenodo_metadata.json").write_text(json.dumps(zen_md, indent=2))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "zenodo_deposit"),
                    help="Output deposit root directory")
    ap.add_argument("--cohorts", default=None,
                    help="Comma-separated cohort filter (default: all registered)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would be written; don't write files")
    args = ap.parse_args()

    deposit_root = Path(args.out)
    if not args.dry_run:
        deposit_root.mkdir(parents=True, exist_ok=True)

    targets = (args.cohorts.split(",") if args.cohorts
               else list(COHORT_REGISTRY.keys()))

    print(f"Deposit root: {deposit_root}", flush=True)
    print(f"Cohorts to bundle: {len(targets)}", flush=True)
    print()

    # Substrate bundle
    print("=== substrate ===", flush=True)
    assemble_substrate_dir(deposit_root, dry_run=args.dry_run)

    # Per-cohort bundles
    state = None
    bundle_records = []
    for cohort in targets:
        if cohort not in COHORT_REGISTRY:
            print(f"⚠ {cohort} not in registry; skipping", flush=True); continue
        meta = COHORT_REGISTRY[cohort]
        print(f"\n=== {cohort} ===", flush=True)
        if args.dry_run:
            print(f"  [dry] disease: {meta['disease']}", flush=True)
            print(f"  [dry] F path: {find_F_path(cohort)}", flush=True)
            continue
        if state is None:
            print("  loading substrate (one-time)…", flush=True)
            state = load_substrate_for_decomposition()
        rec = assemble_cohort_bundle(cohort, meta, state, deposit_root, dry_run=False)
        if rec is not None:
            bundle_records.append(rec)

    if args.dry_run:
        print("\n(dry-run; no files written)", flush=True)
        return

    # Top-level manifest + README + Zenodo metadata
    print("\n=== top-level deposit files ===", flush=True)
    write_top_readme(deposit_root, bundle_records)
    write_manifest_and_zenodo_metadata(deposit_root, bundle_records)
    print(f"\nDeposit assembled at {deposit_root}", flush=True)
    print(f"  {len(bundle_records)} cohort bundles", flush=True)
    print(f"  Next step: hand-edit per-cohort README.md case-study sections, "
          f"then upload via Zenodo web UI or zenodo-python.", flush=True)


if __name__ == "__main__":
    main()
