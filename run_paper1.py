"""End-to-end Paper 1 pipeline runner.

For one cohort: load → MAP-solve F → β/α → α-PCA → signed basins.
Writes outputs to ``results/<cohort>/`` and optionally renders the signed-basin
figure for α-PC1.

Usage
-----
    python3 run_paper1.py Crohn                              # default cohort
    python3 run_paper1.py IDH_glioma --n-components 5
    python3 run_paper1.py Crohn --substrate substrate/graph.json
    python3 run_paper1.py Crohn --out results/Crohn --skip-figure
    python3 run_paper1.py --from-F results/Crohn/F.npz \
        --patient-ids results/Crohn/patient_ids.json \
        --out results/Crohn_replay              # replay from cached F matrix

Inputs
------
Cohort raw data is NOT redistributed with this repo. Loaders in
``benchmarks/per_patient_master.py`` raise ``FileNotFoundError`` with the
expected path if data is missing. The Zenodo deposit (DOI on acceptance)
ships per-cohort F matrices that can be replayed via ``--from-F``.

Outputs (under ``--out``)
-------------------------
- ``F.npz`` — (n_patients, n_nodes) state matrix + patient_ids
- ``beta_alpha.tsv`` — per-patient β, ‖α‖₂, α-PC1..PCk scores
- ``alpha_pc_loadings.npz`` — α-PC components (k × n_nodes)
- ``signed_basins.tsv`` — one row per (PC, sign) listing top basin members
- ``diagnostics.json`` — variance explained, mean smoothness, settings
- ``alpha_pc1_basin.png`` — signed-basin figure for α-PC1 (unless --skip-figure)
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("run_paper1")


COHORT_LOADERS = {
    "Crohn": "load_crohn",
    "Su_COVID": "load_su_covid",
    "Gao_RA": "load_gao_ra",
    "IDH_glioma": "load_idh_glioma",
    "TCGA_IDH_glioma": "load_tcga_idh_glioma",
    "Filbin_COVID": "load_filbin_covid",
    "Erawijantari": "load_erawijantari",
    "KMPLOT_BRCA": "load_kmplot_brca",
    "TCGA_LUAD": "load_tcga_luad",
    "GSE89408_RA": "load_gse89408_ra",
    "HMP2_IBD_CD": "load_hmp2_ibd_cd",
}


def _load_cohort(name: str):
    """Dispatch to the per-cohort loader in benchmarks/per_patient_master.py."""
    if name not in COHORT_LOADERS:
        raise SystemExit(
            f"Unknown cohort '{name}'. Known: {sorted(COHORT_LOADERS)}.\n"
            "If your cohort isn't listed, add a loader to "
            "benchmarks/per_patient_master.py and register it in "
            "COHORT_LOADERS at the top of run_paper1.py.")
    from per_patient_master import (  # noqa: WPS433 — local import after sys.path
        __dict__ as ppm_ns,
    )
    fn = ppm_ns[COHORT_LOADERS[name]]
    log.info("Loading cohort %s via %s()", name, fn.__name__)
    out = fn()
    if len(out) != 4:
        raise RuntimeError(
            f"Loader {fn.__name__} returned {len(out)} items; "
            "expected (prot, metab, ylabel, samples)")
    return out


def _build_modality_setups(prot, metab, samples, geometry):
    """Map prot + metab dicts onto substrate node columns, z-score per modality.

    Returns a list of ModalitySetup objects suitable for solve_map().
    """
    from gizmo.evidence.mappers import GeneMapper, MetaboliteMapper
    from gizmo.inference.projection import ModalitySetup

    # NB: mappers consume the full MultiGraph, not the subgraph geometry —
    # but feature_cols must reference geometry.nid_idx.
    setups = []

    if prot:
        gmap = GeneMapper(geometry._mg)  # type: ignore[attr-defined]
        prot_features = sorted({f for s in samples if s in prot for f in prot[s]})
        prot_node = {}
        for f in prot_features:
            node_id, _conf = gmap.map(f)
            if node_id and node_id in geometry.nid_idx:
                prot_node[f] = node_id
        log.info("Proteomics: %d features → %d substrate nodes",
                 len(prot_features), len(prot_node))
        if prot_node:
            data_p = {}
            vals_by_feat = {f: [] for f in prot_node}
            for s in samples:
                if s not in prot:
                    continue
                for f in prot_node:
                    if f in prot[s]:
                        vals_by_feat[f].append(prot[s][f])
            stats = {}
            for f, vs in vals_by_feat.items():
                if vs:
                    arr = np.asarray(vs, dtype=float)
                    stats[f] = (arr.mean(), arr.std() + 1e-9)
            for s in samples:
                if s not in prot:
                    continue
                d = {}
                for f in prot_node:
                    if f in prot[s] and f in stats:
                        mu, sd = stats[f]
                        d[f] = float((prot[s][f] - mu) / sd)
                data_p[s] = d
            feat_cols_p = [(f, geometry.nid_idx[prot_node[f]]) for f in prot_node]
            setups.append(ModalitySetup(
                label="proteomics", sigma=1.0, diffusion_t=0.0,
                feature_cols=feat_cols_p, data=data_p))

    if metab:
        mmap = MetaboliteMapper(geometry._mg)  # type: ignore[attr-defined]
        metab_features = sorted({f for s in samples if s in metab for f in metab[s]})
        metab_node = {}
        for f in metab_features:
            node_id, _conf = mmap.map(f)
            if node_id and node_id in geometry.nid_idx:
                metab_node[f] = node_id
        log.info("Metabolomics: %d features → %d substrate nodes",
                 len(metab_features), len(metab_node))
        if metab_node:
            data_m = {}
            vals_by_feat = {f: [] for f in metab_node}
            for s in samples:
                if s not in metab:
                    continue
                for f in metab_node:
                    if f in metab[s]:
                        vals_by_feat[f].append(metab[s][f])
            stats = {}
            for f, vs in vals_by_feat.items():
                if vs:
                    arr = np.asarray(vs, dtype=float)
                    stats[f] = (arr.mean(), arr.std() + 1e-9)
            for s in samples:
                if s not in metab:
                    continue
                d = {}
                for f in metab_node:
                    if f in metab[s] and f in stats:
                        mu, sd = stats[f]
                        d[f] = float((metab[s][f] - mu) / sd)
                data_m[s] = d
            feat_cols_m = [(f, geometry.nid_idx[metab_node[f]]) for f in metab_node]
            setups.append(ModalitySetup(
                label="metabolomics", sigma=1.0, diffusion_t=0.0,
                feature_cols=feat_cols_m, data=data_m))

    if not setups:
        raise RuntimeError(
            "No modalities had any features mapping to the substrate. "
            "Check feature naming conventions in the loader.")
    return setups


def _write_outputs(result, geometry, out_dir: Path, ylabel: dict | None = None):
    """Persist the per-cohort Paper 1 artifacts to disk."""
    out_dir.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(out_dir / "F.npz",
                        F=result.F,
                        patient_ids=np.array(result.patient_ids, dtype=object))

    # beta_alpha.tsv
    lines = ["patient_id\tlabel\tbeta\talpha_norm" + "".join(
        f"\talpha_pc{k+1}" for k in range(result.alpha_pc_scores.shape[1]))]
    for i, pid in enumerate(result.patient_ids):
        lab = (ylabel or {}).get(pid, "")
        cells = [str(pid), str(lab),
                 f"{result.beta[i]:.6f}",
                 f"{result.alpha_norm[i]:.6f}"]
        cells += [f"{result.alpha_pc_scores[i, k]:.6f}"
                  for k in range(result.alpha_pc_scores.shape[1])]
        lines.append("\t".join(cells))
    (out_dir / "beta_alpha.tsv").write_text("\n".join(lines) + "\n")

    np.savez_compressed(out_dir / "alpha_pc_loadings.npz",
                        components=result.alpha_pc_components,
                        explained_variance_ratio=result.alpha_pc_explained_variance,
                        node_ids=np.array(geometry.nodes, dtype=object))

    # signed_basins.tsv
    sb_lines = ["pc\tsign\trank\tnode_id\tsymbol\tloading"]
    mg = geometry._mg  # type: ignore[attr-defined]
    for basin in result.signed_basins:
        pc = basin["pc_index"]
        loadings = result.alpha_pc_components[pc - 1]
        for sign, key in [("+", "pos_basin"), ("-", "neg_basin")]:
            members = basin.get(key, {}).get("nodes", [])
            for rank, node_id in enumerate(members, 1):
                idx = geometry.nid_idx.get(node_id)
                if idx is None:
                    continue
                attrs = mg.graph.nodes.get(node_id, {})
                sym = attrs.get("symbol") or attrs.get("name") or node_id
                sb_lines.append(
                    f"{pc}\t{sign}\t{rank}\t{node_id}\t{sym}\t{loadings[idx]:+.6f}")
    (out_dir / "signed_basins.tsv").write_text("\n".join(sb_lines) + "\n")

    diag = dict(result.diagnostics)
    diag["alpha_pc_explained_variance_ratio"] = \
        result.alpha_pc_explained_variance.tolist()
    (out_dir / "diagnostics.json").write_text(json.dumps(diag, indent=2))

    log.info("Wrote outputs to %s", out_dir)


def _render_basin_figure(result, geometry, out_path: Path, *, pc_index: int = 1,
                          top_k: int = 12):
    """Render a simple node × loading bar chart for α-PC<pc_index> basins.

    This is the lightweight figure; the manuscript-grade chord co-structure
    figure lives in benchmarks/figures/build_chord_costructure.py.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib unavailable; skipping figure")
        return

    basin = result.signed_basins[pc_index - 1]
    loadings = result.alpha_pc_components[pc_index - 1]
    mg = geometry._mg  # type: ignore[attr-defined]

    def _label(node_id):
        attrs = mg.graph.nodes.get(node_id, {})
        return attrs.get("symbol") or attrs.get("name") or node_id

    fig, ax = plt.subplots(figsize=(8, 0.32 * 2 * top_k + 1.5))
    pos = basin.get("pos_basin", {}).get("nodes", [])[:top_k]
    neg = basin.get("neg_basin", {}).get("nodes", [])[:top_k]

    nodes_ord = list(neg)[::-1] + list(pos)
    vals = []
    for n in nodes_ord:
        idx = geometry.nid_idx.get(n)
        vals.append(loadings[idx] if idx is not None else 0.0)
    labels = [_label(n) for n in nodes_ord]
    colors = ["#c44e52" if v > 0 else "#4878d0" for v in vals]

    ax.barh(range(len(vals)), vals, color=colors)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel(f"α-PC{pc_index} loading")
    ax.set_title(f"α-PC{pc_index} signed basin (top {top_k}/sign)")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Wrote figure %s", out_path)


def _replay_from_F(F_path: Path, patient_ids_path: Path | None,
                    geometry, n_components: int):
    """Replay β/α + signed basins from a pre-computed F matrix (Zenodo deposit)."""
    from gizmo.inference.projection import (Paper1Result, decompose_beta_alpha,
                                              extract_signed_basin)

    with np.load(F_path, allow_pickle=True) as zf:
        F = np.asarray(zf["F"], dtype=np.float32)
        if "patient_ids" in zf:
            patient_ids = list(zf["patient_ids"])
        else:
            patient_ids = None
    if patient_ids is None and patient_ids_path and patient_ids_path.exists():
        patient_ids = json.loads(patient_ids_path.read_text())
    if patient_ids is None:
        patient_ids = [f"sample_{i}" for i in range(F.shape[0])]

    if F.shape[1] != len(geometry.nodes):
        raise SystemExit(
            f"F matrix has {F.shape[1]} columns but substrate has "
            f"{len(geometry.nodes)} nodes — did F come from a different "
            "substrate build?")

    log.info("Replaying β/α from F: %d patients × %d nodes",
             F.shape[0], F.shape[1])

    beta, alpha_norm, alpha_pc_scores, pca = decompose_beta_alpha(
        F, geometry.log_pr, n_components=n_components)
    signed_basins = []
    for k in range(pca.components_.shape[0]):
        basin = extract_signed_basin(pca.components_[k], geometry, top_k=15)
        basin["pc_index"] = k + 1
        basin["explained_variance_ratio"] = float(pca.explained_variance_ratio_[k])
        signed_basins.append(basin)

    return Paper1Result(
        patient_ids=patient_ids,
        F=F,
        beta=beta,
        alpha_norm=alpha_norm,
        alpha_pc_scores=alpha_pc_scores,
        alpha_pc_components=pca.components_,
        alpha_pc_explained_variance=pca.explained_variance_ratio_,
        signed_basins=signed_basins,
        smoothness=np.zeros(F.shape[0], dtype=np.float32),
        diagnostics={
            "n_patients": F.shape[0],
            "n_nodes": F.shape[1],
            "replay": True,
            "alpha_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        },
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort", nargs="?", default="Crohn",
                    help="Cohort name (see COHORT_LOADERS in this script).")
    ap.add_argument("--substrate", default="substrate/graph.json",
                    help="Path to substrate JSON (default: substrate/graph.json)")
    ap.add_argument("--out", default=None,
                    help="Output directory (default: results/<cohort>)")
    ap.add_argument("--n-components", type=int, default=5,
                    help="Number of α-PCs to extract (default: 5)")
    ap.add_argument("--hub-cap", type=int, default=200,
                    help="Max degree per substrate node (default: 200)")
    ap.add_argument("--skip-figure", action="store_true",
                    help="Don't render the α-PC1 basin figure")
    ap.add_argument("--from-F", default=None,
                    help="Path to a precomputed F.npz to replay β/α from. "
                         "Skips cohort load + MAP solve.")
    ap.add_argument("--patient-ids", default=None,
                    help="Optional patient_ids JSON if --from-F's npz lacks them.")
    args = ap.parse_args()

    out_dir = Path(args.out) if args.out else (REPO / "results" / args.cohort)

    from gizmo.export.json_export import read_json
    from gizmo.inference.projection import (build_biochem_subgraph,
                                              run_paper1_pipeline)

    log.info("Loading substrate from %s", args.substrate)
    mg = read_json(args.substrate)
    log.info("Substrate: %d nodes, %d edges",
             mg.graph.number_of_nodes(), mg.graph.number_of_edges())
    geometry = build_biochem_subgraph(mg, hub_cap=args.hub_cap)
    geometry._mg = mg  # type: ignore[attr-defined]  # for mapper construction later
    log.info("Subgraph (post-hub-cap=%d): %d nodes",
             args.hub_cap, len(geometry.nodes))

    if args.from_F:
        result = _replay_from_F(Path(args.from_F),
                                 Path(args.patient_ids) if args.patient_ids else None,
                                 geometry, args.n_components)
        ylabel = None
    else:
        prot, metab, ylabel, samples = _load_cohort(args.cohort)
        log.info("Cohort %s: n_samples=%d", args.cohort, len(samples))
        modality_setups = _build_modality_setups(prot, metab, samples, geometry)
        result = run_paper1_pipeline(
            geometry, modality_setups, samples,
            n_components=args.n_components)

    _write_outputs(result, geometry, out_dir, ylabel=ylabel)
    if not args.skip_figure:
        _render_basin_figure(result, geometry, out_dir / "alpha_pc1_basin.png",
                              pc_index=1)

    print("Pipeline complete.")
    print(f"  outputs:  {out_dir}")
    print(f"  α-PC EVR: {result.alpha_pc_explained_variance.tolist()}")


if __name__ == "__main__":
    main()
