"""Stage 31 v14c — WGCNA vs GIZMO at same gene universe.

v14b had universe confound: WGCNA used measured features only; GIZMO modules
expanded via substrate. Different universes, unfair comparison.

v14c restricts BOTH to the same gene universe (cohort-measured features →
genes), then compares modular pathway coherence:
  - WGCNA: correlation-Louvain modules on cohort features, restricted to genes
    that exist in pathway map
  - GIZMO_restricted: for each v6 GIZMO module, restrict its gene_symbols to
    only genes present in the cohort feature universe, then compute coherence

Tests: given the SAME gene universe, does substrate-Louvain organize more
pathway-coherently than correlation-Louvain? This is the honest comparison.

Output: stage31_v14c_wgcna_same_universe.json
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
import networkx as nx

REPO = Path(__file__).resolve().parent.parent.parent
RESULTS = REPO / "benchmarks" / "results" / "unsupervised"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))
sys.path.insert(0, str(REPO / "benchmarks" / "unsupervised_stratification"))


CORR_THRESH = 0.5
MIN_MODULE_SIZE = 5
LOUVAIN_RES = 1.0


def coherence(gene_set, g2p):
    pw = Counter()
    for g in gene_set:
        for p in g2p.get(g, set()): pw[p] += 1
    if not pw: return None
    total = sum(pw.values())
    return pw.most_common(1)[0][1] / total


def build_gene_pathway_map(mg):
    g2p = defaultdict(set)
    for n, a in mg.graph.nodes(data=True):
        if a.get("node_type") != "reaction": continue
        pws = set(a.get("pathways") or [])
        if not pws: continue
        for s in (a.get("gene_symbols") or []):
            if s: g2p[s.upper()] |= pws
    return g2p


def main():
    from gizmo.export.json_export import read_json
    from gizmo.evidence.mappers import GeneMapper
    from per_patient_master import (load_crohn, load_su_covid, load_erawijantari,
                                       load_gao_ra, load_filbin_covid, load_idh_glioma,
                                       load_tcga_idh_glioma, load_corevitas,
                                       load_kmplot_brca, load_tcga_luad,
                                       load_gse89408_ra, load_hmp2_ibd_cd)
    from stage31_v14b_wgcna_coherence import wgcna_modules

    LOADERS = {
        "Crohn": load_crohn, "Su_COVID": load_su_covid,
        "Erawijantari": load_erawijantari, "Gao_RA": load_gao_ra,
        "Filbin_COVID": load_filbin_covid, "IDH_glioma": load_idh_glioma,
        "TCGA_IDH_glioma": load_tcga_idh_glioma, "CorEvitas_RA": load_corevitas,
        "KMPLOT_BRCA": load_kmplot_brca, "TCGA_LUAD": load_tcga_luad,
        "GSE89408_RA": load_gse89408_ra, "HMP2_IBD_CD": load_hmp2_ibd_cd,
    }

    print("Loading bundle…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    gmap = GeneMapper(mg)
    g2p = build_gene_pathway_map(mg)

    def sym_translate(feat):
        feat_u = str(feat).upper().strip()
        if feat_u in g2p: return feat_u
        try:
            r = gmap.map(feat)
            if r and r[0]:
                a = mg.graph.nodes.get(r[0], {})
                s = (a.get("symbol") or a.get("name") or "").upper()
                if s: return s
        except: pass
        return None

    v6 = json.loads((RESULTS / "stage31_v6_full_modules.json").read_text())
    cells_v6 = v6["by_cohort_design"]

    print("\n" + "=" * 110, flush=True)
    print("STAGE 31 v14c — WGCNA vs GIZMO at SAME gene universe", flush=True)
    print("=" * 110, flush=True)

    results = []
    for cohort, loader in LOADERS.items():
        print(f"\n--- {cohort} ---", flush=True)
        try:
            prot, metab, ylab, common = loader()
        except Exception as exc:
            print(f"  loader failed: {exc}"); continue

        # Build cohort gene universe (genes that ANY cohort feature translates to)
        gene_universe = set()
        for fd in [prot, metab]:
            if fd is None: continue
            features = set().union(*[set(fd[s]) for s in common if s in fd])
            for f in features:
                s = sym_translate(f)
                if s: gene_universe.add(s)
        gene_universe = {g for g in gene_universe if g in g2p}
        if len(gene_universe) < 20:
            print(f"  gene universe too small ({len(gene_universe)})"); continue
        print(f"  Cohort gene universe (mapped to pathway): {len(gene_universe)} genes",
              flush=True)

        # WGCNA modules on cohort features (use prot if available, else metab)
        feat_data = prot if (prot is not None and len(prot) > 0) else metab
        if feat_data is None: continue
        wmods, _ = wgcna_modules(feat_data, common, sym_translate)
        wgcna_coh = []
        for m in wmods:
            genes_in_universe = set(m["genes"]) & gene_universe
            if len(genes_in_universe) < 3: continue
            c = coherence(genes_in_universe, g2p)
            if c is not None: wgcna_coh.append(c)

        # GIZMO modules RESTRICTED to gene universe
        gizmo_coh_restricted = []
        gizmo_coh_full = []
        for cell_key, cd in cells_v6.items():
            if not (cell_key == cohort or cell_key.startswith(cohort + "/")): continue
            for m_data in cd.get("modules", []):
                full_genes = {g.upper() for g in (m_data.get("gene_symbols") or []) if g}
                if len(full_genes) < 5: continue
                # Coherence on FULL gene set (v14b baseline)
                c_full = coherence(full_genes, g2p)
                if c_full is not None: gizmo_coh_full.append(c_full)
                # Coherence RESTRICTED to cohort universe
                restricted = full_genes & gene_universe
                if len(restricted) >= 3:
                    c_res = coherence(restricted, g2p)
                    if c_res is not None: gizmo_coh_restricted.append(c_res)
            break  # one design per cohort

        if not wgcna_coh or not gizmo_coh_restricted:
            print(f"  insufficient modules"); continue

        wgcna_med = np.median(wgcna_coh)
        gizmo_full_med = np.median(gizmo_coh_full) if gizmo_coh_full else None
        gizmo_res_med = np.median(gizmo_coh_restricted)

        print(f"  WGCNA modules (n={len(wgcna_coh)}):              median tpf = {wgcna_med:.3f}",
              flush=True)
        full_str = f"{gizmo_full_med:.3f}" if gizmo_full_med is not None else "n/a"
        print(f"  GIZMO modules FULL (n={len(gizmo_coh_full)}):    median tpf = {full_str}", flush=True)
        print(f"  GIZMO modules RESTRICTED (n={len(gizmo_coh_restricted)}):  "
              f"median tpf = {gizmo_res_med:.3f}", flush=True)
        ratio_restricted = (gizmo_res_med / wgcna_med) if wgcna_med > 0 else None
        print(f"  Ratio (GIZMO-restricted / WGCNA): {ratio_restricted:.2f}", flush=True)

        results.append({
            "cohort": cohort,
            "gene_universe_size": len(gene_universe),
            "wgcna_n": len(wgcna_coh),
            "wgcna_median_tpf": float(wgcna_med),
            "gizmo_full_n": len(gizmo_coh_full),
            "gizmo_full_median_tpf": float(gizmo_full_med) if gizmo_full_med else None,
            "gizmo_restricted_n": len(gizmo_coh_restricted),
            "gizmo_restricted_median_tpf": float(gizmo_res_med),
            "ratio_restricted": float(ratio_restricted) if ratio_restricted else None,
        })

    # Aggregate
    print("\n" + "=" * 110, flush=True)
    ratios = [r["ratio_restricted"] for r in results if r["ratio_restricted"]]
    if ratios:
        wins = sum(1 for r in ratios if r > 1)
        print(f"AGGREGATE — same-gene-universe comparison", flush=True)
        print(f"  GIZMO-restricted/WGCNA median tpf ratio:  mean={np.mean(ratios):.2f}  "
              f"median={np.median(ratios):.2f}", flush=True)
        print(f"  Cohorts where GIZMO-restricted > WGCNA: {wins}/{len(ratios)}", flush=True)

    out_path = RESULTS / "stage31_v14c_wgcna_same_universe.json"
    out_path.write_text(json.dumps({"results": results,
                                       "config": {"corr_thresh": CORR_THRESH,
                                                  "min_module_size": MIN_MODULE_SIZE}},
                                      indent=2, default=str))
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
