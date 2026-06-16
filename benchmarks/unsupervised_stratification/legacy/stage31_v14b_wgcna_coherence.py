"""Stage 31 v14b — WGCNA-style modular coherence comparator (right null for v14).

v14 compared GIZMO modules to TopPR flat-list — a category-error null because
TopPR doesn't produce modules. Right null: another method that ALSO produces
modules from the same data.

WGCNA-equivalent here (simplified):
  1. Compute Spearman correlation between cohort features across patients
  2. Build weighted adjacency: w(i,j) = |corr(i,j)| if > threshold else 0
  3. Louvain community detection on the correlation graph
  4. Each community = WGCNA-style module (data-driven, no substrate prior)

For each cohort:
  - Get GIZMO modules (from v6 partition) → pathway coherence per module
  - Get WGCNA modules (from cohort feature correlations) → pathway coherence per
    module via feat→gene mapping
  - Compare: per-module top_pathway_frac, n_distinct_pathways

Test: are GIZMO modules more pathway-coherent than WGCNA modules ON THE SAME DATA?
If yes (3-5× difference like v14 vs TopPR), the modular-coherence claim survives.
If no, the claim was vs the wrong null.

Output: stage31_v14b_wgcna_coherence.json
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
from scipy.stats import spearmanr
import networkx as nx

REPO = Path(__file__).resolve().parent.parent.parent
RESULTS = REPO / "benchmarks" / "results" / "unsupervised"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))
sys.path.insert(0, str(REPO / "benchmarks" / "unsupervised_stratification"))


CORR_THRESH = 0.5         # |Spearman r| above this counts as an edge in WGCNA graph
MIN_MODULE_SIZE = 5       # min features per WGCNA module
LOUVAIN_RES = 1.0


def coherence_metrics(gene_set, gene_to_pathways):
    pw_counter = Counter()
    for g in gene_set:
        pws = gene_to_pathways.get(g, set())
        for p in pws: pw_counter[p] += 1
    if not pw_counter:
        return {"n_genes": len(gene_set), "n_distinct_pathways": 0,
                  "top_pathway_frac": None, "entropy": None}
    total = sum(pw_counter.values())
    top_count = pw_counter.most_common(1)[0][1]
    probs = np.array([c/total for c in pw_counter.values()])
    entropy = float(-(probs * np.log(probs + 1e-12)).sum())
    return {"n_genes": len(gene_set), "n_distinct_pathways": len(pw_counter),
              "top_pathway_frac": top_count / total, "entropy": entropy,
              "top_pathway_id": pw_counter.most_common(1)[0][0]}


def build_gene_pathway_map(mg):
    g2p = defaultdict(set)
    for n, a in mg.graph.nodes(data=True):
        if a.get("node_type") != "reaction": continue
        pws = set(a.get("pathways") or [])
        if not pws: continue
        for s in (a.get("gene_symbols") or []):
            if s: g2p[s.upper()] |= pws
    return g2p


def wgcna_modules(feature_data, common, sym_translate, min_size=MIN_MODULE_SIZE,
                    thresh=CORR_THRESH):
    """Build WGCNA-style modules from cohort feature correlations + Louvain."""
    # Get feature x sample matrix
    features = sorted(set().union(*[set(feature_data[s]) for s in common if s in feature_data]))
    if len(features) < 20: return [], []
    X = np.array([[feature_data[s].get(f, np.nan) for f in features] for s in common
                    if s in feature_data], dtype=float)
    # Drop features with too many NAs
    nan_frac = np.isnan(X).mean(axis=0)
    keep = nan_frac < 0.5
    X = X[:, keep]; features = [f for f, k in zip(features, keep) if k]
    if X.shape[1] < 20: return [], []
    # Replace NA with col-mean
    col_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_means, inds[1])

    # Spearman correlation matrix (feature × feature)
    # For large feature counts, this is expensive; truncate to top 2000 by variance
    if X.shape[1] > 2000:
        var = np.var(X, axis=0)
        top_idx = np.argsort(-var)[:2000]
        X = X[:, top_idx]
        features = [features[i] for i in top_idx]

    n_feat = X.shape[1]
    # Spearman ranks
    ranks = np.zeros_like(X)
    for j in range(n_feat):
        ranks[:, j] = np.argsort(np.argsort(X[:, j]))
    # Correlation via ranks
    ranks_z = (ranks - ranks.mean(axis=0)) / (ranks.std(axis=0) + 1e-12)
    corr = (ranks_z.T @ ranks_z) / ranks_z.shape[0]
    np.fill_diagonal(corr, 0.0)
    A = np.abs(corr)

    # Louvain on weighted graph (only edges above threshold)
    edges = [(i, j, A[i, j]) for i in range(n_feat) for j in range(i+1, n_feat)
              if A[i, j] >= thresh]
    if len(edges) < 50: return [], features
    g = nx.Graph()
    g.add_weighted_edges_from(edges)
    coms = nx.community.louvain_communities(g, weight="weight", seed=0,
                                                resolution=LOUVAIN_RES)
    modules = []
    for c in coms:
        if len(c) >= min_size:
            module_features = [features[i] for i in c]
            # Translate features → gene symbols
            module_genes = set()
            for f in module_features:
                sym = sym_translate(f)
                if sym: module_genes.add(sym.upper())
            if module_genes:
                modules.append({"n_features": len(module_features),
                                "features": module_features[:10],
                                "genes": sorted(module_genes)})
    return modules, features


def main():
    from gizmo.export.json_export import read_json
    from gizmo.evidence.mappers import GeneMapper, MetaboliteMapper
    from per_patient_master import (load_crohn, load_su_covid, load_erawijantari,
                                       load_gao_ra, load_filbin_covid, load_idh_glioma,
                                       load_tcga_idh_glioma, load_corevitas,
                                       load_kmplot_brca, load_tcga_luad,
                                       load_gse89408_ra, load_hmp2_ibd_cd)

    LOADERS = {
        "Crohn": load_crohn, "Su_COVID": load_su_covid,
        "Erawijantari": load_erawijantari, "Gao_RA": load_gao_ra,
        "Filbin_COVID": load_filbin_covid, "IDH_glioma": load_idh_glioma,
        "TCGA_IDH_glioma": load_tcga_idh_glioma, "CorEvitas_RA": load_corevitas,
        "KMPLOT_BRCA": load_kmplot_brca, "TCGA_LUAD": load_tcga_luad,
        "GSE89408_RA": load_gse89408_ra, "HMP2_IBD_CD": load_hmp2_ibd_cd,
    }

    print("Loading bundle + gene-pathway map…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    gmap = GeneMapper(mg); mmap = MetaboliteMapper(mg)
    g2p = build_gene_pathway_map(mg)
    print(f"  {len(g2p)} genes have pathway annotations", flush=True)

    # Translation: feature name → gene symbol (rough)
    def sym_translate(feat):
        # Try as gene symbol first
        feat_u = str(feat).upper().strip()
        if feat_u in g2p: return feat_u
        # Try gene mapper
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

    print("\n" + "=" * 115, flush=True)
    print("STAGE 31 v14b — WGCNA-style modular coherence comparator", flush=True)
    print("=" * 115, flush=True)

    results = []
    for cohort, loader in LOADERS.items():
        print(f"\n--- {cohort} ---", flush=True)
        try:
            prot, metab, ylab, common = loader()
        except Exception as exc:
            print(f"  loader failed: {exc}"); continue

        # WGCNA on proteomics (if available)
        cohort_results = {"cohort": cohort, "designs": {}}
        for design_name, feat_data in [("prot", prot), ("metab", metab)]:
            if feat_data is None or len(feat_data) < 20: continue
            wmods, feats_used = wgcna_modules(feat_data, common, sym_translate)
            if not wmods: continue
            wgcna_coh = []
            for m in wmods:
                c = coherence_metrics(set(m["genes"]), g2p)
                if c["top_pathway_frac"] is not None:
                    wgcna_coh.append(c["top_pathway_frac"])
            if not wgcna_coh: continue
            wgcna_mean_tpf = float(np.mean(wgcna_coh))
            wgcna_median_tpf = float(np.median(wgcna_coh))
            wgcna_n_mods = len(wgcna_coh)

            # Compare to GIZMO modules for this cohort (use any design)
            gizmo_tpfs = []
            for cell_key, cd in cells_v6.items():
                if not cell_key.startswith(cohort + "/") and cell_key != cohort: continue
                for m_data in cd.get("modules", []):
                    genes = set((m_data.get("gene_symbols") or []))
                    genes = {g.upper() for g in genes if g}
                    if len(genes) < 5: continue
                    c = coherence_metrics(genes, g2p)
                    if c["top_pathway_frac"] is not None:
                        gizmo_tpfs.append(c["top_pathway_frac"])
                break  # one design's modules per cohort for comparison
            gizmo_mean_tpf = float(np.mean(gizmo_tpfs)) if gizmo_tpfs else None
            gizmo_median_tpf = float(np.median(gizmo_tpfs)) if gizmo_tpfs else None
            gizmo_n_mods = len(gizmo_tpfs)

            ratio = (gizmo_median_tpf / wgcna_median_tpf
                       if gizmo_median_tpf and wgcna_median_tpf and wgcna_median_tpf > 0
                       else None)
            print(f"  {design_name:<6}: WGCNA n_mod={wgcna_n_mods:3d} mean_tpf={wgcna_mean_tpf:.3f} "
                  f"median_tpf={wgcna_median_tpf:.3f}  |  "
                  f"GIZMO n_mod={gizmo_n_mods:3d} median_tpf={gizmo_median_tpf or 0:.3f}  "
                  f"ratio(G/W) = {ratio if ratio else 'n/a':.2f}",
                  flush=True)

            cohort_results["designs"][design_name] = {
                "wgcna_n_modules": wgcna_n_mods,
                "wgcna_mean_top_pathway_frac": wgcna_mean_tpf,
                "wgcna_median_top_pathway_frac": wgcna_median_tpf,
                "gizmo_n_modules": gizmo_n_mods,
                "gizmo_median_top_pathway_frac": gizmo_median_tpf,
                "ratio_gizmo_over_wgcna": ratio,
            }
        results.append(cohort_results)

    # Aggregate
    print("\n" + "=" * 115, flush=True)
    ratios = []
    for r in results:
        for design_name, d in r.get("designs", {}).items():
            if d.get("ratio_gizmo_over_wgcna"):
                ratios.append(d["ratio_gizmo_over_wgcna"])
    if ratios:
        print(f"AGGREGATE — GIZMO/WGCNA median top_pathway_frac ratio:", flush=True)
        print(f"  mean = {np.mean(ratios):.2f}", flush=True)
        print(f"  median = {np.median(ratios):.2f}", flush=True)
        print(f"  cohorts where GIZMO > WGCNA: {sum(1 for r in ratios if r > 1)}/{len(ratios)}",
              flush=True)

    out_path = RESULTS / "stage31_v14b_wgcna_coherence.json"
    out_path.write_text(json.dumps({"results": results,
                                       "config": {"corr_thresh": CORR_THRESH,
                                                  "min_module_size": MIN_MODULE_SIZE}},
                                      indent=2, default=str))
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
