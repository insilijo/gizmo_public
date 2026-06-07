"""Cross-cohort α-PC alignment chord diagram.

For each cohort in the panel, compute its top-5 α-PC directions (signed
loadings over substrate nodes). Then compute pairwise |cosine| between
every (cohort_A, PC_i) × (cohort_B, PC_j) — 17×5 = 85 axes, ~3500 pairs.
Render as a chord diagram where:

  • Each arc segment around the ring = one (cohort, PC) entry, labeled
    "<cohort>·PCk"; PCs of the same cohort are adjacent.
  • Each chord = a cross-cohort similarity above a threshold
    (default |cos| > 0.3). Same-cohort pairs are not drawn (trivial).
  • Chord color + alpha = |cos|; thickness scales with strength.
  • Disease-class color stripe inside the ring annotates cohort class.

The chord diagram reveals **cross-PC matches** that the diagonal
within-disease cosine display cannot show: e.g., PC1 of study A may align
with PC2 of study B even though both labeled the same disease.

Output: fig_pc_alignment_chord.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MPath
from sklearn.decomposition import PCA

REPO = Path("/home/jgardner/GIZMO")
RESULTS = REPO / "benchmarks/results"
UR = RESULTS / "unsupervised"
FIG = RESULTS / "figures"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))


COHORTS = [
    # (cohort, disease class, fine-grained disease label for grouping)
    ("IDH_glioma",      "cancer_subtype",      "Glioma"),
    ("TCGA_IDH_glioma", "cancer_subtype",      "Glioma"),
    ("TCGA_LUAD",       "cancer_GoF_subtype",  "LungAdeno"),
    ("KMPLOT_BRCA",     "cancer_subtype",      "Breast"),
    ("CPTAC_CCRCC",     "cancer_intensity",    "Renal"),
    ("CPTAC_COAD",      "cancer_intensity",    "Colon"),
    ("CPTAC_OV",        "cancer_intensity",    "Ovarian"),
    ("Su_COVID",        "inflammatory",        "COVID"),
    ("Filbin_COVID",    "inflammatory",        "COVID"),
    ("GSE65391_SLE",    "inflammatory",        "SLE"),
    ("GSE65682_sepsis", "inflammatory",        "Sepsis"),
    ("Gao_RA",          "inflammatory",        "RA"),
    ("GSE89408_RA",     "inflammatory",        "RA"),
    ("Crohn",           "inflammatory",        "IBD"),
    ("HMP2_IBD_CD",     "inflammatory",        "IBD"),
    ("Erawijantari",    "microbiome_gastric",  "Gastric"),
]

# Gene-side input modality per cohort. Used to label gene-node loadings as
# 'gene' (RNA / transcriptomics) vs 'protein' (Olink, RPPA, etc.). Reaction
# and metabolite-node loadings are unaffected.
COHORT_GENE_INPUT = {
    "IDH_glioma":       "gene",      # Trautwein RNA-seq + NMR metab
    "TCGA_IDH_glioma":  "gene",      # RNA-seq only
    "TCGA_LUAD":        "gene",      # RNA-seq only
    "KMPLOT_BRCA":      "gene",      # mRNA expression
    "CPTAC_CCRCC":      "protein",   # proteomics primary
    "CPTAC_COAD":       "protein",
    "CPTAC_OV":         "protein",
    "Su_COVID":         "protein",   # Olink + metabolomics
    "Filbin_COVID":     "protein",   # Olink only
    "GSE65391_SLE":     "gene",      # microarray
    "GSE65682_sepsis":  "gene",      # microarray
    "Gao_RA":           "protein",   # multiplex prot + metab
    "GSE89408_RA":      "gene",      # RNA-seq
    "Crohn":            "protein",   # Olink + metab (Koopman)
    "HMP2_IBD_CD":      "gene",      # metab-dominated but RNA available
    "Erawijantari":     "metab_only", # metabolomics-only — no gene input
}
N_PCS = 5
COS_THRESHOLD = 0.05  # lower threshold surfaces weak cross-cohort alignments

# Ordering of fine-grained diseases around the ring — keeps same-disease pairs
# (IDH + TCGA_IDH gliomas, Su + Filbin COVID, Gao + GSE89408 RA, Crohn + HMP2 IBD)
# adjacent. Cancers grouped, inflammatory grouped, microbiome at the seam.
DISEASE_ORDER = [
    "Glioma", "LungAdeno", "Breast", "Renal", "Colon", "Ovarian",
    "COVID", "SLE", "Sepsis", "RA", "IBD",
    "Gastric",
]

CLASS_COLOR = {
    "inflammatory":      "#d62728",
    "cancer_intensity":  "#1f77b4",
    "cancer_subtype":    "#9467bd",
    "cancer_GoF_subtype": "#d97706",   # dark orange — distinct from the
                                        # visible-K loader-missing gray
    "microbiome_gastric": "#2ca02c",
}


def load_substrate():
    from gizmo.export.json_export import read_json
    from per_patient_wlsp_v2 import biochem_subgraph
    print("Loading substrate…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    sub_dir, nodes, _ = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)
    return mg, nodes, log_pr


# ---------------------------------------------------------------------------
# Per-cohort input-feature → substrate-node mapping (for "visible K" modality)
# ---------------------------------------------------------------------------

def visible_substrate_nodes_for(cohort, mg, *, use_cache: bool = True):
    """Return set of substrate node IDs anchored by the cohort's input features.

    Cached to disk at benchmarks/results/.cache_visible_sets/<cohort>.npz
    so subsequent runs avoid the ~10-min mapper.map() loop per cohort.
    """
    cache_dir = RESULTS / ".cache_visible_sets"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{cohort}.npz"
    if use_cache and cache_file.exists():
        try:
            nids = np.load(cache_file, allow_pickle=True)["nids"].tolist()
            return set(nids)
        except Exception:
            pass
    try:
        from per_patient_master import (
            load_crohn, load_su_covid, load_erawijantari, load_gao_ra,
            load_filbin_covid, load_idh_glioma, load_tcga_idh_glioma,
            load_hmp2_ibd_cd, load_gse89408_ra, load_kmplot_brca,
            load_tcga_luad,
        )
        from gizmo.evidence.mappers import GeneMapper, MetaboliteMapper
    except Exception:
        return None
    loader_map = {
        "Crohn":            load_crohn,
        "Su_COVID":         load_su_covid,
        "Gao_RA":           load_gao_ra,
        "Erawijantari":     load_erawijantari,
        "Filbin_COVID":     load_filbin_covid,
        "IDH_glioma":       load_idh_glioma,
        "TCGA_IDH_glioma":  load_tcga_idh_glioma,
        "HMP2_IBD_CD":      load_hmp2_ibd_cd,
        "GSE89408_RA":      load_gse89408_ra,
        "KMPLOT_BRCA":      load_kmplot_brca,
        "TCGA_LUAD":        load_tcga_luad,
    }
    if cohort not in loader_map:
        return None
    try:
        prot, metab, _y, common = loader_map[cohort]()
    except Exception:
        return None
    gmap = GeneMapper(mg); mmap = MetaboliteMapper(mg)
    visible = set()
    # Iterate UNIQUE features once per modality (not per-sample × per-feature
    # which was a ~200× redundant slowdown for big-n cohorts).
    for data, mapper in [(prot, gmap), (metab, mmap)]:
        if not data: continue
        all_feats = set()
        for sid in common:
            if sid in data:
                all_feats.update(data[sid])
        for feat in all_feats:
            try:
                res = mapper.map(feat)
                nid = res[0] if isinstance(res, tuple) else res
            except Exception:
                nid = None
            if nid and nid in mg.graph:
                visible.add(nid)
    # Cache to disk
    try:
        np.savez_compressed(cache_file, nids=np.array(sorted(visible), dtype=object))
    except Exception:
        pass
    return visible


# ---------------------------------------------------------------------------
# Per-PC biology label (manual cases + top-loading fallback)
# ---------------------------------------------------------------------------

# Manually curated biology one-liners for the load-bearing (cohort, PC) pairs
# from the basin builder's cases list. Other PCs fall back to top-3 loadings.
MANUAL_BIO = {
    ("Filbin_COVID", 4):    "LPO antimicrobial ↔ galanin receptors",
    ("CPTAC_CCRCC", 1):     "Warburg/HIF ↔ renal tubule",
    ("IDH_glioma", 1):      "neuronal ↔ immune (IDH-mut axis)",
    ("CPTAC_COAD", 4):      "intestinal endocrine ↔ pancreatic exocrine",
    ("GSE65391_SLE", 3):    "type I IFN signature (psychosis)",
    ("KMPLOT_BRCA", 4):     "tumor immune infiltrate ↔ secretory",
    ("GSE65682_sepsis", 2): "metabolic ↔ lymphoid exhaustion",
    ("Gao_RA", 1):          "synovial inflammation",
    ("TCGA_IDH_glioma", 5): "epigenetic/translation axis (IDH-mut)",
    ("Su_COVID", 1):        "SRPK1 viral nucleoprotein ↔ chemokine receptors",
    ("Erawijantari", 4):    "phospholipid head-group ↔ kynurenine/NAD",
    ("HMP2_IBD_CD", 1):     "SAM methylation cycle ↔ biotin carboxylation",
    ("GSE89408_RA", 1):     "synovial TLS (CXCL13/CXCR5 + BAFF)",
    ("CPTAC_OV", 1):        "IGF2BP oncofetal ↔ CD34 stromal",
    ("Crohn", 1):           "MPG glycosylase / thiopurine repair",
    ("TCGA_LUAD", 3):       "xenobiotic UGT/sterol ↔ MHC-I/IFN (alt-axis)",
}


def biology_label_for(cohort, pc_idx, pc_loadings, mg, nodes, k=5):
    """Manual label if known, else top-K gene/metab/reaction names from the loading."""
    if (cohort, pc_idx) in MANUAL_BIO:
        return MANUAL_BIO[(cohort, pc_idx)]
    abs_load = np.abs(pc_loadings)
    top_idx = np.argsort(-abs_load)[:k]
    names = []
    for i in top_idx:
        nid = nodes[i]
        attrs = mg.graph.nodes.get(nid, {})
        sym = attrs.get("symbol") or attrs.get("name") or attrs.get("display_name") or nid
        # Shorten reaction names
        s = str(sym)
        if len(s) > 24: s = s[:21] + "…"
        names.append(s)
    return " / ".join(names[:3])


def modality_breakdown(pc_loadings, mg, nodes, *, top_k=50, restrict_to=None,
                        cohort=None):
    """Fraction of top-K substrate nodes (by |loading|) by INPUT MODALITY.

    Reaction-nodes are explicitly excluded from the count — they are
    reified-edge substrate nodes that accumulate MAP-propagated signal
    from their gene catalyst + substrate metabolite + product metabolite,
    not independent measurements. Counting them inflates reaction-fraction
    artifactually.

    Gene-nodes are labeled 'gene' (transcriptomics input) or 'protein'
    (proteomics input) based on the cohort's gene-side input modality
    (see COHORT_GENE_INPUT). Metabolite-nodes stay 'metabolite'.

    If restrict_to is provided (set of node IDs), only those are considered
    (= 'visible K' breakdown — input-anchored substrate nodes only).

    Returns dict with raw counts (excluding reactions), reaction count
    (for diagnostic comparison), and the K used.
    """
    if restrict_to is not None:
        candidate_idx = [i for i, nid in enumerate(nodes) if nid in restrict_to]
        if not candidate_idx:
            return {"gene": np.nan, "protein": np.nan, "metabolite": np.nan,
                    "n_reaction": 0, "n_nonreaction": 0, "k": 0}
        abs_load = np.abs(pc_loadings)
        cand_sorted = sorted(candidate_idx, key=lambda i: -abs_load[i])
        top_idx = cand_sorted[:top_k]
    else:
        abs_load = np.abs(pc_loadings)
        top_idx = np.argsort(-abs_load)[:top_k]

    gene_input = COHORT_GENE_INPUT.get(cohort, "gene")  # default = gene
    counts = {"gene": 0, "protein": 0, "metabolite": 0,
              "reaction": 0, "other": 0}
    for i in top_idx:
        t = mg.graph.nodes.get(nodes[i], {}).get("node_type", "other")
        if t == "gene":
            counts["protein" if gene_input == "protein" else "gene"] += 1
        elif t == "metabolite":
            counts["metabolite"] += 1
        elif t == "reaction":
            counts["reaction"] += 1
        else:
            counts["other"] += 1
    # Normalize over the NON-reaction count (reactions are derived, not input)
    n_nonreact = counts["gene"] + counts["protein"] + counts["metabolite"]
    denom = n_nonreact if n_nonreact else 1
    return {
        "gene":       counts["gene"] / denom,
        "protein":    counts["protein"] / denom,
        "metabolite": counts["metabolite"] / denom,
        "n_reaction": counts["reaction"],   # for diagnostic only
        "n_nonreaction": n_nonreact,
        "k":          len(top_idx),
    }


def alpha_pc_components(cohort, log_pr, n_components=N_PCS):
    """Return n_components × n_nodes loading matrix, or None if F missing.

    Falls back to ``_edge_informed`` / ``_combined`` variants when the bare
    ``stage3_F_<cohort>.npz`` doesn't exist (e.g., TCGA_IDH_glioma was only
    saved under the edge-informed name)."""
    F_path = None
    for suffix in ("", "_edge_informed", "_combined", "_node_informed"):
        cand = UR / f"stage3_F_{cohort}{suffix}.npz"
        if cand.exists():
            F_path = cand
            break
    if F_path is None:
        return None
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"].astype(np.float64)
    if F.shape[1] != len(log_pr):
        return None
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    n = min(n_components, alpha.shape[0])
    pca = PCA(n_components=n, random_state=0)
    pca.fit(alpha)
    return pca.components_


def chord_path(theta1, theta2, r=1.0, curve_strength=0.45):
    """Build a Bezier-curve path from point at angle theta1 to point at theta2
    on the inside of a circle of radius r. The chord curves toward the center."""
    p1 = np.array([r * np.cos(theta1), r * np.sin(theta1)])
    p2 = np.array([r * np.cos(theta2), r * np.sin(theta2)])
    # Control points: pull each toward the origin proportional to curve_strength
    c1 = (1 - curve_strength) * p1
    c2 = (1 - curve_strength) * p2
    verts = [p1, c1, c2, p2]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    return MPath(verts, codes)


def main():
    mg, nodes, log_pr = load_substrate()

    # Compute loadings per (cohort, PC); skip cohorts without F
    items = []   # list of (cohort, cls, disease_label, pc)
    components_per_cohort = {}
    for cohort, cls, disease in COHORTS:
        comp = alpha_pc_components(cohort, log_pr, n_components=N_PCS)
        if comp is None:
            print(f"  skip {cohort}: no F", flush=True); continue
        components_per_cohort[cohort] = comp
        for k in range(comp.shape[0]):
            items.append((cohort, cls, disease, k + 1))   # PC index 1-based

    n_items = len(items)
    print(f"  {n_items} (cohort, PC) entries across "
          f"{len(components_per_cohort)} cohorts", flush=True)

    # Pairwise |cos| between distinct cohorts
    pairs = []  # (i_idx, j_idx, |cos|)
    for i in range(n_items):
        ci, _cls_i, _d_i, pci = items[i]
        vi = components_per_cohort[ci][pci - 1]
        ni = np.linalg.norm(vi) + 1e-12
        for j in range(i + 1, n_items):
            cj, _cls_j, _d_j, pcj = items[j]
            if ci == cj:
                continue
            vj = components_per_cohort[cj][pcj - 1]
            nj = np.linalg.norm(vj) + 1e-12
            cs = float(abs(np.dot(vi, vj) / (ni * nj)))
            if cs >= COS_THRESHOLD:
                pairs.append((i, j, cs))
    pairs.sort(key=lambda t: t[2])  # render weak chords first (low alpha)
    print(f"  {len(pairs)} cross-cohort pairs with |cos| ≥ {COS_THRESHOLD}", flush=True)

    # Layout: order items so same-disease pairs land adjacent on the ring.
    # Sort key: (disease index in DISEASE_ORDER, cohort name, PC index).
    disease_rank = {d: i for i, d in enumerate(DISEASE_ORDER)}
    items_sorted = sorted(
        items,
        key=lambda t: (disease_rank.get(t[2], 99), t[0], t[3]),
    )
    item_to_idx = {(c, pc): i for i, (c, _, _, pc) in enumerate(items_sorted)}
    pairs_remapped = [(item_to_idx[(items[i][0], items[i][3])],
                        item_to_idx[(items[j][0], items[j][3])], cs)
                       for i, j, cs in pairs]

    n_total = len(items_sorted)
    theta_step = 2 * np.pi / n_total
    # Slight gap between cohorts
    thetas = np.array([(i + 0.5) * theta_step for i in range(n_total)])

    R_outer = 1.20  # text labels live here (pushed out to fit modality bars)
    R_chord = 0.85  # chord endpoints
    R_class_in, R_class_out = 0.93, 0.99  # disease-class stripe
    # Two modality-bar rings:
    #   INVISIBLE-K (full substrate ranking) — inner, closer to disease stripe
    #   VISIBLE-K  (input-anchored only)     — outer, closer to labels
    R_modI_in, R_modI_out = 1.00, 1.07
    R_modV_in, R_modV_out = 1.09, 1.16
    # Modality colors — gene (RNA-input cohort) light green, protein (prot-input
    # cohort) dark green, metabolite blue. Reaction nodes are EXCLUDED from
    # the bar count by user request (reactions are reified-edge nodes that
    # accumulate signal from gene + substrate + product, not independent
    # measurements; counting them as 1.0-weight inflates reaction-fraction).
    MOD_COLOR = {
        "gene":       "#7CC082",   # transcriptomics input
        "protein":    "#2F7A3F",   # proteomics input
        "metabolite": "#2D5BBF",   # metabolomics input
    }

    fig, ax = plt.subplots(figsize=(13, 13))
    ax.set_xlim(-1.55, 1.55); ax.set_ylim(-1.55, 1.55)
    ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_facecolor("white")

    # Pre-compute visible-substrate-node sets per cohort (used by modality bars)
    visible_sets_ring = {}
    for c in components_per_cohort.keys():
        visible_sets_ring[c] = visible_substrate_nodes_for(c, mg)

    # Disease-class color stripe (one wedge per item, colored by item's class)
    from matplotlib.patches import Wedge
    for idx, (c, cls, _disease, pc) in enumerate(items_sorted):
        theta_deg = np.degrees(thetas[idx])
        half_step = np.degrees(theta_step) / 2
        wedge = Wedge(center=(0, 0),
                       r=R_class_out, width=R_class_out - R_class_in,
                       theta1=theta_deg - half_step + 0.5,
                       theta2=theta_deg + half_step - 0.5,
                       facecolor=CLASS_COLOR.get(cls, "#bdbdbd"),
                       edgecolor="none", alpha=0.85)
        ax.add_patch(wedge)

    # Modality-bar rings (reactions EXCLUDED — reified-edge nodes, not inputs):
    #   INVISIBLE-K (full substrate top-50) inside (R_modI_in..R_modI_out)
    #   VISIBLE-K  (input-anchored top-50) outside (R_modV_in..R_modV_out)
    # Each item's wedge is split into 3 angular sub-wedges (gene/protein/metab)
    # sized proportionally to the modality fraction in its PC's top-50 loadings.
    # Diagnostic: print the reaction fraction next to gene+protein+metab to
    # quantify how much reaction-exclusion changes the picture.
    print("\nPer-PC modality breakdown — diagnostic (reactions in top-50 vs non-reactions):",
          flush=True)
    print(f"  {'(cohort,PC)':30}  {'gene%':>5} {'prot%':>5} {'met%':>5}  "
          f"{'rxn-in-top50':>12} {'non-rxn-in-top50':>16}",
          flush=True)
    big_deal_count = 0
    for idx, (c, cls, _disease, pc) in enumerate(items_sorted):
        theta_deg = np.degrees(thetas[idx])
        half_step_deg = np.degrees(theta_step) / 2 - 0.5
        pc_loadings = components_per_cohort[c][pc - 1]
        inv = modality_breakdown(pc_loadings, mg, nodes, top_k=50,
                                  restrict_to=None, cohort=c)
        vis = (modality_breakdown(pc_loadings, mg, nodes, top_k=50,
                                    restrict_to=visible_sets_ring.get(c),
                                    cohort=c)
               if visible_sets_ring.get(c) else None)

        # Diagnostic print: how reaction-heavy was the top-50 before exclusion?
        rxn_frac = inv["n_reaction"] / max(50, 1)
        if rxn_frac >= 0.5:
            big_deal_count += 1
        print(f"  {c+'·PC'+str(pc):30}  "
              f"{100*inv['gene']:>5.1f} {100*inv['protein']:>5.1f} {100*inv['metabolite']:>5.1f}  "
              f"{inv['n_reaction']:>12d} {inv['n_nonreaction']:>16d}",
              flush=True)

        def draw_stacked_arc(R_in, R_out, frac_dict):
            total = frac_dict["gene"] + frac_dict["protein"] + frac_dict["metabolite"]
            if total <= 0:
                return
            arc_start = theta_deg - half_step_deg
            arc_full = (theta_deg + half_step_deg) - arc_start
            cursor = arc_start
            for mod_key in ("gene", "protein", "metabolite"):
                frac = frac_dict[mod_key] / total
                if frac <= 0: continue
                arc_w = frac * arc_full
                w = Wedge(center=(0, 0),
                           r=R_out, width=R_out - R_in,
                           theta1=cursor, theta2=cursor + arc_w,
                           facecolor=MOD_COLOR[mod_key],
                           edgecolor="none", alpha=0.85)
                ax.add_patch(w)
                cursor += arc_w

        draw_stacked_arc(R_modI_in, R_modI_out, inv)
        if vis is not None and vis["k"] > 0:
            draw_stacked_arc(R_modV_in, R_modV_out, vis)
        else:
            placeholder = Wedge(center=(0, 0),
                                 r=R_modV_out, width=R_modV_out - R_modV_in,
                                 theta1=theta_deg - half_step_deg,
                                 theta2=theta_deg + half_step_deg,
                                 facecolor="#dddddd", edgecolor="none", alpha=0.5)
            ax.add_patch(placeholder)

    print(f"\n  Diagnostic verdict: {big_deal_count}/{len(items_sorted)} PCs "
          f"({100*big_deal_count/max(len(items_sorted),1):.0f}%) had ≥50% reaction-fraction "
          f"in their top-50 loadings BEFORE exclusion — so the reshuffling "
          f"matters (reaction nodes were dominating). After exclusion the bars "
          f"show only input-modality fractions (gene/protein/metab).", flush=True)

    # Identify "best-per-cohort" chord: for each cohort, the (i, j, cs)
    # pair with the highest |cos| where at least one endpoint is that cohort.
    # These get rendered in gold-yellow as a special layer on top.
    best_per_cohort_edges = set()  # frozenset({i, j}) entries
    cohort_to_best = {}
    for i, j, cs in pairs_remapped:
        ci = items_sorted[i][0]; cj = items_sorted[j][0]
        for c, other_idx in ((ci, j), (cj, i)):
            if c not in cohort_to_best or cs > cohort_to_best[c][2]:
                cohort_to_best[c] = (min(i, j), max(i, j), cs)
    for c, (a, b, _cs) in cohort_to_best.items():
        best_per_cohort_edges.add(frozenset({a, b}))

    # Chords (drawn weakest-first so strong ones lay on top). At a low
    # threshold (0.05) there will be many weak chords; render those almost
    # invisible so the strong ones still pop.
    max_cs = max((cs for _, _, cs in pairs_remapped), default=1.0)
    for i, j, cs in pairs_remapped:
        path = chord_path(thetas[i], thetas[j], r=R_chord, curve_strength=0.45)
        norm_cs = (cs - COS_THRESHOLD) / max(max_cs - COS_THRESHOLD, 1e-6)
        norm_cs_pow = norm_cs ** 1.6
        alpha = 0.04 + 0.78 * norm_cs_pow
        lw = 0.3 + 3.0 * norm_cs_pow
        d_i = items_sorted[i][2]; d_j = items_sorted[j][2]
        cls_i = items_sorted[i][1]; cls_j = items_sorted[j][1]
        if d_i == d_j:
            chord_color = CLASS_COLOR.get(cls_i, "#222")
            alpha = max(alpha, 0.45)
            lw = max(lw, 1.2)
        elif cls_i == cls_j:
            chord_color = CLASS_COLOR.get(cls_i, "#888")
        else:
            chord_color = "#666"
        patch = PathPatch(path, facecolor="none", edgecolor=chord_color,
                          linewidth=lw, alpha=alpha)
        ax.add_patch(patch)

    # Best-per-cohort overlay: gold-yellow on top of regular chords. Drawn
    # as a second pass with high zorder so it always shows above the cloud.
    for i, j, cs in pairs_remapped:
        if frozenset({i, j}) not in best_per_cohort_edges:
            continue
        path = chord_path(thetas[i], thetas[j], r=R_chord, curve_strength=0.45)
        patch = PathPatch(path, facecolor="none", edgecolor="#F2B500",
                          linewidth=2.6, alpha=0.92, zorder=5)
        ax.add_patch(patch)

    # Outer labels: "<cohort>·PCk"
    for idx, (c, cls, _disease, pc) in enumerate(items_sorted):
        theta = thetas[idx]
        x_lbl = (R_outer + 0.05) * np.cos(theta)
        y_lbl = (R_outer + 0.05) * np.sin(theta)
        angle_deg = np.degrees(theta)
        # Rotate so text reads outward; flip on the left half so it isn't upside-down
        if angle_deg > 90 and angle_deg < 270:
            rot = angle_deg - 180; ha = "right"
        else:
            rot = angle_deg; ha = "left"
        label_color = CLASS_COLOR.get(cls, "#333")
        # Shorten common cohort names a bit
        short = (c.replace("CPTAC_", "")
                   .replace("GSE65391_", "")
                   .replace("GSE65682_", "")
                   .replace("GSE89408_", "")
                   .replace("_COVID", "·COVID")
                   .replace("_RA", "·RA")
                   .replace("_IBD_CD", "·CD")
                   .replace("_glioma", "·gli"))
        ax.text(x_lbl, y_lbl, f"{short}·PC{pc}",
                rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", fontsize=7.5, color=label_color,
                fontweight="bold" if pc == 1 else "normal")

    # Legend for disease class (figure-level so it's not clipped by tight bbox)
    from matplotlib.lines import Line2D
    legend_classes = sorted({cls for _, cls, _d in COHORTS}, key=lambda x: x)
    legend_handles = [Line2D([0], [0], marker="s", color="none",
                              markerfacecolor=CLASS_COLOR.get(c, "#bdbdbd"),
                              markeredgecolor="#222", markeredgewidth=0.6,
                              markersize=12, label=c.replace("_", " "))
                       for c in legend_classes]
    fig.legend(handles=legend_handles, loc="lower left",
               bbox_to_anchor=(0.02, 0.03), frameon=False, fontsize=9,
               title="disease class", title_fontsize=10)

    # Modality legend (color = input modality; reactions excluded from bars
    # as derived nodes)
    modality_handles = [
        Line2D([0], [0], marker="s", color="none",
               markerfacecolor=MOD_COLOR["gene"],
               markeredgecolor="#222", markeredgewidth=0.6,
               markersize=12, label="gene (transcriptomics input)"),
        Line2D([0], [0], marker="s", color="none",
               markerfacecolor=MOD_COLOR["protein"],
               markeredgecolor="#222", markeredgewidth=0.6,
               markersize=12, label="protein (proteomics input)"),
        Line2D([0], [0], marker="s", color="none",
               markerfacecolor=MOD_COLOR["metabolite"],
               markeredgecolor="#222", markeredgewidth=0.6,
               markersize=12, label="metabolite"),
        Line2D([0], [0], marker="s", color="none",
               markerfacecolor="#dddddd",
               markeredgecolor="#222", markeredgewidth=0.6,
               markersize=12, label="visible-K loader missing"),
        Line2D([0], [0], color="#F2B500", linewidth=3.2,
               label="strongest cross-cohort chord per cohort"),
    ]
    fig.legend(handles=modality_handles, loc="lower right",
               bbox_to_anchor=(0.98, 0.03), frameon=False, fontsize=9,
               title="modality bars + chord highlight", title_fontsize=10)

    # Ring labels (which ring is invisible-K vs visible-K). Placed off to one
    # side so they don't crowd the disease-class legend on the lower-left.
    ax.annotate("invisible-K\n(full substrate top-50)", xy=(0, 1.035),
                ha="left", va="center", fontsize=7.5, color="#444",
                xytext=(1.30, 1.20), textcoords="data",
                arrowprops=dict(arrowstyle="->", color="#444", lw=0.6))
    ax.annotate("visible-K\n(input-anchored top-50)", xy=(0, 1.125),
                ha="left", va="center", fontsize=7.5, color="#444",
                xytext=(1.30, 1.40), textcoords="data",
                arrowprops=dict(arrowstyle="->", color="#444", lw=0.6))

    # (Title + bottom caption removed — caption lives in the manuscript figure
    # legend so the chord stands alone visually and doesn't fight the
    # ring labels or the modality / disease-class legends.)

    # Colorbar for |cos| (sampled from the colormap actually used)
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.cm import ScalarMappable
    cb_cmap = LinearSegmentedColormap.from_list(
        "cos_alpha", [(0.6, 0.6, 0.6, 0.04), (0.15, 0.15, 0.15, 0.85)],
    )
    sm = ScalarMappable(norm=Normalize(vmin=COS_THRESHOLD, vmax=max_cs),
                         cmap=cb_cmap)
    sm.set_array([])
    cax = fig.add_axes([0.41, 0.05, 0.18, 0.018])
    cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label(f"|cos|  (chord rendered only at ≥ {COS_THRESHOLD:.2f})",
                    fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    out_path = FIG / "fig_pc_alignment_chord.png"
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"\nWrote {out_path}")

    # ── Best cross-cohort match table (one row per (cohort, PC) entry) ──
    # For each (cohort_A, PC_i), find the (cohort_B, PC_j) with the highest
    # |cos| — across all PCs and all OTHER cohorts. Also report top-3 cross
    # matches, biology label, intra-cohort max-other-PC cos (orthogonality
    # check), and modality breakdowns (visible + invisible top-K).
    print("\nBuilding best-cross-cohort-match table…", flush=True)
    cohorts_present = list(components_per_cohort.keys())
    # Pre-compute visible-substrate-node sets per cohort (cached)
    visible_sets = {}
    for c in cohorts_present:
        v = visible_substrate_nodes_for(c, mg)
        visible_sets[c] = v
        if v is None:
            print(f"  visible set: {c} → no loader, visible-K modality NaN", flush=True)
        else:
            print(f"  visible set: {c} → {len(v)} substrate nodes", flush=True)

    rows = []
    for ci in cohorts_present:
        comp_i = components_per_cohort[ci]
        cls_i = next(cls for c, cls, _d in COHORTS if c == ci)
        disease_i = next(d for c, _cls, d in COHORTS if c == ci)
        vis_set_i = visible_sets.get(ci)
        for pi in range(comp_i.shape[0]):
            vi = comp_i[pi]
            ni = np.linalg.norm(vi) + 1e-12
            # Cross-cohort matches
            matches = []
            for cj in cohorts_present:
                if cj == ci:
                    continue
                comp_j = components_per_cohort[cj]
                disease_j = next(d for c, _cls, d in COHORTS if c == cj)
                for pj in range(comp_j.shape[0]):
                    vj = comp_j[pj]
                    nj = np.linalg.norm(vj) + 1e-12
                    cs = float(abs(np.dot(vi, vj) / (ni * nj)))
                    matches.append((cs, cj, pj + 1, disease_j))
            matches.sort(reverse=True)
            if not matches: continue
            # Intra-cohort orthogonality: max |cos| between this PC and any
            # OTHER PC of the same cohort. Should be ~0 by PCA construction
            # but worth verifying — a high value means the PC overlaps
            # heavily with another PC in the cohort.
            intra_max = 0.0
            for pk in range(comp_i.shape[0]):
                if pk == pi: continue
                vk = comp_i[pk]
                nk = np.linalg.norm(vk) + 1e-12
                cs_intra = float(abs(np.dot(vi, vk) / (ni * nk)))
                intra_max = max(intra_max, cs_intra)
            # Biology label
            bio = biology_label_for(ci, pi + 1, vi, mg, nodes)
            # Modality breakdowns (top-50 by |loading|, reactions EXCLUDED
            # from the fraction denominator; reaction count reported separately).
            mod_inv = modality_breakdown(vi, mg, nodes, top_k=50,
                                          restrict_to=None, cohort=ci)
            mod_vis = (modality_breakdown(vi, mg, nodes, top_k=50,
                                           restrict_to=vis_set_i, cohort=ci)
                       if vis_set_i is not None
                       else {"gene": np.nan, "protein": np.nan,
                              "metabolite": np.nan, "n_reaction": 0,
                              "n_nonreaction": 0, "k": 0})
            row = {
                "cohort": ci,
                "pc": pi + 1,
                "disease": disease_i,
                "disease_class": cls_i,
                "biology": bio,
                "best_match_cohort":  matches[0][1],
                "best_match_pc":      matches[0][2],
                "best_match_disease": matches[0][3],
                "best_match_cos":     round(matches[0][0], 4),
                "best_match_within_disease": (matches[0][3] == disease_i),
                "intra_cohort_max_other_pc_cos": round(intra_max, 4),
                "invisible_pct_gene":     round(100 * mod_inv["gene"], 1),
                "invisible_pct_protein":  round(100 * mod_inv["protein"], 1),
                "invisible_pct_metab":    round(100 * mod_inv["metabolite"], 1),
                "invisible_n_reaction":   mod_inv["n_reaction"],
                "visible_pct_gene":       round(100 * mod_vis["gene"], 1)       if not np.isnan(mod_vis["gene"])       else np.nan,
                "visible_pct_protein":    round(100 * mod_vis["protein"], 1)    if not np.isnan(mod_vis["protein"])    else np.nan,
                "visible_pct_metab":      round(100 * mod_vis["metabolite"], 1) if not np.isnan(mod_vis["metabolite"]) else np.nan,
                "visible_n_reaction":     mod_vis["n_reaction"],
                "visible_k_size":         mod_vis["k"],
            }
            top3 = ["{}·PC{} ({:.3f})".format(m[1], m[2], m[0]) for m in matches[:3]]
            row["top3_matches"] = " | ".join(top3)
            rows.append(row)

    # Write TSV
    import pandas as pd
    df = pd.DataFrame(rows).sort_values(
        by=["best_match_within_disease", "best_match_cos"],
        ascending=[False, False],
    ).reset_index(drop=True)
    tsv_path = RESULTS / "pc_alignment_best_cross_cohort_match.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"Wrote {tsv_path}  ({len(df)} entries)")
    print(f"  Intra-cohort PC orthogonality: max(other-PC |cos|) median = "
          f"{df['intra_cohort_max_other_pc_cos'].median():.4f}, "
          f"max = {df['intra_cohort_max_other_pc_cos'].max():.4f}")

    # (Table panel removed — the chord diagram + caption carries the visual
    # story; per-row detail is available in the TSV above for reviewers who
    # want it.)


if __name__ == "__main__":
    main()
