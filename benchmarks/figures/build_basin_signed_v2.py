"""Polished signed-basin figures — back to the original two-pole layout
(fig_signed_basin_* style) with pathway boxes + strip + boxplots added.

Layout: single combined axes, + basin nodes anchored at x=+1.6, − basin at x=−1.6,
bridges in the middle. Spring layout with bias.
"""
from __future__ import annotations

import sys
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA

REPO = Path("/home/jgardner/GIZMO")
RESULTS = REPO / "benchmarks/results"
UR = RESULTS / "unsupervised"
FIG = RESULTS / "figures"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))


def load_substrate():
    from gizmo.export.json_export import read_json
    from per_patient_wlsp_v2 import biochem_subgraph
    print("Loading substrate…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    sub_dir, nodes, nid_idx = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)
    return mg, sub, nodes, nid_idx, log_pr


def find_F_path(cohort):
    for cand in [UR / f"stage3_F_{cohort}.npz",
                 UR / f"stage3_F_{cohort}_combined.npz",
                 UR / f"stage3_F_{cohort}_edge_informed.npz",
                 UR / f"stage3_F_{cohort}_node_informed.npz"]:
        if cand.exists():
            return cand
    return None


def decompose_unit_norm(F, log_pr):
    F_norm = np.linalg.norm(F, axis=1, keepdims=True) + 1e-12
    F = F / F_norm
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F.mean(axis=1, keepdims=True)
    cov = ((F - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F - F_mean - beta[:, None] * (x - x_mean)[None, :]
    return beta, alpha


def short_name(nid, mg, max_len=24):
    attrs = mg.graph.nodes.get(nid, {})
    name = attrs.get("name") or attrs.get("symbol") or nid
    return str(name)[:max_len], attrs.get("node_type", "?")


def steiner_expand(seeds, basin_sub, max_extra=40):
    """Pairwise shortest paths within basin_sub. Returns (all_nodes, bridges)."""
    seeds = list(seeds)
    if len(seeds) < 2:
        return set(seeds), set()
    bridges = set()
    for i, s1 in enumerate(seeds):
        if s1 not in basin_sub: continue
        for s2 in seeds[i+1:]:
            if s2 not in basin_sub: continue
            try:
                path = nx.shortest_path(basin_sub, s1, s2)
                for p in path:
                    if p not in seeds:
                        bridges.add(p)
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                continue
        if len(bridges) >= max_extra: break
    return set(seeds) | bridges, bridges


def basin_main_cc_nodes(pc, sub, nodes):
    """Return the set of substrate node IDs that would appear in this PC's
    rendered basin layout (largest positive CC + largest negative CC of the
    same-sign induced subgraphs). Used to filter halo highlight sets to
    nodes that are guaranteed to be drawn."""
    pos_nids = [nodes[i] for i in range(len(nodes)) if pc[i] > 0]
    neg_nids = [nodes[i] for i in range(len(nodes)) if pc[i] < 0]
    pos_sub = sub.subgraph(pos_nids)
    neg_sub = sub.subgraph(neg_nids)
    pos_ccs = sorted(nx.connected_components(pos_sub), key=len, reverse=True)
    neg_ccs = sorted(nx.connected_components(neg_sub), key=len, reverse=True)
    out = set()
    if pos_ccs: out |= set(pos_ccs[0])
    if neg_ccs: out |= set(neg_ccs[0])
    return out


def pick_basin_nodes(pc, sub, nodes, nid_idx, mg, n_top_per_basin=10,
                       force_include=None):
    """Pick top-n NON-REACTION seeds per basin and Steiner-expand WITHIN basin.

    force_include: optional set of node IDs that should be added to whichever
    basin their PC-sign places them in, regardless of |loading| rank. Used to
    guarantee that shared-with-partner substrate nodes appear in the basin
    image when this PC is drawn as part of a cross-cohort similarity pair.
    """
    pos_nids = [nodes[i] for i in range(len(nodes)) if pc[i] > 0]
    neg_nids = [nodes[i] for i in range(len(nodes)) if pc[i] < 0]
    pos_sub = sub.subgraph(pos_nids)
    neg_sub = sub.subgraph(neg_nids)

    def cc_l2sq(cc):
        return float((pc[[nid_idx[n] for n in cc]] ** 2).sum())

    pos_ccs = sorted(nx.connected_components(pos_sub), key=cc_l2sq, reverse=True)
    neg_ccs = sorted(nx.connected_components(neg_sub), key=cc_l2sq, reverse=True)
    if not pos_ccs or not neg_ccs:
        return [], [], set(), set(), 0.0, 0.0
    pos_cc = pos_ccs[0]; neg_cc = neg_ccs[0]
    total_mass = float((pc ** 2).sum())
    pos_mass = cc_l2sq(pos_cc) / total_mass
    neg_mass = cc_l2sq(neg_cc) / total_mass

    force_include = set(force_include or [])

    # Top-N non-reaction seeds per basin (+ force-include nodes that live in
    # this basin's main CC so they're guaranteed to appear in the layout).
    # Nodes outside the main CC are NOT force-included — adding them as
    # isolated points puts the halo somewhere visually disconnected from
    # the basin sub-graph, which is worse than dropping them.
    def top_non_reaction(cc, k):
        non_rx = [n for n in cc
                  if mg.graph.nodes.get(n, {}).get("node_type") != "reaction"]
        top = sorted(non_rx, key=lambda n: -abs(pc[nid_idx[n]]))[:k]
        forced_here = [n for n in force_include
                       if n in cc and n not in top
                       and mg.graph.nodes.get(n, {}).get("node_type") != "reaction"]
        return top + forced_here

    pos_seeds = top_non_reaction(pos_cc, n_top_per_basin)
    neg_seeds = top_non_reaction(neg_cc, n_top_per_basin)

    # Steiner expand WITHIN basin (same-sign nodes), then fallback to full sub if needed
    pos_set, pos_bridges = steiner_expand(pos_seeds, pos_sub.subgraph(pos_cc), max_extra=30)
    neg_set, neg_bridges = steiner_expand(neg_seeds, neg_sub.subgraph(neg_cc), max_extra=30)
    # Also walk through full substrate to catch cross-sign reaction bridges
    pos_set_full, pos_bridges_full = steiner_expand(pos_seeds, sub, max_extra=20)
    neg_set_full, neg_bridges_full = steiner_expand(neg_seeds, sub, max_extra=20)
    pos_bridges |= (pos_bridges_full - set(neg_seeds) - set(pos_seeds))
    neg_bridges |= (neg_bridges_full - set(neg_seeds) - set(pos_seeds))
    pos_set = set(pos_seeds) | pos_bridges
    neg_set = set(neg_seeds) | neg_bridges

    return list(pos_set), list(neg_set), pos_bridges, neg_bridges, pos_mass, neg_mass


# ---------------- Reactome ----------------

def load_pathway_names():
    out = {}
    path = REPO / "data/raw/reactome/ReactomePathways.txt"
    if not path.exists(): return out
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == "Homo sapiens":
                out[parts[0]] = parts[1]
    return out


def load_leaf_pathways():
    relations = REPO / "data/raw/reactome/ReactomePathwaysRelation.txt"
    if not relations.exists(): return set(), {}
    all_paths, has_children = set(), set()
    parent_of = {}
    with open(relations) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                parent, child = parts[0], parts[1]
                has_children.add(parent)
                all_paths.add(parent); all_paths.add(child)
                parent_of[child] = parent
    return all_paths - has_children, parent_of


# ---------------- Metadata ----------------

CLINICAL_OUTCOME_FIELDS = {
    # Disease-state / clinical-outcome fields to exclude from the demographic strip
    "active_vs_control", "tumor-vs-normal", "SLE-vs-healthy", "sepsis-vs-control",
    "COVID-vs-healthy", "Crohn-vs-control", "RA-vs-OA", "IDH-mut-vs-WT",
    "BRCA-vs-other", "LUAD-subtype", "case-vs-control",
    "psychosis", "seizure", "cva", "pleurisy", "pericarditis", "lupus_headache",
    "myositis", "organic_brain_syndrome", "mucosal_ulcers", "vasculitis",
    "trop_72h", "creat_0_cat", "ddimer_0_cat", "ldh_0_cat", "crp_0_cat",
    "abs_neut_0_cat", "abs_lymph_0_cat", "abs_mono_0_cat",
    "kidney", "heart", "lung", "immuno", "diabetes", "htn",
    "acuity_0", "acuity_3", "acuity_7", "acuity_28", "acuity_max",
    "resp_symp", "fever_sympt", "gi_symp",
    "diagnosis", "idh_mut", "mgmt_methylated", "stage", "grade",
    "vital.status", "tumor_size_cm", "tumor_site", "stromalscore",
    "transcriptome_subtype", "original_subtype", "polyps_present", "cin",
    "vascular_invasion", "thiopurine", "treatment_received",
    "anticoagulant", "analgesic", "diabetesmed", "dyslipidemia",
    "hypertension", "study_group", "surgery_type",
    "lymphocyte_percent", "neutrophil_percent", "neutrophil_count",
    "thrombocytopenia", "visual_disturbance", "urinary_casts",
    "tumorpurity", "n_dose", "visit", "visit_count", "mortality_event_28days",
}
DEMOGRAPHIC_FIELDS = {
    "age", "age_cat", "sex", "gender", "bmi", "bmi_cat",
    "tobacco_smoking_hist", "ethnicity", "race", "weight", "height",
}


def load_metadata_atlas():
    df = pd.read_csv(RESULTS / "axis_metadata_discrimination_extended.tsv", sep="\t")
    df = df.copy()
    case_labels = {
        "CPTAC_CCRCC": "tumor-vs-normal", "CPTAC_COAD": "tumor-vs-normal",
        "GSE65391_SLE": "SLE-vs-healthy", "Filbin_COVID": "COVID-vs-healthy",
        "Gao_RA": "RA-vs-OA", "IDH_glioma": "IDH-mut-vs-WT",
        "GSE65682_sepsis": "sepsis-vs-control", "KMPLOT_BRCA": "BRCA-vs-other",
        "TCGA_IDH_glioma": "IDH-mut-vs-WT",
    }
    df.loc[df["metadata"] == "active_vs_control", "metadata"] = (
        df.loc[df["metadata"] == "active_vs_control", "cohort"]
          .map(case_labels).fillna("case-vs-control"))
    df["strength"] = np.where(df["metric"] == "auc",
                                2 * (df["value"] - 0.5).abs(),
                                df["value"].abs())
    df["metric_pref"] = np.where(df["metric"] == "auc", 0, 1)
    top = (df.sort_values(["strength", "metric_pref"], ascending=[False, True])
              .groupby(["cohort", "axis"]).head(1))
    return top, df  # return full df too so we can find non-clinical top


def top_non_clinical(df_full, cohort, pc_idx):
    """Return the top demographic/non-clinical metadata field for cohort × PC."""
    pc_key = f"α-PC{pc_idx+1}"
    sub = df_full[(df_full["cohort"] == cohort) & (df_full["axis"] == pc_key)].copy()
    if sub.empty:
        return None
    # Exclude clinical outcomes; prefer demographic fields
    sub_demo = sub[sub["metadata"].isin(DEMOGRAPHIC_FIELDS)]
    if sub_demo.empty:
        sub_demo = sub[~sub["metadata"].isin(CLINICAL_OUTCOME_FIELDS)]
    if sub_demo.empty:
        return None
    sub_demo = sub_demo.sort_values(["strength", "metric_pref"],
                                      ascending=[False, True])
    if sub_demo.iloc[0]["strength"] < 0.01:
        return None
    return sub_demo.iloc[0].to_dict()


def metadata_label(cohort, pc_idx, atlas):
    pc_key = f"α-PC{pc_idx+1}"
    row = atlas[(atlas["cohort"] == cohort) & (atlas["axis"] == pc_key)]
    if row.empty: return None, None, None
    r = row.iloc[0]
    return str(r["metadata"]), float(r["value"]), str(r["metric"])


_METADATA_CACHE = {}


def get_cohort_metadata(cohort):
    """Use existing cohort-specific loaders from axis_metadata_extended.py."""
    if cohort in _METADATA_CACHE:
        return _METADATA_CACHE[cohort]
    md = None
    try:
        from benchmarks.diagnostics.axis_metadata_extended import (
            load_kmplot_metadata, load_cptac_metadata, load_tcga_idh_glioma_metadata,
            load_idh_glioma_trautwein_metadata, load_gao_ra_metadata,
            load_crohn_metadata, load_filbin_metadata, load_erawijantari_metadata,
            load_hmp2_metadata, load_tcga_luad_metadata, load_su_covid_metadata,
            load_gse_series_metadata,
        )
        loader_map = {
            "KMPLOT_BRCA": load_kmplot_metadata,
            "CPTAC_CCRCC": lambda: load_cptac_metadata("CPTAC_CCRCC"),
            "CPTAC_COAD": lambda: load_cptac_metadata("CPTAC_COAD"),
            "CPTAC_OV": lambda: load_cptac_metadata("CPTAC_OV"),
            "TCGA_IDH_glioma": load_tcga_idh_glioma_metadata,
            "IDH_glioma": load_idh_glioma_trautwein_metadata,
            "Gao_RA": load_gao_ra_metadata,
            "Crohn": load_crohn_metadata,
            "Filbin_COVID": load_filbin_metadata,
            "Erawijantari": load_erawijantari_metadata,
            "HMP2_IBD_CD": load_hmp2_metadata,
            "TCGA_LUAD": load_tcga_luad_metadata,
            "Su_COVID": load_su_covid_metadata,
        }
        if cohort in loader_map:
            md = loader_map[cohort]()
        elif cohort.startswith("GSE"):
            # GSE series — try standard cache path
            geo_id = cohort.split("_")[0]
            geo_file = Path.home() / f".cache/geo/{geo_id}/{geo_id}_series_matrix.txt.gz"
            if geo_file.exists():
                md = load_gse_series_metadata(geo_file)
    except Exception as e:
        print(f"  [warn] metadata load for {cohort} failed: {e}", flush=True)
    _METADATA_CACHE[cohort] = md
    return md


def _normalize_pid(pid, cohort):
    """Normalize patient ID for cross-source matching."""
    p = str(pid).lower()
    if cohort.startswith("CPTAC"):
        # CPTAC F-file: c3l-00004_t → c3l-00004
        import re
        p = re.sub(r"_[tn]$", "", p)
    return p


def demographic_labels_for(cohort, patient_ids, field):
    md = get_cohort_metadata(cohort)
    if md is None or field not in md.columns:
        return None
    id_cols = ("patient_id", "sample_id", "id", "patient")
    id_col = None
    for c in id_cols:
        if c in md.columns:
            id_col = c; break
    if id_col is None:
        return None
    md = md.copy()
    md["pid_norm"] = md[id_col].astype(str).str.lower().apply(
        lambda x: _normalize_pid(x, cohort))
    pids = [_normalize_pid(p, cohort) for p in patient_ids]
    vals_series = pd.to_numeric(md[field], errors="coerce")
    vals_dict = dict(zip(md["pid_norm"], vals_series))
    vals = np.array([vals_dict.get(p, np.nan) for p in pids])
    mask = ~np.isnan(vals)
    if mask.sum() < 10:
        return None
    field_lower = field.lower()
    if field_lower in ("sex", "gender"):
        labels = np.where(np.isnan(vals), -1, (vals == 1).astype(int))
        return labels, "male", "female"
    thr = float(np.median(vals[mask]))
    labels = np.where(np.isnan(vals), -1, (vals > thr).astype(int))
    return labels, f"high {field}", f"low {field}"


def patient_labels_for(cohort, patient_ids):
    pids = [str(p) for p in patient_ids]
    if cohort.startswith("CPTAC"):
        labels = np.array([1 if pid.upper().endswith("_T") else
                            0 if pid.upper().endswith("_N") else -1
                            for pid in pids])
        if (labels >= 0).sum() < 5: return None
        return labels, "tumor", "adjacent-normal"
    if cohort == "Filbin_COVID":
        try:
            path = Path.home() / ".cache/filbin_mgh_covid/Clinical_Metadata.xlsx"
            if not path.exists(): return None
            md = pd.read_excel(path)
            md["pid_lower"] = md["Public ID"].astype(int).astype(str).str.lower()
            trop = pd.to_numeric(md["Trop_72h"], errors="coerce")
            pid_to_trop = dict(zip(md["pid_lower"], trop))
            vals = np.array([pid_to_trop.get(str(p).lower(), np.nan) for p in patient_ids])
            mask = ~np.isnan(vals)
            if mask.sum() < 10: return None
            thr = float(np.median(vals[mask]))
            labels = np.where(np.isnan(vals), -1, np.where(vals > thr, 1, 0))
            return labels, "high Trop_72h", "low Trop_72h"
        except Exception:
            return None
    return None


# ---------------- Pathway overlay ----------------

def find_pathway_groups(all_module_nodes, mg, sub, pathway_names, leaf_pathways,
                         pathway_parent, min_members=2):
    """Find leaf pathways with ≥ min_members nodes in module; fall back to parent
    pathways if no leaf qualifies."""
    node_pathways = {}
    for n in all_module_nodes:
        n_attrs = mg.graph.nodes.get(n, {})
        pws = set()
        if n_attrs.get("node_type") == "reaction":
            for p in (n_attrs.get("pathways") or []):
                if p in leaf_pathways and p in pathway_names:
                    pws.add((p, pathway_names[p]))
        else:
            for nb in sub.neighbors(n):
                nb_attrs = mg.graph.nodes.get(nb, {})
                if nb_attrs.get("node_type") == "reaction":
                    for p in (nb_attrs.get("pathways") or []):
                        if p in leaf_pathways and p in pathway_names:
                            pws.add((p, pathway_names[p]))
        node_pathways[n] = pws

    pw_to_nodes = defaultdict(set)
    for n, pws in node_pathways.items():
        for pid, pname in pws:
            pw_to_nodes[(pid, pname)].add(n)
    groups = [(k, v) for k, v in pw_to_nodes.items() if len(v) >= min_members]
    groups.sort(key=lambda x: -len(x[1]))
    if groups:
        return groups[:4]
    # Parent fallback
    parent_to_nodes = defaultdict(set)
    for n in all_module_nodes:
        n_attrs = mg.graph.nodes.get(n, {})
        if n_attrs.get("node_type") != "reaction":
            for nb in sub.neighbors(n):
                nb_attrs = mg.graph.nodes.get(nb, {})
                if nb_attrs.get("node_type") == "reaction":
                    for p in (nb_attrs.get("pathways") or []):
                        parent = pathway_parent.get(p)
                        if parent and parent in pathway_names:
                            parent_to_nodes[(parent, pathway_names[parent])].add(n)
    groups = [(k, v) for k, v in parent_to_nodes.items() if len(v) >= 3]
    groups.sort(key=lambda x: -len(x[1]))
    return groups[:4]


# ---------------- Drawing ----------------

def draw_basin_single(ax, pos_nodes, neg_nodes, pos_bridges, neg_bridges,
                      sub, nid_idx, pc, mg, pos_mass, neg_mass,
                      pathway_names, leaf_pathways, pathway_parent,
                      highlight_nodes=None):
    keep = list(set(pos_nodes) | set(neg_nodes))
    G = sub.subgraph(keep).copy()
    pos_set, neg_set = set(pos_nodes), set(neg_nodes)
    seed_pos = pos_set - pos_bridges
    seed_neg = neg_set - neg_bridges

    # Layout each basin INDEPENDENTLY on its own subgraph, then place side-by-side.
    # This lets nodes within each basin spread freely instead of being pulled to
    # the centerline by cross-basin forces.
    G_pos = G.subgraph(set(pos_nodes)).copy()
    G_neg = G.subgraph(set(neg_nodes)).copy()
    rng = np.random.default_rng(42)
    # Wider spring k + larger scale → nodes spread further horizontally
    # to match the wider figsize (24×14 instead of 18×16).
    if G_pos.number_of_edges() > 0:
        pos_layout = nx.spring_layout(G_pos, k=3.0, iterations=400, seed=42, scale=3.5)
    else:
        pos_layout = {n: np.array([rng.uniform(-2.0, 2.0), rng.uniform(-2.5, 2.5)])
                      for n in G_pos.nodes()}
    if G_neg.number_of_edges() > 0:
        neg_layout = nx.spring_layout(G_neg, k=3.0, iterations=400, seed=42, scale=3.5)
    else:
        neg_layout = {n: np.array([rng.uniform(-2.0, 2.0), rng.uniform(-2.5, 2.5)])
                      for n in G_neg.nodes()}
    # Shift basins apart horizontally (matched to larger scale above)
    layout = {}
    for n, p in neg_layout.items():
        layout[n] = np.array([p[0] - 5.5, p[1]])
    for n, p in pos_layout.items():
        layout[n] = np.array([p[0] + 5.5, p[1]])
    # If any node was in both (cross-sign bridge), pos wins
    for n in G_pos.nodes():
        if n in pos_layout:
            layout[n] = np.array([pos_layout[n][0] + 5.5, pos_layout[n][1]])

    # ---- Pathway overlay (combined seeds + bridges for both basins) ----
    pos_groups = find_pathway_groups(list(pos_set), mg, sub, pathway_names,
                                       leaf_pathways, pathway_parent)
    neg_groups = find_pathway_groups(list(neg_set), mg, sub, pathway_names,
                                       leaf_pathways, pathway_parent)
    pw_cmap_pos = plt.cm.Reds(np.linspace(0.25, 0.55, max(len(pos_groups), 3)))
    pw_cmap_neg = plt.cm.Blues(np.linspace(0.25, 0.55, max(len(neg_groups), 3)))

    for groups, cmap, basin_color in [(pos_groups, pw_cmap_pos, "pos"),
                                        (neg_groups, pw_cmap_neg, "neg")]:
        for pi, ((pid, pname), gnodes) in enumerate(groups):
            pts = np.array([layout[n] for n in gnodes if n in layout])
            if len(pts) < 2: continue
            color = cmap[pi % len(cmap)]
            edge_color = (color[0]*0.6, color[1]*0.6, color[2]*0.6, 1.0)
            centroid = pts.mean(axis=0)
            # Cytoscape-style rounded rectangle container (axis-aligned bounding box)
            pad = 0.45
            x_min, x_max = pts[:, 0].min() - pad, pts[:, 0].max() + pad
            y_min, y_max = pts[:, 1].min() - pad, pts[:, 1].max() + pad
            w, h = x_max - x_min, y_max - y_min
            rect = mpatches.FancyBboxPatch(
                (x_min, y_min), w, h,
                boxstyle="round,pad=0.1,rounding_size=0.4",
                facecolor=color, alpha=0.30,
                edgecolor=edge_color, linewidth=2.0, zorder=0)
            ax.add_patch(rect)
            # Pathway label — wrap long names; position ABOVE the rectangle with
            # a generous clearance + a white background so it survives any overlap.
            pname_wrapped = "\n".join(textwrap.wrap(pname, width=28))
            ax.text(centroid[0], y_max + 0.85, pname_wrapped, fontsize=8.5,
                     color=edge_color, ha="center", va="bottom",
                     fontweight="bold", style="italic", zorder=15,
                     linespacing=1.1,
                     bbox=dict(facecolor="white", alpha=0.95,
                               edgecolor=edge_color, linewidth=0.6, pad=2))

    # ---- Edges — prominent dark grey ----
    nx.draw_networkx_edges(G, layout, edge_color="#444", alpha=0.7, width=2.5, ax=ax)

    # ---- Highlight halos (drawn ABOVE seed/bridge nodes; zorder 6) ----
    # Gold rings around any node in highlight_nodes that landed in the layout.
    # Used to mark substrate nodes shared with a paired basin in the cross-
    # cohort similarity figure — visually links the two basins. Drawn in
    # front of nodes so they're never obscured.
    if highlight_nodes:
        highlight_color = "#F2B500"
        for n in highlight_nodes:
            if n not in layout: continue
            x, y = layout[n]
            # Bright inner ring (thick, just outside the node circle)
            halo = mpatches.Circle((x, y), 0.95, fill=False,
                                     edgecolor=highlight_color,
                                     linewidth=6.5, alpha=1.0, zorder=6)
            ax.add_patch(halo)
            # Soft outer ring for additional pop
            outer = mpatches.Circle((x, y), 1.30, fill=False,
                                      edgecolor=highlight_color,
                                      linewidth=2.5, alpha=0.6, zorder=6)
            ax.add_patch(outer)

    # ---- Bridge nodes ----
    bridge_color_pos = to_rgba("#D9351F", alpha=0.4)
    bridge_color_neg = to_rgba("#2D5BBF", alpha=0.4)
    abs_loads = np.array([abs(pc[nid_idx[n]]) for n in G.nodes()])
    max_load = abs_loads.max() + 1e-9

    shape_for_type = {"gene": "o", "reaction": "s", "metabolite": "^"}
    for tname, shape in shape_for_type.items():
        # Bridges — zorder 3 (above pathway box + labels)
        bnodes = [n for n in pos_bridges if mg.graph.nodes.get(n, {}).get("node_type") == tname]
        bnodes_neg = [n for n in neg_bridges if mg.graph.nodes.get(n, {}).get("node_type") == tname]
        for nl, color in [(bnodes, bridge_color_pos), (bnodes_neg, bridge_color_neg)]:
            if nl:
                sizes = [120 + 400 * (abs(pc[nid_idx[n]]) / max_load) ** 0.7 for n in nl]
                colors = [color] * len(nl)
                coll = nx.draw_networkx_nodes(G, layout, nodelist=nl, node_size=sizes,
                                                node_color=colors, node_shape=shape,
                                                edgecolors="none", ax=ax)
                coll.set_zorder(3)

    # ---- Seed nodes — zorder 5 (over bridges, edges, pathway boxes & labels) ----
    seed_color_pos = to_rgba("#D9351F", alpha=0.95)
    seed_color_neg = to_rgba("#2D5BBF", alpha=0.95)
    for tname, shape in shape_for_type.items():
        spos = [n for n in seed_pos if mg.graph.nodes.get(n, {}).get("node_type") == tname]
        sneg = [n for n in seed_neg if mg.graph.nodes.get(n, {}).get("node_type") == tname]
        for nl, color in [(spos, seed_color_pos), (sneg, seed_color_neg)]:
            if nl:
                sizes = [400 + 2200 * (abs(pc[nid_idx[n]]) / max_load) ** 0.7 for n in nl]
                colors = [color] * len(nl)
                coll = nx.draw_networkx_nodes(G, layout, nodelist=nl, node_size=sizes,
                                                node_color=colors, node_shape=shape,
                                                edgecolors="white", linewidths=2.0, ax=ax)
                coll.set_zorder(5)

    # ---- Labels on seeds — BELOW the node, wrapped to multiple lines so full
    # names show without truncation. No bbox = doesn't obscure pathway labels.
    for n in (seed_pos | seed_neg):
        if n not in layout: continue
        x, y = layout[n]
        name = short_name(n, mg)[0]
        # Wrap long names (multi-line) instead of truncating
        wrapped = "\n".join(textwrap.wrap(name, width=16))
        ax.text(x, y - 0.30, wrapped, ha="center", va="top",
                 fontsize=8, color="#0a0a0a", fontweight="bold",
                 zorder=20, linespacing=1.0)

    ax.set_axis_off()


def draw_patient_strip(ax, pc_scores, patient_labels, pos_name, neg_name,
                        pos_color="#D9351F", neg_color="#2D5BBF"):
    """Strip plot with overlaid boxplots per metadata group."""
    ax.set_facecolor("#FAFBFC")
    if pc_scores is None or len(pc_scores) < 5:
        ax.text(0.5, 0.5, "(no patient scores)", ha="center", va="center",
                 fontsize=9, color="#888", transform=ax.transAxes)
        ax.set_axis_off(); return
    n = len(pc_scores)
    rng = np.random.default_rng(7)

    def _box_at(scores, y_center, color, w=0.5):
        ax.boxplot([scores], positions=[y_center], widths=[w],
                    vert=False, patch_artist=True, showfliers=False,
                    whiskerprops=dict(color=color, linewidth=1.2),
                    capprops=dict(color=color, linewidth=1.2),
                    medianprops=dict(color="black", linewidth=1.5),
                    boxprops=dict(facecolor=color, alpha=0.45,
                                  edgecolor=color, linewidth=1.5))

    if patient_labels is not None:
        pos_mask = patient_labels == 1
        neg_mask = patient_labels == 0
        other_mask = patient_labels == -1
        if pos_mask.any():
            y = 0.8 + 0.3 * rng.standard_normal(int(pos_mask.sum()))
            ax.scatter(pc_scores[pos_mask], y, c=pos_color, s=22, alpha=0.75,
                        edgecolors="white", linewidths=0.4, zorder=3,
                        label=f"{pos_name} (n={int(pos_mask.sum())})")
            _box_at(pc_scores[pos_mask], 2.0, pos_color)
        if neg_mask.any():
            y = -0.8 + 0.3 * rng.standard_normal(int(neg_mask.sum()))
            ax.scatter(pc_scores[neg_mask], y, c=neg_color, s=22, alpha=0.75,
                        edgecolors="white", linewidths=0.4, zorder=3,
                        label=f"{neg_name} (n={int(neg_mask.sum())})")
            _box_at(pc_scores[neg_mask], -2.0, neg_color)
        if other_mask.any():
            y = 0.3 * rng.standard_normal(int(other_mask.sum()))
            ax.scatter(pc_scores[other_mask], y, c="#cccccc", s=12, alpha=0.4,
                        edgecolors="none", zorder=2,
                        label=f"unlabeled (n={int(other_mask.sum())})")
        ax.legend(loc="upper right", fontsize=9, frameon=False, markerscale=1.2)
        ax.set_ylim(-3, 3); ax.set_yticks([])
    else:
        y = rng.uniform(-0.6, 0.6, n)
        ax.scatter(pc_scores, y, c="#666", s=18, alpha=0.6,
                    edgecolors="white", linewidths=0.3, zorder=3,
                    label=f"all patients (n={n})")
        _box_at(pc_scores, 1.4, "#666", w=0.6)
        ax.legend(loc="upper right", fontsize=9, frameon=False)
        ax.set_ylim(-1.5, 2.2); ax.set_yticks([])

    ax.set_xlabel("α-PC score  (low ← → high)", fontsize=11)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.axvline(0, color="#999", lw=0.7, linestyle="--", alpha=0.5)


_PALETTE = {"bg": "#F8F9FB", "accent": "#6E4FB0"}


def render_one_basin(state, cohort, pc_idx, label, *,
                       highlight_nodes=None, output_suffix=""):
    """Render a single basin figure to ``fig_basin_v6_<cohort>_pc<k>{suffix}.png``.

    state: dict produced by ``load_render_state()`` (substrate + metadata).
    highlight_nodes: optional set of substrate node IDs to gold-halo.
    output_suffix: e.g. ``"_pair3_highlighted"`` to avoid clobbering the
        canonical basin PNG; default empty preserves the canonical name.
    """
    mg, sub, nodes, nid_idx, log_pr = (
        state["mg"], state["sub"], state["nodes"],
        state["nid_idx"], state["log_pr"],
    )
    atlas, atlas_full = state["atlas"], state["atlas_full"]
    pathway_names = state["pathway_names"]
    leaf_pathways = state["leaf_pathways"]
    pathway_parent = state["pathway_parent"]
    PALETTE = _PALETTE

    F_path = find_F_path(cohort)
    if F_path is None:
        print(f"  skip {cohort}: no F", flush=True)
        return None
    print(f"\n=== {cohort} α-PC{pc_idx+1}"
          + (f"  (suffix={output_suffix}, "
             f"highlight={len(highlight_nodes or set())} nodes)"
             if highlight_nodes else "") + " ===", flush=True)
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"].astype(np.float64)
    if F.shape[1] != len(nodes):
        print(f"  skip {cohort}: F shape mismatch", flush=True)
        return None
    _, alpha = decompose_unit_norm(F, log_pr)
    pca = PCA(n_components=pc_idx + 1, random_state=0)
    scores = pca.fit_transform(alpha)
    pc = pca.components_[pc_idx]
    pc_scores = scores[:, pc_idx]

    md_field, md_val, md_metric = metadata_label(cohort, pc_idx, atlas)
    if md_val is not None:
        if md_metric == "auc":
            md_str = f"{md_field}  (AUC = {max(md_val, 1-md_val):.3f})"
        else:
            md_str = f"{md_field}  (|ρ| = {abs(md_val):.3f})"
    else:
        md_str = ""

    patient_labels, pos_name, neg_name = None, None, None
    if "patient_ids" in fd.files:
        lab = patient_labels_for(cohort, fd["patient_ids"])
        if lab is not None:
            patient_labels, pos_name, neg_name = lab

    # Force-include highlight nodes in basin pick so the gold halos land on
    # nodes that actually exist in the rendered layout.
    pos_nodes, neg_nodes, pos_bridges, neg_bridges, pos_mass, neg_mass = pick_basin_nodes(
        pc, sub, nodes, nid_idx, mg, n_top_per_basin=10,
        force_include=highlight_nodes,
    )
    print(f"  + basin: {len(pos_nodes)} nodes ({len(pos_bridges)} bridges); "
          f"− basin: {len(neg_nodes)} nodes ({len(neg_bridges)} bridges)", flush=True)
    if highlight_nodes:
        in_layout = sum(1 for n in highlight_nodes
                        if n in set(pos_nodes) | set(neg_nodes))
        print(f"  highlight coverage: {in_layout}/{len(highlight_nodes)} shared nodes in basin layout",
              flush=True)

    fig = plt.figure(figsize=(24, 14))
    fig.patch.set_facecolor(PALETTE["bg"])
    gs = fig.add_gridspec(4, 1, height_ratios=[0.06, 1.8, 0.12, 0.12], hspace=0.05)
    ax_header = fig.add_subplot(gs[0])
    ax_net = fig.add_subplot(gs[1])
    ax_strip = fig.add_subplot(gs[2])
    ax_strip_demo = fig.add_subplot(gs[3])
    for a in (ax_header, ax_net, ax_strip, ax_strip_demo):
        a.set_facecolor(PALETTE["bg"])

    ax_header.set_axis_off()
    ax_header.text(0.0, 0.75,
                    f"{cohort.replace('_', ' ')}    ·    α-PC{pc_idx+1}",
                    fontsize=24, fontweight="bold", color="#101010",
                    transform=ax_header.transAxes)
    ax_header.text(0.0, 0.20, label, fontsize=14, color="#444",
                    transform=ax_header.transAxes)
    if md_str:
        ax_header.text(1.0, 0.75, "Top metadata explanation",
                        fontsize=11, color="#888", ha="right",
                        transform=ax_header.transAxes)
        ax_header.text(1.0, 0.20, md_str, fontsize=18, fontweight="bold",
                        color=PALETTE["accent"], ha="right",
                        transform=ax_header.transAxes)

    draw_basin_single(ax_net, pos_nodes, neg_nodes, pos_bridges, neg_bridges,
                       sub, nid_idx, pc, mg, pos_mass, neg_mass,
                       pathway_names, leaf_pathways, pathway_parent,
                       highlight_nodes=highlight_nodes)
    ax_net.relim(); ax_net.autoscale_view()
    cur_xlim = ax_net.get_xlim(); cur_ylim = ax_net.get_ylim()
    pad = 0.5
    ax_net.set_xlim(cur_xlim[0] - pad, cur_xlim[1] + pad)
    ax_net.set_ylim(cur_ylim[0] - pad, cur_ylim[1] + pad)
    ax_net.text(0.02, 0.97, "−  basin (low-score patient end)",
                 fontsize=12, color="#2D5BBF", fontweight="bold",
                 transform=ax_net.transAxes, va="top")
    ax_net.text(0.02, 0.93, f"{neg_mass*100:.1f}% of |PC|² mass",
                 fontsize=9, color="#666", style="italic",
                 transform=ax_net.transAxes, va="top")
    ax_net.text(0.98, 0.97, "+  basin (high-score patient end)",
                 fontsize=12, color="#D9351F", fontweight="bold",
                 transform=ax_net.transAxes, va="top", ha="right")
    ax_net.text(0.98, 0.93, f"{pos_mass*100:.1f}% of |PC|² mass",
                 fontsize=9, color="#666", style="italic",
                 transform=ax_net.transAxes, va="top", ha="right")

    draw_patient_strip(ax_strip, pc_scores, patient_labels, pos_name, neg_name)

    demo_row = top_non_clinical(atlas_full, cohort, pc_idx)
    if demo_row is not None and "patient_ids" in fd.files:
        demo_field = demo_row["metadata"]
        demo_val = demo_row["value"]
        demo_metric = demo_row["metric"]
        demo_lab = demographic_labels_for(cohort, fd["patient_ids"], demo_field)
        if demo_lab is not None:
            d_labels, d_pos, d_neg = demo_lab
            draw_patient_strip(ax_strip_demo, pc_scores, d_labels, d_pos, d_neg,
                                pos_color="#7B3F8C", neg_color="#3F8C7B")
            if demo_metric == "auc":
                demo_str = f"{demo_field} (AUC = {max(demo_val, 1-demo_val):.3f})"
            else:
                demo_str = f"{demo_field} (|ρ| = {abs(demo_val):.3f})"
            ax_strip_demo.set_xlabel(f"top non-clinical metadata: {demo_str}",
                                      fontsize=10, color="#555", style="italic")
        else:
            ax_strip_demo.text(0.5, 0.5,
                                f"(no demographic data for {demo_field})",
                                ha="center", va="center", fontsize=9, color="#aaa",
                                transform=ax_strip_demo.transAxes)
            ax_strip_demo.set_axis_off()
    else:
        ax_strip_demo.text(0.5, 0.5, "(no non-clinical metadata available)",
                            ha="center", va="center", fontsize=9, color="#aaa",
                            transform=ax_strip_demo.transAxes)
        ax_strip_demo.set_axis_off()

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
                markeredgecolor="white", markeredgewidth=1.5, markersize=10,
                label="gene"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#888",
                markeredgecolor="white", markeredgewidth=1.5, markersize=10,
                label="reaction"),
        Line2D([0], [0], marker="h", color="w", markerfacecolor="#888",
                markeredgecolor="white", markeredgewidth=1.5, markersize=11,
                label="metabolite"),
    ]
    if highlight_nodes:
        handles.append(Line2D([0], [0], marker="o", color="w",
                               markerfacecolor="none", markeredgecolor="#F2B500",
                               markeredgewidth=2.5, markersize=12,
                               label="shared with paired basin"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
                frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.0))

    out_path = FIG / f"fig_basin_v6_{cohort}_pc{pc_idx+1}{output_suffix}.png"
    plt.savefig(out_path, dpi=140, bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close()
    print(f"  wrote {out_path}", flush=True)
    return out_path


def load_render_state():
    """Load substrate + metadata once; reused across many render_one_basin calls."""
    mg, sub, nodes, nid_idx, log_pr = load_substrate()
    atlas, atlas_full = load_metadata_atlas()
    pathway_names = load_pathway_names()
    leaf_pathways, pathway_parent = load_leaf_pathways()
    print(f"  pathway names: {len(pathway_names)}; leaf pathways: {len(leaf_pathways)}",
          flush=True)
    return {
        "mg": mg, "sub": sub, "nodes": nodes, "nid_idx": nid_idx,
        "log_pr": log_pr, "atlas": atlas, "atlas_full": atlas_full,
        "pathway_names": pathway_names, "leaf_pathways": leaf_pathways,
        "pathway_parent": pathway_parent,
    }


def main():
    state = load_render_state()

    cases = [
        # Main figure picks (Section 2 canonical examples)
        ("Filbin_COVID", 3, "cardiorenal — LPO antimicrobial ↔ galanin receptors"),
        ("CPTAC_CCRCC", 0, "tumor-vs-normal — renal tubule ↔ Warburg/HIF"),
        ("IDH_glioma", 0, "IDH-mut — neuronal ↔ immune"),
        ("CPTAC_COAD", 3, "tumor — intestinal endocrine ↔ pancreatic exocrine"),
        ("GSE65391_SLE", 2, "psychosis — type I IFN signature"),
        ("KMPLOT_BRCA", 3, "grade — tumor immune infiltrate ↔ secretory"),
        ("GSE65682_sepsis", 1, "sepsis — metabolic ↔ lymphoid exhaustion"),
        ("Gao_RA", 0, "RA diagnosis — synovial inflammation"),
        ("TCGA_IDH_glioma", 4, "IDH-mut — epigenetic/translation axis"),
        # Supplementary S10 — full-coverage basins for the remaining panel cohorts
        ("Su_COVID",     0, "COVID-vs-healthy — α-PC1 (substrate-friendly inflammatory)"),
        ("Erawijantari", 3, "anticoagulant — α-PC4 (gastric microbiome modulation)"),
        ("HMP2_IBD_CD",  0, "IBD-CD — α-PC1 (microbiome-driven)"),
        ("GSE89408_RA",  0, "RA-vs-OA — α-PC1 (synovial inflammation, RNA-only)"),
        ("CPTAC_OV",     0, "tumor-vs-normal — α-PC1 (at-ceiling for any method)"),
        ("Crohn",        0, "thiopurine — α-PC1 (small-n scope-limit demonstration)"),
        ("TCGA_LUAD",    2, "LUAD-subtype — α-PC3 (GoF scope-limit; substrate-routed signal weak)"),
        # Additional basins for the top-10 cross-cohort similarity panel
        # (paired against the corresponding entries already in the cases list)
        ("IDH_glioma",       1, "α-PC2 (cross-cohort match TCGA_IDH α-PC2, |cos|=0.77)"),
        ("IDH_glioma",       2, "α-PC3 (cross-cohort match TCGA_IDH α-PC4, |cos|=0.66)"),
        ("IDH_glioma",       3, "α-PC4 (cross-cohort match TCGA_IDH α-PC3, |cos|=0.47)"),
        ("TCGA_IDH_glioma",  0, "α-PC1 (IDH-mut axis paired with Trautwein α-PC1)"),
        ("TCGA_IDH_glioma",  1, "α-PC2 (IDH-mut axis paired with Trautwein α-PC2)"),
        ("TCGA_IDH_glioma",  2, "α-PC3 (matches Trautwein α-PC4 / GoF-comparator)"),
        ("TCGA_IDH_glioma",  3, "α-PC4 (matches Trautwein α-PC3 — cross-PC rotation)"),
        ("KMPLOT_BRCA",      2, "α-PC3 (matches GSE89408_RA α-PC1, cross-disease |cos|=0.50)"),
        ("CPTAC_COAD",       0, "α-PC1 (matches CPTAC_OV α-PC1, |cos|=0.69 — shared cancer-intensity hub axis)"),
        ("TCGA_LUAD",        0, "α-PC1 (matches IDH_glioma α-PC2 — cross-disease cancer axis)"),
        ("GSE65391_SLE",     0, "α-PC1 (matches GSE65682_sepsis α-PC2, |cos|=0.51 — autoimmune/inflammation)"),
        ("GSE65682_sepsis",  4, "α-PC5 (matches GSE65391_SLE α-PC3, |cos|=0.52 — shared autoimmune axis)"),
    ]

    for cohort, pc_idx, label in cases:
        render_one_basin(state, cohort, pc_idx, label)


if __name__ == "__main__":
    main()
