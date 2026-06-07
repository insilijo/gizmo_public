"""Per-chord shared sub-graph view: the substrate co-structure that
*undergirds* each top-N cross-cohort α-PC similarity.

For each of the top-N chords from the alignment diagram, render the
substrate sub-graph induced by the **top-K contributors to the cosine**
(nodes ranked by |PC_A[i] · PC_B[i]| descending). The cosine of two
unit-norm PCs is literally the sum of these products, so these are the
nodes that *make* |cos| high.

Each node carries two pieces of agreement information at once:

  • **Fill color** = sign in PC_A (red = +basin, blue = −basin in cohort A)
  • **Edge color** = sign in PC_B (red = +basin, blue = −basin in cohort B)

  → Same-color fill+edge = sign agrees (constructive contribution to cos).
    Different fill+edge   = sign flip   (destructive contribution).

Node size scales with |PC_A · PC_B| — the bigger the dot, the more it
contributes to the chord's cosine. Edges drawn between any two selected
nodes that share a substrate edge.

Inputs:
  pc_alignment_best_cross_cohort_match.tsv

Output:
  fig_chord_costructure.png       (5 × 2 grid of per-pair sub-graph panels)
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
from sklearn.decomposition import PCA

REPO = Path("/home/jgardner/GIZMO")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))
sys.path.insert(0, str(REPO / "benchmarks/figures"))

# Reuse Reactome pathway tooling from the basin builder
import build_basin_signed_v2 as basin_mod  # noqa: E402

RESULTS = REPO / "benchmarks/results"
UR = RESULTS / "unsupervised"
FIG = RESULTS / "figures"
TOP_N = 10
# Per-pair seed selection uses a coverage target instead of a fixed top-K
# (the chord cosine is heavy-tailed; an arbitrary K either over- or
# under-shoots the dominant biology). Greedy-connected growth from the
# strongest contributor within each sign class, stopping when cumulative
# |contribution| reaches COVERAGE_TARGET of that sign's total.
COVERAGE_TARGET = 0.30    # 30% of the sign-class contribution to cos.
MARGINAL_DROPOFF = 0.20   # elbow stop: stop when next candidate's
                           # |contribution| drops below 20% of the strongest
                           # candidate's |contribution|. 0.20 catches the
                           # actual elbow without truncating concentrated
                           # pairs to a single seed.
                           # Stop at WHICHEVER comes first (coverage or elbow).
MAX_BRIDGES = 30          # max Steiner bridge nodes added to connect seeds
POS_COLOR = "#D9351F"     # red — +basin / positive PC loading
NEG_COLOR = "#2D5BBF"     # blue — −basin / negative PC loading
BRIDGE_COLOR = "#bbbbbb"  # gray — Steiner bridge, not a contributor
EDGE_GRAY = "#666"


def load_substrate():
    from gizmo.export.json_export import read_json
    from per_patient_wlsp_v2 import biochem_subgraph
    print("Loading substrate…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    sub_dir, nodes, _ = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)
    # Also return the FULL substrate (no hub cap) for Steiner connection.
    # The hub-capped `sub` is appropriate for MAP / basin work, but for
    # connecting top-cosine-contributor seeds across the substrate we need
    # the unrestricted graph — otherwise high-degree intermediate hubs
    # (transcription factors, signaling kinases, etc.) get masked out and
    # otherwise-connected seeds appear disconnected.
    full_g = mg.graph.to_undirected() if mg.graph.is_directed() else mg.graph
    return mg, sub, full_g, nodes, log_pr


def alpha_pc_components(cohort, log_pr, n_components=5):
    F_path = None
    for suffix in ("", "_edge_informed", "_combined", "_node_informed"):
        cand = UR / f"stage3_F_{cohort}{suffix}.npz"
        if cand.exists():
            F_path = cand; break
    if F_path is None: return None
    fd = np.load(F_path, allow_pickle=True)
    F = fd["F"].astype(np.float64)
    if F.shape[1] != len(log_pr): return None
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    n = min(n_components, alpha.shape[0])
    pca = PCA(n_components=n, random_state=0)
    pca.fit(alpha)
    # Return UNIT-NORM components (already returned that way by sklearn)
    return pca.components_


def greedy_connected_seeds(abs_contrib, sign_mask, nodes, sub,
                              coverage_target=COVERAGE_TARGET,
                              marginal_dropoff=MARGINAL_DROPOFF,
                              max_hops=2):
    """Greedy-connected seed selection within one sign-class.

    Stops at the FIRST of:
      (a) cumulative |contribution| ≥ coverage_target × sign_class_total
      (b) next candidate's |contribution| < marginal_dropoff × max contribution
          (elbow — beyond this point new seeds add negligible biology)
      (c) no reachable candidate remains

    Pure non-K stopping; (a) bounds coverage, (b) bounds tail length.
    """
    cand_idx = np.where(sign_mask)[0]
    if len(cand_idx) == 0:
        return [], 0.0
    total_class = float(abs_contrib[cand_idx].sum())
    if total_class <= 0:
        return [], 0.0
    sorted_idx = cand_idx[np.argsort(-abs_contrib[cand_idx])]
    ordered = [(nodes[i], abs_contrib[i]) for i in sorted_idx if nodes[i] in sub]
    if not ordered:
        return [], 0.0
    max_contrib_class = float(ordered[0][1])
    cutoff_marginal = marginal_dropoff * max_contrib_class
    selected = [ordered[0][0]]
    cumulative = float(ordered[0][1])
    reachable = set(
        nx.single_source_shortest_path_length(sub, ordered[0][0],
                                                cutoff=max_hops).keys()
    )
    for cand, contrib_val in ordered[1:]:
        if cumulative / total_class >= coverage_target:
            break
        if float(contrib_val) < cutoff_marginal:
            break
        if cand in selected: continue
        if cand in reachable:
            selected.append(cand)
            cumulative += float(contrib_val)
            try:
                reachable |= set(
                    nx.single_source_shortest_path_length(sub, cand,
                                                            cutoff=max_hops).keys()
                )
            except (nx.NetworkXError, nx.NodeNotFound):
                continue
    return selected, cumulative / total_class


def steiner_bridge(seeds, sub, max_extra=MAX_BRIDGES):
    """Connect seed nodes into a single subgraph using an approximate Steiner tree.

    Reduces the seed terminals to those that lie in the substrate's largest
    connected component (cross-component seeds cannot be connected at all),
    then computes an approximate Steiner tree spanning them. Caps the bridge
    count by truncating the seed list per CC if needed.
    Returns (all_nodes, bridges, seed_set, dropped_seeds).
    """
    from networkx.algorithms.approximation import steiner_tree
    seeds = [s for s in seeds if s in sub]
    if len(seeds) < 2:
        return set(seeds), set(), set(seeds), []
    # Keep only seeds in the substrate's largest CC (Steiner tree requires
    # connectivity); report the rest as dropped.
    ccs = list(nx.connected_components(sub))
    biggest = max(ccs, key=len) if ccs else set()
    in_main = [s for s in seeds if s in biggest]
    dropped = [s for s in seeds if s not in biggest]
    if len(in_main) < 2:
        return set(in_main), set(), set(in_main), dropped
    try:
        tree = steiner_tree(sub, in_main, method="kou")
    except Exception:
        # Fallback: union of pairwise shortest paths
        tree_nodes = set(in_main)
        for i, s1 in enumerate(in_main):
            for s2 in in_main[i+1:]:
                try:
                    p = nx.shortest_path(sub, s1, s2)
                    tree_nodes.update(p)
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    continue
        all_nodes = tree_nodes
        bridges = all_nodes - set(in_main)
        return all_nodes, bridges, set(in_main), dropped
    all_nodes = set(tree.nodes())
    bridges = all_nodes - set(in_main)
    # If the Steiner tree dragged in too many bridges, trim by removing the
    # lowest-degree bridge nodes until we hit max_extra.
    if len(bridges) > max_extra:
        bridge_deg = sorted(bridges, key=lambda n: tree.degree(n))
        for b in bridge_deg:
            if len(bridges) <= max_extra: break
            # Only safely removable if degree-1 (leaf) — leaf bridges can be
            # pruned without disconnecting the seed terminals.
            if tree.degree(b) == 1:
                bridges.discard(b)
                all_nodes.discard(b)
    return all_nodes, bridges, set(in_main), dropped


def label_for(nid, mg):
    attrs = mg.graph.nodes.get(nid, {})
    sym = attrs.get("symbol") or attrs.get("name") or attrs.get("display_name") or nid
    s = str(sym)
    if len(s) > 18: s = s[:17] + "…"
    return s


def node_type_first(nid, mg):
    return (mg.graph.nodes.get(nid, {}).get("node_type", "?") or "?")[0]


def render_back_to_back_dotplot(ax, mg, nodes, vA, vB, seed_set,
                                   pathway_groups, cA_label, cB_label):
    """Back-to-back dot plot of signed loadings, grouped vertically by
    Reactome pathway.

    Y axis = seed nodes, ordered by pathway membership; pathway names are
    drawn as italic header rows between groups.
    Left dot at x = −\|vA[i]\|  — its color encodes sign of vA[i].
    Right dot at x = +\|vB[i]\| — its color encodes sign of vB[i].
    Same color on both sides ⇒ constructive contribution; different ⇒ sign-flip.
    """
    idx_of = {n: i for i, n in enumerate(nodes)}

    seed_to_first_pw = {}
    for (_pid, pname), gnodes in (pathway_groups or []):
        for n in gnodes:
            seed_to_first_pw.setdefault(n, pname)
    grouped = defaultdict(list)
    unassigned = []
    for n in seed_set:
        if n in seed_to_first_pw:
            grouped[seed_to_first_pw[n]].append(n)
        else:
            unassigned.append(n)

    rows = []   # list of (label, kind, node_id|None) — pathway grouping is
                # implicit (sort order) only; explicit pathway header rows
                # have been removed because they crowded the gene labels.
                # Pathway info now lives in the hypergeometric column.
    for pname in sorted(grouped.keys()):
        nodes_in_pw = sorted(grouped[pname],
                               key=lambda n: -(abs(vA[idx_of[n]])
                                                + abs(vB[idx_of[n]])))
        for n in nodes_in_pw:
            rows.append((label_for(n, mg), "node", n))
    if unassigned:
        unassigned_sorted = sorted(unassigned,
                                     key=lambda n: -(abs(vA[idx_of[n]])
                                                       + abs(vB[idx_of[n]])))
        for n in unassigned_sorted:
            rows.append((label_for(n, mg), "node", n))

    # Determine x-range from max abs loading among seeds
    max_abs = max((max(abs(vA[idx_of[n]]), abs(vB[idx_of[n]]))
                    for n in seed_set if n in idx_of), default=0.001)

    for ri, (label, _kind, n) in enumerate(rows):
        y = -ri
        i = idx_of[n]
        la = vA[i]; lb = vB[i]
        ax.scatter([-abs(la)], [y], s=180,
                    c=[POS_COLOR if la > 0 else NEG_COLOR],
                    edgecolor="black", linewidth=0.6, zorder=3)
        ax.scatter([+abs(lb)], [y], s=180,
                    c=[POS_COLOR if lb > 0 else NEG_COLOR],
                    edgecolor="black", linewidth=0.6, zorder=3)
        # Center label with white background
        ax.text(0, y, label, ha="center", va="center", fontsize=7.5,
                 color="#101010", zorder=10,
                 bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                            edgecolor="#ccc", linewidth=0.4, alpha=0.95))

    ax.axvline(0, color="#444", linewidth=0.7, alpha=0.7)
    ax.set_xlim(-max_abs * 1.25, max_abs * 1.25)
    ax.set_ylim(-len(rows) - 0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel("|loading|", fontsize=9)
    ax.set_title(f"←  {cA_label}                {cB_label}  →",
                  fontsize=10, fontweight="bold")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    # Sign legend (small) in lower-left corner
    ax.text(-max_abs * 1.20, -len(rows) + 0.0,
             "● red = +basin   ● blue = −basin",
             fontsize=7, color="#555", va="bottom")


def _pathway_membership(mg, sub_full, nodes, pathway_names, leaf_pathways):
    """Return {pathway_id: frozenset(member_node_ids)} for every Reactome
    leaf pathway with ≥3 substrate members. Members are reactions tagged
    with the pathway, plus their direct gene/metabolite neighbors in the
    substrate (matches the basin-builder's pathway-overlay convention).
    """
    pw_members = defaultdict(set)
    for n in nodes:
        attrs = mg.graph.nodes.get(n, {})
        if attrs.get("node_type") != "reaction":
            continue
        for pid in (attrs.get("pathways") or []):
            if pid not in leaf_pathways or pid not in pathway_names:
                continue
            pw_members[pid].add(n)
            for nb in sub_full.neighbors(n):
                nb_t = mg.graph.nodes.get(nb, {}).get("node_type")
                if nb_t in ("gene", "metabolite", "reaction"):
                    pw_members[pid].add(nb)
    return {pid: frozenset(m) for pid, m in pw_members.items() if len(m) >= 3}


def hypergeom_top_pathways(vec, nodes, pw_members, pathway_names,
                             foreground_k=100, max_paths=12):
    """Hypergeometric enrichment of a PC's top-|loading| nodes against
    Reactome leaf pathways. Returns list of (pathway_id, pname, neglog10p,
    n_in_fg, n_in_path) sorted by neglog10p descending. Bonferroni adjusted.
    """
    from scipy.stats import hypergeom
    N = len(nodes)
    abs_load = np.abs(vec)
    fg_idx = set(np.argsort(-abs_load)[:foreground_k].tolist())
    fg_ids = {nodes[i] for i in fg_idx}
    K = len(fg_ids)
    rows = []
    for pid, members in pw_members.items():
        n = len(members)
        k = len(fg_ids & members)
        if k < 2:
            continue
        p = float(hypergeom.sf(k - 1, N, n, K))
        rows.append((pid, pathway_names.get(pid, pid), p, k, n))
    # Bonferroni
    if not rows:
        return []
    n_tests = len(rows)
    rows = [(pid, name, min(1.0, p * n_tests), k, n) for pid, name, p, k, n in rows]
    rows.sort(key=lambda r: r[2])
    out = []
    for pid, name, p_adj, k, n in rows[:max_paths]:
        neglogp = -np.log10(max(p_adj, 1e-300))
        out.append((pid, name, neglogp, k, n))
    return out


def render_hypergeom_back_to_back(ax, mg, nodes, vA, vB, pw_members,
                                    pathway_names, cA_label, cB_label,
                                    foreground_k=100, max_paths=12):
    """Back-to-back hypergeometric pathway-enrichment bars.

    For each PC (vA, vB) separately, run hypergeometric enrichment of the
    top-K nodes by |loading| against Reactome leaf pathways. Union the two
    cohorts' top-significant pathways, sort by combined -log10(p),
    then plot back-to-back bars: PCₐ -log10(padj) on the left axis,
    PCᵦ -log10(padj) on the right.
    """
    top_a = hypergeom_top_pathways(vA, nodes, pw_members, pathway_names,
                                      foreground_k=foreground_k,
                                      max_paths=max_paths)
    top_b = hypergeom_top_pathways(vB, nodes, pw_members, pathway_names,
                                      foreground_k=foreground_k,
                                      max_paths=max_paths)
    a_by_id = {r[0]: r for r in top_a}
    b_by_id = {r[0]: r for r in top_b}
    union_ids = list(dict.fromkeys([r[0] for r in top_a] + [r[0] for r in top_b]))
    rows = []
    for pid in union_ids:
        ra = a_by_id.get(pid, (pid, pathway_names.get(pid, pid), 0.0, 0, 0))
        rb = b_by_id.get(pid, (pid, pathway_names.get(pid, pid), 0.0, 0, 0))
        combined = ra[2] + rb[2]
        rows.append((pid, ra[1], ra[2], rb[2], combined))
    rows.sort(key=lambda r: -r[4])
    rows = rows[:max_paths]
    if not rows:
        ax.text(0.5, 0.5, "(no enriched pathways)", ha="center", va="center",
                 fontsize=9, color="#888", transform=ax.transAxes)
        ax.set_axis_off()
        return
    max_x = max(max(r[2], r[3]) for r in rows) * 1.05
    # Reserve ~half the panel width on the right for pathway labels so bars
    # are unobstructed (was: labels overlaid on bars, hid the smaller side).
    label_x = max_x * 1.15
    bar_height = 0.7
    for ri, (pid, name, neglog_a, neglog_b, _) in enumerate(rows):
        y = -ri
        ax.barh(y, -neglog_a, height=bar_height, color=POS_COLOR,
                 alpha=0.8, edgecolor="black", linewidth=0.4)
        ax.barh(y, +neglog_b, height=bar_height, color="#2D5BBF",
                 alpha=0.8, edgecolor="black", linewidth=0.4)
        # Pathway name on the right margin, left-justified, with a leader
        # line so the row is unambiguous even when many bars are tiny.
        wrapped = "\n".join(textwrap.wrap(name, width=42))
        ax.text(label_x, y, wrapped, ha="left", va="center", fontsize=7.5,
                 color="#101010", linespacing=1.05, zorder=10)
    ax.axvline(0, color="#444", linewidth=0.7, alpha=0.7)
    # Significance gridlines at p_adj = 0.05 and 0.01 (-log10 = 1.30, 2.0)
    for tick in (1.30, 2.0, 3.0):
        if tick <= max_x:
            ax.axvline(-tick, color="#aaa", linewidth=0.4, linestyle=":")
            ax.axvline(+tick, color="#aaa", linewidth=0.4, linestyle=":")
    # Extend xlim well past the right bar so the label band fits on canvas
    ax.set_xlim(-max_x, max_x * 2.5)
    ax.set_ylim(-len(rows) - 0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel("−log₁₀(p_adj)   (Bonferroni)", fontsize=9)
    ax.set_title(f"←  {cA_label}                {cB_label}  →",
                  fontsize=10, fontweight="bold")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    ax.text(-max_x * 0.95, 0.6,
             f"top-{foreground_k} |loading| nodes per PC; dotted = p_adj 0.05/0.01/0.001",
             fontsize=7, color="#555", va="bottom")


def render_pair_panel(ax, mg, sub_full, nodes, vA, vB, cA, pA, cB, pB, cs,
                       same_disease,
                       coverage_target=COVERAGE_TARGET,
                       global_layout=None,
                       pathway_names=None,
                       leaf_pathways=None,
                       pathway_parent=None,
                       show_bridges=True):
    """Render one shared-substructure panel.

    vA, vB: unit-norm PC loadings (length n_nodes) for cohort A and cohort B.
    Seed selection is *coverage-based greedy-connected* (not top-K). Two
    passes — constructive (vA·vB > 0) and destructive (vA·vB < 0) — each
    starts from its single highest-|contribution| node and grows by adding
    the next-strongest candidate that is within 2 substrate hops of any
    already-selected seed, stopping when cumulative |contribution| reaches
    COVERAGE_TARGET of that sign class's total — no K cap.
    This produces representative, connected biology instead of a scatter of
    top-individual contributors.
    Seeds carry the sign-pair coloring (fill = PCₐ sign, edge = PCᵦ sign);
    bridges drawn small and gray after a Steiner solve through the substrate.
    """
    # ── PC sign canonicalization ─────────────────────────────────────────
    # PCA component signs are arbitrary, so two cohorts measuring the same
    # biology may end up with PCs pointing in opposite directions. The
    # reported |cos| is invariant to this, but the per-node sign comparison
    # is NOT — every constructive feature would appear "sign-flipped" if
    # cohort B's PC happened to be flipped. Canonicalize by flipping vB
    # so that vA @ vB ≥ 0, then run all downstream selection + coloring on
    # the aligned PCs.
    if float(vA @ vB) < 0:
        vB = -vB
    contrib = vA * vB                # signed contribution per node to cos
    abs_contrib = np.abs(contrib)
    idx_of = {nid: i for i, nid in enumerate(nodes)}

    # Single-pool selection by signed contribution. After canonicalization,
    # the cosine is non-negative, so the top contributors are the biggest
    # constructive contributors (the carriers of the chord). Sign-flipped
    # nodes only appear if they have comparable mass to the carriers, not
    # just because |contribution| ranks them.
    constructive_mask = contrib > 0
    sel_seeds, sel_cov = greedy_connected_seeds(
        contrib, constructive_mask, nodes, sub_full,
        coverage_target=COVERAGE_TARGET)
    seed_nids = list(sel_seeds)
    top_idx = [idx_of[n] for n in seed_nids if n in idx_of]

    all_nids, bridges, seed_set, dropped_seeds = steiner_bridge(
        seed_nids, sub_full, max_extra=MAX_BRIDGES)
    if dropped_seeds:
        print(f"    {cA}·PC{pA} ↔ {cB}·PC{pB}: dropped {len(dropped_seeds)} "
              f"seeds outside main CC", flush=True)
    sub_g = sub_full.subgraph(all_nids).copy()
    # Aspect-stretch the layout so the panel's horizontal real estate is
    # actually used — panels are landscape (wider than tall), so multiply
    # x coordinates by ASPECT_STRETCH after the spring solve.
    ASPECT_STRETCH = 1.7
    if global_layout is not None:
        layout = {}
        backup = nx.spring_layout(sub_g, k=2.5, iterations=400, seed=42,
                                   scale=1.0) if sub_g.number_of_edges() > 0 else {}
        for n in all_nids:
            if n in global_layout:
                layout[n] = global_layout[n]
            elif n in backup:
                layout[n] = backup[n]
            else:
                layout[n] = np.array([0.0, 0.0])
    elif sub_g.number_of_edges() > 0:
        # True per-panel force-directed: larger k pushes nodes apart, more
        # iterations let the layout converge to a separated arrangement.
        layout = nx.spring_layout(sub_g, k=2.5, iterations=400, seed=42,
                                    scale=1.0)
    else:
        n = len(all_nids)
        layout = {nid: np.array([np.cos(2*np.pi*i/n), np.sin(2*np.pi*i/n)])
                   for i, nid in enumerate(all_nids)}
    # Stretch x so nodes spread across the landscape panel
    layout = {nid: np.array([p[0] * ASPECT_STRETCH, p[1]])
              for nid, p in layout.items()}

    # Reactome pathway grouping for THIS pair's seeds. Drawn BEFORE edges/nodes
    # so the boxes sit in the background. Uses ≥2-member rule (like basin).
    panel_pathway_groups = []
    if pathway_names and leaf_pathways is not None and pathway_parent is not None:
        panel_pathway_groups = basin_mod.find_pathway_groups(
            list(seed_set), mg, sub_full, pathway_names,
            leaf_pathways, pathway_parent, min_members=2)
        if panel_pathway_groups:
            _draw_pathway_groups(ax, panel_pathway_groups, layout,
                                   base_cmap="tab10")

    # When bridges are hidden, drop bridge↔seed and bridge↔bridge edges —
    # otherwise they'd dangle into empty space and look broken.
    if show_bridges:
        nx.draw_networkx_edges(sub_g, layout, edge_color=EDGE_GRAY,
                                 alpha=0.55, width=1.2, ax=ax)
    else:
        seed_only_edges = [(u, v) for u, v in sub_g.edges()
                             if u in seed_set and v in seed_set]
        if seed_only_edges:
            nx.draw_networkx_edges(sub_g, layout, edgelist=seed_only_edges,
                                     edge_color=EDGE_GRAY,
                                     alpha=0.55, width=1.2, ax=ax)

    max_contrib = abs_contrib[top_idx].max() + 1e-12
    pos_contrib = sum(1 for i in top_idx if contrib[i] > 0)
    neg_contrib = len(top_idx) - pos_contrib

    # Cosine coverage: how much of the dot product PCₐ·PCᵦ do the displayed
    # seeds actually account for? Users intuitively expect the top contributors
    # to dominate, but the cosine is typically distributed across many small
    # contributors — so this fraction is usually well under 100%.
    cos_total = float(vA @ vB)           # signed; equals cs * sign-of-vA·vB
    seed_sum = float(sum(contrib[i] for i in top_idx))
    cov_pct = 100.0 * seed_sum / cos_total if abs(cos_total) > 1e-12 else 0.0
    # Also report absolute-value coverage: |seed contributions| / sum of |all
    # contributions| — how much of the total magnitude is shown.
    abs_total = float(np.sum(abs_contrib))
    abs_seed = float(np.sum(abs_contrib[top_idx]))
    abs_cov_pct = 100.0 * abs_seed / abs_total if abs_total > 1e-12 else 0.0

    shape_for_type = {"gene": "o", "metabolite": "h", "reaction": "s"}
    # Bridges first (smaller, behind), then seeds (bigger, in front).
    # Per-pair panels hide bridges by default to reduce visual crowding —
    # the bridges still participate in Steiner connectivity + layout, just
    # not drawn. Set show_bridges=True to enable the small gray dots.
    if show_bridges:
        for tname, shape in shape_for_type.items():
            bridge_nids = [n for n in bridges
                            if mg.graph.nodes.get(n, {}).get("node_type") == tname]
            if bridge_nids:
                nx.draw_networkx_nodes(sub_g, layout, nodelist=bridge_nids,
                                        node_size=80, node_color=BRIDGE_COLOR,
                                        edgecolors="#888", linewidths=0.8,
                                        node_shape=shape, ax=ax)
    for tname, shape in shape_for_type.items():
        seed_nids_t = [n for n in seed_set
                        if mg.graph.nodes.get(n, {}).get("node_type") == tname]
        if not seed_nids_t: continue
        sizes, face, edge = [], [], []
        for n in seed_nids_t:
            i_global = idx_of[n]
            size_rel = (abs_contrib[i_global] / max_contrib) ** 0.7
            sizes.append(220 + 1600 * size_rel)
            face.append(POS_COLOR if vA[i_global] > 0 else NEG_COLOR)
            edge.append(POS_COLOR if vB[i_global] > 0 else NEG_COLOR)
        nx.draw_networkx_nodes(sub_g, layout, nodelist=seed_nids_t,
                                node_size=sizes, node_color=face,
                                edgecolors=edge, linewidths=3.5,
                                node_shape=shape, ax=ax)

    # Labels only on SEEDS — bridges are unlabeled to keep the panel readable
    for n in seed_set:
        if n not in layout: continue
        x, y = layout[n]
        label = label_for(n, mg)
        ax.text(x, y - 0.08, label, ha="center", va="top",
                 fontsize=7.0, color="#0a0a0a", fontweight="bold",
                 zorder=20,
                 bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                           edgecolor="none", alpha=0.85))

    # Frame + title (account for the horizontal stretch)
    ax.set_xlim(-1.3 * ASPECT_STRETCH, 1.3 * ASPECT_STRETCH)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal"); ax.set_axis_off()

    short_a = (cA.replace("CPTAC_", "").replace("GSE65391_", "")
                  .replace("GSE65682_", "").replace("GSE89408_", "")
                  .replace("_glioma", "·gli").replace("_COVID", "·COVID")
                  .replace("_IBD_CD", "·CD"))
    short_b = (cB.replace("CPTAC_", "").replace("GSE65391_", "")
                  .replace("GSE65682_", "").replace("GSE89408_", "")
                  .replace("_glioma", "·gli").replace("_COVID", "·COVID")
                  .replace("_IBD_CD", "·CD"))
    title = f"{short_a}·PC{pA}  ↔  {short_b}·PC{pB}"
    sub_title = (f"|cos|={cs:.3f}    {pos_contrib} constructive · "
                 f"{neg_contrib} sign-flipped    ·    "
                 f"shown seeds = {cov_pct:.0f}% of net cos "
                 f"({abs_cov_pct:.0f}% of |contrib|)")
    badge_color = "#10b981" if same_disease else "#6b7280"
    ax.set_title(title, fontsize=11, fontweight="bold", pad=12)
    ax.text(0.5, 1.06, sub_title, transform=ax.transAxes,
             ha="center", va="bottom", fontsize=8.5, color=badge_color)

    # Hand back per-pair state the caller can use to render a companion
    # back-to-back dot plot (loadings per seed, sorted by pathway).
    # vA_canon / vB_canon are the sign-canonicalized PCs — the caller MUST
    # use these for the dot plot so the dots match the network coloring.
    return {
        "seed_set":         seed_set,
        "pathway_groups":   panel_pathway_groups,
        "short_a":          short_a,
        "short_b":          short_b,
        "pc_a_label":       f"{short_a}·PC{pA}",
        "pc_b_label":       f"{short_b}·PC{pB}",
        "cov_pct":          cov_pct,
        "abs_cov_pct":      abs_cov_pct,
        "vA_canon":         vA,
        "vB_canon":         vB,
    }


def _draw_pathway_groups(ax, groups, layout, base_cmap="tab10"):
    """Draw Cytoscape-style rounded-box outlines around pathway-grouped node
    sets. `groups` is the list returned by basin_mod.find_pathway_groups —
    each item is ((pathway_id, pathway_name), set_of_node_ids)."""
    if not groups: return
    cmap = plt.get_cmap(base_cmap)
    for gi, ((_pid, pname), gnodes) in enumerate(groups):
        pts = np.array([layout[n] for n in gnodes if n in layout])
        if len(pts) < 2: continue
        color = cmap(gi % cmap.N)
        edge_color = (color[0]*0.55, color[1]*0.55, color[2]*0.55, 1.0)
        pad = 0.05
        x_min, x_max = pts[:, 0].min() - pad, pts[:, 0].max() + pad
        y_min, y_max = pts[:, 1].min() - pad, pts[:, 1].max() + pad
        w, h = x_max - x_min, y_max - y_min
        rect = mpatches.FancyBboxPatch(
            (x_min, y_min), w, h,
            boxstyle="round,pad=0.03,rounding_size=0.06",
            facecolor=color, alpha=0.18,
            edgecolor=edge_color, linewidth=1.8, zorder=0)
        ax.add_patch(rect)
        # Pathway label above the box
        pname_wrapped = "\n".join(textwrap.wrap(pname, width=26))
        ax.text((x_min + x_max) / 2, y_max + 0.04, pname_wrapped,
                 fontsize=8, color=edge_color, ha="center", va="bottom",
                 fontweight="bold", style="italic", zorder=15,
                 linespacing=1.05,
                 bbox=dict(facecolor="white", alpha=0.9,
                           edgecolor=edge_color, linewidth=0.6, pad=2))


def render_cross_disease_overlap(out_path, mg, sub_full, nodes, comps, top,
                                   pathway_names=None, leaf_pathways=None,
                                   pathway_parent=None):
    """One figure showing which substrate nodes recur as top contributors
    across multiple top-N chord pairs — and across distinct diseases.

    For every top-N pair, collect the constructive + destructive seed
    nodes (same pools as the per-pair panels). For each unique substrate
    node, count (a) how many top-N pairs it appeared in, and (b) how
    many distinct diseases those appearances span. Recurrent nodes
    (≥ 2 pairs) get drawn as a single Steiner-connected sub-graph.
    Node size ∝ pair appearances; fill color encodes disease span
    (1 = within-disease repeat, 2 = mild cross-disease, 3+ = broad
    cross-disease hub).
    """
    idx_of = {nid: i for i, nid in enumerate(nodes)}
    node_pairs = defaultdict(list)   # nid -> list of (pair_idx, dA, dB)
    for i, row in top.iterrows():
        cA, pA, dA = row["cohort"], int(row["pc"]), row["disease"]
        cB, pB, dB = (row["best_match_cohort"], int(row["best_match_pc"]),
                       row["best_match_disease"])
        comp_a = comps.get(cA); comp_b = comps.get(cB)
        if comp_a is None or comp_b is None: continue
        if comp_a.shape[0] < pA or comp_b.shape[0] < pB: continue
        vA = comp_a[pA - 1]; vB = comp_b[pB - 1]
        # Canonicalize cohort-B PC sign so the per-node sign comparison is
        # not flipped by arbitrary PCA orientation; then single-pool greedy-
        # connected selection by signed contribution.
        if float(vA @ vB) < 0:
            vB = -vB
        contrib = vA * vB
        constructive_mask = contrib > 0
        sel_seeds, _ = greedy_connected_seeds(
            contrib, constructive_mask, nodes, sub_full,
            coverage_target=COVERAGE_TARGET)
        for nid in sel_seeds:
            node_pairs[nid].append((i + 1, dA, dB))

    recurrent = {nid: aps for nid, aps in node_pairs.items() if len(aps) >= 2}
    print(f"  cross-disease overlap: {len(recurrent)} recurrent substrate nodes "
          f"(out of {len(node_pairs)} unique seeds across {len(top)} pairs)",
          flush=True)
    if not recurrent:
        print("  no recurrent nodes — skipping figure", flush=True)
        return None

    disease_span = {
        nid: len({d for _, dA, dB in aps for d in (dA, dB)})
        for nid, aps in recurrent.items()
    }
    pair_count = {nid: len(aps) for nid, aps in recurrent.items()}

    # Steiner-connect the recurrent seeds (keep bridges modest)
    all_nids, bridges, seed_set, dropped = steiner_bridge(
        list(recurrent.keys()), sub_full, max_extra=40)
    if dropped:
        print(f"  dropped {len(dropped)} recurrent nodes outside main CC",
              flush=True)
    sub_g = sub_full.subgraph(all_nids).copy()

    # Spring layout — larger figure → larger k for spacing
    if sub_g.number_of_edges() > 0:
        layout = nx.spring_layout(sub_g, k=1.6, iterations=600, seed=42,
                                    scale=1.0)
    else:
        n = len(all_nids)
        layout = {nid: np.array([np.cos(2*np.pi*i/n), np.sin(2*np.pi*i/n)])
                   for i, nid in enumerate(all_nids)}

    fig, ax = plt.subplots(figsize=(20, 16))
    fig.patch.set_facecolor("white")

    # Pathway grouping — draw BEFORE edges/nodes so the boxes sit in the back
    pathway_groups = []
    if pathway_names and leaf_pathways is not None and pathway_parent is not None:
        pathway_groups = basin_mod.find_pathway_groups(
            list(seed_set), mg, sub_full, pathway_names,
            leaf_pathways, pathway_parent, min_members=2)
        print(f"  {len(pathway_groups)} pathway groups span ≥2 recurrent nodes",
              flush=True)
        _draw_pathway_groups(ax, pathway_groups, layout, base_cmap="tab10")

    # Edges
    nx.draw_networkx_edges(sub_g, layout, edge_color=EDGE_GRAY,
                             alpha=0.45, width=1.0, ax=ax)

    # Bridges (small gray, behind)
    shape_for_type = {"gene": "o", "metabolite": "h", "reaction": "s"}
    for tname, shape in shape_for_type.items():
        bnids = [n for n in bridges
                 if mg.graph.nodes.get(n, {}).get("node_type") == tname]
        if bnids:
            nx.draw_networkx_nodes(sub_g, layout, nodelist=bnids,
                                    node_size=70, node_color=BRIDGE_COLOR,
                                    edgecolors="#888", linewidths=0.6,
                                    node_shape=shape, ax=ax)

    # Recurrent seeds — color by disease span, size by pair count
    # Colormap: 1 disease = light blue (within-disease repeat),
    #           2 diseases = orange, 3+ diseases = bright red
    color_for_span = {1: "#74add1", 2: "#fdae61", 3: "#d73027"}
    def span_color(s):
        return color_for_span.get(min(s, 3), "#d73027")

    max_count = max(pair_count.values())
    for tname, shape in shape_for_type.items():
        rnids = [n for n in seed_set
                 if mg.graph.nodes.get(n, {}).get("node_type") == tname]
        if not rnids: continue
        sizes, face = [], []
        for n in rnids:
            pc_n = pair_count[n]
            sizes.append(250 + 1400 * (pc_n / max_count))
            face.append(span_color(disease_span[n]))
        nx.draw_networkx_nodes(sub_g, layout, nodelist=rnids,
                                node_size=sizes, node_color=face,
                                edgecolors="#222", linewidths=1.5,
                                node_shape=shape, ax=ax)

    # Labels on the most recurrent nodes (top 25 by pair_count)
    label_pool = sorted(seed_set, key=lambda n: -pair_count[n])[:25]
    for n in label_pool:
        if n not in layout: continue
        x, y = layout[n]
        ax.text(x, y - 0.06, f"{label_for(n, mg)}\n({pair_count[n]}p / "
                              f"{disease_span[n]}d)",
                 ha="center", va="top", fontsize=9, color="#101010",
                 fontweight="bold", zorder=20, linespacing=1.0,
                 bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                           edgecolor="none", alpha=0.85))

    ax.set_axis_off()
    ax.set_title("Cross-chord, cross-disease overlap — substrate nodes that "
                  "show up as top contributors in ≥ 2 of the top-10 chord pairs.\n"
                  "Color = how many DISTINCT diseases the node appears across "
                  "(blue 1, orange 2, red 3+).  Size = # of top-10 pairs the "
                  "node appears in.  Gray dots = Steiner bridges connecting the "
                  "recurrent nodes through the full substrate.",
                  fontsize=13, pad=14, loc="left")

    # Color legend
    from matplotlib.lines import Line2D
    legend_handles = [
        mpatches.Patch(facecolor=color_for_span[1], edgecolor="#222",
                        label="1 disease (within-disease repeat)"),
        mpatches.Patch(facecolor=color_for_span[2], edgecolor="#222",
                        label="2 diseases (cross-disease, modest)"),
        mpatches.Patch(facecolor=color_for_span[3], edgecolor="#222",
                        label="3+ diseases (cross-disease hub)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
                markeredgecolor="#888", markersize=10, label="gene"),
        Line2D([0], [0], marker="h", color="w", markerfacecolor="#888",
                markeredgecolor="#888", markersize=11, label="metabolite"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#888",
                markeredgecolor="#888", markersize=10, label="reaction"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=6,
                frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.0))

    plt.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  wrote {out_path}", flush=True)

    # Companion TSV: ranked list of recurrent nodes with their biology +
    # which pairs / diseases they appear across.
    tsv_path = out_path.with_suffix(".tsv")
    rows = []
    for nid, aps in sorted(recurrent.items(),
                             key=lambda kv: (-len(kv[1]),
                                              -disease_span[kv[0]],
                                              kv[0])):
        attrs = mg.graph.nodes.get(nid, {})
        sym = attrs.get("symbol") or attrs.get("name") \
              or attrs.get("display_name") or nid
        ntype = attrs.get("node_type", "?")
        diseases = sorted({d for _, dA, dB in aps for d in (dA, dB)})
        pair_idx_list = sorted({i for i, _, _ in aps})
        rows.append({
            "node_id":      nid,
            "symbol":       sym,
            "node_type":    ntype,
            "n_pairs":      len(aps),
            "n_diseases":   disease_span[nid],
            "diseases":     "|".join(diseases),
            "pair_indices": ",".join(str(i) for i in pair_idx_list),
        })
    df_overlap = pd.DataFrame(rows)
    df_overlap.to_csv(tsv_path, sep="\t", index=False)
    print(f"  wrote {tsv_path}  ({len(df_overlap)} recurrent nodes)", flush=True)

    # Pathway-level rollup: for each Reactome pathway that contains ≥2 of the
    # recurrent nodes, record its members + the pair/disease coverage.
    if pathway_names and leaf_pathways is not None and pathway_parent is not None:
        pw_rows = []
        for (pid, pname), gnodes in pathway_groups:
            members = sorted(gnodes, key=lambda n: -pair_count.get(n, 0))
            n_pairs_pw = len({i for n in members for i, _, _ in recurrent.get(n, [])})
            diseases_pw = sorted({d for n in members
                                    for _, dA, dB in recurrent.get(n, [])
                                    for d in (dA, dB)})
            member_syms = [
                (mg.graph.nodes.get(n, {}).get("symbol")
                  or mg.graph.nodes.get(n, {}).get("name")
                  or mg.graph.nodes.get(n, {}).get("display_name")
                  or n)
                for n in members
            ]
            pw_rows.append({
                "pathway_id":    pid,
                "pathway_name":  pname,
                "n_members":     len(members),
                "n_pairs":       n_pairs_pw,
                "n_diseases":    len(diseases_pw),
                "diseases":      "|".join(diseases_pw),
                "members":       " | ".join(str(s) for s in member_syms),
            })
        if pw_rows:
            pw_tsv = out_path.parent / (out_path.stem + "_pathways.tsv")
            pd.DataFrame(pw_rows).to_csv(pw_tsv, sep="\t", index=False)
            print(f"  wrote {pw_tsv}  ({len(pw_rows)} conserved pathways)",
                  flush=True)
    return out_path


def main():
    # NOTE: the data-poor-pair filter (below, after substrate load) will
    # replace this initial `top` selection. We still build it here to give
    # the substrate + alpha-pc-components loop the right cohort set.
    df = pd.read_csv(RESULTS / "pc_alignment_best_cross_cohort_match.tsv", sep="\t")
    df["edge"] = df.apply(
        lambda r: tuple(sorted([(r["cohort"], r["pc"]),
                                  (r["best_match_cohort"], r["best_match_pc"])])),
        axis=1,
    )
    top = (df.drop_duplicates("edge")
             .sort_values("best_match_cos", ascending=False)
             .reset_index(drop=True))   # full list — filter later
    print(f"Candidate cross-cohort pairs: {len(top)}")
    print(top.head(TOP_N + 5)[
            ["cohort", "pc", "best_match_cohort", "best_match_pc",
              "best_match_cos", "best_match_within_disease"]]
          .to_string(index=False))

    mg, sub, sub_full, nodes, log_pr = load_substrate()

    needed = set(top["cohort"]) | set(top["best_match_cohort"])
    comps = {}
    for c in needed:
        comps[c] = alpha_pc_components(c, log_pr, n_components=5)
        if comps[c] is None:
            print(f"  ⚠ no F for {c}", flush=True)

    # ---- Reactome pathway tables (loaded once, used by per-pair + overlap) ----
    pathway_names = basin_mod.load_pathway_names()
    leaf_pathways, pathway_parent = basin_mod.load_leaf_pathways()
    print(f"  loaded {len(pathway_names)} pathway names, "
          f"{len(leaf_pathways)} leaf pathways", flush=True)
    print("  building pathway → substrate-member map for hypergeometric tests…",
          flush=True)
    pw_members = _pathway_membership(mg, sub_full, nodes, pathway_names,
                                       leaf_pathways)
    print(f"  pathway membership map: {len(pw_members)} pathways with ≥3 substrate members",
          flush=True)

    # ---- Filter data-poor chord pairs from the top-N display ----
    # A "data-poor" pair is one whose greedy-connected seed selection
    # produces fewer than MIN_SEEDS_PER_PAIR seeds — they show as a
    # near-empty panel that doesn't add information. Drop those and pull
    # in the next-ranked pair from the full TSV.
    MIN_SEEDS_PER_PAIR = 5
    df_all = pd.read_csv(RESULTS / "pc_alignment_best_cross_cohort_match.tsv",
                          sep="\t")
    df_all["edge"] = df_all.apply(
        lambda r: tuple(sorted([(r["cohort"], r["pc"]),
                                  (r["best_match_cohort"], r["best_match_pc"])])),
        axis=1,
    )
    df_all = (df_all.drop_duplicates("edge")
                     .sort_values("best_match_cos", ascending=False)
                     .reset_index(drop=True))
    print(f"\nScreening {len(df_all)} candidate pairs for ≥{MIN_SEEDS_PER_PAIR} seeds…",
          flush=True)
    kept_rows = []
    for _, row in df_all.iterrows():
        if len(kept_rows) >= TOP_N:
            break
        cA, pA = row["cohort"], int(row["pc"])
        cB, pB = row["best_match_cohort"], int(row["best_match_pc"])
        comp_a = comps.get(cA); comp_b = comps.get(cB)
        if comp_a is None or comp_b is None: continue
        if comp_a.shape[0] < pA or comp_b.shape[0] < pB: continue
        vA_p = comp_a[pA - 1]; vB_p = comp_b[pB - 1]
        if float(vA_p @ vB_p) < 0:
            vB_p = -vB_p
        contrib_p = vA_p * vB_p
        seeds, _ = greedy_connected_seeds(
            contrib_p, contrib_p > 0, nodes, sub_full,
            coverage_target=COVERAGE_TARGET)
        n_seeds = len(seeds)
        verdict = "kept" if n_seeds >= MIN_SEEDS_PER_PAIR else "skipped"
        print(f"  {cA}·PC{pA} ↔ {cB}·PC{pB}  |cos|={row['best_match_cos']:.3f}  "
              f"seeds={n_seeds}  → {verdict}", flush=True)
        if n_seeds >= MIN_SEEDS_PER_PAIR:
            kept_rows.append(row)
    top = pd.DataFrame(kept_rows).reset_index(drop=True)
    print(f"\nUsing {len(top)} kept pairs for figure", flush=True)

    # (Union-layout pinning removed — per-pair panels use their own
    # force-directed layout for clearer separation. The Fig 2h cross-
    # disease overlap handles inter-panel recurrence tracking.)
    global_layout = None

    from matplotlib.lines import Line2D

    def build_legend_handles():
        items = [
            ("PC_A fill = +basin",   POS_COLOR, "fill"),
            ("PC_A fill = −basin",   NEG_COLOR, "fill"),
            ("PC_B edge = +basin",   POS_COLOR, "edge"),
            ("PC_B edge = −basin",   NEG_COLOR, "edge"),
        ]
        out = []
        for label, color, kind in items:
            if kind == "fill":
                out.append(mpatches.Patch(facecolor=color, edgecolor="black",
                                            linewidth=0.8, label=label))
            else:
                out.append(mpatches.Patch(facecolor="white", edgecolor=color,
                                            linewidth=3, label=label))
        out += [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
                    markeredgecolor="#888", markersize=10, label="gene"),
            Line2D([0], [0], marker="h", color="w", markerfacecolor="#888",
                    markeredgecolor="#888", markersize=11, label="metabolite"),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#888",
                    markeredgecolor="#888", markersize=10, label="reaction"),
        ]
        return out

    # --- Per-pair individual PNGs (readable labels) ---
    print(f"\nRendering {len(top)} individual per-pair panels…", flush=True)
    individual_paths = []
    for i, row in top.iterrows():
        cA, pA = row["cohort"], int(row["pc"])
        cB, pB = row["best_match_cohort"], int(row["best_match_pc"])
        cs = float(row["best_match_cos"])
        same_disease = bool(row["best_match_within_disease"])
        comp_a = comps.get(cA); comp_b = comps.get(cB)
        if comp_a is None or comp_b is None:
            continue
        vA = comp_a[pA - 1]
        vB = comp_b[pB - 1]

        # Per-pair figure: THREE-column layout
        #   col 1 — network sub-graph (left)
        #   col 2 — back-to-back loading dot plot (middle)
        #   col 3 — back-to-back hypergeometric pathway enrichment (right)
        fig_i = plt.figure(figsize=(30, 11))
        gs_i = fig_i.add_gridspec(1, 3, width_ratios=[1.55, 1.0, 1.0],
                                    wspace=0.10)
        ax_net = fig_i.add_subplot(gs_i[0, 0])
        ax_dot = fig_i.add_subplot(gs_i[0, 1])
        ax_hyp = fig_i.add_subplot(gs_i[0, 2])
        panel_info = render_pair_panel(
            ax_net, mg, sub_full, nodes, vA, vB,
            cA, pA, cB, pB, cs, same_disease,
            global_layout=None,
            pathway_names=pathway_names,
            leaf_pathways=leaf_pathways,
            pathway_parent=pathway_parent)
        render_back_to_back_dotplot(
            ax_dot, mg, nodes,
            panel_info["vA_canon"], panel_info["vB_canon"],
            panel_info["seed_set"], panel_info["pathway_groups"],
            cA_label=panel_info["pc_a_label"],
            cB_label=panel_info["pc_b_label"])
        render_hypergeom_back_to_back(
            ax_hyp, mg, nodes,
            panel_info["vA_canon"], panel_info["vB_canon"],
            pw_members, pathway_names,
            cA_label=panel_info["pc_a_label"],
            cB_label=panel_info["pc_b_label"])
        fig_i.suptitle(
            f"Chord #{i+1} co-structure — greedy-connected seeds covering "
            f"{int(COVERAGE_TARGET*100)}% of each sign's contribution to cos "
            "(gray dots = Steiner bridges)",
            fontsize=12, y=0.97)
        fig_i.legend(handles=build_legend_handles(), loc="lower center",
                      ncol=7, frameon=False, fontsize=10,
                      bbox_to_anchor=(0.5, 0.0))
        out_i = FIG / f"fig_chord_costructure_pair{i+1}.png"
        fig_i.savefig(out_i, dpi=160, bbox_inches="tight", facecolor="white")
        plt.close(fig_i)
        individual_paths.append(out_i)
        print(f"  wrote {out_i}", flush=True)

    # --- Grid overview: 10 rows × 3 columns
    #     col 1 — network sub-graph; col 2 — loading dot plot;
    #     col 3 — hypergeometric pathway enrichment back-to-back.
    nrows = len(top)
    fig = plt.figure(figsize=(26, 70))
    gs_main = fig.add_gridspec(nrows, 3, width_ratios=[1.55, 1.0, 1.0],
                                 hspace=0.22, wspace=0.10)
    fig.suptitle(
        f"Top cross-cohort α-PC similarities — greedy-connected seed pools "
        f"covering {int(COVERAGE_TARGET*100)}% of each sign's contribution "
        f"to cos. Left = network sub-graph; middle = back-to-back loading "
        f"dot plot; right = hypergeometric Reactome pathway enrichment "
        f"back-to-back per cohort.",
        fontsize=14, fontweight="bold", y=0.992)
    for ax_idx, (i, row) in enumerate(top.iterrows()):
        cA, pA = row["cohort"], int(row["pc"])
        cB, pB = row["best_match_cohort"], int(row["best_match_pc"])
        cs = float(row["best_match_cos"])
        same_disease = bool(row["best_match_within_disease"])
        comp_a = comps.get(cA); comp_b = comps.get(cB)
        ax_net = fig.add_subplot(gs_main[ax_idx, 0])
        ax_dot = fig.add_subplot(gs_main[ax_idx, 1])
        ax_hyp = fig.add_subplot(gs_main[ax_idx, 2])
        if comp_a is None or comp_b is None:
            ax_net.text(0.5, 0.5, f"missing F for {cA} or {cB}",
                          ha="center", va="center",
                          transform=ax_net.transAxes, fontsize=10, color="#888")
            ax_net.set_axis_off()
            ax_dot.set_axis_off()
            ax_hyp.set_axis_off()
            continue
        vA = comp_a[pA - 1]
        vB = comp_b[pB - 1]
        panel_info = render_pair_panel(
            ax_net, mg, sub_full, nodes, vA, vB,
            cA, pA, cB, pB, cs, same_disease,
            global_layout=None,
            pathway_names=pathway_names,
            leaf_pathways=leaf_pathways,
            pathway_parent=pathway_parent)
        render_back_to_back_dotplot(
            ax_dot, mg, nodes,
            panel_info["vA_canon"], panel_info["vB_canon"],
            panel_info["seed_set"], panel_info["pathway_groups"],
            cA_label=panel_info["pc_a_label"],
            cB_label=panel_info["pc_b_label"])
        render_hypergeom_back_to_back(
            ax_hyp, mg, nodes,
            panel_info["vA_canon"], panel_info["vB_canon"],
            pw_members, pathway_names,
            cA_label=panel_info["pc_a_label"],
            cB_label=panel_info["pc_b_label"])
    fig.legend(handles=build_legend_handles(), loc="lower center",
                ncol=7, frameon=False, fontsize=9,
                bbox_to_anchor=(0.5, 0.0))
    out_path = FIG / "fig_chord_costructure.png"
    plt.savefig(out_path, dpi=100, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\nWrote {out_path}  (grid overview)")
    print(f"Wrote {len(individual_paths)} per-pair PNGs as fig_chord_costructure_pair*.png")

    # Cross-disease overlap figure (one PNG, recurrent nodes across all pairs)
    print("\nBuilding cross-disease overlap figure…", flush=True)
    overlap_path = FIG / "fig_chord_costructure_cross_disease_overlap.png"
    render_cross_disease_overlap(overlap_path, mg, sub_full, nodes, comps, top,
                                   pathway_names=pathway_names,
                                   leaf_pathways=leaf_pathways,
                                   pathway_parent=pathway_parent)


if __name__ == "__main__":
    main()
