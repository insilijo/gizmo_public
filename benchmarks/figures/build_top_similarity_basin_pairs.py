"""Side-by-side basin pair figure for the top-N cross-cohort α-PC similarities.

For each of the top-N highest |cos| cross-cohort pairs (deduped), render
the two basins FRESH with a gold halo around every substrate node that
appears in BOTH basins' top-15 |loadings|. The halos visually link the
two basins so the reader can find the agreement at a glance — not just
read a list of shared gene/metabolite names.

Inputs:
  pc_alignment_best_cross_cohort_match.tsv  (one row per (cohort, PC))
  + substrate (via build_basin_signed_v2.load_render_state)

Output:
  fig_top_similarity_basin_pairs.png   (PIL composite of all pair PNGs)
  fig_basin_v6_<cohort>_pc<k>_pair<i>_highlighted.png  (per-pair, per-side basins)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from sklearn.decomposition import PCA

REPO = Path("/home/jgardner/GIZMO")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))
sys.path.insert(0, str(REPO / "benchmarks/figures"))

import build_basin_signed_v2 as basin_mod  # noqa: E402

RESULTS = REPO / "benchmarks/results"
FIG = RESULTS / "figures"
TOP_N = 10


def basin_path(cohort: str, pc: int, pair_idx: int) -> Path:
    return FIG / f"fig_basin_v6_{cohort}_pc{pc}_pair{pair_idx}_highlighted.png"


def _alpha_pc_components(cohort, log_pr, n_components=5):
    UR = RESULTS / "unsupervised"
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
    return pca.components_


def shared_top_nodes(comp_a, pc_a, comp_b, pc_b, mg, nodes, k=15,
                      visible_a=None, visible_b=None):
    """Return shared node IDs + labels + counts.

    If visible_a / visible_b are provided (sets of node IDs that will
    actually appear in each basin's rendered layout), the returned shared
    set is further filtered down to nodes present in BOTH layouts. This
    guarantees every halo appears in both panels — no orphan halos.
    """
    if comp_a is None or comp_b is None:
        return set(), [], 0, 0
    abs_a = np.abs(comp_a[pc_a - 1])
    abs_b = np.abs(comp_b[pc_b - 1])
    top_a_idx = set(np.argsort(-abs_a)[:k].tolist())
    top_b_idx = set(np.argsort(-abs_b)[:k].tolist())
    shared_idx = top_a_idx & top_b_idx
    union_idx  = top_a_idx | top_b_idx
    if visible_a is not None and visible_b is not None:
        shared_idx = {i for i in shared_idx
                       if nodes[i] in visible_a and nodes[i] in visible_b}
    shared_ids = {nodes[i] for i in shared_idx}
    shared_labels = []
    for i in sorted(shared_idx, key=lambda j: -abs_a[j] - abs_b[j]):
        attrs = mg.graph.nodes.get(nodes[i], {})
        sym = attrs.get("symbol") or attrs.get("name") or attrs.get("display_name") or nodes[i]
        ntype = attrs.get("node_type", "?")[0]
        shared_labels.append(f"{sym} [{ntype}]")
    return shared_ids, shared_labels, len(shared_idx), len(union_idx)


def main():
    df = pd.read_csv(RESULTS / "pc_alignment_best_cross_cohort_match.tsv", sep="\t")
    df["edge"] = df.apply(
        lambda r: tuple(sorted([(r["cohort"], r["pc"]),
                                  (r["best_match_cohort"], r["best_match_pc"])])),
        axis=1,
    )
    top = (df.drop_duplicates("edge")
             .sort_values("best_match_cos", ascending=False)
             .head(TOP_N)
             .reset_index(drop=True))
    print(f"Top {len(top)} cross-cohort pairs:")
    print(top[["cohort", "pc", "best_match_cohort", "best_match_pc",
                "best_match_cos", "best_match_within_disease"]]
          .to_string(index=False))

    # Load substrate + metadata once via the basin module
    print("\nLoading substrate + atlas (one-time)…", flush=True)
    state = basin_mod.load_render_state()
    mg, sub, sub_nodes, log_pr = (state["mg"], state["sub"],
                                     state["nodes"], state["log_pr"])

    # Pre-compute α-PC components for every cohort referenced in the top pairs
    needed_cohorts = set(top["cohort"]) | set(top["best_match_cohort"])
    comps = {}
    for c in needed_cohorts:
        comps[c] = _alpha_pc_components(c, log_pr, n_components=5)
        if comps[c] is None:
            print(f"  ⚠ {c}: no F or shape mismatch — pair will be skipped",
                  flush=True)

    # Pre-compute basin main-CC node sets for every (cohort, PC) entry that
    # appears in the top-N pairs. Used to filter shared halos down to nodes
    # that are guaranteed to be drawn in BOTH basins' rendered layouts.
    # (Crossovers — same node, opposite PC sign in the two cohorts — are
    # preserved because the node is still in each cohort's main basin CC,
    # just on opposite sides; only nodes that fall outside the main CC on
    # at least one side get dropped.)
    print("\nPre-computing basin main-CC node sets…", flush=True)
    basin_cc_cache = {}   # (cohort, pc_idx) -> set of node IDs in basin layout
    for _, row in top.iterrows():
        for cohort, pc in ((row["cohort"], int(row["pc"])),
                           (row["best_match_cohort"], int(row["best_match_pc"]))):
            key = (cohort, pc)
            if key in basin_cc_cache: continue
            comp = comps.get(cohort)
            if comp is None or comp.shape[0] < pc:
                basin_cc_cache[key] = set(); continue
            basin_cc_cache[key] = basin_mod.basin_main_cc_nodes(
                comp[pc - 1], sub, sub_nodes)
            print(f"    {cohort}·PC{pc}: {len(basin_cc_cache[key])} nodes in basin CC",
                  flush=True)

    # Compute shared-node sets + render highlighted basins per pair
    print(f"\nRendering {len(top)} pairs × 2 basins each with gold halos on "
          f"shared top-15 substrate nodes (filtered to nodes that appear in "
          f"BOTH basins' main CCs)…", flush=True)
    pair_renders = []   # one entry per pair: dict with side PNG paths + meta
    for i, row in top.iterrows():
        cA, pA = row["cohort"], int(row["pc"])
        cB, pB = row["best_match_cohort"], int(row["best_match_pc"])
        cs = float(row["best_match_cos"])
        same_disease = bool(row["best_match_within_disease"])

        shared_ids, shared_labels, n_shared, n_union = shared_top_nodes(
            comps.get(cA), pA, comps.get(cB), pB, mg, sub_nodes, k=15,
            visible_a=basin_cc_cache.get((cA, pA)),
            visible_b=basin_cc_cache.get((cB, pB)))
        shared_jaccard = n_shared / max(n_union, 1)

        # Render fresh basins for this pair with the SAME highlight set on both
        # sides so the gold rings on the left match the right.
        suffix = f"_pair{i+1}_highlighted"
        label_A = (f"paired with {cB} α-PC{pB}  (|cos|={cs:.3f}, "
                   f"{n_shared} shared / {n_union} union)")
        label_B = (f"paired with {cA} α-PC{pA}  (|cos|={cs:.3f}, "
                   f"{n_shared} shared / {n_union} union)")
        path_A = basin_mod.render_one_basin(
            state, cA, pA - 1, label_A,
            highlight_nodes=shared_ids, output_suffix=suffix,
        )
        path_B = basin_mod.render_one_basin(
            state, cB, pB - 1, label_B,
            highlight_nodes=shared_ids, output_suffix=suffix,
        )
        pair_renders.append({
            "idx": i + 1, "cA": cA, "pA": pA, "cB": cB, "pB": pB,
            "cs": cs, "same_disease": same_disease,
            "n_shared": n_shared, "n_union": n_union,
            "jaccard": shared_jaccard,
            "shared_labels": shared_labels,
            "path_A": path_A, "path_B": path_B,
        })

    # ---- PIL composite (side-by-side basins + slim text band) ----
    cell_w, cell_h = 1200, 700
    header_h      = 60
    shared_h      = 110   # slim — bulk of overlap-comm now visible IN basins
    gutter        = 8
    side_pad      = 20
    top_pad       = 110
    bottom_pad    = 20
    n_rows = len(pair_renders)
    total_w = side_pad * 2 + 2 * cell_w + gutter
    total_h = top_pad + n_rows * (cell_h + header_h + shared_h + gutter) + bottom_pad

    composite = Image.new("RGB", (total_w, total_h), color="white")
    draw = ImageDraw.Draw(composite)
    try:
        font_title  = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 28)
        font_hdr    = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 22)
        font_label  = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 18)
    except (OSError, IOError):
        font_title = font_hdr = font_label = ImageFont.load_default()

    draw.text((side_pad, 30),
              "Top cross-cohort α-PC similarities — gold halos mark substrate nodes "
              "shared between the paired basins (top-15 |loadings| intersection)",
              fill="black", font=font_title)
    draw.text((side_pad, 70),
              "Identical gold rings on the left and right basins are the visual "
              "link the shared-node list alone couldn't show.",
              fill="#666", font=font_label)

    missing = []
    row_cursor = top_pad
    for pr in pair_renders:
        i = pr["idx"]
        row_y = row_cursor
        header_color = "#10b981" if pr["same_disease"] else "#6b7280"
        header_text = (f"#{i}    |cos| = {pr['cs']:.3f}    "
                       + ("same disease ✓" if pr["same_disease"] else "cross-disease")
                       + f"    ·    shared top-15 = {pr['n_shared']}   "
                       + f"(Jaccard = {pr['jaccard']:.2f})")
        draw.rectangle([0, row_y, total_w, row_y + header_h], fill="#f3f4f6")
        draw.text((side_pad, row_y + 15), header_text,
                  fill=header_color, font=font_hdr)

        cell_y = row_y + header_h
        for col_idx, (cohort, pc, path) in enumerate(
                ((pr["cA"], pr["pA"], pr["path_A"]),
                 (pr["cB"], pr["pB"], pr["path_B"]))):
            cell_x = side_pad + col_idx * (cell_w + gutter)
            if path is not None and path.exists():
                with Image.open(path) as img:
                    img_resized = img.resize((cell_w, cell_h),
                                              Image.Resampling.LANCZOS)
                    composite.paste(img_resized.convert("RGB"),
                                    (cell_x, cell_y))
            else:
                draw.rectangle([cell_x, cell_y,
                                 cell_x + cell_w, cell_y + cell_h],
                                outline="#cccccc", width=1)
                draw.text((cell_x + cell_w // 2 - 200, cell_y + cell_h // 2),
                          f"missing render", fill="#888", font=font_label)
                missing.append(f"pair{i}_{cohort}_pc{pc}")
            draw.text((cell_x + 10, cell_y + cell_h - 28),
                      f"{cohort} · α-PC{pc}", fill="white", font=font_label,
                      stroke_width=2, stroke_fill="black")

        # Slim shared-node text band (now supplementary — primary link is the halos)
        sh_y = cell_y + cell_h
        sh_color = "#fef3c7" if pr["shared_labels"] else "#f9fafb"
        draw.rectangle([0, sh_y, total_w, sh_y + shared_h], fill=sh_color)
        draw.text((side_pad + 4, sh_y + 8),
                  f"shared top-15 nodes (gold halos, n = {pr['n_shared']}):",
                  fill="#5b3f00" if pr["shared_labels"] else "#666",
                  font=font_label)
        if pr["shared_labels"]:
            per_line = 8
            line_y = sh_y + 38
            for chunk_start in range(0, len(pr["shared_labels"]), per_line):
                chunk = pr["shared_labels"][chunk_start:chunk_start + per_line]
                draw.text((side_pad + 4, line_y),
                          "  ·  ".join(chunk),
                          fill="#2d2200", font=font_label)
                line_y += 26
                if line_y > sh_y + shared_h - 26:
                    break
        else:
            draw.text((side_pad + 4, sh_y + 50),
                      "(no overlap in top-15 — basins share the SUBSPACE not "
                      "the same top nodes)",
                      fill="#888", font=font_label)

        row_cursor = sh_y + shared_h + gutter

    out_path = FIG / "fig_top_similarity_basin_pairs.png"
    composite.save(out_path, optimize=True)
    print(f"\nWrote {out_path}  ({total_w}×{total_h} px)")
    if missing:
        print(f"\n⚠ {len(missing)} basin renders missing:")
        for m in missing:
            print(f"   {m}")


if __name__ == "__main__":
    main()
