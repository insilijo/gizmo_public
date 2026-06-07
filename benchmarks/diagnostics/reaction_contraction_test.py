"""Test whether reaction nodes are load-bearing for MAP integration.

Hypothesis being tested: reactions are reified-edge nodes that integrate
signal from {catalyzing gene, substrate metabolite, product metabolite}.
Counter-hypothesis: reactions are 2-step pass-throughs that just attenuate
substrate↔product propagation without adding integration.

Test design
-----------
For ONE cohort (IDH_glioma — has paired NMR + RNA, clean α-PC1 basin biology):

  1. Build STANDARD substrate: reactions as nodes, standard topology.
  2. Build CONTRACTED substrate: remove every reaction R; for its (G_R, S_R, P_R) =
     (catalyzing genes, substrate metabolites, product metabolites), add edges
       S × P  (substrate↔product direct, 1-hop instead of 2 via R)
       G × S  (gene↔substrate direct)
       G × P  (gene↔product direct)
     so the integration that used to happen AT R is now distributed across edges.
  3. Run MAP on both substrates. Compute α/β decomposition, α-PCA.
  4. Compare:
       • Top-50 loadings overlap (Jaccard) for α-PC1
       • Signed-basin gene/metabolite sets (Jaccard)
       • |F| Spearman rho between standard and contracted at gene+metab nodes
       • α-PC1 score AUC against the cohort's binary label

Verdict
-------
H1 (integration): contracted basin biology degrades; AUC drops; top-K Jaccard low.
H2 (pass-through): contracted basin biology preserved; AUC ≈ same; high Jaccard.

Output: reaction_contraction_test_<cohort>.tsv + side-by-side basin PNG.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, cg, expm_multiply
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

REPO = Path("/home/jgardner/GIZMO")
RESULTS = REPO / "benchmarks/results"
UR = RESULTS / "unsupervised"
FIG = RESULTS / "figures"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))

_PROPAGATION_TYPES = frozenset({"reaction", "metabolite", "gene"})
COHORT = "IDH_glioma"   # paired NMR + RNA, clean α-PC1 basin biology


def standard_subgraph(mg, hub_cap: int = 200):
    g = mg.graph
    keep = [n for n, a in g.nodes(data=True)
            if a.get("node_type") in _PROPAGATION_TYPES
            and not a.get("is_currency")
            and g.degree(n) <= hub_cap]
    sub = g.subgraph(keep).copy()
    if sub.is_directed():
        sub = sub.to_undirected()
    nodes = sorted(sub.nodes())
    nid_idx = {n: i for i, n in enumerate(nodes)}
    return sub, nodes, nid_idx


def contracted_subgraph(mg, hub_cap: int = 200):
    """Build a substrate with reactions removed + bypass edges added.

    Strategy: take the standard subgraph, then for each reaction R, find
    {catalysts, substrates, products}, add bypass edges (S↔P, G↔S, G↔P),
    remove R. Result: same gene/metabolite node set, no reaction nodes,
    densified edges around the formerly-reaction-mediated transformations.
    """
    sub, _, _ = standard_subgraph(mg, hub_cap=hub_cap)
    sub = sub.copy()

    reaction_nodes = [n for n, a in sub.nodes(data=True)
                       if a.get("node_type") == "reaction"]
    new_edges = []
    for R in reaction_nodes:
        neighbors = list(sub.neighbors(R))
        catalysts = []
        substrates = []
        products = []
        for nbr in neighbors:
            nbr_type = sub.nodes[nbr].get("node_type", "")
            if nbr_type == "gene":
                catalysts.append(nbr)
            elif nbr_type == "metabolite":
                # Use edge data role if present, else treat as both substrate+product
                edata = sub.get_edge_data(R, nbr, default={})
                role = (edata.get("role") or "").lower()
                if "substrate" in role:
                    substrates.append(nbr)
                elif "product" in role:
                    products.append(nbr)
                else:
                    substrates.append(nbr); products.append(nbr)
        # Bypass edges
        for s in substrates:
            for p in products:
                if s != p:
                    new_edges.append((s, p))
        for g in catalysts:
            for s in substrates: new_edges.append((g, s))
            for p in products:   new_edges.append((g, p))
        sub.remove_node(R)
    # Add new edges (deduped)
    sub.add_edges_from(set(new_edges))
    nodes = sorted(sub.nodes())
    nid_idx = {n: i for i, n in enumerate(nodes)}
    print(f"  Contracted substrate: removed {len(reaction_nodes)} reactions; "
          f"added {len(set(new_edges))} bypass edges. "
          f"Remaining nodes: {len(nodes)}", flush=True)
    return sub, nodes, nid_idx


def laplacian_from(sub, nodes):
    A = nx.adjacency_matrix(sub, nodelist=nodes).astype(float)
    deg = np.asarray(A.sum(axis=1)).ravel()
    L = (sp.diags(deg) - A).tocsr()
    # Compute PageRank ONCE (not per node — that was 38k× slower)
    pr = nx.pagerank(sub)
    log_pr = np.log10(
        np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15
    )
    return L, log_pr


def load_cohort_data(mg):
    """Return (modality_setups, patient_ids, labels) for IDH_glioma."""
    from per_patient_master import load_idh_glioma
    from gizmo.evidence.mappers import GeneMapper, MetaboliteMapper
    prot, metab, y, common = load_idh_glioma()
    patient_ids = list(common)
    labels = [y.get(s, "unknown") for s in patient_ids]
    label_bin = np.array([1 if l == "active" else (0 if l == "control" else -1)
                          for l in labels])
    gmap = GeneMapper(mg); mmap = MetaboliteMapper(mg)

    def setup(data, mapper, label):
        if not data: return None
        feats = sorted(set().union(*[set(data[s]) for s in patient_ids if s in data]))
        feat_to_node = {}
        for f in feats:
            try:
                res = mapper.map(f)
                nid = res[0] if isinstance(res, tuple) else res
            except Exception: nid = None
            if nid: feat_to_node[f] = nid
        valid_feats = list(feat_to_node.keys())
        if not valid_feats: return None
        X = np.zeros((len(patient_ids), len(valid_feats)))
        for i, s in enumerate(patient_ids):
            if s not in data: continue
            for j, f in enumerate(valid_feats):
                X[i, j] = data[s].get(f, 0.0)
        X_log = np.sign(X) * np.log1p(np.abs(X) + 1e-12)
        gs = float(np.nanstd(X_log)) or 1.0
        x_norm = np.nan_to_num(X_log / gs, nan=0.0)
        sigma = float(max(np.std(x_norm[x_norm != 0]) * 0.25, 0.05))
        zdata = {patient_ids[i]: {valid_feats[j]: float(x_norm[i, j])
                                    for j in range(len(valid_feats))}
                 for i in range(len(patient_ids))}
        return {"label": label, "sigma": sigma, "t": 0.5,
                "feat_to_node": feat_to_node, "data": zdata}

    setups = []
    for data, mapper, lbl in [(prot, gmap, "prot"), (metab, mmap, "metab")]:
        s = setup(data, mapper, lbl)
        if s is not None: setups.append(s)
    return setups, patient_ids, label_bin


def solve_map(L, nid_idx, n_nodes, modality_setups, patient_ids,
              lambda_smooth=None, lambda_calibration=0.1, ridge_alpha=1e-6):
    if lambda_smooth is None:
        mean_dw = float(np.mean([1.0 / (s["sigma"] ** 2) for s in modality_setups]))
        lambda_smooth = lambda_calibration * mean_dw / (0.05 * n_nodes)
    print(f"  λ_smooth = {lambda_smooth:.5f}", flush=True)
    base_M = lambda_smooth * L + sp.diags(np.full(n_nodes, ridge_alpha))
    F = np.zeros((len(patient_ids), n_nodes), dtype=np.float32)
    for i, sid in enumerate(patient_ids):
        rhs = np.zeros(n_nodes); per_anchor = np.zeros(n_nodes)
        for s in modality_setups:
            if sid not in s["data"]: continue
            obs = s["data"][sid]
            sigma_inv = 1.0 / (s["sigma"] ** 2)
            x_sparse = np.zeros(n_nodes)
            for f, nid in s["feat_to_node"].items():
                if nid in nid_idx and f in obs:
                    x_sparse[nid_idx[nid]] += float(obs[f])
                    per_anchor[nid_idx[nid]] += sigma_inv
            x_smoothed = expm_multiply(-s["t"] * L, x_sparse) if s["t"] > 0 else x_sparse
            rhs += x_smoothed * sigma_inv
        M_i = (base_M + sp.diags(per_anchor)).tocsr()
        M_diag = M_i.diagonal()
        M_inv = 1.0 / np.where(M_diag > 1e-12, M_diag, 1.0)
        Minv_op = LinearOperator(M_i.shape, matvec=lambda v, w=M_inv: v * w,
                                  dtype=M_i.dtype)
        try:
            f_sol, _ = cg(M_i, rhs, M=Minv_op, rtol=1e-8, maxiter=2000)
        except Exception:
            f_sol = np.zeros(n_nodes)
        F[i] = f_sol.astype(np.float32)
        if (i + 1) % 20 == 0:
            print(f"    patient {i+1}/{len(patient_ids)}", flush=True)
    return F


def alpha_pca(F, log_pr, n_components=5):
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    pca = PCA(n_components=n_components, random_state=0)
    scores = pca.fit_transform(alpha)
    return beta, alpha, scores, pca


def top_k_at_nodes(F, nodes, mg, k=50):
    mean_abs_F = np.abs(F).mean(axis=0)
    top_idx = np.argsort(-mean_abs_F)[:k]
    return [(nodes[i], mg.graph.nodes.get(nodes[i], {}).get("symbol")
              or mg.graph.nodes.get(nodes[i], {}).get("name")
              or nodes[i], float(mean_abs_F[i])) for i in top_idx]


def main():
    from gizmo.export.json_export import read_json
    print(f"Loading substrate…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")

    # ── STANDARD substrate ────────────────────────────────────────────────
    print("\n=== Standard substrate (reactions as nodes) ===", flush=True)
    sub_std, nodes_std, idx_std = standard_subgraph(mg, hub_cap=200)
    print(f"  {len(nodes_std)} nodes, {sub_std.number_of_edges()} edges", flush=True)
    L_std, log_pr_std = laplacian_from(sub_std, nodes_std)

    setups, patient_ids, label_bin = load_cohort_data(mg)
    print(f"  Cohort {COHORT}: {len(patient_ids)} patients, {len(setups)} modalities", flush=True)
    F_std = solve_map(L_std, idx_std, len(nodes_std), setups, patient_ids)
    beta_std, alpha_std, scores_std, pca_std = alpha_pca(F_std, log_pr_std)
    top_std = top_k_at_nodes(F_std, nodes_std, mg, k=50)

    # ── CONTRACTED substrate ──────────────────────────────────────────────
    print("\n=== Contracted substrate (reactions collapsed to bypass edges) ===",
          flush=True)
    sub_con, nodes_con, idx_con = contracted_subgraph(mg, hub_cap=200)
    print(f"  {len(nodes_con)} nodes, {sub_con.number_of_edges()} edges", flush=True)
    L_con, log_pr_con = laplacian_from(sub_con, nodes_con)
    F_con = solve_map(L_con, idx_con, len(nodes_con), setups, patient_ids)
    beta_con, alpha_con, scores_con, pca_con = alpha_pca(F_con, log_pr_con)
    top_con = top_k_at_nodes(F_con, nodes_con, mg, k=50)

    # ── Comparison ────────────────────────────────────────────────────────
    print("\n=== COMPARISON ===", flush=True)
    set_std = {nid for nid, _, _ in top_std}
    set_con = {nid for nid, _, _ in top_con}
    # Restrict to common node set (excluding reactions)
    common_nodes = set(nodes_std) & set(nodes_con)
    set_std_gm = set_std & common_nodes
    set_con_gm = set_con & common_nodes
    jacc_top50 = len(set_std & set_con) / max(len(set_std | set_con), 1)
    jacc_top50_gm = (len(set_std_gm & set_con_gm) /
                     max(len(set_std_gm | set_con_gm), 1)
                     if set_std_gm | set_con_gm else 0.0)
    print(f"  Top-50 Jaccard (all node types):           {jacc_top50:.3f}")
    print(f"  Top-50 Jaccard (gene+metab only, common):  {jacc_top50_gm:.3f}")

    # |F| Spearman ρ at gene+metab nodes (common to both substrates)
    from scipy.stats import spearmanr
    common_list = sorted(common_nodes)
    f_std_vec = np.abs(F_std).mean(axis=0)[[idx_std[n] for n in common_list]]
    f_con_vec = np.abs(F_con).mean(axis=0)[[idx_con[n] for n in common_list]]
    rho, p = spearmanr(f_std_vec, f_con_vec)
    print(f"  |F̄| Spearman ρ (gene+metab nodes, n={len(common_list)}): "
          f"{rho:.3f} (p={p:.2e})")

    # α-PC1 AUC vs label (binary)
    mask = label_bin >= 0
    if mask.sum() >= 6 and len(np.unique(label_bin[mask])) >= 2:
        try:
            auc_std = roc_auc_score(label_bin[mask], scores_std[mask, 0])
            auc_con = roc_auc_score(label_bin[mask], scores_con[mask, 0])
            auc_std = max(auc_std, 1 - auc_std)
            auc_con = max(auc_con, 1 - auc_con)
            print(f"  α-PC1 AUC vs label:  standard = {auc_std:.3f}   "
                  f"contracted = {auc_con:.3f}   Δ = {auc_con - auc_std:+.3f}")
        except Exception as e:
            print(f"  α-PC1 AUC failed: {e}")

    # Top-10 hits side-by-side
    print(f"\n  Top-10 |F̄| nodes side-by-side:")
    print(f"  {'STANDARD (with reactions)':50}  {'CONTRACTED (reactions removed)':50}")
    for ((n1, l1, v1), (n2, l2, v2)) in zip(top_std[:15], top_con[:15]):
        t1 = mg.graph.nodes.get(n1, {}).get("node_type", "?")[0]
        t2 = mg.graph.nodes.get(n2, {}).get("node_type", "?")[0]
        s1 = f"{v1:.4f}  {l1} [{t1}]"[:48]
        s2 = f"{v2:.4f}  {l2} [{t2}]"[:48]
        print(f"  {s1:50}  {s2:50}")

    # Save summary TSV
    import pandas as pd
    df_out = pd.DataFrame([{
        "cohort": COHORT, "n_patients": len(patient_ids),
        "n_nodes_std": len(nodes_std), "n_nodes_con": len(nodes_con),
        "n_edges_std": sub_std.number_of_edges(),
        "n_edges_con": sub_con.number_of_edges(),
        "top50_jaccard_all": jacc_top50,
        "top50_jaccard_gene_metab_common": jacc_top50_gm,
        "F_spearman_gene_metab": float(rho),
        "auc_std": float(auc_std) if 'auc_std' in dir() else np.nan,
        "auc_con": float(auc_con) if 'auc_con' in dir() else np.nan,
    }])
    out_tsv = RESULTS / f"reaction_contraction_test_{COHORT}.tsv"
    df_out.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nWrote {out_tsv}")

    # Save F matrices for further analysis if needed
    np.savez_compressed(UR / f"stage3_F_{COHORT}_reaction_contracted.npz",
                        F=F_con, patient_ids=patient_ids,
                        nodes=np.array(nodes_con, dtype=object))


if __name__ == "__main__":
    main()
