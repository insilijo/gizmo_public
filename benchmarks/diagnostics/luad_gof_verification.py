"""Verify whether canonical substrate can/can't represent LUAD GoF biology.

Step 1: Does the raw RNA contain the KRAS-mut vs EGFR-mut MAPK-target signal?
   — If no: GoF claim fails at the data level; substrate isn't to blame.

Step 2: Does the canonical-substrate α-PCA find that signal when restricted
        to single-driver KRAS-mut vs EGFR-mut LUAD?
   — If no: substrate can't route the signal even when it exists in data.
   — If yes: substrate CAN route it; the full-cohort null is a mixture issue.

Outcome:
   Step1 PASS + Step2 FAIL  → GoF explanation supported (Paper 4 makes sense)
   Step1 PASS + Step2 PASS  → mixture/heterogeneity, not GoF; Paper 4 wrong fix
   Step1 FAIL               → signal not in data; deeper representational gap
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

REPO = Path("/home/jgardner/GIZMO")
RESULTS = REPO / "benchmarks/results"
UR = RESULTS / "unsupervised"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))

ENTREZ_TO_GENE = {
    3845: "KRAS", 1956: "EGFR", 673: "BRAF", 9817: "KEAP1", 6794: "STK11",
    4780: "NFE2L2", 7157: "TP53", 4233: "MET", 238: "ALK", 6098: "ROS1",
    5728: "PTEN", 55193: "PBRM1", 8085: "KMT2D", 5290: "PIK3CA", 5295: "PIK3R1",
}

MAPK_TARGETS = ["DUSP6", "DUSP4", "SPRY2", "SPRY4", "ETV4", "ETV5",
                 "FOS", "EGR1", "PHLDA1", "CCND1"]
EGFR_TARGETS = ["HBEGF", "AREG", "EREG", "BTC", "EPGN", "TGFA", "LRIG1"]


def load_mutations():
    raw = json.load(open(REPO / "data/cohorts/TCGA_LUAD_mutations_cbio.json"))
    rows = []
    for r in raw:
        gene = ENTREZ_TO_GENE.get(r["entrezGeneId"])
        if gene is None:
            continue
        # cBioPortal sampleId format: TCGA-XX-XXXX-01
        # F file format: tcga-xx-xxxx (lowercase, 3-segment)
        sid = r["sampleId"]
        m = re.match(r"^(TCGA-\w{2}-\w{4})", sid)
        if not m:
            continue
        pid = m.group(1).lower()
        rows.append({"patient_id": pid, "gene": gene,
                      "protein_change": r.get("proteinChange"),
                      "mutation_type": r.get("mutationType")})
    df = pd.DataFrame(rows)
    print(f"Mutation table: {len(df)} records across {df.patient_id.nunique()} patients", flush=True)
    return df


def build_mutation_status(mut_df, all_patients):
    """Per-patient binary status for each driver gene."""
    status = pd.DataFrame(index=all_patients)
    for gene in ENTREZ_TO_GENE.values():
        muts = mut_df[mut_df.gene == gene]
        # Keep only nonsynonymous coding mutations as functional
        # (cBioPortal "Missense_Mutation" / "Nonsense_Mutation" / "Frame_Shift_*" / "Splice_Site")
        funct = muts[muts.mutation_type.isin(["Missense_Mutation", "Nonsense_Mutation",
                                                "Frame_Shift_Del", "Frame_Shift_Ins",
                                                "Splice_Site", "In_Frame_Del",
                                                "In_Frame_Ins"])]
        status[gene] = status.index.isin(set(funct.patient_id)).astype(int)
    return status


def load_F_and_decompose():
    import networkx as nx
    from gizmo.export.json_export import read_json
    from per_patient_wlsp_v2 import biochem_subgraph
    print("Loading substrate…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    sub_dir, nodes, nid_idx = biochem_subgraph(mg, hub_cap=200)
    sub = sub_dir.to_undirected() if sub_dir.is_directed() else sub_dir
    pr = nx.pagerank(sub)
    log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)
    print("Loading TCGA_LUAD F…", flush=True)
    fd = np.load(UR / "stage3_F_TCGA_LUAD.npz", allow_pickle=True)
    F = fd["F"].astype(np.float64)
    pids = [str(p).lower() for p in fd["patient_ids"]]
    # Unit-norm decomposition
    F_norm = np.linalg.norm(F, axis=1, keepdims=True) + 1e-12
    F_unit = F / F_norm
    x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
    F_mean = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
    beta = (cov / x_var).ravel()
    alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
    return alpha, np.array(pids), nodes, mg, log_pr


def load_raw_rna():
    """Load TCGA_LUAD RSEM-normalized RNA matrix; columns = sample IDs (truncate to 3-seg)."""
    expr_path = (Path.home() / ".cache/tcga_luad/"
                  "gdac.broadinstitute.org_LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__"
                  "Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/"
                  "LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt")
    print(f"Loading raw RNA from {expr_path}…", flush=True)
    df = pd.read_csv(expr_path, sep="\t", skiprows=[1], low_memory=False)
    df = df.rename(columns={"Hybridization REF": "gene"})
    df["gene"] = df["gene"].astype(str).str.split("|").str[0]
    df = df[df["gene"].notna() & (df["gene"] != "?")]
    df = df.drop_duplicates(subset="gene", keep="first").set_index("gene")
    # Truncate sample IDs to 3-segment lowercase
    new_cols = []
    for c in df.columns:
        m = re.match(r"^(tcga-\w{2}-\w{4})", str(c).lower())
        new_cols.append(m.group(1) if m else str(c).lower())
    df.columns = new_cols
    # Drop columns with -10 to -19 vial codes (normals) — but only have 3-seg now,
    # so this isn't relevant. Just keep numeric.
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def main():
    print("== Loading mutation data ==", flush=True)
    mut_df = load_mutations()

    print("\n== Loading α (decomposed F) ==", flush=True)
    alpha, pids_F, nodes, mg, log_pr = load_F_and_decompose()

    print("\n== Building per-patient mutation status ==", flush=True)
    status = build_mutation_status(mut_df, pids_F)
    n_mut_total = status.sum(axis=0)
    print("Per-gene mutation counts in our F cohort (n=508):", flush=True)
    for g, c in n_mut_total.sort_values(ascending=False).items():
        print(f"  {g}: {c}", flush=True)

    # Build single-driver subsets (exclude major co-mutants of OTHER GoF drivers)
    # Major drivers: KRAS, EGFR, BRAF, MET, ALK, ROS1
    gof_drivers = ["KRAS", "EGFR", "BRAF", "MET", "ALK", "ROS1"]
    n_gof = status[gof_drivers].sum(axis=1)
    pure_kras = (status["KRAS"] == 1) & (n_gof == 1)
    pure_egfr = (status["EGFR"] == 1) & (n_gof == 1)
    print(f"\nSingle-driver subsets (no other major GoF):", flush=True)
    print(f"  KRAS-only mut: {pure_kras.sum()}", flush=True)
    print(f"  EGFR-only mut: {pure_egfr.sum()}", flush=True)

    if pure_kras.sum() < 10 or pure_egfr.sum() < 10:
        print("FATAL: subsets too small. Aborting.", flush=True)
        return

    # ========================
    # STEP 1 — raw RNA signal
    # ========================
    print("\n========== STEP 1: raw RNA differential expression ==========", flush=True)
    expr = load_raw_rna()
    # Restrict expression to patients we have
    expr_pids = [p for p in expr.columns.tolist() if p in set(pids_F)]
    print(f"  RNA matrix: {expr.shape[0]} genes × {expr.shape[1]} samples", flush=True)
    print(f"  Patients overlapping with F: {len(expr_pids)}", flush=True)

    kras_pids = [p for p in pids_F[pure_kras.values] if p in set(expr.columns)]
    egfr_pids = [p for p in pids_F[pure_egfr.values] if p in set(expr.columns)]
    print(f"  KRAS-only in RNA: {len(kras_pids)}, EGFR-only in RNA: {len(egfr_pids)}", flush=True)

    from scipy.stats import mannwhitneyu, ranksums

    def de_test(genes, label):
        results = []
        for g in genes:
            if g not in expr.index:
                continue
            k = expr.loc[g, kras_pids].dropna()
            e = expr.loc[g, egfr_pids].dropna()
            if len(k) < 5 or len(e) < 5:
                continue
            # log-transform
            k_log = np.log2(k.values.astype(float) + 1)
            e_log = np.log2(e.values.astype(float) + 1)
            u, p = mannwhitneyu(k_log, e_log, alternative="two-sided")
            results.append({"gene": g, "mean_kras_log2": float(k_log.mean()),
                             "mean_egfr_log2": float(e_log.mean()),
                             "diff_kras_minus_egfr": float(k_log.mean() - e_log.mean()),
                             "p": float(p)})
        df = pd.DataFrame(results)
        print(f"\n  {label}: (KRAS-mut - EGFR-mut, log2 RSEM):", flush=True)
        for _, r in df.iterrows():
            sign_marker = "🠕" if r.diff_kras_minus_egfr > 0 else "🠗"
            print(f"    {r.gene:<10s} ΔlogE = {r.diff_kras_minus_egfr:+.3f}  p = {r.p:.2e}", flush=True)
        return df

    mapk_df = de_test(MAPK_TARGETS, "MAPK targets (should be HIGHER in KRAS-mut)")
    egfr_df = de_test(EGFR_TARGETS, "EGFR targets (should be HIGHER in EGFR-mut)")

    # Step 1 verdict
    mapk_up_in_kras = (mapk_df.diff_kras_minus_egfr > 0).mean() if len(mapk_df) else 0.0
    egfr_up_in_egfr = (egfr_df.diff_kras_minus_egfr < 0).mean() if len(egfr_df) else 0.0
    n_mapk_sig = (mapk_df.p < 0.1).sum() if len(mapk_df) else 0
    n_egfr_sig = (egfr_df.p < 0.1).sum() if len(egfr_df) else 0
    print(f"\n  Step 1 summary:", flush=True)
    print(f"    MAPK targets higher in KRAS-mut: {mapk_up_in_kras*100:.0f}% ({n_mapk_sig}/{len(mapk_df)} at p<0.1)", flush=True)
    print(f"    EGFR targets higher in EGFR-mut: {egfr_up_in_egfr*100:.0f}% ({n_egfr_sig}/{len(egfr_df)} at p<0.1)", flush=True)

    step1_pass = (mapk_up_in_kras > 0.5 and n_mapk_sig >= 2) or (egfr_up_in_egfr > 0.5 and n_egfr_sig >= 2)
    print(f"  → Step 1 {'PASS' if step1_pass else 'FAIL'} — signal {'IS' if step1_pass else 'IS NOT'} in raw RNA", flush=True)

    # ========================
    # STEP 2 — substrate α-PCA on single-driver subset
    # ========================
    print("\n========== STEP 2: substrate α-PCA on KRAS-only vs EGFR-only subset ==========", flush=True)
    # Subset α and labels
    mask = pure_kras.values | pure_egfr.values
    alpha_sub = alpha[mask]
    pids_sub = pids_F[mask]
    label = pure_kras.values[mask].astype(int)  # 1 = KRAS-mut, 0 = EGFR-mut
    print(f"  Subset: {alpha_sub.shape[0]} patients ({label.sum()} KRAS / {(1-label).sum()} EGFR)", flush=True)

    n_pcs = 5
    pca = PCA(n_components=n_pcs, random_state=0)
    scores = pca.fit_transform(alpha_sub)
    print(f"\n  α-PC discrimination (KRAS-mut vs EGFR-mut, AUC; sign-agnostic):", flush=True)
    auc_results = []
    for k in range(n_pcs):
        s = scores[:, k]
        try:
            auc = roc_auc_score(label, s)
            auc_flip = roc_auc_score(label, -s)
            best_auc = max(auc, auc_flip)
            polarity = "+→KRAS" if auc >= auc_flip else "+→EGFR"
            print(f"    α-PC{k+1}  AUC = {best_auc:.3f}  ({polarity})  "
                  f"variance explained = {pca.explained_variance_ratio_[k]:.3f}", flush=True)
            auc_results.append({"pc": k+1, "auc": best_auc,
                                 "polarity": polarity,
                                 "var_explained": float(pca.explained_variance_ratio_[k])})
        except Exception as e:
            print(f"    α-PC{k+1}  AUC failed: {e}", flush=True)

    best_pc_auc = max(r["auc"] for r in auc_results)
    print(f"\n  Best α-PC AUC: {best_pc_auc:.3f}", flush=True)
    step2_pass = best_pc_auc >= 0.70
    print(f"  → Step 2 {'PASS' if step2_pass else 'FAIL'} — substrate {'CAN' if step2_pass else 'CANNOT'} route the signal", flush=True)

    # ========================
    # VERDICT
    # ========================
    print("\n========== VERDICT ==========", flush=True)
    if step1_pass and not step2_pass:
        verdict = "GoF explanation SUPPORTED — signal in data, substrate can't route → Paper 4 (edge injection) is the right fix"
    elif step1_pass and step2_pass:
        verdict = "GoF explanation FALSIFIED — substrate can route at subset scale; full-cohort null is mixture-heterogeneity, not GoF"
    elif not step1_pass:
        verdict = "GoF explanation FALSIFIED — signal not in raw RNA at the cohort/sample size; deeper representational issue or modality-insufficient"
    else:
        verdict = "Indeterminate"
    print(f"  {verdict}", flush=True)

    # Save
    out = {
        "step1_signal_in_rna": {
            "mapk_df": mapk_df.to_dict(orient="records") if len(mapk_df) else [],
            "egfr_df": egfr_df.to_dict(orient="records") if len(egfr_df) else [],
            "pass": bool(step1_pass),
        },
        "step2_substrate_routing": {
            "n_kras": int(label.sum()), "n_egfr": int((1-label).sum()),
            "per_pc_auc": auc_results,
            "best_auc": best_pc_auc,
            "pass": bool(step2_pass),
        },
        "verdict": verdict,
    }
    out_path = RESULTS / "luad_gof_verification.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
