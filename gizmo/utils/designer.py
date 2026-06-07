"""Method designer / sample-size planner.

Translates the framework's empirical findings into actionable
recommendations for a planned cohort study. Given:
  - planned cohort size + class balance
  - available omics + panel sizes
  - tissue / sample-site
  - analysis goal
returns:
  - feasible analyses (group σ / discovery / per-patient classification / response prediction)
  - recommended evidence model (welch / mod_t / plsda / diablo)
  - recommended propagation method (heat / signed Laplacian)
  - recommended kernel (WL / SP) when applicable
  - expected AUC range based on similar cohorts in our 6-cohort benchmark
  - caveats and alternative-method suggestions

Empirical decision rules derived from PAPER_OUTLINE.md results
across IDH, Crohn, Su COVID, RA-vs-OA, Statin, Erawijantari.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal, Optional


Goal = Literal["discovery", "group_sigma", "per_patient_classification",
                "response_prediction"]
Tissue = Literal["disease_site", "peripheral", "in_vitro"]
Regime = Literal["A_focal", "B_multi_faceted", "C_divergent", "unknown"]


@dataclass
class CohortSpec:
    n_total: int
    n_case: int
    n_control: int
    omics: list[str]                # subset of {"trans", "prot", "metab"}
    panel_sizes: dict[str, int]     # e.g., {"trans": 22000, "metab": 1050}
    tissue: Tissue
    goal: Goal
    regime: Regime = "unknown"
    paired_omics: bool = False      # are samples paired across omics?
    longitudinal: bool = False      # paired pre/post visits available?


@dataclass
class Recommendation:
    feasible: dict[str, bool]
    evidence_model: str
    propagation: str
    kernel: Optional[str]
    expected_loo_auc: tuple[float, float]
    caveats: list[str]
    alternatives: list[str]
    rationale: list[str]


# ---------------------------------------------------------------------------
# Decision rules (empirical, derived from 6-cohort benchmark)
# ---------------------------------------------------------------------------

def _class_imbalance_ratio(spec: CohortSpec) -> float:
    """Returns max(n_case, n_control) / min(...) — higher = more imbalanced."""
    if min(spec.n_case, spec.n_control) == 0:
        return float("inf")
    return max(spec.n_case, spec.n_control) / min(spec.n_case, spec.n_control)


def _has_dense_metab(spec: CohortSpec) -> bool:
    return "metab" in spec.omics and spec.panel_sizes.get("metab", 0) >= 500


def _has_sparse_prot_trans(spec: CohortSpec) -> bool:
    for omic in ("prot", "trans"):
        if omic in spec.omics and spec.panel_sizes.get(omic, 0) < 1000:
            return True
    return False


def recommend(spec: CohortSpec) -> Recommendation:
    feasible = {
        "discovery": True,                                  # always feasible
        "group_sigma": min(spec.n_case, spec.n_control) >= 5,
        "per_patient_classification": False,
        "response_prediction": False,
    }
    caveats: list[str] = []
    rationale: list[str] = []
    alternatives: list[str] = []

    n_min_arm = min(spec.n_case, spec.n_control)
    imbalance = _class_imbalance_ratio(spec)

    # --- Per-patient classification feasibility ---
    if spec.n_total >= 100 and n_min_arm >= 10:
        feasible["per_patient_classification"] = True
        rationale.append(
            f"n_total={spec.n_total} ≥ 100 and minority arm ≥ 10 → "
            "kernel-SVM stability achievable"
        )
    elif spec.n_total >= 50:
        feasible["per_patient_classification"] = True
        caveats.append(
            f"n_total={spec.n_total} (50-100) is marginal for kernel-SVM; "
            "expect high seed variance and wide hold-out CIs"
        )
    else:
        caveats.append(
            f"n_total={spec.n_total} < 50: kernel-SVM unstable (Crohn n=33, "
            "Gao RA n=28 both gave random-or-flipped LOO AUC). Per-patient "
            "classification unreliable; group-level σ analysis still works."
        )

    # --- Class imbalance warnings ---
    if imbalance >= 5:
        caveats.append(
            f"Class ratio {imbalance:.1f}:1 — stratified CV / hold-out give "
            "flipped predictions with class_weight='balanced'. Use LOO with "
            "pooled predictions (most reliable at imbalance ≥ 5:1)."
        )
    if min(spec.n_case, spec.n_control) < 10:
        caveats.append(
            f"Minority arm = {n_min_arm} samples — too few for stratified "
            "5-fold CV (folds get 2-3 minority samples). Expect kernel "
            "collapse on minority class."
        )

    # --- Response prediction (longitudinal Δ) ---
    if spec.longitudinal and feasible["per_patient_classification"]:
        feasible["response_prediction"] = True
        rationale.append("Longitudinal pairing enables V_post − V_pre Δ as "
                          "evidence for response prediction.")
        if "metab" not in spec.omics:
            caveats.append(
                "Longitudinal trans-only Δ alone failed on CorEvitas RA "
                "(n=159, Good vs Poor random AUC). Response prediction "
                "may need disease-site metabolomics (per RA cohort-transfer "
                "finding).")

    # --- Evidence model selection ---
    if spec.regime == "C_divergent" or (spec.paired_omics and _has_dense_metab(spec)):
        evidence_model = "diablo"
        rationale.append(
            "DIABLO multi-block sPLS-DA recommended for cross-block coupling "
            "(empirical: Su COVID both_directional mean_pert -1.2σ → +1.7σ "
            "with DIABLO, +2.9σ lift)."
        )
    elif n_min_arm >= 5:
        evidence_model = "moderated_t"
        rationale.append(
            "Moderated-t (Smyth 2004) variance shrinkage stable at small n; "
            "default for Regime A/B cohorts."
        )
    else:
        evidence_model = "welch_t"
        caveats.append(
            f"n_min_arm={n_min_arm} < 5 — moderated-t shrinkage unstable; "
            "Welch's t-test + log2_fc as effect_size."
        )

    # --- Propagation method (density-driven) ---
    has_dense_metab = _has_dense_metab(spec)
    has_sparse = _has_sparse_prot_trans(spec)
    if has_dense_metab and not has_sparse:
        propagation = "signed_laplacian"
        rationale.append(
            "Signed Laplacian for dense metabolomics (≥500 features) — "
            "preserves mass-action direction. Empirical: Su COVID node-only "
            "SP+SVC LOO AUC 0.957 with signed Lap; heat diffusion failed at "
            "0.31 (smears direction)."
        )
    elif has_sparse and not has_dense_metab:
        propagation = "heat_diffusion"
        rationale.append(
            "Heat diffusion for sparse trans/prot evidence — broadcasts "
            "magnitude to fill density gaps. Empirical: Su COVID edge prot "
            "heat+SP LOO AUC 1.000; signed Lap+SP failed at 0.42."
        )
    else:
        propagation = "both_for_comparison"
        rationale.append(
            "Mixed density profile — run both signed Laplacian and heat "
            "diffusion; method-disagreement reactions form a paradigm-"
            "specific subnetwork (see PAPER_OUTLINE.md)."
        )

    # --- Kernel choice ---
    kernel: Optional[str] = None
    if feasible["per_patient_classification"]:
        kernel = "SP"
        rationale.append(
            "Shortest-Path kernel for sparse multi-module submodular "
            "subgraphs (median 14 modules per patient). WL kernel collapses "
            "on disconnected components."
        )

    # --- Tissue / disease-site warning ---
    if spec.tissue == "peripheral" and "metab" not in spec.omics:
        caveats.append(
            "Peripheral transcriptomics alone is far from disease site for "
            "many diseases (e.g., RA whole-blood failed; synovial metab "
            "succeeded at +3.4σ in same framework). Consider adding "
            "metabolomics or moving to disease-site biopsy."
        )

    if spec.tissue == "in_vitro":
        caveats.append(
            "In vitro pharmacologic experiments (e.g., Statin GSE57071, "
            "n~12) typically have very small n. Group-level σ feasible; "
            "per-patient kernel SVM not."
        )

    # --- Expected LOO AUC range ---
    auc_low, auc_high = 0.50, 0.65   # default: random
    if feasible["per_patient_classification"]:
        if has_dense_metab and spec.paired_omics and spec.n_total >= 200:
            auc_low, auc_high = 0.85, 0.97
        elif has_sparse and spec.n_total >= 200:
            auc_low, auc_high = 0.85, 1.00
        elif spec.n_total >= 100:
            auc_low, auc_high = 0.65, 0.85
        else:
            auc_low, auc_high = 0.55, 0.75
    if imbalance >= 10:
        auc_low *= 0.85   # imbalance penalty

    # --- Alternative methods to compare ---
    if "trans" in spec.omics or "prot" in spec.omics:
        alternatives.append(
            "ssGSEA + SVC: pathway-level baseline. Empirically tied with "
            "graph methods on Su COVID prot (LOO 1.000); cannot handle "
            "metab features — gene-set databases lack metabolite pathways."
        )
    if spec.paired_omics:
        alternatives.append(
            "MOFA+: matrix-factorization joint embedding. Replaces graph "
            "with linear factor decomposition. Useful comparison baseline."
        )
    if has_sparse:
        alternatives.append(
            "Heat diffusion + WL/SP graph kernel: same downstream pipeline, "
            "unsigned propagation. Use for sparse evidence."
        )

    return Recommendation(
        feasible=feasible,
        evidence_model=evidence_model,
        propagation=propagation,
        kernel=kernel,
        expected_loo_auc=(auc_low, auc_high),
        caveats=caveats,
        alternatives=alternatives,
        rationale=rationale,
    )


def format_recommendation(rec: Recommendation) -> str:
    """Pretty-printed recommendation suitable for paper Methods or planning docs."""
    lines = ["### GIZMO method-design recommendation", ""]
    lines.append("**Feasible analyses:**")
    for k, v in rec.feasible.items():
        mark = "✓" if v else "✗"
        lines.append(f"  - {mark} {k}")
    lines.append("")
    lines.append(f"**Recommended evidence model**: `{rec.evidence_model}`")
    lines.append(f"**Recommended propagation**: `{rec.propagation}`")
    if rec.kernel:
        lines.append(f"**Recommended graph kernel**: `{rec.kernel}` "
                       "(submodular subgraph + class-balanced SVC)")
    lo, hi = rec.expected_loo_auc
    lines.append(f"**Expected LOO AUC range**: {lo:.2f} – {hi:.2f}")
    lines.append("")
    if rec.rationale:
        lines.append("**Rationale**:")
        for r in rec.rationale:
            lines.append(f"  - {r}")
        lines.append("")
    if rec.caveats:
        lines.append("**Caveats**:")
        for c in rec.caveats:
            lines.append(f"  - ⚠ {c}")
        lines.append("")
    if rec.alternatives:
        lines.append("**Alternative methods to compare against**:")
        for a in rec.alternatives:
            lines.append(f"  - {a}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    examples = [
        ("Crohn-like (small paired multi-omic)", CohortSpec(
            n_total=33, n_case=16, n_control=17,
            omics=["prot", "metab"],
            panel_sizes={"prot": 140, "metab": 50},
            tissue="peripheral", goal="per_patient_classification",
            paired_omics=True,
        )),
        ("Su-like (large paired plasma multi-omic)", CohortSpec(
            n_total=270, n_case=252, n_control=18,
            omics=["prot", "metab"],
            panel_sizes={"prot": 382, "metab": 1050},
            tissue="peripheral", goal="per_patient_classification",
            paired_omics=True,
        )),
        ("Hypothetical: large balanced metab-only", CohortSpec(
            n_total=500, n_case=250, n_control=250,
            omics=["metab"],
            panel_sizes={"metab": 800},
            tissue="disease_site", goal="per_patient_classification",
        )),
        ("RA-CorEvitas-like (longitudinal trans-only)", CohortSpec(
            n_total=159, n_case=99, n_control=60,
            omics=["trans"],
            panel_sizes={"trans": 22000},
            tissue="peripheral", goal="response_prediction",
            longitudinal=True,
        )),
    ]
    for name, spec in examples:
        print("=" * 70)
        print(f"## {name}")
        print(f"   spec: n={spec.n_total} ({spec.n_case}/{spec.n_control}), "
               f"omics={spec.omics}, tissue={spec.tissue}, goal={spec.goal}")
        print()
        rec = recommend(spec)
        print(format_recommendation(rec))
        print()
