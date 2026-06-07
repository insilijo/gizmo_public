"""
GIZMO ↔ GrAndMA integration helpers.

Converts GrAndMA's ``AnalysisResult.table_data`` format directly into a
``SampleContext`` and runs reaction scoring — no CSV files, no CLI, no Dash
dependency required.

GrAndMA table_data format
-------------------------
T-test / ANOVA results are stored as::

    {
        "columns": ["feature", "log2_fc", "p_value", "p_adj", ...],
        "rows":    [["ATP",     1.23,      0.001,     0.01,  ...], ...]
    }

GIZMO maps:
    feature  → FEATURE_ID  (metabolite name, ChEBI ID, InChIKey, gene symbol…)
    log2_fc  → EFFECT_SIZE  (also accepts "fold_change", "log2fc", "effect_size")
    p_value  → P_VALUE      (also "pvalue", "p.value")
    p_adj    → FDR          (also "fdr", "q_value", "padj", "p_adj")

Usage (inside a Django view or background task)::

    from gizmo.integration import score_from_table, build_context_from_table

    # Quick path — score in one call
    scores = score_from_table(mg, result.table_data, assay_type="metabolomics")

    # Fine-grained path — inspect / merge multiple analyses
    ctx = build_context_from_table(mg, result.table_data,
                                   assay_type="metabolomics",
                                   sample_id=str(result.pk))
    ctx2 = build_context_from_table(mg, rna_result.table_data,
                                    assay_type="transcriptomics",
                                    sample_id=str(rna_result.pk))
    ctx.merge(ctx2)
    scores = score_reactions(mg, ctx)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from gizmo.analysis.currency import (
    compute_conditional_currency_edges, flag_currency_metabolites,
)
from gizmo.evidence.mappers import GeneMapper, MetaboliteMapper
from gizmo.evidence.model import EvidenceRecord, SampleContext
from gizmo.scoring.reaction_scorer import ReactionScore, score_reactions

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Graph-state setup: ensures graph has gene nodes, currency flags, and
# conditional-currency edge annotations so BP / heuristic scorers can run
# without each caller re-flagging from scratch. Idempotent and cached on
# the graph itself via a sentinel attribute.
# ---------------------------------------------------------------------------

def ensure_scoring_ready(mg) -> set[tuple[str, str]]:
    """
    Idempotent setup:
      1. Ensure reaction-derived gene nodes are present  (enrich once).
      2. Hard-currency-flag metabolite nodes.
      3. Compute (and cache) conditional-currency edges to skip in BP.

    Returns the set of conditional-currency edges. The caller passes this
    to ``BPConfig(skip_edges=...)`` (the heuristic scorer ignores it).

    Cached on the graph as ``mg.graph.graph["_gizmo_scoring_ready"]`` so
    subsequent calls within the same process are no-ops.
    """
    g = mg.graph
    state = g.graph.get("_gizmo_scoring_ready")
    if state is not None:
        return state

    # 1. Gene-node enrichment (idempotent; emits 0 new nodes if already done)
    try:
        from gizmo.sources.gene_enrichment import enrich_graph_genes_from_reactions
        counts = enrich_graph_genes_from_reactions(mg)
        if counts.get("new_genes", 0):
            log.info("integration: enriched +%d gene nodes, +%d edges",
                     counts.get("new_genes", 0), counts.get("new_edges", 0))
    except ImportError:
        log.debug("integration: gene_enrichment not available; skipping")

    # 2. Hard currency flagging (no borderline; conditional edges handle those)
    flag_currency_metabolites(mg, degree_threshold_k=None,
                               include_borderline=False)

    # 3. Conditional-currency edges (α-KG / SAM / Glu / Gln / etc. as
    #    cofactors in HDM / KMT / transaminase / α-KG-dioxygenase contexts,
    #    preserved as features in TCA / Glu metab / methionine cycle).
    skip_edges = compute_conditional_currency_edges(mg)
    log.info("integration: computed %d conditional-currency skip edges",
             len(skip_edges))

    g.graph["_gizmo_scoring_ready"] = skip_edges
    return skip_edges

# ---------------------------------------------------------------------------
# Column alias resolution
# ---------------------------------------------------------------------------

_EFFECT_ALIASES = {"log2_fc", "log2fc", "fold_change", "effect_size", "log2foldchange", "lfc"}
_PVALUE_ALIASES = {"p_value", "pvalue", "p.value", "pval"}
_FDR_ALIASES    = {"p_adj", "padj", "fdr", "q_value", "qvalue", "adjusted_p"}


def _resolve_columns(columns: list[str]) -> dict[str, int]:
    """
    Return a dict mapping logical field → column index for the columns present.

    Raises ValueError if 'feature' or an effect-size column is not found.
    """
    lower = [c.lower() for c in columns]
    idx: dict[str, int] = {}

    # feature (required)
    for i, c in enumerate(lower):
        if c in ("feature", "feature_id", "biochemical", "gene", "gene_id", "metabolite",
                 "hmdb", "hmdb_id", "compound", "compound_id"):
            idx["feature"] = i
            break
    if "feature" not in idx:
        raise ValueError(
            f"No feature-ID column found in {columns}. "
            "Expected one of: feature, feature_id, biochemical, gene."
        )

    # effect size (required)
    for i, c in enumerate(lower):
        if c in _EFFECT_ALIASES:
            idx["effect_size"] = i
            break
    if "effect_size" not in idx:
        raise ValueError(
            f"No effect-size column found in {columns}. "
            f"Expected one of: {sorted(_EFFECT_ALIASES)}."
        )

    # p-value (optional)
    for i, c in enumerate(lower):
        if c in _PVALUE_ALIASES:
            idx["p_value"] = i
            break

    # FDR (optional)
    for i, c in enumerate(lower):
        if c in _FDR_ALIASES:
            idx["fdr"] = i
            break

    return idx


# ---------------------------------------------------------------------------
# Context builder
# ---------------------------------------------------------------------------

def build_context_from_table(
    mg,
    table_data: dict[str, Any],
    assay_type: str = "metabolomics",
    sample_id: str = "",
    cohort_id: str = "",
    p_value_threshold: float = 1.0,
    fdr_threshold: float = 1.0,
) -> SampleContext:
    """
    Convert a GrAndMA ``table_data`` dict into a ``SampleContext``.

    Parameters
    ----------
    mg            : GizmoGraph
    table_data    : ``{"columns": [...], "rows": [...]}``
    assay_type    : "metabolomics" | "transcriptomics" | "proteomics"
    sample_id     : optional label (e.g. AnalysisResult pk as string)
    cohort_id     : optional cohort label
    p_value_threshold : rows with p_value > this are skipped (default: keep all)
    fdr_threshold     : rows with fdr > this are skipped (default: keep all)

    Returns
    -------
    SampleContext populated with mapped EvidenceRecords.
    """
    columns: list[str] = table_data.get("columns", [])
    rows:    list[list] = table_data.get("rows",    [])

    if not columns:
        raise ValueError("table_data['columns'] is empty")

    col_idx = _resolve_columns(columns)
    fi_i  = col_idx["feature"]
    eff_i = col_idx["effect_size"]
    pv_i  = col_idx.get("p_value")
    fdr_i = col_idx.get("fdr")

    # Dispatch by assay type. Metabolite-class assays (anything observing
    # small-molecule abundance) route to MetaboliteMapper; gene-class
    # assays (mRNA, protein, phospho-protein abundance) route to
    # GeneMapper. Per-assay observation noise is handled separately in
    # gizmo.inference.model.DEFAULT_ASSAY_SIGMAS.
    METABOLITE_ASSAYS = {"metabolomics", "lipidomics", "small_molecule"}
    GENE_ASSAYS = {"transcriptomics", "proteomics", "phosphoproteomics",
                    "translatomics", "ribosomics", "expression"}
    mapper: MetaboliteMapper | GeneMapper
    if assay_type in METABOLITE_ASSAYS:
        mapper = MetaboliteMapper(mg)
        node_type = "metabolite"
    elif assay_type in GENE_ASSAYS:
        mapper = GeneMapper(mg)
        node_type = "gene"
    else:
        # Default: any non-metabolite-class assay falls through to gene
        # mapping. Log a warning so unrecognised assay types are visible.
        log.warning(
            "build_context_from_table: assay_type=%r not recognised; "
            "defaulting to gene-mapper. Known: %s",
            assay_type, sorted(METABOLITE_ASSAYS | GENE_ASSAYS),
        )
        mapper = GeneMapper(mg)
        node_type = "gene"

    ctx = SampleContext(sample_id=sample_id, cohort_id=cohort_id)

    n_mapped = n_filtered = n_unmapped = 0

    for row in rows:
        try:
            feature_id  = str(row[fi_i]).strip()
            effect_size = float(row[eff_i])
        except (IndexError, TypeError, ValueError):
            continue

        p_value = float(row[pv_i]) if pv_i is not None and row[pv_i] is not None else None
        fdr     = float(row[fdr_i]) if fdr_i is not None and row[fdr_i] is not None else None

        # Threshold filtering
        if p_value is not None and p_value > p_value_threshold:
            n_filtered += 1
            continue
        if fdr is not None and fdr > fdr_threshold:
            n_filtered += 1
            continue

        node_id, confidence = mapper.map(feature_id)
        if node_id is None:
            n_unmapped += 1
            continue

        ctx.add(EvidenceRecord(
            feature_id  = feature_id,
            node_id     = node_id,
            node_type   = node_type,
            effect_size = effect_size,
            direction   = 1 if effect_size >= 0 else -1,
            p_value     = p_value,
            fdr         = fdr,
            confidence  = confidence,
            assay_type  = assay_type,
            sample_id   = sample_id,
        ))
        n_mapped += 1

    log.info(
        "build_context_from_table [%s]: %d rows → %d mapped, %d unmapped, %d filtered",
        assay_type, len(rows), n_mapped, n_unmapped, n_filtered,
    )
    return ctx


# ---------------------------------------------------------------------------
# Convenience: score in one call
# ---------------------------------------------------------------------------

def score_from_table(
    mg,
    table_data: dict[str, Any],
    assay_type: str = "metabolomics",
    sample_id: str = "",
    cohort_id: str = "",
    p_value_threshold: float = 1.0,
    fdr_threshold: float = 1.0,
    min_evidence: int = 1,
    engine: str = "heuristic",
    bp_kwargs: dict | None = None,
    **score_kwargs,
) -> list[ReactionScore]:
    """
    Build a SampleContext from a GrAndMA table_data dict and score reactions.

    Parameters
    ----------
    mg, table_data, assay_type, sample_id, cohort_id,
    p_value_threshold, fdr_threshold
        — forwarded to :func:`build_context_from_table`
    min_evidence
        — minimum mapped features per reaction (default 1; heuristic only)
    engine
        — "heuristic" (default, ``score_reactions``) or "bp"
        (``run_bayesian_inference`` + reaction-level posteriors).
    bp_kwargs
        — passed to ``BPConfig`` when ``engine="bp"`` (e.g.
        ``{"restrict_to_observed_hops": 2, "damping": 0.7, "max_iter": 80}``).
        Conditional-currency skip edges are auto-populated from
        ``ensure_scoring_ready(mg)`` and need not be passed explicitly.
    **score_kwargs
        — forwarded to :func:`score_reactions` when engine="heuristic"
        (e.g. ``hub_penalty=False``).

    Returns
    -------
    List of ReactionScore sorted by |score| descending.
    """
    # Run scoring-readiness setup BEFORE context-building, so gene-node
    # enrichment lands before GeneMapper indexes the graph.
    ensure_scoring_ready(mg)
    ctx = build_context_from_table(
        mg, table_data,
        assay_type        = assay_type,
        sample_id         = sample_id,
        cohort_id         = cohort_id,
        p_value_threshold = p_value_threshold,
        fdr_threshold     = fdr_threshold,
    )
    scores = _score_with_engine(
        mg, ctx, engine, min_evidence=min_evidence,
        bp_kwargs=bp_kwargs, **score_kwargs,
    )
    return scores


def _score_with_engine(
    mg,
    ctx: SampleContext,
    engine: str,
    *,
    min_evidence: int = 1,
    bp_kwargs: dict | None = None,
    **score_kwargs,
) -> list[ReactionScore]:
    """Dispatch heuristic vs BP engine and return ReactionScore list."""
    skip_edges = ensure_scoring_ready(mg)

    if engine == "heuristic":
        scores = score_reactions(mg, ctx, min_evidence=min_evidence,
                                  **score_kwargs)
        scores.sort(key=lambda r: abs(r.score), reverse=True)
        return scores

    if engine == "bp":
        from gizmo.inference import run_bayesian_inference, BPConfig
        cfg_kwargs = {
            "restrict_to_observed_hops": 2,
            "damping": 0.7,
            "max_iter": 80,
            "skip_currency": True,
            "skip_edges": skip_edges,
        }
        cfg_kwargs.update(bp_kwargs or {})
        cfg = BPConfig(**cfg_kwargs)
        res = run_bayesian_inference(mg, ctx, cfg)
        return _bp_posteriors_to_scores(mg, res)

    raise ValueError(f"unknown engine={engine!r}; expected 'heuristic' or 'bp'")


def _bp_posteriors_to_scores(mg, res) -> list[ReactionScore]:
    """Convert BPResult posteriors → ReactionScore list, reactions only."""
    g = mg.graph
    scores: list[ReactionScore] = []
    for nid, post in res.posteriors.items():
        attrs = g.nodes.get(nid, {})
        if attrs.get("node_type") != "reaction":
            continue
        p_down = float(post[0])
        p_normal = float(post[1])
        p_up = float(post[2])
        # signed magnitude in [-1, +1]; abs(score) = total perturbation
        signed = p_up - p_down
        scores.append(ReactionScore(
            reaction_id   = nid,
            score         = signed,             # signed; |score| = perturbation
            direction     = signed,
            reversible    = bool(attrs.get("reversible", False)),
            evidence_count = 0,                 # BP does not track per-rxn evidence
            pathway_ids   = list(attrs.get("pathways") or []),
            ec_numbers    = list(attrs.get("ec_numbers") or []),
            confidence    = 1.0 - p_normal,
            notes         = (f"BP posterior: P(down)={p_down:.3f} "
                              f"P(normal)={p_normal:.3f} P(up)={p_up:.3f}"),
        ))
    scores.sort(key=lambda r: abs(r.score), reverse=True)
    return scores


# ---------------------------------------------------------------------------
# Multi-assay helper: merge metabolomics + transcriptomics in one call
# ---------------------------------------------------------------------------

@dataclass
class PerSampleResult:
    """Sample × reaction posterior matrix from per-sample BP.

    Attributes
    ----------
    sample_ids   : list[str]    — order matches axis 0 of ``posteriors``
    reaction_ids : list[str]    — order matches axis 1 of ``posteriors``
    posteriors   : np.ndarray   — shape (n_samples, n_reactions, 3),
                                  3 columns are P(DOWN), P(NORMAL), P(UP)
    bp_meta      : dict         — convergence flags / iteration counts per sample

    Convenience properties
    ----------------------
    signed_scores        : (n_samples, n_reactions) = P(UP) - P(DOWN)
    perturbation_scores  : (n_samples, n_reactions) = 1 - P(NORMAL)
    """
    sample_ids:   list
    reaction_ids: list
    posteriors:   "np.ndarray"
    bp_meta:      dict = field(default_factory=dict)

    @property
    def signed_scores(self):
        return self.posteriors[:, :, 2] - self.posteriors[:, :, 0]

    @property
    def perturbation_scores(self):
        return 1.0 - self.posteriors[:, :, 1]

    def reaction_score_dict(self, sample_id: str,
                              metric: str = "signed") -> dict[str, float]:
        """Return {reaction_id: score} for one sample.
        ``metric`` ∈ {'signed', 'perturbation'}."""
        if sample_id not in self.sample_ids:
            raise KeyError(sample_id)
        idx = self.sample_ids.index(sample_id)
        if metric == "signed":
            scores = self.signed_scores[idx]
        elif metric == "perturbation":
            scores = self.perturbation_scores[idx]
        else:
            raise ValueError(f"unknown metric={metric!r}")
        return dict(zip(self.reaction_ids, scores))

    def reaction_phenotype_assoc(
        self,
        phenotype_values: dict[str, "str | float"],
        method: str = "auto",
    ) -> dict[str, float]:
        """For each reaction, test whether its perturbation score across
        samples differs by phenotype.

        Parameters
        ----------
        phenotype_values : {sample_id: value}
            Categorical (e.g. "RA"/"OA") → Mann-Whitney U.
            Continuous (e.g. age, DAS28) → Spearman correlation.
            ``method`` can override.

        Returns
        -------
        {reaction_id: p_value}  — uncorrected p; FDR adjustment is the
        caller's responsibility.
        """
        import numpy as np
        # Filter to samples we have phenotype for
        keep_idx = [i for i, sid in enumerate(self.sample_ids)
                    if sid in phenotype_values]
        if len(keep_idx) < 4:
            return {}
        scores = self.perturbation_scores[keep_idx]   # (n_keep, n_reactions)
        kept_values = [phenotype_values[self.sample_ids[i]] for i in keep_idx]

        # Detect categorical vs continuous
        is_cat = (method == "categorical" or
                   (method == "auto" and not all(isinstance(v, (int, float))
                                                  for v in kept_values)))

        out: dict[str, float] = {}
        if is_cat:
            from scipy.stats import mannwhitneyu
            cats = sorted({str(v) for v in kept_values})
            if len(cats) != 2:
                # Use Kruskal-Wallis for >2 groups
                from scipy.stats import kruskal
                groups = {c: [scores[i] for i, v in enumerate(kept_values)
                              if str(v) == c] for c in cats}
                for j, rxn_id in enumerate(self.reaction_ids):
                    arrs = [np.array([g[j] for g in groups[c]]) for c in cats]
                    if any(np.var(a) == 0 for a in arrs):
                        out[rxn_id] = 1.0
                        continue
                    try:
                        _, p = kruskal(*arrs)
                        out[rxn_id] = float(p)
                    except Exception:
                        out[rxn_id] = 1.0
            else:
                a_idx = [i for i, v in enumerate(kept_values) if str(v) == cats[0]]
                b_idx = [i for i, v in enumerate(kept_values) if str(v) == cats[1]]
                for j, rxn_id in enumerate(self.reaction_ids):
                    a = scores[a_idx, j]
                    b = scores[b_idx, j]
                    if np.var(a) + np.var(b) == 0:
                        out[rxn_id] = 1.0
                        continue
                    try:
                        _, p = mannwhitneyu(a, b, alternative="two-sided")
                        out[rxn_id] = float(p)
                    except Exception:
                        out[rxn_id] = 1.0
        else:
            from scipy.stats import spearmanr
            cont_vals = np.asarray(kept_values, dtype=float)
            for j, rxn_id in enumerate(self.reaction_ids):
                col = scores[:, j]
                if np.var(col) == 0:
                    out[rxn_id] = 1.0
                    continue
                try:
                    _, p = spearmanr(col, cont_vals)
                    out[rxn_id] = float(p) if not np.isnan(p) else 1.0
                except Exception:
                    out[rxn_id] = 1.0
        return out


# --------------------------------------------------------------------------
# Per-sample BP worker (module-scope so multiprocessing.Pool can pickle it)
# --------------------------------------------------------------------------

_PERSAMPLE_GLOBAL: dict = {}


def _per_sample_init(mg, cfg, feature_node, references, feature_assay,
                      node_type, apply_gate=False):
    """Pool worker initializer: cache graph + cfg + mapper-resolved features
    in worker globals so each task only carries the per-sample data."""
    _PERSAMPLE_GLOBAL["mg"] = mg
    _PERSAMPLE_GLOBAL["cfg"] = cfg
    _PERSAMPLE_GLOBAL["feature_node"] = feature_node
    _PERSAMPLE_GLOBAL["references"] = references
    _PERSAMPLE_GLOBAL["feature_assay"] = feature_assay
    _PERSAMPLE_GLOBAL["node_type"] = node_type
    _PERSAMPLE_GLOBAL["apply_gate"] = apply_gate


def _per_sample_run(args):
    """Pool worker task: run BP for one (sample_id, sample_features) pair."""
    sample_id, sample_features = args
    mg            = _PERSAMPLE_GLOBAL["mg"]
    cfg           = _PERSAMPLE_GLOBAL["cfg"]
    feature_node  = _PERSAMPLE_GLOBAL["feature_node"]
    references    = _PERSAMPLE_GLOBAL["references"]
    feature_assay = _PERSAMPLE_GLOBAL["feature_assay"]
    node_type     = _PERSAMPLE_GLOBAL["node_type"]
    apply_gate    = _PERSAMPLE_GLOBAL.get("apply_gate", False)
    from gizmo.inference import run_bayesian_inference
    ctx = SampleContext(sample_id=str(sample_id))
    # Per-sample gene evidence map (used for the asymmetric enzyme gate
    # if the assay routed via GeneMapper). Keyed by gene_node_id rather
    # than feature_id so the gate code can look up directly.
    sample_gene_evidence: dict[str, float] = {}
    for f, (nid, conf) in feature_node.items():
        if f not in sample_features or f not in references:
            continue
        v = sample_features[f]
        if v is None:
            continue
        effect = float(v) - references[f]
        ctx.add(EvidenceRecord(
            feature_id=f, node_id=nid, node_type=node_type,
            effect_size=effect,
            direction=1 if effect >= 0 else -1,
            confidence=conf, assay_type=feature_assay,
            sample_id=str(sample_id),
        ))
        if node_type == "gene":
            sample_gene_evidence[nid] = effect
    res = run_bayesian_inference(mg, ctx, cfg)
    posteriors_dict = {nid: tuple(post) for nid, post in res.posteriors.items()}

    # Asymmetric enzyme-availability gate. Applied multiplicatively to
    # each reaction's posterior perturbation magnitude (1-P(NORMAL)),
    # not to the directional sign — gene evidence is capacity-limiting,
    # not state-driving. See gizmo.scoring.enzyme_gate for the rationale.
    gate_dict: dict[str, float] = {}
    if apply_gate and sample_gene_evidence:
        from gizmo.scoring.enzyme_gate import compute_enzyme_gate
        gate_dict = compute_enzyme_gate(mg, sample_gene_evidence)
        if gate_dict:
            # Rescale posteriors: keep direction (P(UP)/P(DOWN) ratio),
            # but dampen the perturbed mass by gate; redistribute to NORMAL.
            for rxn_id, gate_val in gate_dict.items():
                post = posteriors_dict.get(rxn_id)
                if post is None or gate_val >= 1.0:
                    continue
                p_down, p_normal, p_up = post
                # Move (1 - gate) * (P(UP) + P(DOWN)) into NORMAL
                pert_mass = p_up + p_down
                kept_pert = pert_mass * gate_val
                if pert_mass > 1e-12:
                    scale = kept_pert / pert_mass
                    new_p_up = p_up * scale
                    new_p_down = p_down * scale
                else:
                    new_p_up = p_up
                    new_p_down = p_down
                new_p_normal = 1.0 - new_p_up - new_p_down
                posteriors_dict[rxn_id] = (new_p_down, new_p_normal, new_p_up)

    meta = {"iterations": res.iterations,
            "converged":  res.converged,
            "n_evidence": sum(1 for _ in ctx.records(mapped_only=True)),
            "n_gated":    sum(1 for v in gate_dict.values() if v < 1.0)}
    return sample_id, posteriors_dict, meta


def score_per_sample_bp(
    mg,
    feature_matrix: dict[str, dict[str, float]],
    feature_assay: str = "transcriptomics",
    reference_method: str = "cohort_mean",
    control_samples: list[str] | None = None,
    bp_kwargs: dict | None = None,
    progress: bool = False,
    n_workers: int | None = None,
    apply_enzyme_gate: bool = True,
) -> PerSampleResult:
    """
    Run BP per sample and return a sample × reaction posterior matrix.

    For each sample S, the per-feature effect_size used as evidence is the
    deviation of S's value from a reference value (cohort or control mean).

    Parameters
    ----------
    mg
        GizmoGraph
    feature_matrix
        ``{sample_id: {feature_id: value}}`` — wide-format expression /
        abundance matrix on the original (typically log2) scale.
    feature_assay
        ``"transcriptomics"`` (default) routes via ``GeneMapper``,
        ``"metabolomics"`` via ``MetaboliteMapper``. Other values follow
        the same dispatch rules as ``build_context_from_table``.
    reference_method
        ``"cohort_mean"`` (default) — reference = mean across all samples.
        ``"control_mean"`` — mean across ``control_samples`` only.
        ``"median"`` — median across all samples.
    control_samples
        List of sample ids to use as controls (only relevant when
        ``reference_method="control_mean"``).
    bp_kwargs
        Forwarded to ``BPConfig`` (e.g. ``{"max_iter": 80}``).
        Conditional-currency skip edges are auto-populated from
        ``ensure_scoring_ready(mg)``.
    progress
        If True, log per-sample progress.
    n_workers
        Number of process workers to use for parallel per-sample BP. None
        or <=1 = serial (default; lowest overhead for small cohorts).
        On Linux, workers are fork()ed and inherit the graph + mapper
        cache via copy-on-write so memory cost is roughly one graph
        regardless of pool size.
    apply_enzyme_gate
        When True (default) and ``feature_assay`` is gene-class, applies
        the asymmetric enzyme-availability gate per sample: reactions
        whose catalyzing genes are strongly DOWN in this sample have
        their perturbation magnitude dampened. Restores the metabolite-
        side AUC that symmetric gene-edge coupling collapses (verified
        on GSE190504 IDH cohort).

    Returns
    -------
    PerSampleResult
    """
    import numpy as np
    from gizmo.inference import run_bayesian_inference, BPConfig

    skip_edges = ensure_scoring_ready(mg)
    g = mg.graph

    # ---- 1. Compute reference per feature ------------------------------
    sample_ids = sorted(feature_matrix.keys())
    if not sample_ids:
        raise ValueError("empty feature_matrix")
    all_features: set[str] = set()
    for sid in sample_ids:
        all_features.update(feature_matrix[sid].keys())
    feature_list = sorted(all_features)

    if reference_method == "cohort_mean":
        ref_samples = sample_ids
    elif reference_method == "control_mean":
        if not control_samples:
            raise ValueError("control_samples required when "
                             "reference_method='control_mean'")
        ref_samples = [s for s in control_samples if s in feature_matrix]
        if len(ref_samples) < 3:
            raise ValueError(f"only {len(ref_samples)} control samples "
                             "available; need >=3")
    elif reference_method == "median":
        ref_samples = sample_ids
    else:
        raise ValueError(f"unknown reference_method={reference_method!r}")

    references: dict[str, float] = {}
    for f in feature_list:
        vals = [feature_matrix[s].get(f) for s in ref_samples
                if f in feature_matrix[s]]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        references[f] = (sum(vals) / len(vals)
                          if reference_method != "median"
                          else sorted(vals)[len(vals) // 2])

    # ---- 2. Identify reaction node ordering once -----------------------
    rxn_ids = sorted(nid for nid, attrs in g.nodes(data=True)
                      if attrs.get("node_type") == "reaction")
    rxn_idx = {rid: j for j, rid in enumerate(rxn_ids)}
    n_samples = len(sample_ids)
    n_rxns = len(rxn_ids)
    posteriors = np.zeros((n_samples, n_rxns, 3), dtype=np.float64)
    posteriors[:, :, 1] = 1.0   # default = NORMAL prior
    bp_meta: dict[str, dict] = {}

    # ---- 3. Mapper for the feature assay -------------------------------
    METABOLITE_ASSAYS = {"metabolomics", "lipidomics", "small_molecule"}
    GENE_ASSAYS = {"transcriptomics", "proteomics", "phosphoproteomics",
                    "translatomics", "ribosomics", "expression"}
    if feature_assay in METABOLITE_ASSAYS:
        from gizmo.evidence.mappers import MetaboliteMapper
        mapper = MetaboliteMapper(mg)
        node_type = "metabolite"
    elif feature_assay in GENE_ASSAYS:
        from gizmo.evidence.mappers import GeneMapper
        mapper = GeneMapper(mg)
        node_type = "gene"
    else:
        from gizmo.evidence.mappers import GeneMapper
        mapper = GeneMapper(mg)
        node_type = "gene"
        log.warning("score_per_sample_bp: assay_type=%r unknown; default to gene",
                    feature_assay)

    # Pre-compute feature → node_id, conf
    feature_node: dict[str, tuple[str, float]] = {}
    for f in feature_list:
        nid, conf = mapper.map(f)
        if nid:
            feature_node[f] = (nid, conf)

    # ---- 4. Run BP per sample -----------------------------------------
    cfg_kwargs = {
        "restrict_to_observed_hops": 2,
        "damping": 0.7,
        "max_iter": 80,
        "skip_currency": True,
        "skip_edges": skip_edges,
    }
    cfg_kwargs.update(bp_kwargs or {})
    cfg = BPConfig(**cfg_kwargs)

    # Each task carries (sample_id, {feature: value}) — small. Worker init
    # does heavy graph + cfg + mapper sharing.
    tasks = [(sid, feature_matrix[sid]) for sid in sample_ids]

    if n_workers is None or n_workers <= 1:
        # Serial path
        _per_sample_init(mg, cfg, feature_node, references, feature_assay,
                          node_type, apply_gate=apply_enzyme_gate)
        for i, task in enumerate(tasks):
            sid = task[0]
            _, posts, meta = _per_sample_run(task)
            for nid, post in posts.items():
                j = rxn_idx.get(nid)
                if j is None:
                    continue
                posteriors[i, j, :] = post
            bp_meta[sid] = meta
            if progress:
                log.info("score_per_sample_bp: %d/%d %s done (%d evidence, %d iter)",
                         i + 1, n_samples, sid, meta["n_evidence"],
                         meta["iterations"])
    else:
        # Parallel path. ``fork()`` (Linux default) shares the graph via
        # copy-on-write so memory cost is one graph regardless of pool size.
        from multiprocessing import Pool, get_context
        ctx_mp = get_context("fork")   # explicit; spawn would re-pickle mg
        idx_for_sid = {sid: i for i, sid in enumerate(sample_ids)}
        with ctx_mp.Pool(
            processes=n_workers,
            initializer=_per_sample_init,
            initargs=(mg, cfg, feature_node, references, feature_assay,
                       node_type, apply_enzyme_gate),
        ) as pool:
            for k, (sid, posts, meta) in enumerate(
                pool.imap_unordered(_per_sample_run, tasks, chunksize=1)
            ):
                i = idx_for_sid[sid]
                for nid, post in posts.items():
                    j = rxn_idx.get(nid)
                    if j is None:
                        continue
                    posteriors[i, j, :] = post
                bp_meta[sid] = meta
                if progress:
                    log.info("score_per_sample_bp: %d/%d %s done (%d evidence, %d iter)",
                             k + 1, n_samples, sid, meta["n_evidence"],
                             meta["iterations"])

    return PerSampleResult(
        sample_ids   = sample_ids,
        reaction_ids = rxn_ids,
        posteriors   = posteriors,
        bp_meta      = bp_meta,
    )


def score_per_sample_laplacian(
    mg,
    feature_matrix: dict[str, dict[str, float]],
    feature_assay: str = "transcriptomics",
    reference_method: str = "cohort_mean",
    control_samples: list[str] | None = None,
    alpha: float = 0.5,
    progress: bool = False,
    n_workers: int | None = None,
    laplacian_kwargs: dict | None = None,
) -> PerSampleResult:
    """
    Drop-in Laplacian replacement for ``score_per_sample_bp``.

    Same input/output API as ``score_per_sample_bp``: returns a
    ``PerSampleResult`` with sample × reaction × 3 posterior tensor.

    Why prefer this over BP:
      - Direction-aware (signed evidence + signed mass-action edges)
      - Topology-controlled (D^{-1/2} normalization)
      - No iterative convergence: single closed-form linear solve per sample
      - ~200× faster than BP on the same graph
      - Validated on 6 cohorts including 2 truly individually-paired
        (Crohn's, Su 2020 COVID); see benchmarks/results/METHOD_HISTORY.md

    Parameters mirror ``score_per_sample_bp``. ``alpha`` controls smoothing
    (α=0.5 was best across all cohorts in our sweep). ``laplacian_kwargs``
    are forwarded to ``LaplacianConfig`` (e.g. ``substrate_weight``,
    ``modifier_weight``).

    Drops the asymmetric enzyme-gate path: the gate produces output
    identical to metab-only on RA/IDH but inflates σ via perm narrowing
    rather than signal extraction, so we don't carry it forward.
    """
    import numpy as np
    from gizmo.inference.laplacian import (
        run_laplacian_inference, LaplacianConfig,
    )

    skip_edges = ensure_scoring_ready(mg)
    g = mg.graph

    sample_ids = sorted(feature_matrix.keys())
    if not sample_ids:
        raise ValueError("empty feature_matrix")
    all_features: set[str] = set()
    for sid in sample_ids:
        all_features.update(feature_matrix[sid].keys())
    feature_list = sorted(all_features)

    if reference_method == "cohort_mean":
        ref_samples = sample_ids
    elif reference_method == "control_mean":
        if not control_samples:
            raise ValueError("control_samples required when "
                             "reference_method='control_mean'")
        ref_samples = [s for s in control_samples if s in feature_matrix]
        if len(ref_samples) < 3:
            raise ValueError(f"only {len(ref_samples)} control samples "
                             "available; need >=3")
    elif reference_method == "median":
        ref_samples = sample_ids
    else:
        raise ValueError(f"unknown reference_method={reference_method!r}")

    references: dict[str, float] = {}
    for f in feature_list:
        vals = [feature_matrix[s].get(f) for s in ref_samples
                if f in feature_matrix[s]]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        references[f] = (sum(vals) / len(vals)
                          if reference_method != "median"
                          else sorted(vals)[len(vals) // 2])

    rxn_ids = sorted(nid for nid, attrs in g.nodes(data=True)
                      if attrs.get("node_type") == "reaction")
    rxn_idx = {rid: j for j, rid in enumerate(rxn_ids)}
    n_samples = len(sample_ids)
    n_rxns = len(rxn_ids)
    posteriors = np.zeros((n_samples, n_rxns, 3), dtype=np.float64)
    posteriors[:, :, 1] = 1.0
    sample_meta: dict[str, dict] = {}

    METABOLITE_ASSAYS = {"metabolomics", "lipidomics", "small_molecule"}
    GENE_ASSAYS = {"transcriptomics", "proteomics", "phosphoproteomics",
                    "translatomics", "ribosomics", "expression"}
    if feature_assay in METABOLITE_ASSAYS:
        from gizmo.evidence.mappers import MetaboliteMapper
        mapper = MetaboliteMapper(mg)
        node_type = "metabolite"
    else:
        from gizmo.evidence.mappers import GeneMapper
        mapper = GeneMapper(mg)
        node_type = "gene"

    feature_node: dict[str, tuple[str, float]] = {}
    if feature_assay in METABOLITE_ASSAYS:
        # Cohort-aware metab name handling: suffixes, KEGG-IDs, Metabolon
        # lipid notation, optional PubChem name → InChIKey fallback.
        from gizmo.evidence.feature_normalize import map_with_fallback
        pubchem_cache = (laplacian_kwargs or {}).get("pubchem_name_cache")
        for f in feature_list:
            nid, conf = map_with_fallback(mapper, f, pubchem_cache)
            if nid:
                feature_node[f] = (nid, conf)
    else:
        for f in feature_list:
            nid, conf = mapper.map(f)
            if nid:
                feature_node[f] = (nid, conf)

    cfg_kwargs = {"alpha": alpha, "skip_edges": skip_edges}
    # Strip the cache key (it's just used by us for mapping, not LaplacianConfig)
    cfg_kwargs.update({k: v for k, v in (laplacian_kwargs or {}).items()
                        if k != "pubchem_name_cache"})
    cfg = LaplacianConfig(**cfg_kwargs)

    for i, sid in enumerate(sample_ids):
        ctx = SampleContext(sample_id=str(sid))
        n_evidence = 0
        for f, (nid, conf) in feature_node.items():
            if f not in feature_matrix[sid] or f not in references:
                continue
            v = feature_matrix[sid][f]
            if v is None:
                continue
            effect = float(v) - references[f]
            ctx.add(EvidenceRecord(
                feature_id=f, node_id=nid, node_type=node_type,
                effect_size=effect,
                direction=1 if effect >= 0 else -1,
                confidence=conf, assay_type=feature_assay,
                sample_id=str(sid),
            ))
            n_evidence += 1
        res = run_laplacian_inference(mg, ctx, cfg)
        for nid, post in res.posteriors.items():
            j = rxn_idx.get(nid)
            if j is None:
                continue
            posteriors[i, j, :] = post
        sample_meta[sid] = {"n_evidence": n_evidence, "alpha": alpha}
        if progress:
            log.info("score_per_sample_laplacian: %d/%d %s done (%d evidence)",
                     i + 1, n_samples, sid, n_evidence)

    return PerSampleResult(
        sample_ids   = sample_ids,
        reaction_ids = rxn_ids,
        posteriors   = posteriors,
        bp_meta      = sample_meta,
    )


def score_from_tables(
    mg,
    tables: list[dict[str, Any]],
    assay_types: list[str],
    sample_ids: list[str] | None = None,
    p_value_threshold: float = 1.0,
    fdr_threshold: float = 1.0,
    min_evidence: int = 1,
    engine: str = "heuristic",
    bp_kwargs: dict | None = None,
    **score_kwargs,
) -> list[ReactionScore]:
    """
    Score reactions from multiple GrAndMA table_data dicts simultaneously.

    Typical use: combine a metabolomics t-test result with a transcriptomics
    differential expression result for the same comparison.

    Parameters
    ----------
    mg            : GizmoGraph
    tables        : list of table_data dicts (one per assay)
    assay_types   : list of assay type strings, same length as tables
    sample_ids    : optional list of sample/analysis IDs, same length as tables
    engine        : "heuristic" (default) or "bp"
    bp_kwargs     : forwarded to ``BPConfig`` when ``engine="bp"``
    """
    if len(tables) != len(assay_types):
        raise ValueError("tables and assay_types must have the same length")
    if sample_ids is None:
        sample_ids = [""] * len(tables)

    # Run scoring-readiness setup BEFORE context-building (see score_from_table)
    ensure_scoring_ready(mg)

    # Build a shared context by merging all tables
    merged_ctx = SampleContext(sample_id="merged")
    for table, atype, sid in zip(tables, assay_types, sample_ids):
        ctx = build_context_from_table(
            mg, table,
            assay_type        = atype,
            sample_id         = sid,
            p_value_threshold = p_value_threshold,
            fdr_threshold     = fdr_threshold,
        )
        for rec in ctx.records():
            merged_ctx.add(rec)

    return _score_with_engine(
        mg, merged_ctx, engine, min_evidence=min_evidence,
        bp_kwargs=bp_kwargs, **score_kwargs,
    )
