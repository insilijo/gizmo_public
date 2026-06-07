"""
CSV ingestion for multiomic evidence.

Standard column conventions
----------------------------
Metabolomics / proteomics CSV:

    FEATURE_ID   required  ChEBI, HMDB, InChIKey, PubChem CID, Metabolon name, or graph node ID
    EFFECT_SIZE  required  signed log2FC or z-score
    P_VALUE      optional
    FDR          optional
    SAMPLE_ID    optional  overrides the sample_id argument
    NOTES        optional

Transcriptomics CSV:

    FEATURE_ID   required  gene symbol, ENSG ID, or graph node ID
    EFFECT_SIZE  required
    P_VALUE      optional
    FDR          optional
    SAMPLE_ID    optional
    NOTES        optional

Column names are case-insensitive.  Common aliases are also accepted:
  BIOCHEMICAL, COMPOUND, NAME            → FEATURE_ID
  LOG2FC, LOG2_FC, FC, LOGFC            → EFFECT_SIZE
  PVAL, P.VALUE, P_VAL                  → P_VALUE
  Q_VALUE, QVALUE, PADJ, ADJ_P         → FDR
  GENE, GENE_SYMBOL, SYMBOL, GENE_ID   → FEATURE_ID (transcriptomics)

Usage::

    from gizmo.evidence.ingest import load_metabolomics_csv, load_transcriptomics_csv

    ctx = SampleContext(sample_id="PAT_001")
    n_met = load_metabolomics_csv("data/metabolomics.csv", mg, ctx)
    n_tx  = load_transcriptomics_csv("data/rna.csv", mg, ctx)
    ctx.print_summary()
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from gizmo.evidence.mappers import GeneMapper, MetaboliteMapper
from gizmo.evidence.model import EvidenceRecord, SampleContext

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Column name aliases
# ---------------------------------------------------------------------------

_FEATURE_ALIASES = {
    "biochemical", "compound", "compound_name", "name",
    "feature", "feature_id", "id",
    "gene", "gene_symbol", "symbol", "gene_id",
    "metabolite",
}
_EFFECT_ALIASES = {
    "effect_size", "log2fc", "log2_fc", "logfc", "log2foldchange",
    "fc", "foldchange", "zscore", "z_score", "beta", "coef", "coefficient",
}
_PVAL_ALIASES  = {"p_value", "pval", "p.value", "p_val", "pvalue"}
_FDR_ALIASES   = {"fdr", "q_value", "qvalue", "padj", "adj_p", "adj.p.val", "bh"}
_SAMPLE_ALIASES = {"sample_id", "sample", "patient", "patient_id", "subject"}
_NOTES_ALIASES  = {"notes", "note", "annotation"}


def _find_col(df: pd.DataFrame, aliases: set[str]) -> Optional[str]:
    """Return the first column name (case-insensitive) that matches an alias."""
    lower = {c.lower(): c for c in df.columns}
    for alias in aliases:
        if alias in lower:
            return lower[alias]
    return None


def _rename_df(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise column names to standard keys."""
    rename = {}
    lower  = {c.lower(): c for c in df.columns}

    for std, aliases in [
        ("FEATURE_ID",  _FEATURE_ALIASES),
        ("EFFECT_SIZE", _EFFECT_ALIASES),
        ("P_VALUE",     _PVAL_ALIASES),
        ("FDR",         _FDR_ALIASES),
        ("SAMPLE_ID",   _SAMPLE_ALIASES),
        ("NOTES",       _NOTES_ALIASES),
    ]:
        for alias in aliases:
            if alias in lower and std not in df.columns:
                rename[lower[alias]] = std
                break

    return df.rename(columns=rename)


# ---------------------------------------------------------------------------
# Generic loader
# ---------------------------------------------------------------------------

def _load_csv(
    path: str | Path,
    mg,
    ctx: SampleContext,
    mapper,
    assay_type: str,
    node_type: str,
    sample_id: Optional[str] = None,
    min_confidence: float = 0.0,
) -> int:
    """
    Core CSV loader.  Returns number of records added to ctx.
    """
    df = pd.read_csv(path)
    df = _rename_df(df)

    if "FEATURE_ID" not in df.columns:
        raise ValueError(
            f"Cannot find feature ID column in {path}.  "
            f"Expected one of: {sorted(_FEATURE_ALIASES)}"
        )
    if "EFFECT_SIZE" not in df.columns:
        raise ValueError(
            f"Cannot find effect size column in {path}.  "
            f"Expected one of: {sorted(_EFFECT_ALIASES)}"
        )

    n_added   = 0
    n_unmapped = 0

    for _, row in df.iterrows():
        fid = str(row["FEATURE_ID"]).strip()
        if not fid or fid.lower() in ("nan", "na", ""):
            continue

        try:
            effect = float(row["EFFECT_SIZE"])
        except (ValueError, TypeError):
            continue

        node_id, conf = mapper.map(fid)
        if conf < min_confidence:
            n_unmapped += 1
            if node_id is None:
                node_id = None  # keep unmapped records for traceability

        p_val = _safe_float(row.get("P_VALUE"))
        fdr   = _safe_float(row.get("FDR"))
        sid   = str(row["SAMPLE_ID"]).strip() if "SAMPLE_ID" in row and str(row.get("SAMPLE_ID", "")) not in ("nan", "") else (sample_id or ctx.sample_id)
        notes = str(row.get("NOTES", "")).strip() if "NOTES" in row else ""

        ctx.add(EvidenceRecord(
            feature_id  = fid,
            node_id     = node_id,
            node_type   = node_type,
            effect_size = effect,
            p_value     = p_val,
            fdr         = fdr,
            confidence  = conf,
            assay_type  = assay_type,
            sample_id   = sid,
            cohort_id   = ctx.cohort_id,
            notes       = notes,
        ))
        n_added += 1

    n_mapped = n_added - n_unmapped
    log.info(
        "%s: %d rows → %d records (%d mapped, %d unmapped)",
        assay_type, len(df), n_added, n_mapped, n_unmapped,
    )
    return n_added


# ---------------------------------------------------------------------------
# Public loaders
# ---------------------------------------------------------------------------

def load_metabolomics_csv(
    path: str | Path,
    mg,
    ctx: SampleContext,
    sample_id: Optional[str] = None,
    min_confidence: float = 0.0,
) -> int:
    """
    Load a metabolomics CSV into a SampleContext.

    FEATURE_ID column should contain ChEBI IDs, HMDB IDs (e.g. HMDB0000001),
    InChIKeys, PubChem CIDs, Metabolon biochemical names, or existing graph node IDs.

    Returns number of EvidenceRecords added.
    """
    mapper = MetaboliteMapper(mg)
    return _load_csv(
        path, mg, ctx, mapper,
        assay_type="metabolomics",
        node_type="metabolite",
        sample_id=sample_id,
        min_confidence=min_confidence,
    )


def load_transcriptomics_csv(
    path: str | Path,
    mg,
    ctx: SampleContext,
    sample_id: Optional[str] = None,
    min_confidence: float = 0.0,
) -> int:
    """
    Load a transcriptomics (RNA-seq / microarray) CSV into a SampleContext.

    FEATURE_ID column should contain HGNC gene symbols, Ensembl IDs,
    or existing gene node IDs.

    Returns number of EvidenceRecords added.
    """
    mapper = GeneMapper(mg)
    return _load_csv(
        path, mg, ctx, mapper,
        assay_type="transcriptomics",
        node_type="gene",
        sample_id=sample_id,
        min_confidence=min_confidence,
    )


def load_proteomics_csv(
    path: str | Path,
    mg,
    ctx: SampleContext,
    sample_id: Optional[str] = None,
    min_confidence: float = 0.0,
) -> int:
    """
    Load a proteomics CSV into a SampleContext.

    FEATURE_ID column should contain gene symbols or Ensembl IDs
    (protein → gene node mapping).

    Returns number of EvidenceRecords added.
    """
    mapper = GeneMapper(mg)
    return _load_csv(
        path, mg, ctx, mapper,
        assay_type="proteomics",
        node_type="gene",
        sample_id=sample_id,
        min_confidence=min_confidence,
    )


def merge_contexts(*contexts: SampleContext, sample_id: str = "merged") -> SampleContext:
    """
    Merge multiple SampleContexts into one (e.g. to aggregate a cohort).

    All records are copied; sample_id in each record is preserved.
    """
    merged = SampleContext(sample_id=sample_id)
    for ctx in contexts:
        merged.add_all(ctx.records(mapped_only=False))
    return merged


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val) -> Optional[float]:
    try:
        f = float(val)
        return f if f == f else None   # NaN check
    except (TypeError, ValueError):
        return None
