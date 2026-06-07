"""
Core evidence data model for GIZMO.

EvidenceRecord
    One mapped measurement from any omic assay.
    Carries the original feature ID, the resolved graph node,
    the effect size, direction, and confidence.

SampleContext
    A lightweight, non-mutating overlay of EvidenceRecords for a
    single sample or sample group.  Applied on top of the canonical
    graph without modifying it.

Usage::

    from gizmo.evidence.model import EvidenceRecord, SampleContext

    ctx = SampleContext(sample_id="PAT_001", cohort_id="cohort_A")
    ctx.add(EvidenceRecord(
        feature_id="CHEBI:15422",
        node_id="CHEBI:15422",
        node_type="metabolite",
        effect_size=1.85,
        direction=1,
        p_value=0.001,
        fdr=0.01,
        confidence=1.0,
        assay_type="metabolomics",
        sample_id="PAT_001",
    ))
    ctx.summary()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class EvidenceRecord:
    """
    One mapped omic measurement.

    Fields
    ------
    feature_id   : original identifier from the assay (column header, ChEBI ID, gene symbol…)
    node_id      : resolved graph node ID (None if unmapped)
    node_type    : "metabolite" | "gene" | "protein" | "reaction"
    effect_size  : signed effect magnitude (log2FC, z-score, …)
    direction    : +1 up / -1 down / 0 unclear — derived from effect_size sign when not provided
    p_value      : raw p-value (None if not available)
    fdr          : FDR-adjusted q-value (None if not available)
    confidence   : mapping confidence [0, 1] — 1.0 for exact matches, lower for fuzzy/inferred
    assay_type   : "metabolomics" | "transcriptomics" | "proteomics" | "phosphoproteomics"
    sample_id    : sample or patient identifier
    cohort_id    : cohort or study identifier
    notes        : free-text annotation
    """

    feature_id:  str
    node_id:     Optional[str]
    node_type:   str
    effect_size: float
    direction:   int              = 0      # set automatically if 0
    p_value:     Optional[float]  = None
    fdr:         Optional[float]  = None
    confidence:  float            = 1.0
    assay_type:  str              = "metabolomics"
    sample_id:   Optional[str]   = None
    cohort_id:   Optional[str]   = None
    notes:       str              = ""

    def __post_init__(self) -> None:
        if self.direction == 0 and self.effect_size != 0:
            self.direction = 1 if self.effect_size > 0 else -1
        # Clamp confidence
        self.confidence = max(0.0, min(1.0, self.confidence))


class SampleContext:
    """
    Non-mutating overlay of multiomic evidence for one sample or group.

    Accumulates EvidenceRecords and provides fast lookup by node ID,
    assay type, or direction.  The canonical graph is never modified.

    Parameters
    ----------
    sample_id  : sample / patient identifier
    cohort_id  : cohort or study identifier
    """

    def __init__(
        self,
        sample_id: str = "sample",
        cohort_id: Optional[str] = None,
    ) -> None:
        self.sample_id = sample_id
        self.cohort_id = cohort_id
        self._records: list[EvidenceRecord] = []
        self._by_node: dict[str, list[EvidenceRecord]] = {}

    # ------------------------------------------------------------------
    # Building
    # ------------------------------------------------------------------

    def add(self, record: EvidenceRecord) -> None:
        """Add one EvidenceRecord."""
        self._records.append(record)
        if record.node_id:
            self._by_node.setdefault(record.node_id, []).append(record)

    def add_all(self, records: list[EvidenceRecord]) -> None:
        for r in records:
            self.add(r)

    # ------------------------------------------------------------------
    # Lookup
    # ------------------------------------------------------------------

    def by_node(self, node_id: str) -> list[EvidenceRecord]:
        """All records mapped to a specific graph node."""
        return self._by_node.get(node_id, [])

    def mapped_nodes(self, assay_type: Optional[str] = None) -> set[str]:
        """Set of node IDs that have at least one evidence record."""
        if assay_type is None:
            return set(self._by_node)
        return {
            r.node_id for r in self._records
            if r.node_id and r.assay_type == assay_type
        }

    def records(
        self,
        assay_type: Optional[str] = None,
        node_type: Optional[str] = None,
        mapped_only: bool = True,
    ) -> list[EvidenceRecord]:
        """Filtered view of all evidence records."""
        out = self._records
        if mapped_only:
            out = [r for r in out if r.node_id]
        if assay_type:
            out = [r for r in out if r.assay_type == assay_type]
        if node_type:
            out = [r for r in out if r.node_type == node_type]
        return out

    def net_effect(self, node_id: str, assay_type: Optional[str] = None) -> Optional[float]:
        """
        Mean signed effect size across all records for a node.
        Returns None if no evidence.
        """
        recs = self.by_node(node_id)
        if assay_type:
            recs = [r for r in recs if r.assay_type == assay_type]
        if not recs:
            return None
        return sum(r.effect_size * r.confidence for r in recs) / len(recs)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> dict:
        total     = len(self._records)
        mapped    = sum(1 for r in self._records if r.node_id)
        met_count = sum(1 for r in self._records if r.assay_type == "metabolomics" and r.node_id)
        tx_count  = sum(1 for r in self._records if r.assay_type == "transcriptomics" and r.node_id)
        pr_count  = sum(1 for r in self._records if r.assay_type == "proteomics" and r.node_id)
        return {
            "sample_id":        self.sample_id,
            "cohort_id":        self.cohort_id,
            "n_total":          total,
            "n_mapped":         mapped,
            "n_unmapped":       total - mapped,
            "n_mapped_nodes":   len(self._by_node),
            "n_metabolomics":   met_count,
            "n_transcriptomics": tx_count,
            "n_proteomics":     pr_count,
        }

    def print_summary(self) -> None:
        s = self.summary()
        print(f"SampleContext  sample={s['sample_id']}  cohort={s['cohort_id']}")
        print(f"  Total records : {s['n_total']}  (mapped {s['n_mapped']}, "
              f"unmapped {s['n_unmapped']})")
        print(f"  Unique nodes  : {s['n_mapped_nodes']}")
        print(f"  Metabolomics  : {s['n_metabolomics']}")
        print(f"  Transcriptomics: {s['n_transcriptomics']}")
        print(f"  Proteomics    : {s['n_proteomics']}")

    def __len__(self) -> int:
        return len(self._records)

    def __repr__(self) -> str:
        s = self.summary()
        return (f"SampleContext(sample_id={self.sample_id!r}, "
                f"n_mapped={s['n_mapped']}, n_nodes={s['n_mapped_nodes']})")
