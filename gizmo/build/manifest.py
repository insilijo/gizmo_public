"""
GraphManifest — provenance record for a GIZMO graph bundle.

Saved alongside each graph as ``manifest.json`` so that any downstream
consumer can trace exactly how the graph was built.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


@dataclass
class SourceRecord:
    """
    Provenance for one data source used during a build.
    """
    name: str                             # e.g. "reactome", "metanetx"
    version: Optional[str] = None         # release tag / date if known
    url: Optional[str] = None             # canonical homepage or download URL
    license: str = "unknown"
    accessed: str = field(
        default_factory=lambda: datetime.now(timezone.utc).date().isoformat()
    )
    n_records: Optional[int] = None       # rows / entries ingested


@dataclass
class GraphManifest:
    """
    Full provenance record for a built GIZMO graph.

    Stored as ``manifest.json`` inside the graph bundle directory.
    """

    # Identity
    graph_name: str
    built_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )
    gizmo_version: str = ""

    # Build configuration
    build_params: dict[str, Any] = field(default_factory=dict)

    # Sources
    sources: list[SourceRecord] = field(default_factory=list)

    # Graph statistics (populated after build)
    node_counts: dict[str, int] = field(default_factory=dict)
    edge_count: int = 0
    pathway_count: int = 0

    # QC summary (subset of ReadinessReport fields)
    qc: dict[str, Any] = field(default_factory=dict)

    # Free-form notes
    notes: str = ""

    # ------------------------------------------------------------------
    # Convenience constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_graph(
        cls,
        mg,
        graph_name: str,
        build_params: dict | None = None,
        sources: list[SourceRecord] | None = None,
        notes: str = "",
    ) -> "GraphManifest":
        """
        Build a manifest from a completed GizmoGraph.
        Runs `assess_readiness` to populate the QC section.
        """
        from gizmo import __version__

        g = mg.graph
        node_counts = {}
        for nid, attrs in g.nodes(data=True):
            t = attrs.get("node_type", "unknown")
            node_counts[t] = node_counts.get(t, 0) + 1

        # Count unique pathway stIDs across all reaction nodes
        pw_set: set[str] = set()
        for nid, attrs in g.nodes(data=True):
            pw_set.update(attrs.get("pathways") or [])

        # QC summary
        qc: dict[str, Any] = {}
        try:
            from gizmo.analysis.qc import assess_readiness
            r = assess_readiness(mg)
            qc = {
                "is_fba_ready":                r.is_fba_ready,
                "n_dead_end_metabolites":       r.n_dead_end_metabolites,
                "n_orphan_reactions":           r.n_orphan_reactions,
                "n_weakly_connected_components": r.n_weakly_connected_components,
                "reactions_with_ec_fraction":   round(r.reactions_with_ec_fraction, 3),
                "reactions_with_gene_fraction": round(r.reactions_with_gene_fraction, 3),
                "metabolites_with_chebi_fraction": round(r.metabolites_with_chebi_fraction, 3),
                "metabolon_chebi_coverage":     round(r.metabolon_chebi_coverage, 3),
            }
        except Exception:
            pass

        return cls(
            graph_name=graph_name,
            gizmo_version=__version__,
            build_params=build_params or {},
            sources=sources or [],
            node_counts=node_counts,
            edge_count=g.number_of_edges(),
            pathway_count=len(pw_set),
            qc=qc,
            notes=notes,
        )

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def save(self, path: str | Path) -> None:
        """Write manifest to ``path`` as pretty-printed JSON."""
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        d = asdict(self)
        p.write_text(json.dumps(d, indent=2, default=str))

    @classmethod
    def load(cls, path: str | Path) -> "GraphManifest":
        """Load a manifest from a JSON file."""
        raw = json.loads(Path(path).read_text())
        sources = [SourceRecord(**s) for s in raw.pop("sources", [])]
        return cls(sources=sources, **raw)

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def print_summary(self) -> None:
        lines = [
            f"Graph:   {self.graph_name}",
            f"Built:   {self.built_at}  (gizmo {self.gizmo_version})",
            f"Nodes:   {sum(self.node_counts.values())}  "
            + "  ".join(f"{t}={n}" for t, n in sorted(self.node_counts.items())),
            f"Edges:   {self.edge_count}",
            f"Pathways:{self.pathway_count}",
            "Sources:",
        ]
        for s in self.sources:
            lines.append(f"  {s.name:<18} {s.license:<15} accessed {s.accessed}"
                         + (f"  ({s.n_records} records)" if s.n_records else ""))
        if self.qc:
            lines.append("QC:")
            for k, v in self.qc.items():
                lines.append(f"  {k}: {v}")
        print("\n".join(lines))
