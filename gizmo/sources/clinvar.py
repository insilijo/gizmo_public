"""
ClinVar pathogenic variant loader.

License: NCBI ClinVar data is in the public domain (US government).
Source:  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

Downloads the ClinVar variant_summary.txt.gz (tab-delimited, ~200 MB compressed)
and extracts Pathogenic / Likely pathogenic variants, creating VariantNode and
VariantEdge (variant → gene) objects.

Only human (GRCh38/GRCh37) variants with a gene assignment and an rsID or
ClinVar variation ID are included.  Star-1 review status or higher is required
by default to limit noise.
"""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path
from urllib.request import urlretrieve

from gizmo.schema import VariantEdge, VariantNode

log = logging.getLogger(__name__)

_CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
    "variant_summary.txt.gz"
)

# ClinVar clinical significance values to keep
_PATHOGENIC = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
}

# Minimum review status to accept (number of stars ≥ this)
_MIN_STARS = 1

_STAR_MAP = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
}

# Variant consequence aliases → short labels
_CONSEQUENCE_MAP = {
    "single nucleotide variant":   "snv",
    "deletion":                    "deletion",
    "duplication":                 "duplication",
    "insertion":                   "insertion",
    "indel":                       "indel",
    "inversion":                   "inversion",
    "copy number gain":            "cnv_gain",
    "copy number loss":            "cnv_loss",
    "microsatellite":              "microsatellite",
    "protein only":                "protein_only",
}


class ClinVarClient:
    """
    Download and parse ClinVar variant_summary.txt.gz.

    Usage::

        client = ClinVarClient(cache_dir="data/raw/clinvar")
        client.download()
        nodes, edges = client.load_pathogenic(min_stars=1)
    """

    def __init__(self, cache_dir: str | Path = "data/raw/clinvar") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    @property
    def summary_path(self) -> Path:
        return self.cache_dir / "variant_summary.txt.gz"

    def download(self, force: bool = False) -> None:
        if self.summary_path.exists() and not force:
            log.info("ClinVar cache hit: %s", self.summary_path.name)
            return
        log.info("Downloading ClinVar variant_summary (~200 MB) → %s", self.summary_path)
        urlretrieve(_CLINVAR_URL, self.summary_path)

    def load_pathogenic(
        self,
        min_stars: int = _MIN_STARS,
        assembly: str = "GRCh38",
    ) -> tuple[list[VariantNode], list[VariantEdge]]:
        """
        Parse ClinVar and return (VariantNode list, VariantEdge list).

        Parameters
        ----------
        min_stars  : minimum ClinVar review-star rating (0–4)
        assembly   : genome assembly filter ("GRCh38" | "GRCh37" | None = all)

        Returns
        -------
        (nodes, edges) where each edge is variant → gene.
        """
        if not self.summary_path.exists():
            raise FileNotFoundError(
                f"ClinVar summary not found at {self.summary_path}. Run .download() first."
            )

        nodes: list[VariantNode] = []
        edges: list[VariantEdge] = []
        seen: set[str] = set()

        log.info("ClinVar: parsing variant_summary (this takes a minute) …")
        n_rows = 0
        with gzip.open(self.summary_path, "rt", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                n_rows += 1
                if n_rows % 200_000 == 0:
                    log.info("  ClinVar: %dk rows scanned, %d variants kept so far …",
                             n_rows // 1000, len(nodes))
                # Assembly filter
                if assembly and row.get("Assembly", "") != assembly:
                    continue

                # Clinical significance filter
                clnsig = row.get("ClinicalSignificance", "")
                if not any(p in clnsig for p in _PATHOGENIC):
                    continue

                # Review status filter
                review = row.get("ReviewStatus", "").lower()
                stars = _STAR_MAP.get(review, 0)
                if stars < min_stars:
                    continue

                # Gene assignment
                gene_sym = row.get("GeneSymbol", "").strip()
                if not gene_sym or gene_sym in {"-", "na", "N/A"}:
                    continue

                # Node ID: prefer rsID, fall back to VariationID
                rsid_raw = row.get("RS# (dbSNP)", "").strip()
                var_id   = row.get("VariationID", "").strip()
                if rsid_raw and rsid_raw not in {"-1", "-", ""}:
                    node_id = f"rs:{rsid_raw}"
                    rsid    = f"rs{rsid_raw}"
                else:
                    node_id = f"ClinVar:{var_id}"
                    rsid    = None

                if node_id in seen:
                    continue
                seen.add(node_id)

                consequence = _CONSEQUENCE_MAP.get(
                    row.get("Type", "").lower(), row.get("Type", "").lower()
                )

                node = VariantNode(
                    node_id=node_id,
                    clinvar_id=var_id or None,
                    rsid=rsid,
                    gene_symbol=gene_sym,
                    gene_id=f"symbol:{gene_sym}",
                    consequence=consequence or None,
                    clinical_significance=clnsig,
                    condition=row.get("PhenotypeList", "").strip() or None,
                    review_status=review or None,
                )
                nodes.append(node)

                # variant → gene edge
                edges.append(VariantEdge(
                    source=node_id,
                    target=f"symbol:{gene_sym}",
                    edge_type="variant_gene",
                    consequence=consequence or None,
                    source_db="clinvar",
                ))

        log.info(
            "ClinVar: loaded %d pathogenic variant nodes, %d variant-gene edges.",
            len(nodes), len(edges),
        )
        return nodes, edges
