"""
GTEx tissue expression enrichment.

License: GTEx Portal data is available under dbGaP accession phs000424.
         Summary-level (median TPM per tissue) data is freely downloadable.

Downloads the GTEx v10 gene median TPM table and enriches GeneNode
attributes with {tissue: median_tpm} in the ``tissue_expression`` field.

URL: https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/
     GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz

The GCT file is ~8 MB compressed.  It is parsed lazily and cached.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from urllib.request import urlretrieve

log = logging.getLogger(__name__)

_GTEX_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
    "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
)

# Tissues to retain (short labels → GCT column name substrings).
# None = keep all ~54 tissues (larger memory footprint).
DEFAULT_TISSUES = None


class GTExClient:
    """
    Download and parse the GTEx median TPM table.

    Usage::

        client = GTExClient(cache_dir="data/raw/gtex")
        client.download()
        # Enrich the GizmoGraph in-place:
        n = client.enrich_graph(mg, min_tpm=1.0)
        print(f"Enriched {n} gene nodes with tissue expression data.")
    """

    def __init__(self, cache_dir: str | Path = "data/raw/gtex") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._tpm: dict[str, dict[str, float]] | None = None  # {ensg_id: {tissue: tpm}}

    @property
    def gct_path(self) -> Path:
        return self.cache_dir / "GTEx_gene_median_tpm.gct.gz"

    def download(self, force: bool = False) -> None:
        if self.gct_path.exists() and not force:
            log.info("GTEx cache hit: %s", self.gct_path.name)
            return
        log.info("Downloading GTEx median TPM table (~50 MB) → %s", self.gct_path)
        urlretrieve(_GTEX_URL, self.gct_path)

    def _load(self) -> dict[str, dict[str, float]]:
        if self._tpm is not None:
            return self._tpm

        if not self.gct_path.exists():
            raise FileNotFoundError(
                f"GTEx GCT not found at {self.gct_path}. Run .download() first."
            )

        log.info("Parsing GTEx median TPM table …")
        tpm: dict[str, dict[str, float]] = {}

        with gzip.open(self.gct_path, "rt", encoding="utf-8") as fh:
            # GCT format: 2 header lines, then column header, then data
            next(fh)  # version line
            next(fh)  # dimensions line
            header = next(fh).rstrip("\n").split("\t")
            # header[0] = "Name" (ENSG ID), header[1] = "Description" (symbol)
            # header[2:] = tissue names
            tissues = [_clean_tissue(t) for t in header[2:]]

            for line in fh:
                parts = line.rstrip("\n").split("\t")
                ensg_raw = parts[0]  # e.g. "ENSG00000223972.5"
                ensg = ensg_raw.split(".")[0]   # strip version suffix
                vals = parts[2:]
                if len(vals) != len(tissues):
                    continue
                expr: dict[str, float] = {}
                for tissue, v in zip(tissues, vals):
                    try:
                        f = float(v)
                    except ValueError:
                        continue
                    if f > 0:
                        expr[tissue] = round(f, 3)
                if expr:
                    tpm[ensg] = expr

        log.info("GTEx: loaded expression for %d genes across %d tissues.",
                 len(tpm), len(tissues))
        self._tpm = tpm
        return tpm

    def enrich_graph(self, mg, min_tpm: float = 1.0) -> int:
        """
        Enrich gene nodes in ``mg`` in-place with GTEx tissue expression.

        Only tissues where median TPM ≥ ``min_tpm`` are stored to limit
        storage footprint (0 = store all non-zero tissues).

        Returns the number of gene nodes enriched.
        """
        tpm_index = self._load()
        g = mg.graph
        n_enriched = 0

        for nid, attrs in list(g.nodes(data=True)):
            if attrs.get("node_type") != "gene":
                continue

            ensg = attrs.get("ensembl_id") or ""
            # Strip version suffix from stored ID if needed
            ensg_clean = ensg.split(".")[0]

            expr = tpm_index.get(ensg_clean) or tpm_index.get(ensg)
            if not expr:
                continue

            filtered = {t: v for t, v in expr.items() if v >= min_tpm}
            if filtered:
                g.nodes[nid]["tissue_expression"] = filtered
                n_enriched += 1

        log.info("GTEx: enriched %d gene nodes.", n_enriched)
        return n_enriched


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _clean_tissue(raw: str) -> str:
    """Normalise GTEx tissue column names to snake_case short labels."""
    return (
        raw.strip()
        .lower()
        .replace(" - ", "_")
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("-", "_")
    )
