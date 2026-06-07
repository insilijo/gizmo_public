"""
NCBI Gene Ortholog mapper.

License: NCBI data is in the public domain (US government).
Source:  https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz

The file maps (tax_id, gene_id) → (tax_id, gene_id) pairs for orthologous genes.
We combine it with gene_info.gz to resolve Entrez IDs to HGNC symbols/Ensembl IDs.

NCBI taxonomy IDs for common species:
  9606   Homo sapiens
  10090  Mus musculus
  10116  Rattus norvegicus
  7955   Danio rerio (zebrafish)
  9913   Bos taurus
  9823   Sus scrofa
"""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path
from urllib.request import urlretrieve

log = logging.getLogger(__name__)

_ORTHOLOGS_URL  = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz"
_GENE_INFO_URL  = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"

# Common species taxonomy IDs
SPECIES_TAXIDS = {
    "Homo sapiens":       "9606",
    "Mus musculus":       "10090",
    "Rattus norvegicus":  "10116",
    "Danio rerio":        "7955",
    "Bos taurus":         "9913",
    "Sus scrofa":         "9823",
    "Caenorhabditis elegans": "6239",
    "Drosophila melanogaster": "7227",
}


class OrthologMapper:
    """
    Map gene identifiers between species using NCBI ortholog data.

    Usage::

        mapper = OrthologMapper(cache_dir="data/raw/ncbi")
        mapper.download()
        mapper.build(from_taxid="9606", to_taxid="10090")

        # Map a human gene symbol to its mouse ortholog symbol
        mouse_sym = mapper.map_symbol("BRCA2")   # → "Brca2" (or None)
    """

    def __init__(self, cache_dir: str | Path = "data/raw/ncbi") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        # {from_entrez: to_entrez}
        self._entrez_map: dict[str, str] = {}
        # {taxid: {entrez_id: symbol}}
        self._id_to_sym: dict[str, dict[str, str]] = {}
        # {taxid: {symbol_lower: entrez_id}}
        self._sym_to_id: dict[str, dict[str, str]] = {}
        self._from_taxid: str | None = None
        self._to_taxid: str | None = None

    @property
    def orthologs_path(self) -> Path:
        return self.cache_dir / "gene_orthologs.gz"

    @property
    def gene_info_path(self) -> Path:
        return self.cache_dir / "gene_info.gz"

    def download(self, force: bool = False) -> None:
        for url, path in [
            (_ORTHOLOGS_URL, self.orthologs_path),
            (_GENE_INFO_URL, self.gene_info_path),
        ]:
            if path.exists() and not force:
                log.info("NCBI cache hit: %s", path.name)
                continue
            log.info("Downloading %s (~%s) → %s",
                     path.name,
                     "50 MB" if "orthologs" in path.name else "400 MB",
                     path)
            urlretrieve(url, path)

    def build(
        self,
        from_species: str = "Homo sapiens",
        to_species: str = "Mus musculus",
    ) -> "OrthologMapper":
        """
        Build ortholog index for from_species → to_species.

        Parameters
        ----------
        from_species : source species name or NCBI taxonomy ID (string)
        to_species   : target species name or NCBI taxonomy ID (string)
        """
        from_taxid = SPECIES_TAXIDS.get(from_species, from_species)
        to_taxid   = SPECIES_TAXIDS.get(to_species,   to_species)
        self._from_taxid = from_taxid
        self._to_taxid   = to_taxid

        # Build symbol ↔ entrez index for both species
        for taxid in (from_taxid, to_taxid):
            self._id_to_sym[taxid] = {}
            self._sym_to_id[taxid] = {}

        log.info("Building gene symbol index from gene_info.gz …")
        with gzip.open(self.gene_info_path, "rt", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                taxid = row.get("#tax_id", "").strip()
                if taxid not in (from_taxid, to_taxid):
                    continue
                gene_id = row.get("GeneID", "").strip()
                symbol  = row.get("Symbol", "").strip()
                if gene_id and symbol:
                    self._id_to_sym[taxid][gene_id] = symbol
                    self._sym_to_id[taxid][symbol.lower()] = gene_id

        log.info(
            "Gene info: %d %s genes, %d %s genes.",
            len(self._id_to_sym[from_taxid]),
            from_species,
            len(self._id_to_sym[to_taxid]),
            to_species,
        )

        # Build entrez → entrez ortholog map
        log.info("Building ortholog index from gene_orthologs.gz …")
        n = 0
        with gzip.open(self.orthologs_path, "rt", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                tax1 = row.get("tax_id", "").strip()
                gid1 = row.get("GeneID", "").strip()
                tax2 = row.get("Other_tax_id", "").strip()
                gid2 = row.get("Other_GeneID", "").strip()

                if tax1 == from_taxid and tax2 == to_taxid:
                    self._entrez_map[gid1] = gid2
                    n += 1
                elif tax1 == to_taxid and tax2 == from_taxid:
                    self._entrez_map[gid2] = gid1
                    n += 1

        log.info("Ortholog map built: %d pairs.", n)
        return self

    def map_symbol(self, symbol: str) -> str | None:
        """
        Map a gene symbol from the from_species to the to_species.

        Returns the ortholog symbol, or None if not found.
        """
        if not self._from_taxid or not self._to_taxid:
            raise RuntimeError("Call .build() before mapping.")

        sym_lower = symbol.lower()
        from_id = self._sym_to_id.get(self._from_taxid, {}).get(sym_lower)
        if not from_id:
            return None
        to_id = self._entrez_map.get(from_id)
        if not to_id:
            return None
        return self._id_to_sym.get(self._to_taxid, {}).get(to_id)

    def map_batch(self, symbols: list[str]) -> dict[str, str | None]:
        """Map a list of gene symbols in bulk.  Returns {symbol: ortholog_symbol}."""
        return {s: self.map_symbol(s) for s in symbols}
