"""
IUPHAR / Guide to Pharmacology (GtoPdb) source (CC BY-SA 4.0).

https://www.guidetopharmacology.org/

Loads bulk-download TSVs (interactions, ligands, targets) and produces
small-molecule -> sensor-protein binding edges for the GIZMO substrate.

Coverage:
  - GPCR endogenous ligands + drug agonists/antagonists/modulators
  - Nuclear receptor ligands (VDR, RAR, RXR, PPAR, ESR, AR, GR, MR, TR)
  - Xenobiotic sensors (AhR, PXR/NR1I2, CAR/NR1I3, FXR/NR1H4)
  - Adenosine receptors (ADORA1/2A/2B/3) — required for MTX-class anti-inflammatory bridge
  - Ion channels (NMDA, AMPA, GABA, voltage-gated, K+, TRP)

Usage::

    from gizmo.sources.gtopdb import GtoPdbClient
    client = GtoPdbClient()                                    # auto-download
    client = GtoPdbClient(data_dir='data/raw/gtopdb')          # use staged files

    edges = client.binding_edges(
        gene_symbols={'IL6R', 'ADORA2A', 'VDR'},        # filter by target gene
        endogenous_only=False,                           # include drug interactions
        approved_only=False,                             # include non-approved
        species='human',
    )
    # -> list[BindingEdge(ligand_name, ligand_pubchem_sid, target_gene_symbol,
    #                     action, endogenous, approved, affinity_pki)]

The class also exposes resolve_to_substrate(geom, mg) to map ligand and
target identifiers into the substrate's node id space; this is what the
bundle builder uses.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
import requests

log = logging.getLogger(__name__)

_BASE = 'https://www.guidetopharmacology.org/DATA'


@dataclass
class BindingEdge:
    """Single small-molecule -> sensor-protein binding annotation."""

    ligand_name: str
    ligand_pubchem_sid: Optional[str]
    target_gene_symbol: str
    target_uniprot: Optional[str]
    action: Optional[str]         # Agonist / Antagonist / Inhibitor / etc.
    endogenous: bool              # True if ligand is endogenous
    approved: bool                # True if ligand is an approved drug
    affinity_pki_median: Optional[float]  # median pKi from GtoPdb
    species: str = 'human'


class GtoPdbClient:
    """
    Thin loader for GtoPdb bulk-download TSVs.

    Parameters
    ----------
    data_dir : if given, loads from staged files; otherwise auto-downloads
               to a temp directory and caches in memory
    """

    def __init__(self, data_dir: Optional[str | Path] = None) -> None:
        self._data_dir = Path(data_dir) if data_dir else None
        self._interactions: Optional[pd.DataFrame] = None
        self._ligands: Optional[pd.DataFrame] = None
        self._targets: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    # Loaders
    # ------------------------------------------------------------------

    def _load_tsv(self, filename: str) -> pd.DataFrame:
        """Load a GtoPdb TSV, handling the comment-prefixed first line."""
        if self._data_dir and (self._data_dir / filename).exists():
            path = self._data_dir / filename
        else:
            path = self._download(filename)
        return pd.read_csv(
            path, sep='\t', skiprows=1, quotechar='"', low_memory=False
        )

    def _download(self, filename: str) -> Path:
        url = f'{_BASE}/{filename}'
        log.info('downloading %s', url)
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        # write to a temp file under data_dir (or system temp)
        from tempfile import gettempdir
        out_dir = self._data_dir or Path(gettempdir()) / 'gtopdb'
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / filename
        path.write_bytes(resp.content)
        return path

    @property
    def interactions(self) -> pd.DataFrame:
        if self._interactions is None:
            self._interactions = self._load_tsv('interactions.tsv')
        return self._interactions

    @property
    def ligands(self) -> pd.DataFrame:
        if self._ligands is None:
            self._ligands = self._load_tsv('ligands.tsv')
        return self._ligands

    @property
    def targets(self) -> pd.DataFrame:
        if self._targets is None:
            self._targets = self._load_tsv('targets_and_families.tsv')
        return self._targets

    # ------------------------------------------------------------------
    # Query
    # ------------------------------------------------------------------

    def binding_edges(
        self,
        gene_symbols: Optional[Iterable[str]] = None,
        endogenous_only: bool = False,
        approved_only: bool = False,
        species: str = 'human',
    ) -> list[BindingEdge]:
        """
        Return BindingEdge records for the given filter.

        Parameters
        ----------
        gene_symbols : if given, restrict to interactions targeting these
                        gene symbols (None = all human-targets)
        endogenous_only : keep only interactions with Endogenous=True
        approved_only : keep only interactions with Approved=True
        species : species filter for target (default human)
        """
        df = self.interactions
        df = df[df['Target Species'].astype(str).str.lower() == species.lower()]
        if endogenous_only:
            df = df[df['Endogenous'].astype(str).str.lower() == 'true']
        if approved_only:
            df = df[df['Approved'].astype(str).str.lower() == 'true']
        if gene_symbols is not None:
            gene_set = {str(g).upper() for g in gene_symbols}
            df = df[df['Target Gene Symbol'].astype(str)
                     .str.upper().isin(gene_set)]
        edges = []
        for _, row in df.iterrows():
            target_sym = str(row.get('Target Gene Symbol', '') or '').strip()
            if not target_sym: continue
            # GtoPdb concatenates multi-subunit targets with '|'; split
            for t in target_sym.split('|'):
                t = t.strip()
                if not t: continue
                edges.append(BindingEdge(
                    ligand_name=str(row.get('Ligand', '') or '').strip(),
                    ligand_pubchem_sid=(str(row.get('Ligand PubChem SID', ''))
                                          .strip() or None),
                    target_gene_symbol=t,
                    target_uniprot=(str(row.get('Target UniProt ID', ''))
                                      .strip() or None),
                    action=(str(row.get('Action', '')).strip() or None),
                    endogenous=str(row.get('Endogenous', '')).lower() == 'true',
                    approved=str(row.get('Approved', '')).lower() == 'true',
                    affinity_pki_median=row.get('Affinity Median')
                                          if not pd.isna(row.get('Affinity Median'))
                                          else None,
                    species=species.lower(),
                ))
                break  # one row per interaction, not per subunit
        return edges

    # ------------------------------------------------------------------
    # Substrate resolution
    # ------------------------------------------------------------------

    def resolve_to_substrate(
        self,
        geom,
        mg,
        edges: Optional[list[BindingEdge]] = None,
    ) -> list[tuple[str, str, BindingEdge]]:
        """
        Map each BindingEdge to a (ligand_nid, target_nid, edge) triple,
        using substrate gene-symbol nodes for targets and metabolite-name
        match for ligands.

        Returns only edges where BOTH ends resolve to a substrate node id
        present in geom.nid_idx.
        """
        if edges is None:
            edges = self.binding_edges()

        # Build name -> metabolite-nid map (lowercased)
        name_to_nid = {}
        for nid, attrs in mg.graph.nodes.items():
            if attrs.get('node_type') != 'metabolite': continue
            name = (attrs.get('name') or '').lower()
            if name and len(name) > 2:
                name_to_nid[name] = nid
            for syn in (attrs.get('synonyms') or []):
                if isinstance(syn, str) and len(syn) > 2:
                    name_to_nid[syn.lower()] = nid

        resolved = []
        for edge in edges:
            # target: prefer 'symbol:XYZ' form
            target_nid = None
            for cand in [f'symbol:{edge.target_gene_symbol}',
                          edge.target_gene_symbol]:
                if cand in geom.nid_idx:
                    target_nid = cand; break
            if target_nid is None: continue

            ligand_nid = name_to_nid.get(edge.ligand_name.lower())
            if ligand_nid is None: continue
            if ligand_nid not in geom.nid_idx: continue

            resolved.append((ligand_nid, target_nid, edge))
        return resolved


__all__ = ['GtoPdbClient', 'BindingEdge']
