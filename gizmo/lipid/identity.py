"""LipidIdentity — structured representation of a lipid species annotation.

Maps onto the LIPID MAPS resolution hierarchy (Levels 1–5). Used by the
Refmet parser, the lipid sub-graph builder, and the variable-resolution
observation model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, IntEnum
from typing import Optional


class LipidCategory(str, Enum):
    """LIPID MAPS top-level category (Level 1)."""

    GP = "GP"   # Glycerophospholipids
    GL = "GL"   # Glycerolipids
    SP = "SP"   # Sphingolipids
    ST = "ST"   # Sterol lipids
    FA = "FA"   # Fatty acyls
    PR = "PR"   # Prenols
    SL = "SL"   # Saccharolipids
    PK = "PK"   # Polyketides
    UNKNOWN = "UNKNOWN"


class ResolutionLevel(IntEnum):
    """LIPID MAPS resolution hierarchy.

    Higher = more specific.
    """

    CATEGORY = 1        # e.g., Glycerophospholipid
    MAIN_CLASS = 2      # e.g., PC
    SUBCLASS = 3        # e.g., diacyl-PC vs ether-PC vs lyso-PC
    SUM_COMPOSITION = 4  # "direct length" — e.g., PC 36:2
    MOLECULAR_SPECIES = 5  # "specific" — e.g., PC 18:1/18:1


@dataclass(frozen=True)
class ChainSpec:
    """A single acyl/alkyl chain.

    For ether-linkage and other non-acyl bonds, set ``linkage``.
    """

    carbons: int
    double_bonds: int
    linkage: str = "acyl"  # "acyl", "ether" (O-), "vinyl-ether" (P-, plasmalogen)
    hydroxyls: int = 0


@dataclass(frozen=True)
class LipidIdentity:
    """Structured identity of a lipid annotation.

    Built by parsers (e.g. Refmet); consumed by graph builders and the
    observation model.

    Fields:
        raw_name           — original annotation string
        sub_class          — short sub-class identifier (e.g. "PC", "Cer", "GD1")
        category           — LIPID MAPS top-level category
        lmid_prefix        — LIPID MAPS sub-class LMID (8-char + 4 zeros), if known
        resolution_level   — how specific the annotation is (3/4/5)

        # Sum-composition fields (Level 4+)
        total_carbons      — total acyl/alkyl carbons (None for named compounds)
        total_double_bonds — total acyl/alkyl double bonds

        # Molecular-species fields (Level 5)
        chains             — explicit per-chain composition (None at Level 4)

        # Modifications
        ether_linked       — at least one chain is ether-linked (O-)
        plasmalogen        — at least one chain is vinyl-ether (P-)
        is_lyso            — single-chain glycerolipid (LPC, LPE, ...)

        # Sphingolipid-specific (sphingoid base + N-acyl)
        sphingoid_carbons  — sphingoid base carbon count
        sphingoid_db       — sphingoid base double bonds
        sphingoid_oxygens  — number of OHs on sphingoid (1, 2, 3 — phyto)
        deoxy              — 1-deoxy variant (1-DeoxyCer)
        dihydro            — sphinganine-based (dh prefix: dhCer, dhSM, ...)

        # Other modifications
        hydroxyl_count     — explicit ;OH count on whole species (e.g., CAR ...;OH)
        named_compound     — for things like "Cholesterol", "Sphingosine"
                             (no chain composition; the name IS the identity)

        # Resolution metadata
        is_resolved        — at least one of (lmid_prefix, named_compound) is set
    """

    raw_name: str
    sub_class: str
    category: LipidCategory = LipidCategory.UNKNOWN
    lmid_prefix: Optional[str] = None
    resolution_level: ResolutionLevel = ResolutionLevel.SUBCLASS

    # Sum-composition
    total_carbons: Optional[int] = None
    total_double_bonds: Optional[int] = None

    # Molecular-species
    chains: tuple[ChainSpec, ...] = field(default_factory=tuple)

    # Modifications
    ether_linked: bool = False
    plasmalogen: bool = False
    is_lyso: bool = False

    # Sphingolipid-specific
    sphingoid_carbons: Optional[int] = None
    sphingoid_db: Optional[int] = None
    sphingoid_oxygens: Optional[int] = None
    deoxy: bool = False
    dihydro: bool = False

    # Other
    hydroxyl_count: int = 0
    named_compound: bool = False

    # Annotation uncertainty (Refmet "P-X:Y or O-X:(Y+1)" mass-equivalent
    # ambiguity). When True, the species could be either plasmalogen or
    # ether-linked at the reported abundance — the model should treat
    # this as linkage-type uncertainty.
    ambiguous_linkage: bool = False

    @property
    def is_resolved(self) -> bool:
        """True if we mapped this annotation to a known sub-class or compound."""
        return self.lmid_prefix is not None or self.named_compound

    @property
    def is_sphingolipid(self) -> bool:
        return self.category == LipidCategory.SP

    def saturation_class(self) -> str:
        """Rough saturation pool assignment from total double bonds.

        Returns one of: SFA, MUFA, PUFA, UNKNOWN. ω-3/ω-6 split requires
        molecular-species resolution and is deferred to a separate helper.
        """
        if self.total_double_bonds is None:
            return "UNKNOWN"
        if self.total_double_bonds == 0:
            return "SFA"
        if self.total_double_bonds == 1:
            return "MUFA"
        return "PUFA"

    def chain_length_bin(self) -> str:
        """Rough chain-length pool assignment from total carbons.

        Bins follow common lipidomics convention. Sphingolipids include
        sphingoid + N-acyl carbons in total_carbons.

        Returns: short (<14), medium (14-18), long (19-22),
        very-long (≥23), unknown (None).
        """
        if self.total_carbons is None:
            return "UNKNOWN"
        c = self.total_carbons
        if c < 14:
            return "short"
        if c <= 18:
            return "medium"
        if c <= 22:
            return "long"
        return "very-long"
