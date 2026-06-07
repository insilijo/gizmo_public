"""LIPID MAPS sub-class taxonomy.

Maps short sub-class identifiers (as used by Refmet / BioPAN / common
lipidomics nomenclature) to:
- LIPID MAPS LMID prefix (8-char sub-class anchor)
- LIPID MAPS top-level category
- lyso-form flag

Coverage strategy (per ``docs/biopan_benchmark_inventory.md``):
- **Floor:** the 39 sub-classes in BioPAN's ``lipid_nodes.csv``.
- **Extension:** sub-classes BioPAN omits but the Ando PR000713
  benchmark cohort measures at species level (HexCer, Hex2Cer,
  SHexCer, gangliosides, 1-DeoxyCer, CE, Cholesterol, CAR, named FFAs,
  acyl-CoAs, CoQ homologs, dolichols).

LMID prefixes are sourced from the LIPID MAPS Structure Database (LMSD)
sub-class anchors. They use the form ``LM<CAT><MAIN><SUB>0000`` where
the trailing four zeros indicate sub-class level (Level 3).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from gizmo.lipid.identity import LipidCategory


@dataclass(frozen=True)
class SubclassEntry:
    """One row in the sub-class taxonomy table."""

    sub_class: str
    category: LipidCategory
    lmid_prefix: Optional[str]
    is_lyso: bool = False
    in_biopan: bool = False
    notes: str = ""


# ---------------------------------------------------------------------------
# Glycerophospholipids (GP) — BioPAN floor + lyso forms
# ---------------------------------------------------------------------------
_GP: list[SubclassEntry] = [
    SubclassEntry("PC",   LipidCategory.GP, "LMGP01010000", in_biopan=True),
    SubclassEntry("LPC",  LipidCategory.GP, "LMGP01050000", is_lyso=True, in_biopan=True),
    SubclassEntry("PE",   LipidCategory.GP, "LMGP02010000", in_biopan=True),
    SubclassEntry("LPE",  LipidCategory.GP, "LMGP02050000", is_lyso=True, in_biopan=True),
    SubclassEntry("PS",   LipidCategory.GP, "LMGP03010000", in_biopan=True),
    SubclassEntry("LPS",  LipidCategory.GP, "LMGP03050000", is_lyso=True, in_biopan=True),
    SubclassEntry("PG",   LipidCategory.GP, "LMGP04010000", in_biopan=True),
    SubclassEntry("LPG",  LipidCategory.GP, "LMGP04050000", is_lyso=True, in_biopan=True),
    SubclassEntry("PI",   LipidCategory.GP, "LMGP06010000", in_biopan=True),
    SubclassEntry("LPI",  LipidCategory.GP, "LMGP06050000", is_lyso=True, in_biopan=True),
    SubclassEntry("PIP",  LipidCategory.GP, "LMGP06020000", in_biopan=True,
                  notes="Phosphatidylinositol monophosphate"),
    SubclassEntry("PIP2", LipidCategory.GP, "LMGP06030000", in_biopan=True,
                  notes="Phosphatidylinositol bisphosphate"),
    SubclassEntry("PIP3", LipidCategory.GP, "LMGP06040000", in_biopan=True,
                  notes="Phosphatidylinositol trisphosphate"),
    SubclassEntry("PA",   LipidCategory.GP, "LMGP10010000", in_biopan=True),
    SubclassEntry("LPA",  LipidCategory.GP, "LMGP10050000", is_lyso=True, in_biopan=True),
    SubclassEntry("CL",   LipidCategory.GP, "LMGP12010000", in_biopan=True,
                  notes="Cardiolipin"),
]

# ---------------------------------------------------------------------------
# Glycerolipids (GL)
# ---------------------------------------------------------------------------
_GL: list[SubclassEntry] = [
    SubclassEntry("MG",   LipidCategory.GL, "LMGL01010000", in_biopan=True),
    SubclassEntry("DG",   LipidCategory.GL, "LMGL02010000", in_biopan=True),
    SubclassEntry("TG",   LipidCategory.GL, "LMGL03010000", in_biopan=True),
]

# ---------------------------------------------------------------------------
# Sphingolipids (SP) — BioPAN floor + Ando-extension breadth
# ---------------------------------------------------------------------------
_SP: list[SubclassEntry] = [
    # Sphingoid bases
    SubclassEntry("SPB",    LipidCategory.SP, "LMSP01010000", in_biopan=True,
                  notes="Sphingosine"),
    SubclassEntry("dhSPB",  LipidCategory.SP, "LMSP01020000", in_biopan=True,
                  notes="Sphinganine"),
    SubclassEntry("SPBP",   LipidCategory.SP, "LMSP01050000", in_biopan=True,
                  notes="Sphingosine-1-phosphate"),
    SubclassEntry("dhSPBP", LipidCategory.SP, "LMSP01050000", in_biopan=True,
                  notes="Sphinganine-1-phosphate"),

    # Ceramides + dihydro + 1-deoxy
    SubclassEntry("Cer",        LipidCategory.SP, "LMSP02010000", in_biopan=True),
    SubclassEntry("dhCer",      LipidCategory.SP, "LMSP02020000", in_biopan=True,
                  notes="Dihydroceramide"),
    SubclassEntry("Cer1P",      LipidCategory.SP, "LMSP02050000", in_biopan=True),
    SubclassEntry("CerP",       LipidCategory.SP, "LMSP02050000", in_biopan=False,
                  notes="Synonym for Cer1P (some Refmet outputs)"),
    SubclassEntry("1-DeoxyCer", LipidCategory.SP, "LMSP02080000", in_biopan=False,
                  notes="EXTENSION — Ando 2019 headline finding; missing from BioPAN"),

    # Sphingomyelin
    SubclassEntry("SM",     LipidCategory.SP, "LMSP03010000", in_biopan=True),
    SubclassEntry("dhSM",   LipidCategory.SP, "LMSP03010000", in_biopan=True,
                  notes="Dihydrosphingomyelin"),

    # Glycosphingolipids — extension beyond BioPAN
    SubclassEntry("HexCer",  LipidCategory.SP, "LMSP05010000", in_biopan=False,
                  notes="EXTENSION — hexosylceramides"),
    SubclassEntry("Hex2Cer", LipidCategory.SP, "LMSP05020000", in_biopan=False,
                  notes="EXTENSION — dihexosylceramides"),
    SubclassEntry("SHexCer", LipidCategory.SP, "LMSP0501AB00", in_biopan=False,
                  notes="EXTENSION — sulfatides (sulfated hexosylceramides)"),

    # Gangliosides — extension beyond BioPAN
    SubclassEntry("GM3",  LipidCategory.SP, "LMSP0601AA00", in_biopan=False,
                  notes="EXTENSION — monosialoganglioside"),
    SubclassEntry("GD1",  LipidCategory.SP, "LMSP0601AC00", in_biopan=False,
                  notes="EXTENSION — disialoganglioside"),
    SubclassEntry("GD2",  LipidCategory.SP, "LMSP0601AC00", in_biopan=False,
                  notes="EXTENSION — disialoganglioside"),
    SubclassEntry("GD3",  LipidCategory.SP, "LMSP0601AC00", in_biopan=False,
                  notes="EXTENSION — disialoganglioside"),
    SubclassEntry("GT1",  LipidCategory.SP, "LMSP0601AD00", in_biopan=False,
                  notes="EXTENSION — trisialoganglioside"),
    SubclassEntry("GQ1",  LipidCategory.SP, "LMSP0601AE00", in_biopan=False,
                  notes="EXTENSION — tetrasialoganglioside"),
]

# ---------------------------------------------------------------------------
# Sterol lipids (ST) — extension beyond BioPAN
# ---------------------------------------------------------------------------
_ST: list[SubclassEntry] = [
    SubclassEntry("Cholesterol", LipidCategory.ST, "LMST01010001", in_biopan=False,
                  notes="EXTENSION — sterol"),
    SubclassEntry("CE",          LipidCategory.ST, "LMST01020000", in_biopan=False,
                  notes="EXTENSION — cholesteryl esters"),
]

# ---------------------------------------------------------------------------
# Fatty acyls (FA) — extension beyond BioPAN
# ---------------------------------------------------------------------------
_FA: list[SubclassEntry] = [
    SubclassEntry("CAR",     LipidCategory.FA, "LMFA0707",   in_biopan=False,
                  notes="EXTENSION — acylcarnitines"),
    SubclassEntry("FFA",     LipidCategory.FA, "LMFA01",     in_biopan=False,
                  notes="EXTENSION — free fatty acids (named compounds)"),
    SubclassEntry("Acyl-CoA", LipidCategory.FA, "LMFA0707",  in_biopan=False,
                  notes="EXTENSION — acyl-CoA species"),
]

# ---------------------------------------------------------------------------
# Prenols (PR) — extension beyond BioPAN
# ---------------------------------------------------------------------------
_PR: list[SubclassEntry] = [
    SubclassEntry("Coenzyme Q", LipidCategory.PR, "LMPR02010001", in_biopan=False,
                  notes="EXTENSION — ubiquinone homologs"),
    SubclassEntry("Dolichol",   LipidCategory.PR, "LMPR03020001", in_biopan=False,
                  notes="EXTENSION — polyprenols"),
]


# ---------------------------------------------------------------------------
# Ether-linked variants (BioPAN's O- prefix nodes)
#
# Refmet renders these as "<class> O-<C>:<DB>" rather than as a separate
# sub-class string — so the parser handles ether-linkage as a flag on the
# base sub-class entry. These BioPAN-floor entries are kept here for
# reference / parity testing only; the parser does not look them up by
# string.
# ---------------------------------------------------------------------------
_ETHER_REFERENCE: list[SubclassEntry] = [
    SubclassEntry("O-PC",  LipidCategory.GP, "LMGP01020000", in_biopan=True),
    SubclassEntry("O-PE",  LipidCategory.GP, "LMGP02020000", in_biopan=True),
    SubclassEntry("O-LPC", LipidCategory.GP, "LMGP01060000", is_lyso=True, in_biopan=True),
    SubclassEntry("O-LPE", LipidCategory.GP, "LMGP02060000", is_lyso=True, in_biopan=True),
    SubclassEntry("O-LPA", LipidCategory.GP, "LMGP10060000", is_lyso=True, in_biopan=True),
    SubclassEntry("O-DG",  LipidCategory.GL, "LMGL02020000", in_biopan=True),
    # Plasmalogens (P-)
    SubclassEntry("P-PC",  LipidCategory.GP, "LMGP01030000", in_biopan=True),
    SubclassEntry("P-PE",  LipidCategory.GP, "LMGP02030000", in_biopan=True),
    SubclassEntry("P-LPC", LipidCategory.GP, "LMGP01070000", is_lyso=True, in_biopan=True),
    SubclassEntry("P-LPE", LipidCategory.GP, "LMGP02070000", is_lyso=True, in_biopan=True),
    SubclassEntry("P-LPA", LipidCategory.GP, "LMGP10070000", is_lyso=True, in_biopan=True),
]


_ALL_ENTRIES: list[SubclassEntry] = _GP + _GL + _SP + _ST + _FA + _PR + _ETHER_REFERENCE
_BY_NAME: dict[str, SubclassEntry] = {e.sub_class: e for e in _ALL_ENTRIES}


def get_subclass(name: str) -> Optional[SubclassEntry]:
    """Look up a sub-class entry by short name. Case-sensitive.

    Returns None if the sub-class is not in the taxonomy.
    """
    return _BY_NAME.get(name)


def iter_subclasses(
    *,
    category: Optional[LipidCategory] = None,
    biopan_only: bool = False,
) -> list[SubclassEntry]:
    """Iterate sub-class entries, optionally filtered."""
    out = list(_ALL_ENTRIES)
    if category is not None:
        out = [e for e in out if e.category == category]
    if biopan_only:
        out = [e for e in out if e.in_biopan]
    return out
