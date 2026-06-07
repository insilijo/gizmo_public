"""Refmet name parser.

Refmet is the standardized lipid nomenclature served by Metabolomics
Workbench and adopted by LIPID MAPS. The patterns this parser
handles cover the species observed in the BioPAN benchmark cohort
(PR000713, ~1,079 species across four mouse tissues) plus common
nomenclature for compounds BioPAN omits.

Pattern families handled:

1. Sum-composition glycerolipid / glycerophospholipid:
   ``PC 36:2``, ``DG 32:0``, ``TG 48:0``, ``LPC 14:0``, ``CL 70:6``, ``CE 18:2``
2. Ether-linked sum-composition:
   ``DG O-32:1``, ``PC O-32:0``
3. Plasmalogen sum-composition:
   ``PC P-32:1``
4. Sphingolipid sum form (sphingoid only):
   ``SM 14:0;O2``, ``Cer 18:1;O2``
5. Sphingolipid molecular form (sphingoid + N-acyl):
   ``Cer 18:1;O2/24:0``, ``HexCer 18:1;O2/16:0``,
   ``1-DeoxyCer 18:0;O/16:0``, ``GD1 18:1;O2/16:0``
6. Acylcarnitine variants:
   ``CAR 10:0``, ``CAR 10:0;OH``, ``CAR 10:1;OH``
7. Coenzyme Q homologs:  ``Coenzyme Q4`` ... ``Coenzyme Q10``
8. Dolichols:           ``Dolichol-16`` ... ``Dolichol-20``
9. Named sterols:       ``Cholesterol``, ``Cholesterol sulfate``
10. Named sphingoid bases: ``Sphingosine``, ``Sphinganine``
11. Named fatty acids (acid-suffix):
    ``Arachidonic acid``, ``Palmitic acid``, ``Hexacosatrienoic acid``
12. Acyl-CoA species:   ``Palmitoyl-CoA``, ``Arachidonyl-CoA``

Unrecognised names return a ``LipidIdentity`` with
``is_resolved == False`` so callers can route to manual curation.
"""

from __future__ import annotations

import re
from typing import Optional

from gizmo.lipid.identity import (
    ChainSpec,
    LipidCategory,
    LipidIdentity,
    ResolutionLevel,
)
from gizmo.lipid.taxonomy import SubclassEntry, get_subclass


# ---------------------------------------------------------------------------
# Regexes (compiled once)
# ---------------------------------------------------------------------------

# Sum-composition: "PC 36:2", "TG 48:0", "CL 70:6", "CE 18:2"
_SUM_COMP_RE = re.compile(
    r"^(?P<sub>[A-Za-z][A-Za-z0-9]*)\s+(?P<c>\d+):(?P<d>\d+)$"
)

# Ether-linked sum-composition: "DG O-32:1", "PC O-32:0"
_ETHER_RE = re.compile(
    r"^(?P<sub>[A-Za-z][A-Za-z0-9]*)\s+O-(?P<c>\d+):(?P<d>\d+)$"
)

# Plasmalogen: "PC P-32:1"
_PLASMA_RE = re.compile(
    r"^(?P<sub>[A-Za-z][A-Za-z0-9]*)\s+P-(?P<c>\d+):(?P<d>\d+)$"
)

# Sphingolipid sub-class anchors, ordered longest-first so the regex
# alternation matches the most specific token (e.g., "1-DeoxyCer" wins
# over "Cer", "Hex2Cer" wins over "HexCer").
_SP_SUBCLASSES: tuple[str, ...] = (
    "1-DeoxyCer",
    "Hex2Cer", "SHexCer", "HexCer",
    "dhCer", "Cer1P", "CerP", "Cer",
    "dhSM", "SM",
    "GM1", "GM2", "GM3",
    "GD1", "GD2", "GD3",
    "GT1", "GT2", "GT3",
    "GQ1",
    "dhSPBP", "dhSPB", "SPBP", "SPB",
)
_SP_ALT = "|".join(re.escape(s) for s in _SP_SUBCLASSES)

# Sphingolipid sum form (sphingoid only): "SM 14:0;O2", "Cer 18:1;O2"
_SP_SUM_RE = re.compile(
    rf"^(?P<sub>{_SP_ALT})\s+(?P<c>\d+):(?P<d>\d+);O(?P<o>\d*)$"
)

# Sphingolipid molecular form: "Cer 18:1;O2/24:0", "1-DeoxyCer 18:0;O/16:0".
# Also handles trailing FA hydroxyl modifier (";O" after the FA chain) seen
# in HexCer hydroxylated species: "HexCer 18:1;O2/20:1;O".
_SP_MOL_RE = re.compile(
    rf"^(?P<sub>{_SP_ALT})\s+(?P<sc>\d+):(?P<sd>\d+);O(?P<o>\d*)/"
    r"(?P<fc>\d+):(?P<fd>\d+)(?P<fa_oh>;O\d*)?$"
)

# HexCer alternative without the slash (sphingoid+FA collapsed into total
# C/DB), with optional trailing extra ";O" indicating one extra hydroxyl:
# "HexCer 24:2;O2;O", "HexCer 25:1;O2;O".
_SP_SUM_OH_RE = re.compile(
    rf"^(?P<sub>{_SP_ALT})\s+(?P<c>\d+):(?P<d>\d+);O(?P<o>\d*);O(?P<o2>\d*)$"
)

# Refmet ambiguous "P or O" annotation — mass-equivalent plasmalogen vs
# ether-linked. Examples:
#   "PC P-34:0 or PC O-34:1"
#   "LPE P-18:0 or LPE O-18:1"
# Parsed as ether-linked with ambiguous_linkage=True.
_AMBIG_P_OR_O_RE = re.compile(
    r"^(?P<sub>[A-Za-z][A-Za-z0-9]*)\s+P-(?P<c1>\d+):(?P<d1>\d+)\s+or\s+"
    r"(?P=sub)\s+O-(?P<c2>\d+):(?P<d2>\d+)$"
)

# Acylcarnitine: "CAR 10:0", "CAR 10:0;OH"
_CAR_RE = re.compile(
    r"^CAR\s+(?P<c>\d+):(?P<d>\d+)(?:;(?P<oh>OH))?$"
)

# Coenzyme Q: "Coenzyme Q10"
_COQ_RE = re.compile(r"^Coenzyme\s+Q(?P<n>\d+)$")

# Dolichols: "Dolichol-16"
_DOL_RE = re.compile(r"^Dolichol-(?P<n>\d+)$")

# Named acid: ends with " acid", or hydroxylated forms with leading "Hydroxy"
_ACID_RE = re.compile(r"^.+\s+acid$", re.IGNORECASE)

# Acyl-CoA: ends with "-CoA"
_COA_RE = re.compile(r"^.+-CoA$")

# Standalone named compounds we recognise explicitly
_NAMED_COMPOUNDS: dict[str, tuple[str, LipidCategory, str]] = {
    # name → (sub_class, category, lmid_prefix)
    "Cholesterol":             ("Cholesterol", LipidCategory.ST, "LMST01010001"),
    "Cholesterol sulfate":     ("Cholesterol", LipidCategory.ST, "LMST01010001"),
    "Sphingosine":             ("SPB",    LipidCategory.SP, "LMSP01010000"),
    "Sphinganine":             ("dhSPB",  LipidCategory.SP, "LMSP01020000"),
    "Sphingosine-1-phosphate": ("SPBP",   LipidCategory.SP, "LMSP01050000"),
    "Sphinganine-1-phosphate": ("dhSPBP", LipidCategory.SP, "LMSP01050000"),
    # Refmet sometimes uses space instead of hyphen
    "Sphingosine 1-phosphate": ("SPBP",   LipidCategory.SP, "LMSP01050000"),
    "Sphinganine 1-phosphate": ("dhSPBP", LipidCategory.SP, "LMSP01050000"),
}


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

class RefmetParser:
    """Parse Refmet lipid annotations into :class:`LipidIdentity`.

    Stateless; the class form is a convenience for adding instance-level
    config later (custom sub-class overrides, strictness flags, etc.).
    """

    def parse(self, name: str) -> LipidIdentity:
        return parse_refmet(name)


def parse_refmet(name: str) -> LipidIdentity:
    """Parse a Refmet name into a :class:`LipidIdentity`.

    Always returns a ``LipidIdentity`` (never raises). Unrecognised names
    yield ``is_resolved == False``.
    """
    raw = name.strip()
    # Refmet/LIPID MAPS marks tentative or unusual species with a trailing "*"
    # (e.g. odd-chain FAs, rare ceramide chains). Strip and parse normally;
    # the asterisk is a curation flag, not a structural difference.
    parse_target = raw.rstrip("*").rstrip()

    # 1. Named compounds (exact match)
    if parse_target in _NAMED_COMPOUNDS:
        sub, cat, lmid = _NAMED_COMPOUNDS[parse_target]
        return LipidIdentity(
            raw_name=raw,
            sub_class=sub,
            category=cat,
            lmid_prefix=lmid,
            resolution_level=ResolutionLevel.MOLECULAR_SPECIES,
            named_compound=True,
        )

    # 2. Ambiguous "P-X:Y or O-X:(Y+1)" annotation (must precede ETHER/PLASMA
    # because it would otherwise greedy-match neither). Refmet uses this
    # when a peak's mass is consistent with both plasmalogen and ether-linked
    # forms differing by one DB.
    if m := _AMBIG_P_OR_O_RE.match(parse_target):
        sub = m.group("sub")
        c = int(m.group("c2"))   # use the O- form's carbons (same as P-)
        d = int(m.group("d2"))   # use the O- form's double bonds
        entry = get_subclass(sub)
        return LipidIdentity(
            raw_name=raw,
            sub_class=sub,
            category=entry.category if entry else LipidCategory.UNKNOWN,
            lmid_prefix=entry.lmid_prefix if entry else None,
            resolution_level=ResolutionLevel.SUM_COMPOSITION,
            total_carbons=c,
            total_double_bonds=d,
            ether_linked=True,
            ambiguous_linkage=True,
            is_lyso=entry.is_lyso if entry else False,
        )

    # 3. Sphingolipid molecular form (must precede SP_SUM)
    if m := _SP_MOL_RE.match(parse_target):
        return _build_sp_molecular(raw, m)

    # 4. Sphingolipid sum form with extra hydroxyl ";O2;O"
    if m := _SP_SUM_OH_RE.match(parse_target):
        return _build_sp_sum(raw, m, extra_hydroxyl=True)

    # 5. Sphingolipid sum form
    if m := _SP_SUM_RE.match(parse_target):
        return _build_sp_sum(raw, m)

    # 6. Acylcarnitines (must precede generic sum-comp because CAR has hydroxyl variant)
    if m := _CAR_RE.match(parse_target):
        return _build_carnitine(raw, m)

    # 7. Ether-linked sum-comp
    if m := _ETHER_RE.match(parse_target):
        return _build_glycerolipid(raw, m, ether=True)

    # 8. Plasmalogen sum-comp
    if m := _PLASMA_RE.match(parse_target):
        return _build_glycerolipid(raw, m, plasmalogen=True)

    # 9. Coenzyme Q
    if m := _COQ_RE.match(parse_target):
        return LipidIdentity(
            raw_name=raw,
            sub_class="Coenzyme Q",
            category=LipidCategory.PR,
            lmid_prefix="LMPR02010001",
            resolution_level=ResolutionLevel.MOLECULAR_SPECIES,
            named_compound=True,
        )

    # 10. Dolichols
    if m := _DOL_RE.match(parse_target):
        return LipidIdentity(
            raw_name=raw,
            sub_class="Dolichol",
            category=LipidCategory.PR,
            lmid_prefix="LMPR03020001",
            resolution_level=ResolutionLevel.MOLECULAR_SPECIES,
            named_compound=True,
        )

    # 11. Acyl-CoA
    if _COA_RE.match(parse_target):
        return LipidIdentity(
            raw_name=raw,
            sub_class="Acyl-CoA",
            category=LipidCategory.FA,
            lmid_prefix="LMFA0707",
            resolution_level=ResolutionLevel.MOLECULAR_SPECIES,
            named_compound=True,
        )

    # 12. Named fatty acid (ends in " acid")
    if _ACID_RE.match(parse_target):
        return LipidIdentity(
            raw_name=raw,
            sub_class="FFA",
            category=LipidCategory.FA,
            lmid_prefix="LMFA01",
            resolution_level=ResolutionLevel.MOLECULAR_SPECIES,
            named_compound=True,
        )

    # 13. Generic sum-comp (last because acid/CoA/CAR have own patterns)
    if m := _SUM_COMP_RE.match(parse_target):
        return _build_glycerolipid(raw, m)

    # Unresolved
    return LipidIdentity(
        raw_name=raw,
        sub_class="UNKNOWN",
        category=LipidCategory.UNKNOWN,
        resolution_level=ResolutionLevel.SUBCLASS,
    )


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def _build_glycerolipid(
    raw: str,
    m: re.Match,
    *,
    ether: bool = False,
    plasmalogen: bool = False,
) -> LipidIdentity:
    sub = m.group("sub")
    c = int(m.group("c"))
    d = int(m.group("d"))
    entry = get_subclass(sub)
    category = entry.category if entry else LipidCategory.UNKNOWN
    return LipidIdentity(
        raw_name=raw,
        sub_class=sub,
        category=category,
        lmid_prefix=entry.lmid_prefix if entry else None,
        resolution_level=ResolutionLevel.SUM_COMPOSITION,
        total_carbons=c,
        total_double_bonds=d,
        ether_linked=ether,
        plasmalogen=plasmalogen,
        is_lyso=entry.is_lyso if entry else False,
    )


def _build_carnitine(raw: str, m: re.Match) -> LipidIdentity:
    c = int(m.group("c"))
    d = int(m.group("d"))
    oh = m.group("oh")
    entry = get_subclass("CAR")
    return LipidIdentity(
        raw_name=raw,
        sub_class="CAR",
        category=LipidCategory.FA,
        lmid_prefix=entry.lmid_prefix if entry else "LMFA0707",
        resolution_level=ResolutionLevel.SUM_COMPOSITION,
        total_carbons=c,
        total_double_bonds=d,
        hydroxyl_count=1 if oh else 0,
    )


def _resolve_sphingolipid_subclass(raw_sub: str) -> tuple[str, SubclassEntry | None, bool]:
    """Map raw sphingolipid sub-class string to canonical sub-class entry.

    Handles ``1-DeoxyCer`` (deoxy variant of Cer) and ``dhCer`` /
    ``dhSM`` (sphinganine-based, the "dh" prefix).

    Returns (canonical_sub_class, entry_or_None, deoxy_flag).
    """
    # dh* prefix (dihydro = sphinganine-based)
    # Note: the parser already treats dhCer / dhSM as their own entries,
    # so this is just a passthrough lookup; deoxy is checked separately.
    deoxy = raw_sub.lower().startswith("1-deoxy")
    entry = get_subclass(raw_sub)
    return raw_sub, entry, deoxy


def _build_sp_sum(
    raw: str,
    m: re.Match,
    *,
    extra_hydroxyl: bool = False,
) -> LipidIdentity:
    sub_raw = m.group("sub")
    sc = int(m.group("c"))
    sd = int(m.group("d"))
    so = int(m.group("o")) if m.group("o") else 1  # bare ";O" = 1 oxygen
    sub, entry, deoxy = _resolve_sphingolipid_subclass(sub_raw)
    dihydro = sub_raw.lower().startswith("dh")

    extra_oh = 0
    if extra_hydroxyl:
        # Second ";O<n>" group adds hydroxyls; default 1 if bare ";O".
        o2 = m.group("o2") if "o2" in m.groupdict() else None
        extra_oh = int(o2) if o2 else 1

    return LipidIdentity(
        raw_name=raw,
        sub_class=sub,
        category=entry.category if entry else LipidCategory.SP,
        lmid_prefix=entry.lmid_prefix if entry else None,
        resolution_level=ResolutionLevel.SUM_COMPOSITION,
        total_carbons=sc,
        total_double_bonds=sd,
        sphingoid_carbons=sc,
        sphingoid_db=sd,
        sphingoid_oxygens=so,
        deoxy=deoxy,
        dihydro=dihydro,
        hydroxyl_count=extra_oh,
    )


def _build_sp_molecular(raw: str, m: re.Match) -> LipidIdentity:
    sub_raw = m.group("sub")
    sc = int(m.group("sc"))
    sd = int(m.group("sd"))
    so = int(m.group("o")) if m.group("o") else 1
    fc = int(m.group("fc"))
    fd = int(m.group("fd"))
    sub, entry, deoxy = _resolve_sphingolipid_subclass(sub_raw)
    dihydro = sub_raw.lower().startswith("dh")

    # Optional ";O<n>" tail on the FA chain (HexCer hydroxylated species).
    fa_oh = m.group("fa_oh") if "fa_oh" in m.groupdict() else None
    fa_hydroxyls = 0
    if fa_oh:
        # fa_oh starts with ";O" plus optional digits
        digits = fa_oh[2:]
        fa_hydroxyls = int(digits) if digits else 1

    chains = (
        ChainSpec(carbons=sc, double_bonds=sd, linkage="sphingoid", hydroxyls=so),
        ChainSpec(carbons=fc, double_bonds=fd, linkage="acyl", hydroxyls=fa_hydroxyls),
    )

    return LipidIdentity(
        raw_name=raw,
        sub_class=sub,
        category=entry.category if entry else LipidCategory.SP,
        lmid_prefix=entry.lmid_prefix if entry else None,
        resolution_level=ResolutionLevel.MOLECULAR_SPECIES,
        total_carbons=sc + fc,
        total_double_bonds=sd + fd,
        chains=chains,
        sphingoid_carbons=sc,
        sphingoid_db=sd,
        sphingoid_oxygens=so,
        deoxy=deoxy,
        dihydro=dihydro,
        hydroxyl_count=fa_hydroxyls,
    )
