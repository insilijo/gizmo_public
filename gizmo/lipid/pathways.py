"""Curated lipid sub-class transformations for Layer A.

Hand-curated to match BioPAN's manual-curation style (see Gaud et al.
F1000Research 2021), covering BioPAN's 39-node floor and extending into
the sub-classes BioPAN omits but the Ando PR000713 benchmark cohort
measures (1-DeoxyCer, HexCer, gangliosides, CE, Cholesterol, CAR).

Each transformation records:
- substrate sub-class (matches a ``SubclassEntry.sub_class``)
- product sub-class
- the catalysing enzyme(s) by HGNC gene symbol
- pathway tag (Kennedy / Lands cycle / sphingolipid_biosynthesis / etc.)
- reversibility (forward only / both directions)
- EC class hint where unambiguous

References:
- BioPAN paper Table 1 (liver active reactions), Table 2 (FA reactions)
- KEGG / Reactome lipid biosynthesis pathways
- Standard lipid biochemistry texts (Vance & Vance, "Biochemistry of Lipids")

Coverage status vs BioPAN's claimed 94 reactions:
- Glycerophospholipid network (Kennedy + Lands): full parity intended
- Glycerolipid storage cycle: full parity
- Sphingolipid de novo + salvage: full parity + extension (1-DeoxyCer, HexCer)
- Cholesterol/CE: extension only (BioPAN omits)
- FA elongation/desaturation: ELOVL/SCD/FADS cascades enumerated
- Acylcarnitine shuttle: extension only (BioPAN omits)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class LipidTransformation:
    """One sub-class-level lipid transformation.

    Substrate / product strings must match a ``SubclassEntry.sub_class``
    in :mod:`gizmo.lipid.taxonomy`. Validity is enforced by tests.
    """

    substrate: str
    product: str
    enzymes: tuple[str, ...]
    pathway: str
    reversible: bool = False
    ec_class: Optional[str] = None
    notes: str = ""


# ---------------------------------------------------------------------------
# Kennedy pathway — glycerophospholipid de novo synthesis
# ---------------------------------------------------------------------------
_KENNEDY: list[LipidTransformation] = [
    # CDP-DAG branch: PA → CDP-DAG (not modeled at sub-class level — implicit)
    LipidTransformation("PA", "DG", ("LPIN1", "LPIN2", "LPIN3"),
                        pathway="Kennedy", ec_class="3.1.3.4",
                        notes="PAP1 — phosphatidic acid phosphatase"),
    LipidTransformation("DG", "PC", ("CHPT1", "CEPT1"),
                        pathway="Kennedy", ec_class="2.7.8.2"),
    LipidTransformation("DG", "PE", ("CEPT1", "SELENOI"),
                        pathway="Kennedy", ec_class="2.7.8.1",
                        notes="SELENOI is canonical EPT1"),
    LipidTransformation("LPA", "PA", ("AGPAT1", "AGPAT2", "AGPAT3", "AGPAT4", "AGPAT5"),
                        pathway="Kennedy", ec_class="2.3.1.51"),
    # PA → PG via PGS1 + PTPMT1
    LipidTransformation("PA", "PG", ("PGS1", "PTPMT1"),
                        pathway="Kennedy", ec_class="2.7.8.5",
                        notes="via CDP-DAG and phosphatidylglycerolphosphate intermediates"),
    LipidTransformation("PG", "CL", ("CRLS1",),
                        pathway="Kennedy", ec_class="2.7.8.41",
                        notes="cardiolipin synthase"),
    # PI synthesis
    LipidTransformation("PA", "PI", ("CDS1", "CDS2", "CDIPT"),
                        pathway="Kennedy",
                        notes="via CDP-DAG; CDIPT is PI synthase"),
    # Headgroup interconversion
    LipidTransformation("PE", "PC", ("PEMT",),
                        pathway="Kennedy", ec_class="2.1.1.17",
                        notes="PE methylation — major hepatic route"),
    LipidTransformation("PS", "PE", ("PISD",),
                        pathway="Kennedy", ec_class="4.1.1.65",
                        notes="PS decarboxylase"),
    LipidTransformation("PE", "PS", ("PTDSS2",),
                        pathway="Kennedy", ec_class="2.7.8.29",
                        notes="PS synthase 2"),
    LipidTransformation("PC", "PS", ("PTDSS1",),
                        pathway="Kennedy", ec_class="2.7.8.29",
                        notes="PS synthase 1"),
]


# ---------------------------------------------------------------------------
# Lands cycle — lyso-phospholipid acyl remodeling
# ---------------------------------------------------------------------------
_LANDS: list[LipidTransformation] = [
    # PLA2 hydrolysis: PX → LPX
    LipidTransformation("PC", "LPC", ("PLA2G1B", "PLA2G2A", "PLA2G2D",
                                       "PLA2G2E", "PLA2G2F", "PLA2G4A", "PLA2G4B",
                                       "PLA2G4C", "PLA2G4D", "PLA2G4E", "PLA2G4F",
                                       "PLA2G6"),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    LipidTransformation("PE", "LPE", ("PLA2G6", "PLA2G4A", "PLA2G4D", "PLA2G4F"),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    LipidTransformation("PS", "LPS", ("PLA2G2E", "PLA2G2A", "PLA2G2F"),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    LipidTransformation("PI", "LPI", ("PLA2G6",),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    LipidTransformation("PG", "LPG", ("PLA2G4D", "PLA2G4F"),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    LipidTransformation("PA", "LPA", ("PLA2G6",),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    # LPCAT-style re-acylation: LPX → PX
    LipidTransformation("LPC", "PC", ("LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4"),
                        pathway="Lands_cycle", ec_class="2.3.1.23"),
    LipidTransformation("LPE", "PE", ("LPEAT1", "MBOAT2", "LPCAT4"),
                        pathway="Lands_cycle", ec_class="2.3.1.n",
                        notes="LPEAT1 = aliases LPCAT3-MBOAT5; MBOAT2 main re-acylator"),
    LipidTransformation("LPS", "PS", ("MBOAT5",),
                        pathway="Lands_cycle"),
    LipidTransformation("LPA", "PA", ("AGPAT1", "AGPAT2", "AGPAT3", "AGPAT4"),
                        pathway="Lands_cycle", ec_class="2.3.1.51",
                        notes="AGPATs also serve as Lands-cycle LPA acylators"),
]


# ---------------------------------------------------------------------------
# Glycerolipid storage cycle — DG/TG/MG/CE
# ---------------------------------------------------------------------------
_GLYCEROLIPIDS: list[LipidTransformation] = [
    LipidTransformation("DG", "TG", ("DGAT1", "DGAT2"),
                        pathway="glycerolipid_storage", ec_class="2.3.1.20"),
    LipidTransformation("TG", "DG", ("PNPLA2", "LIPE", "PNPLA3"),
                        pathway="glycerolipid_storage", ec_class="3.1.1.3",
                        notes="ATGL/HSL — lipolysis"),
    LipidTransformation("DG", "MG", ("PNPLA2", "LIPE", "PNPLA3"),
                        pathway="glycerolipid_storage", ec_class="3.1.1.3"),
    LipidTransformation("MG", "DG", ("MOGAT1", "MOGAT2", "MOGAT3"),
                        pathway="glycerolipid_storage", ec_class="2.3.1.22"),
    # Cholesterol → CE → Cholesterol (sterol-ester cycle)
    LipidTransformation("Cholesterol", "CE", ("SOAT1", "SOAT2"),
                        pathway="sterol_ester_cycle", ec_class="2.3.1.26",
                        notes="ACAT1/ACAT2 — cholesterol acyltransferase"),
    LipidTransformation("CE", "Cholesterol", ("LIPA", "LIPE", "NCEH1"),
                        pathway="sterol_ester_cycle", ec_class="3.1.1.13",
                        notes="LIPA lysosomal; LIPE/NCEH1 cytosolic"),
]


# ---------------------------------------------------------------------------
# Sphingolipid de novo biosynthesis + salvage
# ---------------------------------------------------------------------------
_SPHINGOLIPIDS: list[LipidTransformation] = [
    # de novo: serine + palmitoyl-CoA → 3-ketosphinganine → sphinganine → dhCer
    # We model only sub-class transitions; serine/palmitoyl-CoA are at
    # metabolite-graph layer, not lipid sub-class layer.
    LipidTransformation("dhSPB", "dhCer", ("CERS1", "CERS2", "CERS3", "CERS4",
                                            "CERS5", "CERS6"),
                        pathway="sphingolipid_de_novo", ec_class="2.3.1.24",
                        notes="N-acylation of sphinganine"),
    LipidTransformation("dhCer", "Cer", ("DEGS1", "DEGS2"),
                        pathway="sphingolipid_de_novo", ec_class="1.14.19.17",
                        notes="dihydroceramide desaturase"),
    LipidTransformation("Cer", "SPB", ("ACER1", "ACER2", "ACER3", "ASAH1", "ASAH2"),
                        pathway="sphingolipid_salvage", ec_class="3.5.1.23"),
    LipidTransformation("SPB", "Cer", ("CERS1", "CERS2", "CERS3", "CERS4",
                                        "CERS5", "CERS6"),
                        pathway="sphingolipid_salvage", ec_class="2.3.1.24"),
    LipidTransformation("Cer", "SM", ("SGMS1", "SGMS2"),
                        pathway="sphingolipid_de_novo", ec_class="2.7.8.27",
                        notes="sphingomyelin synthase"),
    LipidTransformation("SM", "Cer", ("SMPD1", "SMPD2", "SMPD3", "SMPD4"),
                        pathway="sphingolipid_salvage", ec_class="3.1.4.12",
                        notes="sphingomyelinase — neutral and acid forms"),
    LipidTransformation("dhCer", "dhSM", ("SGMS1", "SGMS2"),
                        pathway="sphingolipid_de_novo", ec_class="2.7.8.27"),
    LipidTransformation("Cer", "Cer1P", ("CERK",),
                        pathway="sphingolipid_de_novo", ec_class="2.7.1.138"),
    LipidTransformation("Cer1P", "Cer", ("CERS1", "CERS2"),  # via PPM1F-style phosphatases
                        pathway="sphingolipid_salvage",
                        notes="Cer1P phosphatase — exact enzymes still uncertain"),
    LipidTransformation("SPB", "SPBP", ("SPHK1", "SPHK2"),
                        pathway="sphingolipid_signaling", ec_class="2.7.1.91",
                        notes="sphingosine kinase — produces S1P"),
    LipidTransformation("SPBP", "SPB", ("SGPP1", "SGPP2"),
                        pathway="sphingolipid_signaling", ec_class="3.1.3.n",
                        notes="S1P phosphatase"),
    LipidTransformation("dhSPB", "dhSPBP", ("SPHK1", "SPHK2"),
                        pathway="sphingolipid_signaling", ec_class="2.7.1.91"),

    # Glycosphingolipids — extension beyond BioPAN
    LipidTransformation("Cer", "HexCer", ("UGCG", "UGT8"),
                        pathway="glycosphingolipid_biosynthesis", ec_class="2.4.1.80",
                        notes="UGCG = GlcCer synthase, UGT8 = GalCer synthase"),
    LipidTransformation("HexCer", "Hex2Cer", ("B4GALT6", "B4GALT5"),
                        pathway="glycosphingolipid_biosynthesis", ec_class="2.4.1.n"),
    LipidTransformation("HexCer", "SHexCer", ("GAL3ST1",),
                        pathway="glycosphingolipid_biosynthesis", ec_class="2.8.2.11",
                        notes="cerebroside sulfotransferase — sulfatide synthesis"),
    LipidTransformation("HexCer", "Cer", ("GBA", "GALC"),
                        pathway="glycosphingolipid_catabolism", ec_class="3.2.1.45"),
    LipidTransformation("Hex2Cer", "GM3", ("ST3GAL5",),
                        pathway="ganglioside_biosynthesis", ec_class="2.4.99.9",
                        notes="GM3 synthase"),
    LipidTransformation("GM3", "GD3", ("ST8SIA1",),
                        pathway="ganglioside_biosynthesis", ec_class="2.4.99.8"),
    LipidTransformation("GD3", "GD2", ("B4GALNT1",),
                        pathway="ganglioside_biosynthesis", ec_class="2.4.1.92"),
    LipidTransformation("GD2", "GD1", ("B3GALT4",),
                        pathway="ganglioside_biosynthesis", ec_class="2.4.1.62"),

    # 1-DeoxyCer biosynthesis — Ando aging headline; SPT promiscuity using alanine
    LipidTransformation("dhSPB", "1-DeoxyCer", ("CERS1", "CERS2", "CERS4", "CERS5"),
                        pathway="sphingolipid_de_novo",
                        notes="EXTENSION — 1-deoxysphinganine N-acylation; "
                              "produced when SPT uses L-alanine instead of L-serine "
                              "(SPTLC1/SPTLC2/SPTSSA/SPTSSB). Ando 2019 aging finding."),
]


# ---------------------------------------------------------------------------
# FA elongation / desaturation cascades (BioPAN's molecular-species mode)
# ---------------------------------------------------------------------------
# These don't fit the sub-class transformation model cleanly because they
# operate on individual chain lengths. We record canonical chain
# elongation/desaturation steps as (FFA -> FFA) edges with chain-length
# annotation in `notes`. Layer A heat propagation handles the chain ladder
# at the molecular-species level.
_FA_REMODELING: list[LipidTransformation] = [
    LipidTransformation("FFA", "FFA", ("ELOVL1",),
                        pathway="FA_elongation", ec_class="2.3.1.199",
                        notes="C26 saturated/monounsaturated; C24:0 → C26:0"),
    LipidTransformation("FFA", "FFA", ("ELOVL2",),
                        pathway="FA_elongation", ec_class="2.3.1.199",
                        notes="C20-22 PUFA; C20:5 → C22:5, C22:5 → C24:5 (DHA precursor)"),
    LipidTransformation("FFA", "FFA", ("ELOVL3",),
                        pathway="FA_elongation", ec_class="2.3.1.199",
                        notes="C18-22 sat/MUFA"),
    LipidTransformation("FFA", "FFA", ("ELOVL4",),
                        pathway="FA_elongation", ec_class="2.3.1.199",
                        notes="VLCFA C26+; retina/skin specific"),
    LipidTransformation("FFA", "FFA", ("ELOVL5",),
                        pathway="FA_elongation", ec_class="2.3.1.199",
                        notes="C18-20 PUFA; C18:3 → C20:3"),
    LipidTransformation("FFA", "FFA", ("ELOVL6",),
                        pathway="FA_elongation", ec_class="2.3.1.199",
                        notes="C12-16 saturated/MUFA; C16:0 → C18:0"),
    LipidTransformation("FFA", "FFA", ("ELOVL7",),
                        pathway="FA_elongation", ec_class="2.3.1.199",
                        notes="C16-20 sat; C16:0 → C18:0, C18:0 → C20:0"),
    # Desaturases
    LipidTransformation("FFA", "FFA", ("SCD", "SCD5"),
                        pathway="FA_desaturation", ec_class="1.14.19.1",
                        notes="Δ9 desaturase — C16:0 → C16:1, C18:0 → C18:1"),
    LipidTransformation("FFA", "FFA", ("FADS1",),
                        pathway="FA_desaturation", ec_class="1.14.19.44",
                        notes="Δ5 desaturase — DGLA → AA (ω-6), EPA generation"),
    LipidTransformation("FFA", "FFA", ("FADS2",),
                        pathway="FA_desaturation", ec_class="1.14.19.3",
                        notes="Δ6 desaturase — LA → GLA, ALA → SDA"),
]


# ---------------------------------------------------------------------------
# Acylcarnitine shuttle (β-oxidation entry — extension beyond BioPAN)
# ---------------------------------------------------------------------------
_ACYLCARNITINE: list[LipidTransformation] = [
    LipidTransformation("Acyl-CoA", "CAR", ("CPT1A", "CPT1B", "CPT1C"),
                        pathway="acylcarnitine_shuttle", ec_class="2.3.1.21",
                        notes="EXTENSION — outer mitochondrial membrane (CPT1)"),
    LipidTransformation("CAR", "Acyl-CoA", ("CPT2",),
                        pathway="acylcarnitine_shuttle", ec_class="2.3.1.21",
                        notes="EXTENSION — inner mitochondrial membrane (CPT2)"),
    LipidTransformation("FFA", "Acyl-CoA", ("ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACSL6"),
                        pathway="acylcarnitine_shuttle", ec_class="6.2.1.3",
                        notes="EXTENSION — fatty acid activation"),
]


# ---------------------------------------------------------------------------
# Eicosanoid liberation (PC/PE arachidonate release)
# ---------------------------------------------------------------------------
_EICOSANOID: list[LipidTransformation] = [
    LipidTransformation("PC", "FFA", ("PLA2G4A", "PLA2G4B", "PLA2G4C",
                                       "PLA2G4D", "PLA2G4E", "PLA2G4F"),
                        pathway="eicosanoid_liberation", ec_class="3.1.1.4",
                        notes="cPLA2 family — arachidonate liberation from sn-2"),
    LipidTransformation("PE", "FFA", ("PLA2G6",),
                        pathway="eicosanoid_liberation", ec_class="3.1.1.4",
                        notes="iPLA2 family"),
    LipidTransformation("PI", "FFA", ("PLA2G6",),
                        pathway="eicosanoid_liberation", ec_class="3.1.1.4"),
]


# ---------------------------------------------------------------------------
# Ether / plasmalogen biosynthesis (BioPAN explicit)
# ---------------------------------------------------------------------------
# Ether lipid synthesis goes via DHAP → alkyl-DHAP → 1-alkyl-2-acyl-G3P.
# Modeled as edges between O-PC ↔ PC (acyl exchange of ether form to
# diacyl form is not biochemically real — they are distinct biosynthetic
# branches. The edges below only represent further remodeling within
# the ether/plasmalogen series.)
_ETHER_LIPIDS: list[LipidTransformation] = [
    # Ether lipid de novo synthesis happens at the O-DG / O-PA level;
    # downstream conversion to O-PC / O-PE uses the same Kennedy enzymes.
    LipidTransformation("O-DG", "O-PC", ("CHPT1", "CEPT1"),
                        pathway="ether_lipid_biosynthesis", ec_class="2.7.8.2",
                        notes="ether-DG → ether-PC via Kennedy enzymes"),
    LipidTransformation("O-DG", "O-PE", ("CEPT1", "SELENOI"),
                        pathway="ether_lipid_biosynthesis", ec_class="2.7.8.1"),
    # Plasmalogen formation — vinyl-ether bond formation
    LipidTransformation("O-PE", "P-PE", ("TMEM189",),
                        pathway="plasmalogen_biosynthesis", ec_class="1.14.19.77",
                        notes="plasmanylethanolamine desaturase (PEDS1)"),
    LipidTransformation("O-PC", "P-PC", ("TMEM189",),
                        pathway="plasmalogen_biosynthesis",
                        notes="suspected — TMEM189 may also act on ether-PC"),
    # Lyso ether/plasmalogen forms
    LipidTransformation("O-PC", "O-LPC", ("PLA2G6",),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    LipidTransformation("O-PE", "O-LPE", ("PLA2G6",),
                        pathway="Lands_cycle", ec_class="3.1.1.4"),
    LipidTransformation("P-PC", "P-LPC", ("PLA2G6",),
                        pathway="Lands_cycle"),
    LipidTransformation("P-PE", "P-LPE", ("PLA2G6",),
                        pathway="Lands_cycle"),
]


# ---------------------------------------------------------------------------
# Phosphoinositide cycle (BioPAN's PI/PIP/PIP2/PIP3 nodes)
# ---------------------------------------------------------------------------
_PHOSPHOINOSITIDE: list[LipidTransformation] = [
    LipidTransformation("PI", "PIP", ("PI4K2A", "PI4K2B", "PI4KA", "PI4KB",
                                       "PIP4K2A", "PIP4K2B", "PIP4K2C"),
                        pathway="phosphoinositide_cycle", ec_class="2.7.1.67"),
    LipidTransformation("PIP", "PIP2", ("PIP5K1A", "PIP5K1B", "PIP5K1C",
                                          "PIP4K2A", "PIP4K2B", "PIP4K2C"),
                        pathway="phosphoinositide_cycle", ec_class="2.7.1.149"),
    LipidTransformation("PIP2", "PIP3", ("PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"),
                        pathway="phosphoinositide_cycle", ec_class="2.7.1.153"),
    LipidTransformation("PIP3", "PIP2", ("PTEN", "INPP5D", "INPP5E"),
                        pathway="phosphoinositide_cycle", ec_class="3.1.3.67"),
    LipidTransformation("PIP2", "PIP", ("INPP5A", "INPP5B", "INPP5E", "OCRL", "SYNJ1", "SYNJ2"),
                        pathway="phosphoinositide_cycle", ec_class="3.1.3.36"),
    LipidTransformation("PIP", "PI", ("SAC1", "SAC2", "INPP5K"),
                        pathway="phosphoinositide_cycle", ec_class="3.1.3.n"),
]


# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------
_ALL: list[LipidTransformation] = (
    _KENNEDY
    + _LANDS
    + _GLYCEROLIPIDS
    + _SPHINGOLIPIDS
    + _FA_REMODELING
    + _ACYLCARNITINE
    + _EICOSANOID
    + _ETHER_LIPIDS
    + _PHOSPHOINOSITIDE
)


def iter_transformations(
    *,
    pathway: Optional[str] = None,
    substrate: Optional[str] = None,
    product: Optional[str] = None,
) -> list[LipidTransformation]:
    """Iterate the curated transformation set with optional filters."""
    out = list(_ALL)
    if pathway is not None:
        out = [t for t in out if t.pathway == pathway]
    if substrate is not None:
        out = [t for t in out if t.substrate == substrate]
    if product is not None:
        out = [t for t in out if t.product == product]
    return out


def all_pathways() -> set[str]:
    """Set of canonical pathway tags used in the curated set."""
    return {t.pathway for t in _ALL}


def all_subclasses_referenced() -> set[str]:
    """Sub-class strings appearing as substrate or product in any transformation.

    Used by tests to verify referential integrity against ``taxonomy``.
    """
    s: set[str] = set()
    for t in _ALL:
        s.add(t.substrate)
        s.add(t.product)
    return s


def all_genes_referenced() -> set[str]:
    """All gene symbols appearing as enzymes anywhere in the curated set."""
    g: set[str] = set()
    for t in _ALL:
        g.update(t.enzymes)
    return g
