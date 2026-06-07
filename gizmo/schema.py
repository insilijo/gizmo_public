"""
Node and edge schemas for the GIZMO metabolite graph.

Graph topology: directed multi-typed graph.
  - MetaboliteNode  (node_type="metabolite")
  - ReactionNode    (node_type="reaction")
  - DiseaseNode     (node_type="disease")
  - GeneNode        (node_type="gene")
  - PhenotypeNode   (node_type="phenotype")   — HPO terms
  - DrugNode        (node_type="drug")        — ChEMBL/DrugBank compounds
  - VariantNode     (node_type="variant")     — ClinVar pathogenic variants

Core bipartite edges:
  metabolite → reaction  (substrate/modifier)
  reaction   → metabolite (product)

Clinical overlay edges (DiseaseEdge):
  disease → gene        (genetic association, Open Targets)
  disease → reaction    (pathway association, Reactome disease)
  disease → metabolite  (biomarker / causal metabolite, Orphanet)
  gene    → reaction    (enzyme-reaction, via EC/Reactome)

Phenotype overlay edges (PhenotypeEdge):
  phenotype → disease   (HPO annotation)
  phenotype → gene      (HPO–gene association)
  phenotype → metabolite (known metabolic biomarker)

Drug overlay edges (DrugEdge):
  drug → gene           (drug target, ChEMBL)
  drug → metabolite     (pharmacokinetic effect)

Variant overlay edges (VariantEdge):
  variant → gene        (ClinVar gene–variant)

Stoichiometry on ReactionEdge; association scores on DiseaseEdge.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


# ---------------------------------------------------------------------------
# Edge role enums
# ---------------------------------------------------------------------------

class EdgeRole(str, Enum):
    SUBSTRATE = "substrate"
    PRODUCT = "product"
    MODIFIER = "modifier"    # catalytic / inhibitory


class DiseaseEdgeType(str, Enum):
    GENE_ASSOCIATED = "gene_associated"          # disease ↔ gene (Open Targets)
    PATHWAY_ASSOCIATED = "pathway_associated"    # disease → reaction/pathway
    BIOMARKER = "biomarker"                      # disease → metabolite (Orphanet)
    CAUSAL = "causal"                            # inborn error: enzyme defect → substrate accumulation
    GENE_REACTION = "gene_reaction"              # gene → reaction (enzyme catalysis)


# ---------------------------------------------------------------------------
# Metabolite node
# ---------------------------------------------------------------------------

class MetaboliteNode(BaseModel):
    """
    Represents a chemical species (compartment-aware).

    Canonical node_id: "CHEBI:XXXXX" or "CHEBI:XXXXX@compartment".
    Falls back to "reactome:{stId}" or "pubchem:{cid}" when ChEBI is unknown.
    """

    node_id: str
    chebi_id: Optional[str] = None
    hmdb_id: Optional[str] = None          # HMDB0000001 (7-digit, VMH-compatible)
    vmh_id: Optional[str] = None           # VMH abbreviation (e.g. "glc_D", "phe_L")
    metanetx_id: Optional[str] = None      # MNXM_XXXXXX
    reactome_id: Optional[str] = None      # R-ALL-XXXXXX PhysicalEntity stId
    pubchem_cid: Optional[str] = None
    metabolon_name: Optional[str] = None   # Metabolon BIOCHEMICAL field
    name: str
    formula: Optional[str] = None
    charge: Optional[int] = None
    mass: Optional[float] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    smiles: Optional[str] = None
    compartment: Optional[str] = None
    cas_id: Optional[str] = None           # CAS Registry Number (e.g. "50-18-0")
    # Metabolon platform / chromatography context
    platform: Optional[str] = None        # e.g. "lc/ms neg"
    retention_index: Optional[float] = None
    # Sample-type context
    biofluids: list[str] = Field(default_factory=list)          # ["plasma", "urine", "csf"]
    reference_ranges: dict = Field(default_factory=dict)         # {biofluid: {low, high, unit}}
    # Graph analysis flags
    is_currency: bool = False              # ATP, NAD, H2O etc. — flag but don't delete
    manually_reviewed: bool = False        # True when ChEBI/HMDB assigned via MetabolonCurator
    node_type: str = "metabolite"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Reaction node
# ---------------------------------------------------------------------------

class ReactionNode(BaseModel):
    """Represents a biochemical reaction or transport event."""

    node_id: str                           # "reactome:R-HSA-XXXXXX" or "mnxr:MNXR_XXXXX"
    reactome_id: Optional[str] = None
    metanetx_id: Optional[str] = None     # MNXR_XXXXXX
    name: str
    reversible: bool = False
    direction: Optional[str] = None       # "left-to-right" | "right-to-left" | "bidirectional"
    ec_numbers: list[str] = Field(default_factory=list)
    gene_symbols: list[str] = Field(default_factory=list)  # HGNC symbols of catalysing genes
    pathways: list[str] = Field(default_factory=list)      # Reactome pathway stIDs
    species: Optional[str] = None
    node_type: str = "reaction"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Disease node
# ---------------------------------------------------------------------------

class DiseaseNode(BaseModel):
    """
    Represents a disease entity.

    Canonical node_id: "MONDO:XXXXXXX"
    Cross-references are stored as strings; we do not embed OMIM data.
    """

    node_id: str                           # "MONDO:XXXXXXX"
    mondo_id: Optional[str] = None
    name: str
    synonyms: list[str] = Field(default_factory=list)
    definition: Optional[str] = None
    # Cross-references (ID strings only — no embedded restricted data)
    xref_omim: list[str] = Field(default_factory=list)     # "OMIM:XXXXXX"
    xref_orphanet: list[str] = Field(default_factory=list) # "Orphanet:XXXXX"
    xref_doid: list[str] = Field(default_factory=list)     # "DOID:XXXXX"
    xref_icd10: list[str] = Field(default_factory=list)    # "ICD10:XXXXX"
    xref_mesh: list[str] = Field(default_factory=list)     # "MeSH:DXXXXXX"
    # Classification flags
    is_rare: bool = False
    is_inborn_error_of_metabolism: bool = False
    node_type: str = "disease"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Gene node
# ---------------------------------------------------------------------------

class GeneNode(BaseModel):
    """
    Represents a gene (any species).

    Canonical node_id: "ENSG:ENSGXXXXXXXXXXXX" (Ensembl)
    Falls back to "HGNC:XXXXX" or "symbol:{SYMBOL}".
    """

    node_id: str
    ensembl_id: Optional[str] = None
    hgnc_id: Optional[str] = None
    entrez_id: Optional[str] = None        # NCBI Gene ID (for ortholog mapping)
    symbol: str
    name: Optional[str] = None
    species: str = "Homo sapiens"          # e.g. "Mus musculus", "Rattus norvegicus"
    tissue_expression: dict = Field(default_factory=dict)   # {tissue: median_tpm}
    node_type: str = "gene"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Pathway node  (Reactome / KEGG)
# ---------------------------------------------------------------------------

class PathwayNode(BaseModel):
    """
    Represents a biological pathway.

    Canonical node_id: "reactome:{stId}" or "kegg:{PATHWAY_ID}"
    Edges: pathway → reaction (pathway_reaction), reaction → pathway (reaction_pathway)
    """

    node_id: str                            # "reactome:R-HSA-XXXXXX"
    reactome_id: Optional[str] = None
    kegg_id: Optional[str] = None
    name: str
    species: str = "Homo sapiens"
    parent_pathways: list[str] = Field(default_factory=list)  # parent stIds
    level: int = 0                          # hierarchy depth (0 = top-level)
    node_type: str = "pathway"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Phenotype node  (HPO)
# ---------------------------------------------------------------------------

class PhenotypeNode(BaseModel):
    """
    Represents an HPO phenotype term.

    Canonical node_id: "HP:XXXXXXX"
    Edges: phenotype → disease (phenotype_disease),
           phenotype → gene (phenotype_gene),
           phenotype → metabolite (phenotype_metabolite)
    """

    node_id: str                            # "HP:XXXXXXX"
    hpo_id: Optional[str] = None
    name: str
    definition: Optional[str] = None
    synonyms: list[str] = Field(default_factory=list)
    is_metabolic: bool = False              # True for terms under "Abnormality of metabolism"
    node_type: str = "phenotype"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Drug node  (ChEMBL / DrugBank)
# ---------------------------------------------------------------------------

class DrugNode(BaseModel):
    """
    Represents a drug or clinical compound.

    Canonical node_id: "CHEMBL:{id}" or "DB:{drugbank_id}"
    Edges: drug → gene (drug_target), drug → metabolite (drug_metabolite)
    """

    node_id: str                            # "CHEMBL:CHEMBL1234" or "DB:DB01234"
    chembl_id: Optional[str] = None
    drugbank_id: Optional[str] = None
    name: str
    synonyms: list[str] = Field(default_factory=list)
    max_phase: int = 0                      # 4 = approved, 3/2/1 = clinical, 0 = preclinical
    mechanism: Optional[str] = None        # "INHIBITOR" | "ACTIVATOR" | "MODULATOR"
    atc_codes: list[str] = Field(default_factory=list)
    node_type: str = "drug"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Variant node  (ClinVar)
# ---------------------------------------------------------------------------

class VariantNode(BaseModel):
    """
    Represents a pathogenic or likely-pathogenic genetic variant.

    Canonical node_id: "ClinVar:{variation_id}" or "rs:{rsid}"
    Edges: variant → gene (variant_gene)
    """

    node_id: str                             # "ClinVar:12345" or "rs:12345"
    clinvar_id: Optional[str] = None        # ClinVar variation ID
    rsid: Optional[str] = None              # dbSNP rs ID
    gene_id: Optional[str] = None           # graph gene node ID
    gene_symbol: Optional[str] = None       # HGNC symbol
    consequence: Optional[str] = None       # "missense", "nonsense", "frameshift", "splice"
    clinical_significance: str = "Pathogenic"   # "Pathogenic" | "Likely pathogenic"
    condition: Optional[str] = None         # associated condition/disease name
    review_status: Optional[str] = None     # ClinVar review stars label
    node_type: str = "variant"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Edges
# ---------------------------------------------------------------------------

class ReactionEdge(BaseModel):
    """
    Directed edge between a metabolite and a reaction node.

    Provenance fields allow downstream code to weight edges by evidence quality
    and distinguish curated from computationally inferred participation.
    """

    source: str
    target: str
    role: EdgeRole
    stoichiometry: float = 1.0
    compartment: Optional[str] = None

    # Provenance
    source_db: Optional[str] = None        # "reactome" | "metanetx" | "manual"
    evidence_level: Optional[str] = None   # "curated" | "inferred" | "predicted"
    confidence: float = 1.0                # [0, 1]; 1.0 = curated
    manually_reviewed: bool = False        # True when edge was manually verified/added

    model_config = ConfigDict(frozen=True)


class DiseaseEdge(BaseModel):
    """
    Directed edge connecting disease/gene nodes to metabolite/reaction nodes.
    score: Open Targets association score [0, 1] where available.
    """

    source: str
    target: str
    edge_type: DiseaseEdgeType
    score: Optional[float] = None
    evidence_count: Optional[int] = None
    source_db: Optional[str] = None       # "open_targets" | "orphanet" | "reactome"
    manually_reviewed: bool = False        # True when association was manually verified

    model_config = ConfigDict(frozen=True)


class PathwayEdge(BaseModel):
    """
    Directed edge connecting pathway nodes to reactions or sub-pathways.
    pathway → reaction  (pathway_reaction)
    pathway → pathway   (pathway_hierarchy — parent → child)
    """

    source: str
    target: str
    edge_type: str = "pathway_reaction"    # "pathway_reaction" | "pathway_hierarchy"
    source_db: Optional[str] = None        # "reactome" | "kegg"

    model_config = ConfigDict(frozen=True)


class PhenotypeEdge(BaseModel):
    """
    Directed edge from a phenotype (HPO) node to a disease, gene, or metabolite.
    """

    source: str                             # HP:XXXXXXX
    target: str
    edge_type: str                          # "phenotype_disease" | "phenotype_gene" | "phenotype_metabolite"
    score: Optional[float] = None
    source_db: Optional[str] = None        # "hpo" | "omim"

    model_config = ConfigDict(frozen=True)


class DrugEdge(BaseModel):
    """
    Directed edge from a drug node to a gene (target) or metabolite (effect).
    """

    source: str                             # CHEMBL:xxx or DB:xxx
    target: str
    edge_type: str = "drug_target"         # "drug_target" | "drug_metabolite"
    mechanism: Optional[str] = None        # "INHIBITOR" | "ACTIVATOR" | "MODULATOR"
    max_phase: int = 0
    source_db: Optional[str] = None        # "chembl" | "drugbank"

    model_config = ConfigDict(frozen=True)


class VariantEdge(BaseModel):
    """
    Directed edge from a variant node to its gene.
    """

    source: str                             # ClinVar:xxx or rs:xxx
    target: str                             # gene node ID
    edge_type: str = "variant_gene"
    consequence: Optional[str] = None
    source_db: Optional[str] = None        # "clinvar"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Microbial taxon node  (GeMMA / AGORA2 community)
# ---------------------------------------------------------------------------

class MicrobialTaxonNode(BaseModel):
    """
    Represents a gut microbial taxon contributing to community metabolism.

    Canonical node_id: "taxon:{rank}__{name}" (e.g. "taxon:genus__Prevotella").
    Populated from GeMMA's SlimmedCommunity via
    :func:`~gizmo.integration.gemma.attach_community`.
    """

    node_id: str                            # "taxon:genus__Prevotella"
    name: str                               # display name (e.g. "Prevotella")
    rank: str = "genus"                     # "genus" | "species" | "family"
    ncbi_taxid: Optional[str] = None
    agora2_model_ids: list[str] = Field(default_factory=list)
    mean_abundance: float = 0.0             # mean relative abundance across samples
    sample_abundances: dict[str, float] = Field(default_factory=dict)
    node_type: str = "microbial_taxon"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Microbial edge  (taxon ↔ metabolite, via AGORA2 reactions)
# ---------------------------------------------------------------------------

class MicrobialEdgeRole(str, Enum):
    PRODUCER  = "producer"   # taxon produces this metabolite (reaction product)
    CONSUMER  = "consumer"   # taxon consumes this metabolite (reaction substrate)
    BIDIRECTIONAL = "bidirectional"  # reversible reaction — both roles possible


class MicrobialEdge(BaseModel):
    """
    Directed edge from a microbial taxon node to a metabolite node.

    Encodes which gut taxa produce or consume a metabolite according to their
    AGORA2 genome-scale metabolic models.  Mean reaction weight across samples
    is stored as ``abundance_weight`` — higher = more community capacity for
    this reaction.

    source: microbial taxon node_id
    target: GIZMO metabolite node_id (e.g. "CHEBI:17234" or "HMDB:HMDB0000122")
    """

    source: str                             # "taxon:genus__Prevotella"
    target: str                             # GIZMO metabolite node_id
    edge_type: str = "microbial_metabolite"
    role: MicrobialEdgeRole = MicrobialEdgeRole.BIDIRECTIONAL
    vmh_metabolite_id: Optional[str] = None  # original VMH abbreviation
    reaction_ids: list[str] = Field(default_factory=list)
    abundance_weight: float = 0.0           # mean reaction weight across samples
    sample_weights: dict[str, float] = Field(default_factory=dict)
    source_db: str = "agora2"

    model_config = ConfigDict(frozen=True)


class ToxEdgeType(str, Enum):
    TOX_GENE    = "tox_gene"       # chemical → gene (expression/activity change)
    TOX_DISEASE = "tox_disease"    # chemical → disease (causal / marker)


class ToxEdge(BaseModel):
    """
    Directed edge from a metabolite/chemical node to a gene or disease node,
    representing a toxicological association.

    Sources: CTD, CompTox (EPA), T3DB, ChEMBL ADMET assays.
    """

    source: str                                # metabolite node ID (the chemical)
    target: str                                # gene or disease node ID
    edge_type: ToxEdgeType
    effect_type: Optional[str]  = None         # CTD: "increases^expression", etc.
    organism:    Optional[str]  = None         # "Homo sapiens" | "Rattus norvegicus" …
    p_value:     Optional[float] = None
    reference_count: Optional[int] = None
    direct_evidence: Optional[str] = None      # CTD: "marker/mechanism" | "therapeutic"
    assay_endpoint:  Optional[str] = None      # ChEMBL: "hERG IC50", "LD50", "AMES"
    assay_value:     Optional[float] = None    # numeric endpoint value
    assay_units:     Optional[str]  = None     # "uM", "mg/kg", etc.
    source_db: str = "ctd"

    model_config = ConfigDict(frozen=True)


class ProteinInteractionEdge(BaseModel):
    """
    Undirected protein–protein interaction from STRING (CC BY 4.0).

    Scores are normalised to [0, 1] from STRING's integer range [0, 1000].
    Edge is stored bidirectionally (A→B and B→A) in the DiGraph so that
    standard predecessor/successor traversal works in both directions.

    Source: https://string-db.org
    """

    source: str                      # gene node ID
    target: str                      # gene node ID
    edge_type: str = "protein_interaction"
    combined_score: float = 0.0      # overall confidence
    experimental:   float = 0.0      # physical interaction evidence
    coexpression:   float = 0.0      # mRNA co-expression
    database:       float = 0.0      # curated database support
    textmining:     float = 0.0      # literature co-mention
    source_db: str = "stringdb"

    model_config = ConfigDict(frozen=True)
