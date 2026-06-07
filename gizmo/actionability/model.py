"""
Actionability data models.

Three nested layers:

    DruggabilityScore     — per enzyme / gene target
    PerturbabilityScore   — per reaction  (aggregates enzyme druggability)
    ActionabilityScore    — combines both with an explanation
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class DruggabilityScore:
    """
    Druggability assessment for one gene / enzyme target.

    Attributes
    ----------
    gene_id              : graph node ID
    symbol               : HGNC symbol
    target_chembl_id     : ChEMBL target ID (if resolved)
    n_approved_drugs     : drugs at max_phase == 4
    n_clinical_drugs     : drugs at max_phase >= 2
    max_phase            : highest clinical phase across all known compounds
    best_drug_name       : name of highest-phase compound
    best_action_type     : e.g. INHIBITOR, ACTIVATOR, MODULATOR
    tractability_sm      : Open Targets small-molecule tractability bucket label
    tractability_ab      : Open Targets antibody tractability bucket label
    score                : overall druggability score [0, 1]
                             1.0 = approved drug exists
                             0.6 = Phase 2/3 compound
                             0.3 = Phase 1 / preclinical
                             0.0 = no known compounds
    drugs                : list of drug dicts from ChEMBL
    """

    gene_id:          str
    symbol:           str
    target_chembl_id: Optional[str]    = None
    n_approved_drugs: int              = 0
    n_clinical_drugs: int              = 0
    max_phase:        int              = 0
    best_drug_name:   Optional[str]    = None
    best_action_type: Optional[str]    = None
    tractability_sm:  Optional[str]    = None
    tractability_ab:  Optional[str]    = None
    score:            float            = 0.0
    drugs:            list[dict]       = field(default_factory=list)

    def __repr__(self) -> str:
        return (f"DruggabilityScore({self.symbol!r}, "
                f"score={self.score:.2f}, max_phase={self.max_phase}, "
                f"n_approved={self.n_approved_drugs})")


@dataclass
class PerturbabilityScore:
    """
    Perturbability assessment for one reaction.

    A reaction is perturbable if at least one of its catalysing enzymes
    has a known drug or chemical tool compound.

    Attributes
    ----------
    reaction_id          : graph node ID
    n_druggable_enzymes  : enzymes with any known compound
    max_phase            : highest clinical phase across all enzymes
    best_drug_name       : name of highest-phase compound across enzymes
    best_gene_symbol     : symbol of the most-druggable enzyme
    action_types         : set of action types (INHIBITOR, ACTIVATOR, …)
    cofactor_dependent   : True if reaction has modifier metabolites
    n_alternative_routes : number of parallel reactions sharing a substrate
                           (higher = more redundancy, harder to perturb)
    score                : overall perturbability score [0, 1]
    notes                : brief explanation
    """

    reaction_id:           str
    n_druggable_enzymes:   int        = 0
    max_phase:             int        = 0
    best_drug_name:        Optional[str] = None
    best_gene_symbol:      Optional[str] = None
    action_types:          list[str]  = field(default_factory=list)
    cofactor_dependent:    bool       = False
    n_alternative_routes:  int        = 0
    score:                 float      = 0.0
    notes:                 str        = ""

    def __repr__(self) -> str:
        return (f"PerturbabilityScore({self.reaction_id!r}, "
                f"score={self.score:.2f}, max_phase={self.max_phase})")


@dataclass
class ActionabilityScore:
    """
    Combined perturbability + druggability score for a reaction.

    Designed to answer: *is this reaction actionable with known chemistry?*

    Attributes
    ----------
    reaction_id        : graph node ID
    perturbability     : PerturbabilityScore for this reaction
    enzyme_scores      : DruggabilityScore for each catalysing enzyme
    combined_score     : max(enzyme druggability) × perturbability base
    reaction_score     : evidence score from Phase 4 (if provided)
    priority_score     : combined_score × |reaction_score| — prioritises
                         reactions that are both evidence-supported and actionable
    explanation        : one-liner suitable for analyst review
    """

    reaction_id:      str
    perturbability:   PerturbabilityScore
    enzyme_scores:    list[DruggabilityScore] = field(default_factory=list)
    combined_score:   float = 0.0
    reaction_score:   float = 0.0
    priority_score:   float = 0.0
    explanation:      str   = ""

    def __repr__(self) -> str:
        return (f"ActionabilityScore({self.reaction_id!r}, "
                f"combined={self.combined_score:.2f}, "
                f"priority={self.priority_score:.3f})")
