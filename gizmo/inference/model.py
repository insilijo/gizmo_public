"""
Model specification for heterogeneous MRF + BP over the GIZMO graph.

State space per node: 3-valued categorical ``{DOWN, NORMAL, UP}``. We index
these as ``[0, 1, 2]`` in arrays but logically treat them as signed.

Edge potentials: one coupling constant per edge type, optionally
sign-flipped per edge role. The potential matrix for an edge of type *t*
connecting nodes i→j is

    ψ_t(x_i, x_j) = exp( θ_t · s_t · sgn(x_i) · sgn(x_j) )

where ``sgn(0) = 0``. A positive θ couples the two states (they tend to
agree in sign); s_t = -1 flips the direction (e.g. substrate of a perturbed
reaction should move opposite to product).

Unary potentials: node prior × observation likelihood. The prior is a weak
"normal is more likely" bias controlled by ``normal_bias`` in BPConfig.
Diseases and phenotypes get an extra "normal is much more likely" pull to
avoid every disease lighting up when nothing has been confirmed.

Currency metabolites (ATP, water, cofactors) are filtered at graph-build
time — they're high-degree hubs that destroy information flow in BP.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# 3-state indexing
DOWN, NORMAL, UP = 0, 1, 2
SIGN = np.array([-1.0, 0.0, +1.0])   # sgn() for each state index


# ---------------------------------------------------------------------------
# Coupling constants (θ, s) per edge type
# ---------------------------------------------------------------------------

@dataclass
class Couplings:
    """
    Coupling constants for each edge type. ``theta`` is the strength;
    ``sign`` flips the direction of coupling (+1 = same-sign, -1 = opposite).

    The keys are edge ``role`` (preferred) or ``edge_type`` (fallback). Match
    is case-insensitive.
    """

    values: dict[str, tuple[float, float]] = field(default_factory=dict)

    def theta_sign(self, edge_role: str) -> tuple[float, float]:
        """Return (theta, sign) for an edge; default (0, +1) if unknown."""
        if not edge_role:
            return 0.0, 1.0
        key = edge_role.strip().lower()
        return self.values.get(key, (0.0, 1.0))

    def update(self, other: dict[str, tuple[float, float]]) -> None:
        self.values.update({k.lower(): v for k, v in other.items()})


# ---------------------------------------------------------------------------
# Per-assay observation noise
# ---------------------------------------------------------------------------
# Gaussian observation σ for ``gaussian_obs_likelihood``, keyed by
# ``EvidenceRecord.assay_type``. Metabolomics observes reaction substrates and
# products directly and is the baseline. Proteomics observes enzyme abundance
# (one hop from reaction activity). Transcriptomics observes mRNA, which sits
# upstream of protein via translation and post-translational regulation;
# cross-tissue mRNA–protein Pearson correlation is ~0.4 (Liu et al., Nature
# 2016; Buccitelli & Selbach, Nat Rev Genet 2020), implying σ_transcript ≈
# 2 × σ_protein when both are modelled as noisy observations of the same
# latent enzyme-activity state.
DEFAULT_ASSAY_SIGMAS: dict[str, float] = {
    "metabolomics":    1.0,
    "proteomics":      1.0,
    "transcriptomics": 2.0,
}


DEFAULT_COUPLINGS = Couplings(values={
    # metabolite ↔ reaction edges
    "substrate":          (1.0,  -1.0),   # met DOWN when rxn UP (and vice versa)
    "product":            (1.0,  +1.0),   # met UP   when rxn UP
    "modifier":           (0.3,  +1.0),
    "cofactor":           (0.2,  +1.0),
    # gene / protein. ``catalysis`` is the Reactome-emitted edge type;
    # ``gene_reaction`` is the legacy alias for non-Reactome graphs.
    #
    # theta=0 (gene evidence does NOT propagate direction through BP):
    # transcriptomics doesn't tell you whether a reaction is running UP
    # or DOWN — it tells you whether the *enzyme is available*. Pushing
    # a reaction UP because its enzyme is UP, or DOWN because its enzyme
    # is DOWN, conflicts with the substrate-side metabolite evidence and
    # collapses the canonical-IDH AUC from 0.87 → 0.26 (verified on
    # GSE190504 IDH1-mut vs WT).
    #
    # Gene evidence is instead applied as a post-BP enzyme-availability
    # gate via gizmo.scoring.enzyme_gate.apply_enzyme_gate — a
    # multiplicative dampener on reaction perturbation magnitude when the
    # catalyzing gene is strongly DOWN. Treats transcriptomics as a
    # stoplight (capacity-limiting) rather than a driver.
    "catalysis":          (0.0,  +1.0),
    "gene_reaction":      (0.0,  +1.0),
    "gene_associated":    (0.0,  +1.0),
    "protein_interaction":(0.5,  +1.0),
    # disease / phenotype
    "disease_gene":       (1.2,  +1.0),
    "disease_reaction":   (1.0,  +1.0),
    "phenotype_disease":  (1.5,  +1.0),
    "hpo_disease":        (1.5,  +1.0),
    # pathway
    "pathway_associated": (0.8,  +1.0),
    "pathway":            (0.8,  +1.0),
    # variant
    "variant_gene":       (2.0,  +1.0),
    # chemical similarity (if present)
    "similar_to":         (0.4,  +1.0),
})


# ---------------------------------------------------------------------------
# Unary priors by node type
# ---------------------------------------------------------------------------

NORMAL_BIAS_DEFAULT = 1.5   # how much "normal" is preferred a priori

# Disease / phenotype nodes get a stronger pull to NORMAL so they don't light
# up without confirmation. This encodes the weak-prior-on-disease choice.
STRONG_NORMAL_NODE_TYPES = {"disease", "phenotype"}
STRONG_NORMAL_BIAS = 4.0


def unary_prior(node_type: str, normal_bias: float = NORMAL_BIAS_DEFAULT) -> np.ndarray:
    """Return a 3-vector prior P(X_i)."""
    if node_type in STRONG_NORMAL_NODE_TYPES:
        b = STRONG_NORMAL_BIAS
    else:
        b = normal_bias
    # Unnormalized: [1, b, 1] — normal weighted more than up/down
    p = np.array([1.0, b, 1.0])
    return p / p.sum()


# ---------------------------------------------------------------------------
# Observation likelihoods
# ---------------------------------------------------------------------------

def gaussian_obs_likelihood(
    effect_size: float,
    confidence: float = 1.0,
    sigma: float = 1.0,
    mu_up: float = 1.0,
) -> np.ndarray:
    """
    Likelihood ``P(observation | X)`` for a continuous effect size.

    We assume the underlying state ``X`` corresponds to a mean effect of
    ``{-μ, 0, +μ}`` respectively (for DOWN, NORMAL, UP) with observation
    noise ``σ``. ``confidence`` scales the signal-to-noise — a low-confidence
    mapping (0.3) widens the effective σ, diluting the signal.

    Returns
    -------
    3-vector of likelihoods, not normalised.
    """
    # Confidence rescales effective signal strength; low confidence => weak
    # likelihood split across states
    s = max(0.05, float(confidence))
    e = float(effect_size) * s
    sig = max(0.1, sigma / s)
    means = np.array([-mu_up, 0.0, +mu_up])
    # Gaussian likelihood, computed in log-space so we don't lose ratio
    # information when |effect_size| >> mu_up. The previous implementation
    # clamped exp(-z²/2) at 1e-6, which collapsed all three states to the
    # floor whenever |z|>~5, making strong observations indistinguishable
    # from noise. Working in log-space + subtracting the max keeps the
    # ratios numerically stable for any effect_size.
    z = (e - means) / sig
    log_like = -0.5 * z * z
    log_like -= log_like.max()           # max = 0 → exp gives 1 for the best state
    like = np.exp(log_like)
    # Floor only the smallest values (don't kill the dominant state)
    return np.maximum(like, 1e-30)


def hard_confirmation_likelihood(state: int, confidence: float = 0.95) -> np.ndarray:
    """
    Likelihood for an ICD/HPO-style confirmation: we are confident that
    ``X = state`` but not certain. ``state ∈ {DOWN, NORMAL, UP}``.
    """
    like = np.full(3, (1.0 - confidence) / 2.0)
    like[state] = confidence
    return np.maximum(like, 1e-6)


# ---------------------------------------------------------------------------
# Pairwise potential matrix
# ---------------------------------------------------------------------------

def pairwise_potential(theta: float, sign: float) -> np.ndarray:
    """
    Build the 3×3 edge potential matrix

        ψ(x_i, x_j) = exp( theta · sign · sgn(x_i) · sgn(x_j) )

    Rows index x_i, columns index x_j.
    """
    s = SIGN[:, None] * SIGN[None, :]   # 3×3 of sgn products
    return np.exp(theta * sign * s)


__all__ = [
    "DOWN", "NORMAL", "UP", "SIGN",
    "Couplings", "DEFAULT_COUPLINGS", "DEFAULT_ASSAY_SIGMAS",
    "unary_prior",
    "gaussian_obs_likelihood", "hard_confirmation_likelihood",
    "pairwise_potential",
    "NORMAL_BIAS_DEFAULT", "STRONG_NORMAL_NODE_TYPES", "STRONG_NORMAL_BIAS",
]
