"""
Bayesian inference over the GIZMO knowledge graph.

The headline entry point is :func:`run_bayesian_inference` — given a
``GizmoGraph`` and a ``SampleContext``, it builds a 3-state Markov Random
Field over the (compartment-collapsed, currency-filtered) graph and runs
loopy belief propagation to produce a per-node posterior over perturbation
direction: ``{DOWN: p₋, NORMAL: p₀, UP: p₊}``.

Rationale
---------
Existing GIZMO scoring is heuristic and layer-by-layer: evidence at
metabolites drives reaction scores, which drive pathway/disease scores via
aggregation. That works but it doesn't *combine* multiomic evidence
coherently: a gene with a strong effect and a phenotype confirmation are
treated separately rather than as two observations of a joint latent.

Bayesian propagation treats every node's perturbation state as a latent
variable and every observation (metabolite effect size, gene effect size,
ICD code, HPO term) as a soft likelihood. Typed pairwise edge potentials
encode the biochemistry: substrate-of edges couple metabolite and reaction
in opposite direction, product-of in the same direction, gene→reaction with
a positive coupling, etc.

This is inference on a heterogeneous MRF with loopy belief propagation.
Posteriors at unobserved nodes (e.g. disease nodes with no ICD evidence) are
driven entirely by the evidence flowing through the graph.
"""

from gizmo.inference.bp import run_bayesian_inference, BPConfig, BPResult
from gizmo.inference.model import DEFAULT_COUPLINGS, Couplings

__all__ = [
    "run_bayesian_inference", "BPConfig", "BPResult",
    "DEFAULT_COUPLINGS", "Couplings",
]
