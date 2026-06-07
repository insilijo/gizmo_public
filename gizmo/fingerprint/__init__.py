"""GIZMO multi-resolution fingerprint package.

Public API consumed by `discomarker.methods.gizmo_fingerprint` (see
`/home/jgardner/discomarker/docs/INTEGRATION_DESIGN.md`).

Two responsibilities:

  - **Bundle I/O** — persist a Paper1Result + SubstrateGeometry as a
    directory on disk, and reload it without needing to re-run MAP.
    See :func:`save_bundle`, :func:`load_bundle`.
  - **Fingerprint transform** — collapse the (n_patients, n_nodes)
    F matrix into a small per-patient feature table suitable as input
    to a supervised learner (elastic-net, SHAP, etc.).
    See :func:`multi_resolution_fingerprint`.
"""
from __future__ import annotations

from gizmo.fingerprint.bundle import (
    compute_substrate_modules,
    load_bundle,
    save_bundle,
)
from gizmo.fingerprint.transforms import multi_resolution_fingerprint

__all__ = [
    "compute_substrate_modules",
    "load_bundle",
    "save_bundle",
    "multi_resolution_fingerprint",
]
