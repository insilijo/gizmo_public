"""Bundle save/load and substrate Louvain partitioning.

Bundle directory layout (see ``discomarker/docs/INTEGRATION_DESIGN.md``)::

    bundle/
      F.npz                 # 'F' key — (n_patients, n_nodes) float32
      patient_ids.json      # list[str]
      substrate_nodes.json  # list[{nid, node_type, name}] parallel to F cols
      log_pr.npy            # (n_nodes,) float64
      modules.json          # {str(cluster_idx): [node_ids]}
      eigvecs.npz           # optional, 'eigvecs' key — (n_nodes, k) float32
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from gizmo.inference.projection import SubstrateGeometry


# ---------------------------------------------------------------------------
# Louvain partition on substrate
# ---------------------------------------------------------------------------

def compute_substrate_modules(
    geometry: "SubstrateGeometry",
    *,
    resolution: float = 1.0,
    seed: int = 0,
) -> dict[int, list[str]]:
    """Louvain community partition over the substrate subgraph.

    Returns ``{cluster_idx: [node_ids]}`` with stable int keys.
    Singleton communities are dropped.
    """
    communities = nx.community.louvain_communities(
        geometry.sub, resolution=resolution, seed=seed
    )
    modules: dict[int, list[str]] = {}
    idx = 0
    for comm in communities:
        members = list(comm)
        if len(members) < 2:
            continue
        modules[idx] = sorted(members)
        idx += 1
    return modules


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

def save_bundle(
    path: str | Path,
    *,
    F: np.ndarray,
    patient_ids: list[str],
    geometry: "SubstrateGeometry",
    modules: dict[int, list[str]] | None = None,
    eigvecs: np.ndarray | None = None,
) -> Path:
    """Persist a fingerprint bundle to disk.

    Parameters
    ----------
    path : directory to create / write into (mkdir -p semantics).
    F : (n_patients, n_nodes) per-patient state vectors from
        :func:`gizmo.inference.projection.solve_map`.
    patient_ids : ordered list parallel to F rows.
    geometry : :class:`gizmo.inference.projection.SubstrateGeometry` with
        ``nodes``, ``sub``, and ``log_pr``. Node ``node_type`` and ``name``
        attributes are pulled from ``sub`` for the substrate_nodes manifest.
    modules : Louvain partition; computed from ``geometry.sub`` if None.
    eigvecs : optional top-K Laplacian eigenvectors for GFT-side work.
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    n_patients, n_nodes = F.shape
    if n_nodes != len(geometry.nodes):
        raise ValueError(
            f"F has {n_nodes} columns but geometry has {len(geometry.nodes)} nodes"
        )
    if n_patients != len(patient_ids):
        raise ValueError(
            f"F has {n_patients} rows but {len(patient_ids)} patient_ids supplied"
        )

    np.savez_compressed(path / "F.npz", F=F.astype(np.float32))
    np.save(path / "log_pr.npy", geometry.log_pr.astype(np.float64))

    (path / "patient_ids.json").write_text(json.dumps(list(patient_ids)))

    substrate_nodes = [
        {
            "nid": nid,
            "node_type": geometry.sub.nodes[nid].get("node_type", "unknown"),
            "name": geometry.sub.nodes[nid].get("name", nid),
        }
        for nid in geometry.nodes
    ]
    (path / "substrate_nodes.json").write_text(json.dumps(substrate_nodes))

    if modules is None:
        modules = compute_substrate_modules(geometry)
    modules_serializable = {str(k): list(v) for k, v in modules.items()}
    (path / "modules.json").write_text(json.dumps(modules_serializable))

    if eigvecs is not None:
        np.savez_compressed(path / "eigvecs.npz", eigvecs=eigvecs.astype(np.float32))

    return path


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_bundle(path: str | Path) -> dict[str, Any]:
    """Read a fingerprint bundle directory into an in-memory dict.

    Returned keys (matches the contract in
    ``discomarker/discomarker/methods/gizmo_fingerprint.py``):

      - ``F``           : np.ndarray (n_patients, n_nodes)
      - ``patient_ids`` : list[str]
      - ``nodes``       : list[str]  (node IDs, parallel to F columns)
      - ``node_types``  : list[str]  (parallel to ``nodes``)
      - ``node_names``  : list[str]  (parallel to ``nodes``)
      - ``log_pr``      : np.ndarray (n_nodes,)
      - ``modules``     : dict[int, list[str]]
      - ``eigvecs``     : np.ndarray or None
    """
    path = Path(path)
    if not path.is_dir():
        raise FileNotFoundError(f"bundle directory not found: {path}")

    required = ["F.npz", "patient_ids.json", "substrate_nodes.json",
                "log_pr.npy", "modules.json"]
    missing = [f for f in required if not (path / f).exists()]
    if missing:
        raise FileNotFoundError(
            f"bundle at {path} is missing required files: {missing}"
        )

    F = np.load(path / "F.npz")["F"]
    log_pr = np.load(path / "log_pr.npy")
    patient_ids = json.loads((path / "patient_ids.json").read_text())

    substrate_nodes = json.loads((path / "substrate_nodes.json").read_text())
    nodes = [n["nid"] for n in substrate_nodes]
    node_types = [n["node_type"] for n in substrate_nodes]
    node_names = [n["name"] for n in substrate_nodes]

    modules_raw = json.loads((path / "modules.json").read_text())
    modules = {int(k): list(v) for k, v in modules_raw.items()}

    eigvecs = None
    if (path / "eigvecs.npz").exists():
        eigvecs = np.load(path / "eigvecs.npz")["eigvecs"]

    return {
        "F": F,
        "patient_ids": patient_ids,
        "nodes": nodes,
        "node_types": node_types,
        "node_names": node_names,
        "log_pr": log_pr,
        "modules": modules,
        "eigvecs": eigvecs,
    }
