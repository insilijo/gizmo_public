from __future__ import annotations

import os
from pathlib import Path


def resource_root() -> Path:
    base = os.environ.get("GRANDMA_DATA_DIR")
    return Path(base) / "gizmo" if base else Path("data/resources/gizmo")


def first_existing(*paths: str | Path) -> str:
    for path in paths:
        p = Path(path)
        if p.exists():
            return str(p)
    return str(Path(paths[0]))


def metabolon_csv_default() -> str:
    return first_existing(
        resource_root() / "sources" / "metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv",
        Path("gizmo/sources/metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv"),
    )


def metanetx_prop_default() -> str:
    return first_existing(
        resource_root() / "metanetx" / "chem_prop.tsv",
        Path("data/raw/metanetx/chem_prop.tsv"),
    )


def metanetx_xref_default() -> str:
    return first_existing(
        resource_root() / "metanetx" / "chem_xref.tsv",
        Path("data/raw/metanetx/chem_xref.tsv"),
    )


def overrides_default() -> str:
    return first_existing(
        resource_root() / "curation" / "metabolon_overrides.json",
        Path("data/curation/metabolon_overrides.json"),
        Path("data/curation_overrides.json"),
    )
