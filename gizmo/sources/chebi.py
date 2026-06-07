"""
ChEBI REST client for enriching metabolite structural data.

License: ChEBI data is CC BY 4.0 — https://www.ebi.ac.uk/chebi/aboutChebiForward.do
API docs: https://www.ebi.ac.uk/chebi/webServices.do (SOAP) and
          https://www.ebi.ac.uk/ols4/api (OLS4 REST for ontology lookups)

We prefer the lightweight REST-ish endpoint:
  https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:XXXXX  (XML)
and the programmatic OLS4 API for bulk lookups.
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import requests

log = logging.getLogger(__name__)

_OLS_BASE = "https://www.ebi.ac.uk/ols4/api/"


class ChebiClient:
    """Fetch structural metadata for ChEBI compounds via OLS4 / ChEBI web services."""

    def __init__(self, timeout: int = 30) -> None:
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})
        self.timeout = timeout

    def get_entity(self, chebi_id: str) -> Optional[dict]:
        """
        Return a dict with name, formula, inchikey, smiles, charge for a ChEBI ID.
        chebi_id: "CHEBI:15422" or "15422"
        """
        numeric = re.sub(r"^CHEBI:", "", chebi_id, flags=re.IGNORECASE)
        url = f"https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity"
        # Use the lightweight JSON-LD OLS4 endpoint instead
        ols_url = f"{_OLS_BASE}ontologies/chebi/terms"
        iri = f"http://purl.obolibrary.org/obo/CHEBI_{numeric}"
        try:
            resp = self.session.get(
                ols_url,
                params={"iri": iri},
                timeout=self.timeout,
            )
            resp.raise_for_status()
            data = resp.json()
            terms = data.get("_embedded", {}).get("terms", [])
            if not terms:
                return None
            term = terms[0]
            annotations = term.get("annotation", {})
            return {
                "chebi_id": f"CHEBI:{numeric}",
                "name": term.get("label", ""),
                "formula": _first(annotations.get("has_formula")),
                "inchikey": _first(annotations.get("has_inchikey")),
                "smiles": _first(annotations.get("has_smiles")),
                "charge": _int_or_none(_first(annotations.get("has_charge"))),
                "inchi": _first(annotations.get("has_inchi")),
            }
        except requests.RequestException as exc:
            log.warning("ChEBI lookup failed for %s: %s", chebi_id, exc)
            return None

    def enrich_batch(self, chebi_ids: list[str]) -> dict[str, dict]:
        """Return {chebi_id: entity_dict} for a list of ChEBI IDs."""
        return {cid: ent for cid in chebi_ids if (ent := self.get_entity(cid)) is not None}


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def _first(val: Optional[list]) -> Optional[str]:
    if val and isinstance(val, list):
        return val[0]
    return val  # type: ignore[return-value]


def _int_or_none(val: Optional[str]) -> Optional[int]:
    try:
        return int(val) if val is not None else None
    except (ValueError, TypeError):
        return None
