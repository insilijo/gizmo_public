"""Feature-name normalization for cohort-specific metabolomics outputs.

Different metabolomics platforms emit feature names in different formats:

  - GeMMA / Crohn:    "Hippuric acid_plasma", "Adipic acid_fecal"  (suffix)
  - Erawijantari:     "C02703_Serine O-sulfate"                     (KEGG-ID prefix)
  - Metabolon:        "sphingomyelin (d18:2/24:1, d18:1/24:2)*",     (parenthesized
                       "1-palmityl-2-oleoyl-GPC (O-16:0/18:1)*"       lipid + tentative *)
  - MTBLS / clean:    "Myoinositol", "Glucaric acid"                  (canonical)

This module exposes ``normalize_feature_name(name)`` which yields a
list of name candidates to try in MetaboliteMapper.map(), in priority
order. Caller should map each candidate and take the best hit.
"""
from __future__ import annotations
import re

_TISSUE_SUFFIXES = ["_plasma", "_fecal", "_serum", "_urine",
                     "_blood", "_csf", "_synovial", "_stool"]
_KEGG_PREFIX = re.compile(r"^(C\d{5})_(.+)$")
_HMDB_PREFIX = re.compile(r"^(HMDB\d+)_(.+)$")

# Lipid species-level handling intentionally NOT done here — see
# /home/jgardner/papers/lipid_aware_propagation/ROADMAP.md for the
# follow-up paper. Class-collapse loses chain-length information that
# carries discriminative biology (e.g., desaturation index, elongation
# state). The proper extension is a hierarchical lipid graph
# (species → class → modifier enzyme) with class-aware Laplacian.


def normalize_feature_name(name: str) -> list[str]:
    """Return a list of candidate names to try mapping, in priority order.

    The caller should iterate through them and accept the first
    successful mapping.
    """
    if not name:
        return []
    candidates: list[str] = []
    raw = str(name).strip()

    # Strip trailing tentative markers (Metabolon)
    cleaned = raw.rstrip("*").strip()

    # Try the original name first
    candidates.append(cleaned)

    # Strip tissue suffixes (e.g. "_plasma")
    lower = cleaned.lower()
    for suf in _TISSUE_SUFFIXES:
        if lower.endswith(suf):
            stripped = cleaned[: -len(suf)]
            candidates.append(stripped)
            break

    # KEGG-ID prefix: "C02703_Serine O-sulfate" → try "C02703" then the name
    m = _KEGG_PREFIX.match(cleaned)
    if m:
        candidates.insert(0, m.group(1))   # try KEGG ID first
        candidates.append(m.group(2))

    # HMDB prefix
    m = _HMDB_PREFIX.match(cleaned)
    if m:
        candidates.insert(0, m.group(1))
        candidates.append(m.group(2))

    # Metabolon parenthesized lipid: "sphingomyelin (d18:2/24:1)" → "sphingomyelin"
    if "(" in cleaned:
        base = cleaned.split("(")[0].strip()
        if base and base != cleaned:
            candidates.append(base)

    # Trailing-asterisk versions
    if raw != cleaned:
        candidates.append(raw)

    # De-duplicate while preserving order
    seen = set()
    unique = []
    for c in candidates:
        c = c.strip()
        if c and c not in seen:
            seen.add(c)
            unique.append(c)
    return unique


def map_with_fallback(mapper, feature_name: str,
                        pubchem_cache: dict | None = None):
    """Map a feature name through ``mapper`` using normalized candidates.

    Tier 1: try each normalized candidate against the local mapper.
    Tier 2: if all fail and ``pubchem_cache`` is provided, query PubChem
            name → InChIKey and try mapping the resulting InChIKey.

    Returns (node_id, confidence) of the best successful candidate, or
    (None, 0.0) if all attempts fail.
    """
    candidates = normalize_feature_name(feature_name)
    best_node = None
    best_conf = 0.0
    for cand in candidates:
        node, conf = mapper.map(cand)
        if node and conf > best_conf:
            best_node, best_conf = node, conf
    if best_node is not None or pubchem_cache is None:
        return best_node, best_conf

    # Tier 2: PubChem name → InChIKey lookup
    for cand in candidates:
        ik = pubchem_cache.get(cand.lower())
        if ik:
            node, conf = mapper.map(ik)
            if node and conf > best_conf:
                best_node, best_conf = node, max(0.7, conf - 0.1)
                break
    return best_node, best_conf


def build_pubchem_cache(feature_names: list, cache_path,
                          rate_limit_s: float = 0.20,
                          max_api_calls: int = 1500) -> dict:
    """Build / extend a name → InChIKey cache via PubChem REST API.

    Cached as JSON. Subsequent calls reuse the cache for previously
    resolved names; only new/missing names trigger an API call.
    """
    import json, time, urllib.request, urllib.parse
    from pathlib import Path
    cache_path = Path(cache_path)
    cache: dict[str, str] = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except Exception:
            cache = {}

    n_calls = 0
    seen_lowers = set()
    for name in feature_names:
        if not name:
            continue
        # Try the normalized candidates one by one
        for cand in normalize_feature_name(name):
            key = cand.lower()
            if key in cache or key in seen_lowers:
                continue
            seen_lowers.add(key)
            if n_calls >= max_api_calls:
                continue
            url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                   + urllib.parse.quote(cand)
                   + "/property/InChIKey/JSON")
            try:
                time.sleep(rate_limit_s)
                with urllib.request.urlopen(url, timeout=10) as r:
                    data = json.loads(r.read().decode())
                    props = data.get("PropertyTable", {}).get("Properties", [])
                    if props and props[0].get("InChIKey"):
                        cache[key] = props[0]["InChIKey"]
            except Exception:
                pass
            n_calls += 1

    cache_path.write_text(json.dumps(cache, indent=2))
    return cache
