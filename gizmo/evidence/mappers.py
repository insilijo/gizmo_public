"""
Feature → graph node mapping utilities.

MetaboliteMapper   maps metabolite identifiers (ChEBI, InChIKey, name, PubChem)
                   to MetaboliteNode IDs in the graph.

GeneMapper         maps gene identifiers (HGNC symbol, Ensembl ID)
                   to GeneNode IDs in the graph.

ClinicalMapper     maps EHR / clinical feature identifiers:
                     ICD-10 / ICD-9 codes  → disease nodes  (via xref_icd10)
                     LOINC codes           → metabolite nodes (via loinc_id)
                     Lab analyte names     → metabolite nodes (via MetaboliteMapper)

All mappers are built lazily from a GizmoGraph and cache their
indexes so they can handle large batches efficiently.
"""

from __future__ import annotations

import logging
import re
from typing import Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Common-metabolite name → InChIKey14 (first 14 chars of InChIKey).
# Used by MetaboliteMapper as a last-line fallback before fuzzy match,
# bridging free-text metabolomics names ("Carnitine", "GPC") to the
# graph's abbreviated nomenclature ("CAR", "GPCho", "ChoP", …).
#
# Sourced from PubChem / HMDB canonical InChIKeys. All keys are lowercased.
# ---------------------------------------------------------------------------
COMMON_METABOLITE_INCHIKEY14: dict[str, str] = {
    # Carnitine family
    "carnitine":             "PHIQHXFUZVPYII",
    "l-carnitine":           "PHIQHXFUZVPYII",
    # Creatine + phosphocreatine
    "creatine":              "CVSVTCORWBXHQV",
    "creatine phosphate":    "PHWBOXQYWZNQIN",
    "phosphocreatine":       "PHWBOXQYWZNQIN",
    # Glucose
    "glucose":               "WQZGKKKJIJFFOK",
    "d-glucose":             "WQZGKKKJIJFFOK",
    # Glycerophospholipid mediators
    "gpc":                   "SUHOQUVVVLNYQR",
    "glycerophosphocholine": "SUHOQUVVVLNYQR",
    "phosphocholine":        "YHHSONZFOIEMCP",
    "choline phosphate":     "YHHSONZFOIEMCP",
    # Sulfur amino acids
    "homocysteine":          "FFFHZYDWPBMWHY",
    "l-homocysteine":        "FFFHZYDWPBMWHY",
    # Inositol
    "myo-inositol":          "CDAISMWEOUEBRE",
    "inositol":              "CDAISMWEOUEBRE",
    # Acetylcholine
    "acetylcholine":         "OIPILFWXSMYKGL",
    # Sarcosine
    "sarcosine":             "FSYKKLYZXJSNPZ",
    "n-methylglycine":       "FSYKKLYZXJSNPZ",
}


# ---------------------------------------------------------------------------
# Metabolite mapper
# ---------------------------------------------------------------------------

class MetaboliteMapper:
    """
    Maps metabolite feature identifiers to graph node IDs.

    Lookup order (first hit wins):
      1. Exact node ID match  (e.g. "CHEBI:15422")
      2. chebi_id attribute
      3. hmdb_id attribute  (HMDB0000001 or short form — VMH-compatible)
      4. inchikey attribute  (exact 27-char key)
      5. inchikey connectivity layer  (first 14 chars — catches stereoisomers)
      6. pubchem_cid attribute
      7. metabolon_name / name attribute  (case-insensitive)
    """

    def __init__(self, mg) -> None:
        self._g = mg.graph
        self._chebi:        dict[str, str] = {}   # chebi_id          → node_id
        self._hmdb:         dict[str, str] = {}   # normalised hmdb_id → node_id
        self._vmh:          dict[str, str] = {}   # vmh_id            → node_id
        self._inchikey:     dict[str, str] = {}   # full IK            → node_id
        self._connectivity: dict[str, str] = {}   # IK[:14]            → node_id
        self._name:         dict[str, str] = {}   # lower name         → node_id
        self._pubchem:      dict[str, str] = {}   # str(cid)           → node_id
        self._built = False

    def _build(self) -> None:
        if self._built:
            return

        def _iter_ids(val):
            if not val:
                return
            if isinstance(val, list):
                for v in val:
                    if v:
                        yield str(v).strip()
            else:
                raw = str(val)
                # Support comma-, semicolon-, and pipe-separated IDs
                sep = ";" if ";" in raw else ("|" if "|" in raw else ",")
                for v in raw.split(sep):
                    v = v.strip()
                    if v:
                        yield v

        def _index_node(nid, attrs):
            # Support both a single chebi_id string (possibly comma-separated)
            # and a chebi_ids list for nodes with multiple identifiers.
            for chebi in _iter_ids(attrs.get("chebi_id")):
                self._chebi.setdefault(chebi, nid)
            for chebi in _iter_ids(attrs.get("chebi_ids")):
                self._chebi.setdefault(chebi, nid)

            for hmdb in _iter_ids(attrs.get("hmdb_id")):
                self._hmdb.setdefault(_normalise_hmdb(hmdb), nid)
            for hmdb in _iter_ids(attrs.get("hmdb_ids")):
                self._hmdb.setdefault(_normalise_hmdb(hmdb), nid)

            vmh = attrs.get("vmh_id")
            if vmh:
                self._vmh.setdefault(str(vmh).strip().lower(), nid)

            ik = attrs.get("inchikey") or ""
            if ik and len(ik) >= 14:
                self._inchikey.setdefault(ik, nid)
                self._connectivity.setdefault(ik[:14], nid)

            for name_field in ("name", "metabolon_name"):
                n = attrs.get(name_field)
                if n:
                    name_variants = list(_split_name_variants(n))
                    # Reactome names encode compartment as a trailing
                    # "[cytosol]" / "[nucleoplasm]" suffix; strip it so bare
                    # queries ("atp", "hypoxanthine") can resolve to the
                    # reaction-connected node instead of losing to an
                    # orphan CHEBI with the bare name.
                    bare = _COMPARTMENT_SUFFIX_RE.sub("", str(n)).strip()
                    if bare and bare != str(n).strip():
                        name_variants.append(bare)
                    # Reactome / VMH display names are often abbreviations
                    # ("L-Asp [mitochondrial matrix]", "2OG [cytosol]") that
                    # never match assay-side feature names ("Aspartic acid",
                    # "2-Oxoglutaric acid"). Look up the variant directly
                    # first (so "l-cit" → citrulline overrides "cit" →
                    # citrate), then fall back to stereo-stripped stem so
                    # "L-Asp" still expands as aspartate.
                    for variant in list(name_variants):
                        expansions = _abbreviation_expansions(variant)
                        if not expansions:
                            stem = _strip_stereo_prefix(variant)
                            for expanded in _abbreviation_expansions(stem):
                                expansions.append(expanded)
                                if stem != variant:
                                    prefix = variant[: -len(stem)].rstrip()
                                    expansions.append(
                                        f"{prefix} {expanded}".strip()
                                    )
                        for expanded in expansions:
                            name_variants.append(expanded)
                    # Index a charge-stripped alias for cofactor-style
                    # display names ("NADP+" → also "NADP"; "NAD+" → "NAD";
                    # "GMP-" → "GMP") so users can query without the +/-.
                    for variant in list(name_variants):
                        v = variant.strip()
                        if v.endswith(("+", "-")) and len(v) > 1:
                            name_variants.append(v[:-1].strip())
                    # Also index conjugate-base/parent-acid variants at
                    # build time. Many data sources use "kynurenate" while
                    # Reactome stores "kynurenic acid" — without this,
                    # features keyed on the -ate form never reach the
                    # reaction-connected node. setdefault + bucket order
                    # keeps reaction-connected nodes winning ties.
                    for variant in list(name_variants):
                        name_variants.extend(_salt_form_variants(variant))
                    for variant in name_variants:
                        if variant:
                            self._name.setdefault(variant.lower(), nid)

            # PubChem synonyms (populated by the enrich_pubchem_synonyms
            # management command) — treat them as additional names, with
            # salt-form variants indexed alongside.
            syns = attrs.get("synonyms")
            if isinstance(syns, (list, tuple)):
                for s in syns:
                    variants = list(_split_name_variants(s))
                    for v in list(variants):
                        variants.extend(_salt_form_variants(v))
                    for v in variants:
                        if v:
                            self._name.setdefault(v.lower(), nid)

            cid = attrs.get("pubchem_cid")
            if cid:
                self._pubchem.setdefault(str(int(float(cid))), nid)

        # Precompute which metabolite nodes actually participate in a
        # reaction (have at least one substrate/product edge). The mapper
        # uses this to prefer biologically active nodes over orphan
        # metadata-only duplicates — otherwise a feature like "kynurenate"
        # routes to the Metabolon-sourced CHEBI node with no reaction
        # edges, and ORA / pathway membership sees zero overlap.
        connected: set[str] = set()
        for u, v, attrs in self._g.edges(data=True):
            role = attrs.get("edge_type") or attrs.get("role") or ""
            if role not in ("substrate", "product"):
                continue
            for nid in (u, v):
                if self._g.nodes.get(nid, {}).get("node_type") == "metabolite":
                    connected.add(nid)

        # Four-pass build, priority order (highest wins). Reaction
        # connectivity dominates — a compartment-suffixed node that
        # participates in a reaction still beats an orphan CHEBI node
        # with the bare name, because the compartment-stripped name
        # variant is indexed from pass 2 onward.
        #   1. Has reaction edges + non-compartment node_id.
        #   2. Has reaction edges + compartment node_id.
        #   3. Orphan + non-compartment node_id.
        #   4. Orphan + compartment node_id.
        # setdefault guards against later overwrites, so the earliest
        # pass claims each index slot.
        buckets: list[list[tuple[str, dict]]] = [[], [], [], []]
        for nid, attrs in self._g.nodes(data=True):
            if attrs.get("node_type") != "metabolite":
                continue
            is_compartment = bool(_COMPARTMENT_RE.search(nid))
            has_rxn = nid in connected
            idx = (0 if has_rxn and not is_compartment
                   else 1 if has_rxn
                   else 2 if not is_compartment
                   else 3)
            buckets[idx].append((nid, attrs))
        for bucket in buckets:
            for nid, attrs in bucket:
                _index_node(nid, attrs)

        self._built = True
        log.debug(
            "MetaboliteMapper built: %d chebi, %d hmdb, %d vmh, %d inchikey, %d names "
            "(reaction-connected metabolites: %d)",
            len(self._chebi), len(self._hmdb), len(self._vmh), len(self._inchikey),
            len(self._name), len(connected),
        )

    def map(self, feature_id: str) -> tuple[Optional[str], float]:
        """
        Map a feature identifier to a graph node ID.

        Accepts a single identifier or multiple identifiers separated by
        semicolons, pipes, or commas (e.g. "CHEBI:12345; CHEBI:67890").
        Returns the best match across all tokens.

        Returns
        -------
        (node_id, confidence)
            node_id    : matched node, or None
            confidence : 1.0 exact / 0.9 connectivity / 0.8 name / 0.7 pubchem
        """
        self._build()
        raw = feature_id.strip()
        # Multi-ID: try each token and return the best hit
        sep = ";" if ";" in raw else ("|" if "|" in raw else None)
        if sep:
            best_node, best_conf = None, 0.0
            for token in raw.split(sep):
                token = token.strip()
                if token:
                    node, conf = self._map_single(token)
                    if conf > best_conf:
                        best_node, best_conf = node, conf
            return best_node, best_conf
        return self._map_single(raw)

    def _map_single(self, fid: str) -> tuple[Optional[str], float]:
        """Map a single (non-delimited) identifier to a graph node."""
        g = self._g

        # 1. Exact node ID — but redirect to a reaction-connected twin
        # if the literal node is an orphan. The index was built with
        # reaction-connected metabolites claiming each key first, so
        # self._chebi[norm(fid)] already points to the preferred node.
        if fid in g:
            if g.nodes[fid].get("node_type") == "metabolite":
                cnorm = _normalise_chebi(fid)
                if cnorm:
                    preferred = self._chebi.get(cnorm)
                    if preferred and preferred != fid and preferred in g:
                        return preferred, 1.0
            return fid, 1.0

        # 2. ChEBI normalisation (CHEBI:NNNNN or chebi:NNNNN)
        chebi_norm = _normalise_chebi(fid)
        if chebi_norm:
            hit = self._chebi.get(chebi_norm)
            if hit:
                return hit, 1.0
            if chebi_norm in g:
                return chebi_norm, 1.0

        # 3. HMDB (HMDB0000001 or short form HMDB1)
        if fid.upper().startswith("HMDB"):
            hmdb_norm = _normalise_hmdb(fid)
            hit = self._hmdb.get(hmdb_norm)
            if hit:
                return hit, 1.0

        # 3b. VMH abbreviation (e.g. "glc_D", "phe_L")
        hit = self._vmh.get(fid.lower())
        if hit:
            return hit, 1.0

        # 4. InChIKey exact
        if _looks_like_inchikey(fid):
            hit = self._inchikey.get(fid)
            if hit:
                return hit, 1.0
            # 5. Connectivity layer
            hit = self._connectivity.get(fid[:14])
            if hit:
                return hit, 0.9

        # 6. PubChem CID (numeric string or "CID:NNNN")
        cid = _parse_pubchem(fid)
        if cid:
            hit = self._pubchem.get(cid)
            if hit:
                return hit, 0.85

        # 7. Exact name / synonym / variant (case-insensitive) — treated as
        #    a full-confidence hit since the strings match character-for-character.
        hit = self._name.get(fid.lower())
        if hit:
            return hit, 1.0

        # 7b. Conjugate-base / salt-form normalisation. Many metabolomics
        #     assays report the conjugate base (e.g. "kynurenate") while the
        #     graph indexes the parent acid ("kynurenic acid"). Try a small
        #     set of deterministic rewrites before dropping to fuzzy.
        for alt in _salt_form_variants(fid):
            hit = self._name.get(alt.lower())
            if hit:
                # Slightly below 1.0 so a literal-name hit still wins.
                return hit, 0.95

        # 7c. Common-metabolite synonym → InChIKey14 lookup. Reactome stores
        #     many bread-and-butter metabolites under abbreviated names (CAR,
        #     GPCho, ChoP, HCYS, AcCho, SARC, …) that callers don't know.
        #     Translate the free-text name to a canonical InChIKey14 and
        #     hit the connectivity layer.
        ik14 = COMMON_METABOLITE_INCHIKEY14.get(fid.lower())
        if ik14:
            hit = self._connectivity.get(ik14)
            if hit:
                return hit, 0.9

        # 8. Fuzzy name fallback — catches near-misses and partial hits like
        #    "arabitol" vs an indexed variant of "arabitol/xylitol", or minor
        #    spelling differences. Confidence capped at 0.75 so it ranks below
        #    exact matches; tokens must share a ratio >= 0.85.
        fuzzy_hit, fuzzy_conf = self._fuzzy_name_match(fid)
        if fuzzy_hit:
            return fuzzy_hit, fuzzy_conf

        return None, 0.0

    _FUZZY_CUTOFF = 0.85
    _FUZZY_MIN_SHARED = 5  # require a >= 5-char substring shared with candidate

    def _fuzzy_name_match(self, fid: str) -> tuple[Optional[str], float]:
        """Return the closest name-indexed node using SequenceMatcher.

        Two guards keep the fuzzy fallback from producing biologically
        nonsensical hits:

        1. Substring filter (pre-fuzzy): the query and candidate must share
           at least one substring of length ``_FUZZY_MIN_SHARED``. Without
           this, "Aspartic acid" routed to "acetate" via 3-char prefix
           overlap and a 0.75 SequenceMatcher ratio. Sharing a 5-char
           window keeps legitimate near-misses ("arabinitol"/"arabitol",
           "kynurenate"/"kynurenic acid") while killing the false matches.

        2. SequenceMatcher cutoff (``_FUZZY_CUTOFF``) enforces overall
           similarity once a candidate clears the substring filter.
        """
        from difflib import SequenceMatcher

        q = (fid or "").strip().lower()
        # Short queries (4 chars or fewer) are typically cofactor/nucleotide
        # codes (NADP, ATP, GMP, FAD). Fuzzy matching them produces
        # biologically wrong results — NADP fuzzes to NADPH, ATP to GTP — so
        # we require an exact-name or salt-form hit for short codes.
        if len(q) < 5 or not self._name:
            return None, 0.0

        # Pre-filter on a token or prefix overlap to bound the work. "gluc" in
        # "glucose", "glc", "d-glucose" all share a short character window.
        tokens = {t for t in q.replace("-", " ").split() if len(t) >= 3}
        prefix = q[:3]
        candidates = []
        for name in self._name:
            if name == q:
                continue
            if name.startswith(prefix) or prefix in name:
                candidates.append(name)
                continue
            if q in name or name in q:
                candidates.append(name)
                continue
            if tokens and any(t in name for t in tokens):
                candidates.append(name)

        # Cap work for pathological short queries that match too many names.
        if len(candidates) > 1500:
            candidates = candidates[:1500]

        # 5-char shared-substring set for the query (collapsed across all
        # tokens of length >= MIN). Compared via substring containment, not
        # exact-token match, so "arabinitol" / "arabitol" still pass via
        # "arabi". Hyphens are normalised to spaces first.
        q_norm = q.replace("-", " ")
        m = self._FUZZY_MIN_SHARED
        q_subs: set[str] = set()
        for tok in q_norm.split():
            if len(tok) >= m:
                for i in range(len(tok) - m + 1):
                    q_subs.add(tok[i : i + m])

        best_name, best_ratio = None, 0.0
        for name in candidates:
            cand_norm = name.replace("-", " ")
            # Substring filter: candidate must share at least one M-char
            # window with the query. q_subs is non-empty whenever the
            # query has any token of length >= M, which is guaranteed
            # here because q itself is >= M (early return above).
            if q_subs and not any(sub in cand_norm for sub in q_subs):
                continue
            ratio = SequenceMatcher(None, q, name).ratio()
            if ratio > best_ratio:
                best_ratio = ratio
                best_name = name

        if best_name and best_ratio >= self._FUZZY_CUTOFF:
            # Fuzzy tier sits between exact-name (1.0) and not-found. We
            # return the real ratio, capped at 0.95 so an unambiguous exact
            # match in a later query still wins decisively.
            return self._name[best_name], round(min(0.95, best_ratio), 2)
        return None, 0.0

    def map_batch(
        self,
        feature_ids: list[str],
    ) -> dict[str, tuple[Optional[str], float]]:
        """Map a list of feature IDs, returning {feature_id: (node_id, confidence)}."""
        return {fid: self.map(fid) for fid in feature_ids}

    def map_all_compartments(
        self, feature_id: str,
    ) -> list[tuple[str, float]]:
        """Map a feature to ALL graph nodes representing the same molecule.

        A metabolomics assay typically can't distinguish compartments —
        it sees a pooled signal. But Reactome encodes one node per
        ``(metabolite, compartment)`` pair, with reactions edge-connected
        only to their compartment-specific instance. The default
        :meth:`map` picks a single compartment via bucket priority, so
        any reaction in a *different* compartment never receives the
        observation as evidence (e.g. cytosolic α-KG mapped to the
        mitochondrial ISCIT prevents propagation to peroxisomal IDH).

        This method returns every node sharing the matched node's
        first-14-character InChIKey (structural connectivity layer,
        stereochemistry-agnostic). Each retains the original mapping
        confidence — they're the same molecule, the assay just can't
        distinguish compartments.

        Returns
        -------
        list[(node_id, confidence)]
            Empty list if the feature doesn't map at all.
            Single-element list if the matched node has no InChIKey
            (can't expand without a structural identifier).
        """
        self._build()
        nid, conf = self.map(feature_id)
        if not nid:
            return []
        attrs = self._g.nodes.get(nid, {})
        ik = attrs.get("inchikey") or ""
        if len(ik) < 14:
            return [(nid, conf)]
        ik14 = ik[:14]
        matches: list[tuple[str, float]] = []
        for other_nid, other_attrs in self._g.nodes(data=True):
            if other_attrs.get("node_type") != "metabolite":
                continue
            other_ik = other_attrs.get("inchikey") or ""
            if len(other_ik) >= 14 and other_ik[:14] == ik14:
                matches.append((other_nid, conf))
        # Make the original-bucket-winner come first so callers that
        # only consume the head get the priority node.
        matches.sort(key=lambda x: (x[0] != nid, x[0]))
        return matches or [(nid, conf)]

    def hmdb_to_vmh_id(self, hmdb_id: str) -> Optional[str]:
        """Return the VMH abbreviation for a node matched by *hmdb_id*, or None.

        This is a convenience wrapper for the common GeMMA use-case of translating
        Metabolon HMDB IDs into VMH abbreviations via the graph (rather than
        querying VMH directly).  The graph must have been enriched with VMH IDs
        via :func:`~gizmo.sources.vmh.enrich_graph_vmh` first.
        """
        self._build()
        nid = self._hmdb.get(_normalise_hmdb(hmdb_id))
        if nid is None:
            return None
        return self._g.nodes[nid].get("vmh_id")

    def build_hmdb_to_vmh_map(self, hmdb_ids: list[str]) -> dict[str, Optional[str]]:
        """Batch HMDB → VMH abbreviation lookup via the graph.

        Parameters
        ----------
        hmdb_ids:
            List of HMDB IDs (any normalisation accepted).

        Returns
        -------
        ``{hmdb_id: vmh_abbreviation_or_None}``
        """
        return {h: self.hmdb_to_vmh_id(h) for h in hmdb_ids}


# ---------------------------------------------------------------------------
# Gene mapper
# ---------------------------------------------------------------------------

class GeneMapper:
    """
    Maps gene/protein identifiers to GeneNode IDs in the graph.

    Lookup order:
      1. Exact node ID  (e.g. "ENSG:ENSG00000139618", "symbol:BRCA2")
      2. symbol attribute  (exact, case-insensitive)
      3. ensembl_id attribute  (ENSG prefix)
      4. hgnc_id attribute
    """

    def __init__(self, mg) -> None:
        self._g      = mg.graph
        self._symbol: dict[str, str] = {}   # lower symbol → node_id
        self._ensembl: dict[str, str] = {}  # ENSG ID      → node_id
        self._hgnc:    dict[str, str] = {}  # HGNC:NNNNN   → node_id
        self._built = False

    def _build(self) -> None:
        if self._built:
            return
        for nid, attrs in self._g.nodes(data=True):
            if attrs.get("node_type") != "gene":
                continue
            sym = attrs.get("symbol")
            if sym:
                self._symbol.setdefault(sym.lower(), nid)
            ensg = attrs.get("ensembl_id")
            if ensg:
                self._ensembl.setdefault(ensg, nid)
            hgnc = attrs.get("hgnc_id")
            if hgnc:
                self._hgnc.setdefault(str(hgnc), nid)
        self._built = True
        log.debug(
            "GeneMapper built: %d symbols, %d ensembl",
            len(self._symbol), len(self._ensembl),
        )

    def map(self, feature_id: str) -> tuple[Optional[str], float]:
        """Map a gene identifier.  Returns (node_id, confidence).

        Accepts a single identifier or multiple identifiers separated by
        semicolons or pipes (e.g. "BRCA2; TP53"). Returns the best match.
        """
        self._build()
        raw = feature_id.strip()
        sep = ";" if ";" in raw else ("|" if "|" in raw else None)
        if sep:
            best_node, best_conf = None, 0.0
            for token in raw.split(sep):
                token = token.strip()
                if token:
                    node, conf = self._map_single(token)
                    if conf > best_conf:
                        best_node, best_conf = node, conf
            return best_node, best_conf
        return self._map_single(raw)

    def _map_single(self, fid: str) -> tuple[Optional[str], float]:
        """Map a single (non-delimited) gene identifier."""
        # 1. Exact node ID
        if fid in self._g:
            return fid, 1.0

        # 2. symbol:XXX shorthand
        if fid.startswith("symbol:"):
            sym = fid[7:].lower()
            hit = self._symbol.get(sym)
            return (hit, 1.0) if hit else (None, 0.0)

        # 3. Ensembl ENSG
        if fid.startswith("ENSG") or fid.startswith("ensg"):
            hit = self._ensembl.get(fid.upper().replace("ENSG:", "ENSG"))
            if not hit:
                hit = self._ensembl.get(fid)
            return (hit, 1.0) if hit else (None, 0.0)

        # 4. HGNC
        if fid.startswith("HGNC:") or fid.startswith("hgnc:"):
            hit = self._hgnc.get(fid.upper())
            return (hit, 1.0) if hit else (None, 0.0)

        # 5. Try as bare gene symbol
        hit = self._symbol.get(fid.lower())
        if hit:
            return hit, 0.95

        return None, 0.0

    def map_batch(
        self,
        feature_ids: list[str],
    ) -> dict[str, tuple[Optional[str], float]]:
        return {fid: self.map(fid) for fid in feature_ids}


# ---------------------------------------------------------------------------
# Clinical mapper
# ---------------------------------------------------------------------------

_ICD10_RE = re.compile(r"^[A-Z]\d{1,2}(?:\.\d{1,4})?$")
_ICD9_RE  = re.compile(r"^\d{3}(?:\.\d{1,2})?[A-Z]?$")
_LOINC_RE = re.compile(r"^\d{1,5}-\d$")


class ClinicalMapper:
    """
    Maps EHR / clinical feature identifiers to graph nodes.

    Lookup order (first hit wins):
      1. ICD-10 code  (e.g. "E11.9", "K70.3", "ICD10:E11.9")  → disease node
      2. ICD-9  code  (e.g. "250.00")                          → disease node
      3. MONDO / OMIM / Orphanet ID                            → disease node
      4. LOINC code   (e.g. "2345-7")                          → metabolite node
      5. Lab analyte name / synonym                            → metabolite node
         (delegates to MetaboliteMapper for ChEBI/HMDB/name lookup)

    Returns
    -------
    (node_id, node_type, confidence) — node_type is "disease" or "metabolite"
    """

    def __init__(self, mg) -> None:
        self._g = mg.graph
        self._icd10:    dict[str, str] = {}   # normalised ICD-10 → node_id
        self._loinc:    dict[str, str] = {}   # LOINC code        → node_id
        self._dis_id:   dict[str, str] = {}   # MONDO/OMIM/Orphanet → node_id
        self._hpo_id:   dict[str, str] = {}   # HP:XXXXXXX → node_id
        self._met_mapper = MetaboliteMapper(mg)
        self._built = False

    def _build(self) -> None:
        if self._built:
            return
        for nid, attrs in self._g.nodes(data=True):
            ntype = attrs.get("node_type")
            if ntype == "disease":
                # Index xref_icd10 list stored as ["ICD10:E11.9", ...]
                for xref in (attrs.get("xref_icd10") or []):
                    code = xref.upper().replace("ICD10CM:", "").replace("ICD10:", "").strip()
                    self._icd10.setdefault(code, nid)
                    self._icd10.setdefault(code.replace(".", ""), nid)
                # Index MONDO / OMIM / Orphanet IDs for direct lookups
                for did in (attrs.get("xref_omim") or []):
                    self._dis_id.setdefault(did.upper(), nid)
                for did in (attrs.get("xref_orphanet") or []):
                    self._dis_id.setdefault(did.upper(), nid)
                # Node ID itself (MONDO:XXXXXXX)
                self._dis_id.setdefault(nid.upper(), nid)
            elif ntype == "phenotype":
                # HP:XXXXXXX → node_id
                self._hpo_id.setdefault(nid.upper(), nid)
                hpo = attrs.get("hpo_id")
                if hpo:
                    self._hpo_id.setdefault(hpo.upper(), nid)
            elif ntype == "metabolite":
                loinc = attrs.get("loinc_id")
                if loinc:
                    self._loinc.setdefault(str(loinc).strip(), nid)
        self._built = True
        log.debug(
            "ClinicalMapper built: %d ICD-10, %d disease IDs, %d HPO, %d LOINC",
            len(self._icd10), len(self._dis_id), len(self._hpo_id), len(self._loinc),
        )

    def map(self, feature_id: str) -> tuple[Optional[str], Optional[str], float]:
        """
        Returns (node_id, node_type, confidence).
        node_type is "disease", "metabolite", or None (no match).
        """
        self._build()
        fid = feature_id.strip()
        fid_up = fid.upper()

        # Strip common prefixes
        for pfx in ("ICD10CM:", "ICD10:", "ICD9:", "ICD-10:", "ICD-9:"):
            if fid_up.startswith(pfx):
                fid = fid[len(pfx):].strip()
                fid_up = fid.upper()
                break

        # 1. ICD-10
        if _ICD10_RE.match(fid_up):
            hit = self._icd10.get(fid_up) or self._icd10.get(fid_up.replace(".", ""))
            if hit:
                return hit, "disease", 1.0

        # 2. ICD-9
        if _ICD9_RE.match(fid_up):
            hit = self._icd10.get(fid_up)  # some ICD-9s may be xref'd
            if hit:
                return hit, "disease", 0.9

        # 3. HPO phenotype (HP:XXXXXXX)
        if fid_up.startswith("HP:"):
            hit = self._hpo_id.get(fid_up)
            if hit:
                return hit, "phenotype", 1.0

        # 4. Disease ontology ID (MONDO:, OMIM:, Orphanet:)
        hit = self._dis_id.get(fid_up)
        if hit:
            return hit, "disease", 1.0

        # 4. LOINC
        if _LOINC_RE.match(fid):
            hit = self._loinc.get(fid)
            if hit:
                return hit, "metabolite", 0.9

        # 5. Metabolite name / synonym / ChEBI / HMDB
        nid, conf = self._met_mapper.map(fid)
        if nid:
            return nid, "metabolite", conf

        return None, None, 0.0

    def map_batch(
        self,
        feature_ids: list[str],
    ) -> dict[str, tuple[Optional[str], Optional[str], float]]:
        return {fid: self.map(fid) for fid in feature_ids}


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

_CHEBI_RE      = re.compile(r"chebi[:\s_]?(\d+)", re.IGNORECASE)
_COMPARTMENT_RE = re.compile(r"@")  # matches GIZMO's @compartment node ID convention
# Reactome encodes compartment on names as a trailing "[cytosol]" etc.
_COMPARTMENT_SUFFIX_RE = re.compile(r"\s*\[[^\[\]]+\]\s*$")
_INCHIKEY_RE   = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_PUBCHEM_RE    = re.compile(r"(?:pubchem[:\s_]?|cid[:\s_]?)(\d+)", re.IGNORECASE)


def _normalise_chebi(s: str) -> Optional[str]:
    m = _CHEBI_RE.match(s.strip())
    if m:
        return f"CHEBI:{m.group(1)}"
    return None


def _normalise_hmdb(s: str) -> str:
    """Normalise HMDB IDs to 7-digit zero-padded VMH format (HMDB0000001)."""
    s = s.strip().upper()
    if s.startswith("HMDB"):
        return f"HMDB{s[4:].zfill(7)}"
    return s


def _looks_like_inchikey(s: str) -> bool:
    return bool(_INCHIKEY_RE.match(s.strip()))


def _parse_pubchem(s: str) -> Optional[str]:
    m = _PUBCHEM_RE.match(s.strip())
    if m:
        return m.group(1)
    # bare integer
    if s.strip().isdigit():
        return s.strip()
    return None


_NAME_SPLIT_RE = re.compile(r"\s*(?:/|;|\s+or\s+|\s*\|\s*)\s*", re.IGNORECASE)

# Strip leading stereo / configuration prefix so "L-Asp" / "(R)-2HG" /
# "D-glucose" all reduce to a comparable stem for abbrev expansion.
_STEREO_PREFIX_RE = re.compile(
    r"^(?:\(?[RSDL]\)?-|\(?[+-]\)?-|alpha-|beta-|gamma-|delta-|cis-|trans-|α-|β-|γ-|δ-)",
    re.IGNORECASE,
)

# Biochemistry abbreviations that appear as Reactome / VMH display names but
# are unrecognisable to assays that report the canonical full name.
# Keys are lowercase stems (after compartment + stereoprefix stripping).
# Values are alternate names to register in the name index alongside the
# original display name. Salt-form variants are auto-derived from each entry,
# so we list parent-acid forms here and let _salt_form_variants emit the
# conjugate-base ("-ate") forms automatically.
_NAME_ABBREVIATIONS: dict[str, list[str]] = {
    # Standard amino acid 3-letter codes
    "ala": ["alanine"],
    "arg": ["arginine"],
    "asn": ["asparagine"],
    "asp": ["aspartate", "aspartic acid"],
    "cys": ["cysteine"],
    "gln": ["glutamine"],
    "glu": ["glutamate", "glutamic acid"],
    "gly": ["glycine"],
    "his": ["histidine"],
    "ile": ["isoleucine"],
    "leu": ["leucine"],
    "lys": ["lysine"],
    "met": ["methionine"],
    "phe": ["phenylalanine"],
    "pro": ["proline"],
    "ser": ["serine"],
    "thr": ["threonine"],
    "trp": ["tryptophan"],
    "tyr": ["tyrosine"],
    "val": ["valine"],
    # TCA / glycolysis intermediates (Reactome short forms)
    "2og": [
        "2-oxoglutarate", "2-oxoglutaric acid",
        "alpha-ketoglutarate", "alpha-ketoglutaric acid",
        "2-oxopentanedioic acid", "oxoglutaric acid",
    ],
    "akg": [
        "alpha-ketoglutarate", "alpha-ketoglutaric acid",
        "2-oxoglutarate", "2-oxoglutaric acid",
    ],
    "iscit": ["isocitrate", "isocitric acid"],
    "cit": ["citrate", "citric acid"],
    # L-Cit / D-Cit in Reactome are citrulline, NOT citrate. Override the
    # stereoprefix-stripped form so "L-Cit" doesn't collapse to "cit".
    "l-cit": ["L-citrulline", "citrulline"],
    "d-cit": ["D-citrulline"],
    "oaa": ["oxaloacetate", "oxaloacetic acid"],
    "pyr": ["pyruvate", "pyruvic acid"],
    "lact": ["lactate", "lactic acid"],
    "mal": ["malate", "malic acid"],
    "fum": ["fumarate", "fumaric acid"],
    "fuma": ["fumarate", "fumaric acid"],
    "suc": ["succinate", "succinic acid"],
    "succ": ["succinate", "succinic acid"],
    "2hg": [
        "2-hydroxyglutarate", "2-hydroxyglutaric acid",
        "2-hydroxypentanedioic acid",
    ],
    # Cofactor short forms
    "ac-coa": ["acetyl-coa", "acetyl coa", "acetyl-coenzyme a", "acetylcoa"],
    "accoa":  ["acetyl-coa", "acetyl coa", "acetyl-coenzyme a"],
    # Glycolysis / PPP
    "g6p": ["glucose 6-phosphate", "glucose-6-phosphate"],
    "f6p": ["fructose 6-phosphate", "fructose-6-phosphate"],
    "f1,6bp": ["fructose 1,6-bisphosphate"],
    "pep": ["phosphoenolpyruvate", "phosphoenolpyruvic acid"],
    "3pg": ["3-phosphoglycerate", "3-phosphoglyceric acid"],
    "2pg": ["2-phosphoglycerate", "2-phosphoglyceric acid"],
    "r5p": ["ribose 5-phosphate", "ribose-5-phosphate"],
    "x5p": ["xylulose 5-phosphate", "xylulose-5-phosphate"],
    "e4p": ["erythrose 4-phosphate", "erythrose-4-phosphate"],
    "s7p": ["sedoheptulose 7-phosphate", "sedoheptulose-7-phosphate"],
    "6pg": ["6-phosphogluconate", "6-phosphogluconic acid"],
}


def _abbreviation_expansions(stem: str) -> list[str]:
    """Return canonical full-name variants for a known biochem abbreviation.

    Empty list if the stem isn't a recognised abbreviation. Compartment
    suffix and stereo prefix should be stripped from the input first.
    """
    key = stem.strip().lower()
    if not key:
        return []
    return list(_NAME_ABBREVIATIONS.get(key, []))


def _strip_stereo_prefix(name: str) -> str:
    """Strip a leading L-/D-/R-/S-/cis-/alpha- prefix from a display name."""
    s = (name or "").strip()
    while True:
        m = _STEREO_PREFIX_RE.match(s)
        if not m:
            return s
        s = s[m.end():].strip()


def _salt_form_variants(name: str) -> list[str]:
    """Deterministic conjugate-base ↔ parent-acid rewrites.

    Metabolomics assays often report the salt ("kynurenate", "succinate")
    while knowledge graphs index the parent acid ("kynurenic acid",
    "succinic acid"). Attempting a few cheap string rewrites avoids pushing
    every mismatch down to the expensive fuzzy path.
    """
    s = (name or "").strip().lower()
    if not s:
        return []
    out: list[str] = []

    def _push(v: str) -> None:
        v = v.strip()
        if v and v != s and v not in out:
            out.append(v)

    # "-ate" → parent acid forms. "kynurenate" → "kynurenic acid" or
    # "kynurenoic acid"; "succinate" → "succinic acid".
    if s.endswith("ate"):
        stem = s[:-3]
        _push(stem + "ic acid")
        _push(stem + "oic acid")
    # Reverse direction: "kynurenic acid" → "kynurenate" (cheap to try too).
    if s.endswith("ic acid"):
        _push(s[:-len("ic acid")] + "ate")
    if s.endswith("oic acid"):
        _push(s[:-len("oic acid")] + "oate")
    # "-ite" / "-ous acid" (nitrite ↔ nitrous acid). Rare in metabolomics
    # but cheap to cover.
    if s.endswith("ite"):
        _push(s[:-3] + "ous acid")
    if s.endswith("ous acid"):
        _push(s[:-len("ous acid")] + "ite")
    return out


def _split_name_variants(name) -> list[str]:
    """
    Split a compound name that encodes multiple synonyms into individual
    names — e.g. 'arabitol/xylitol' → ['arabitol/xylitol', 'arabitol', 'xylitol'].

    The original string is always returned first so that an exact-name
    lookup still works. Comma splitting is intentionally NOT applied, since
    many legitimate chemical names contain commas (e.g. '2,3-dihydroxybutyric
    acid').
    """
    s = str(name or "").strip()
    if not s:
        return []
    variants: list[str] = [s]
    if _NAME_SPLIT_RE.search(s):
        for part in _NAME_SPLIT_RE.split(s):
            part = part.strip()
            if part and part not in variants:
                variants.append(part)
    return variants
