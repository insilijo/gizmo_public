"""
Currency metabolite detection and flagging.

Currency metabolites are ubiquitous cofactors/co-substrates that participate
in hundreds of reactions. Including them in metabolite-metabolite co-occurrence
graphs creates spurious edges (everything is connected via ATP/NAD+).

Two-tier strategy:

  1. **Hard currency** — metabolites that are ALWAYS cofactors regardless of
     context (water, ATP, NAD+/H, CoA, O2, CO2, Pi, PPi, H+, FAD/H2). Flagged
     globally on the node via ``is_currency=True``; ``run_bayesian_inference``
     drops these nodes from the BP candidate set entirely.

  2. **Conditional currency** — metabolites that act as cofactors in some
     reaction classes but are the focal feature in others. α-ketoglutarate
     is a transamination cosubstrate AND the substrate of α-KG-dependent
     dioxygenases (HDMs, prolyl hydroxylases, JMJC family) — but inside the
     TCA cycle it IS the metabolite of interest. Same for L-glutamate,
     L-glutamine, succinate, SAM, SAH, acetyl-CoA, formaldehyde.

     For these we use ``compute_conditional_currency_edges(mg)`` which scans
     reaction names against EC-class patterns and yields the set of
     ``(metabolite_node, reaction_node)`` edges to skip in BP propagation.
     Pathway overrides un-skip edges in canonical pathways where the
     metabolite IS the focal feature (TCA cycle, glutamate metabolism,
     methionine cycle, etc.).

  3. Optionally flag by degree heuristic: metabolites with reaction-degree
     > mean + k*std (default k=3) are likely unlisted currencies.

ChEBI IDs sourced from Recon3D / published currency lists:
  Brunk et al. 2018, Thiele et al. 2013, Hucka et al. 2019.
"""

from __future__ import annotations

import re
import statistics
from typing import Optional

from gizmo.graph.network import GizmoGraph

# ---------------------------------------------------------------------------
# Canonical currency metabolite ChEBI IDs
# (compartment-agnostic — match by chebi_id attribute, not node_id)
# ---------------------------------------------------------------------------

KNOWN_CURRENCY_CHEBI: dict[str, str] = {
    # Adenine nucleotides
    "CHEBI:30616": "ATP",
    "CHEBI:456216": "ADP",
    "CHEBI:16027": "AMP",
    "CHEBI:17877": "adenosine",
    # Guanine nucleotides
    "CHEBI:37565": "GTP",
    "CHEBI:58189": "GDP",
    "CHEBI:5978":  "GMP",
    # Pyridine nucleotides — oxidised
    "CHEBI:57540": "NAD+",
    "CHEBI:58349": "NADP+",
    # Pyridine nucleotides — reduced
    "CHEBI:57945": "NADH",
    "CHEBI:57783": "NADPH",
    # Flavin nucleotides
    "CHEBI:57692": "FAD",
    "CHEBI:58307": "FADH2",
    # Coenzyme A and acyl carriers
    "CHEBI:57287": "CoA",
    "CHEBI:57288": "acetyl-CoA",
    # Phosphates
    "CHEBI:43474": "orthophosphate (Pi)",
    "CHEBI:33019": "pyrophosphate (PPi)",
    # Small inorganic / proton carriers
    "CHEBI:15377": "water",
    "CHEBI:15378": "H+ (proton)",
    "CHEBI:16526": "CO2",
    "CHEBI:15379": "O2",
    "CHEBI:29033": "H2",
    # Ubiquinones
    "CHEBI:16389": "ubiquinone (CoQ)",
    "CHEBI:17976": "ubiquinol (CoQH2)",
    # Glutamate / 2-OG (transamination hub)
    "CHEBI:16015": "L-glutamate",
    "CHEBI:16810": "2-oxoglutarate",
    # Thioredoxin system (small, but extremely high degree)
    "CHEBI:15422": "glucose-6-phosphate (sometimes)",   # not always currency; flag for review
}

# ChEBI IDs that are borderline — flag but mark as "review"
BORDERLINE_CURRENCY_CHEBI: set[str] = {
    "CHEBI:15422",   # glucose-6-phosphate (central metabolite, not always excluded)
    "CHEBI:16015",   # L-glutamate
    "CHEBI:16810",   # 2-oxoglutarate
    "CHEBI:17877",   # adenosine
}


# ---------------------------------------------------------------------------
# Canonical currency metabolite InChIKeys (first 14 characters)
#
# Why InChIKey rather than ChEBI: ChEBI IDs accumulate obsolete/secondary
# variants over time, and external sources (MetaNetX chem_xref, Reactome,
# HMDB) cross-reference different versions. The first 14 characters of an
# InChIKey are the structural-connectivity layer — they identify a molecule
# regardless of stereochemistry, charge state, or ChEBI numbering. This
# table joins reliably whatever cross-reference history a graph node
# carries.
#
# InChIKeys derived from canonical ChEBI primary IDs in KNOWN_CURRENCY_CHEBI;
# verified against PubChem / EBI ChEBI flat files.
# ---------------------------------------------------------------------------

KNOWN_CURRENCY_INCHIKEY14: dict[str, str] = {
    # Adenine nucleotides
    "ZKHQWZAMYRWXGA": "ATP",        # ZKHQWZAMYRWXGA-KQYNXXCUSA-N
    "XTWYTFMLZFPYCI": "ADP",        # XTWYTFMLZFPYCI-UHFFFAOYSA-N
    "UDMBCSSLTHHNCD": "AMP",        # UDMBCSSLTHHNCD-KQYNXXCUSA-N
    "OIRDTQYFTABQOQ": "adenosine",
    # Guanine nucleotides
    "XKMLYUALXHKNFT": "GTP",        # XKMLYUALXHKNFT-UHFFFAOYSA-N
    "QGWNDRXFNXRZMB": "GDP",
    "RQFCJASXJCIDSX": "GMP",
    # Pyridine nucleotides — oxidised
    "BAWFJGJZGIEFAR": "NAD+",       # BAWFJGJZGIEFAR-NNYOXOHSSA-N
    "XJLXINKUBYWONI": "NADP+",      # XJLXINKUBYWONI-NNYOXOHSSA-N
    # Pyridine nucleotides — reduced
    "BOPGDPNILDQYTO": "NADH",       # BOPGDPNILDQYTO-NNYOXOHSSA-N
    "ACFIXJIJDZMPPO": "NADPH",      # ACFIXJIJDZMPPO-NNYOXOHSSA-N
    # Flavin nucleotides
    "VWWQXMAJTJZDQX": "FAD",        # VWWQXMAJTJZDQX-UYBVJOGSSA-N
    "YPZRHBJKEMOYQH": "FADH2",
    # Coenzyme A and acetyl form
    "RGJOEKWQDUBAIZ": "CoA",        # RGJOEKWQDUBAIZ-IBOSZNHHSA-N
    "ZSLZBFCDCINBPY": "acetyl-CoA", # ZSLZBFCDCINBPY-ZSJPKINUSA-N
    # Phosphates
    "NBIIXXVUZAFLBC": "Pi",         # NBIIXXVUZAFLBC-UHFFFAOYSA-N
    "XPPKVPWEQAFLFU": "PPi",        # XPPKVPWEQAFLFU-UHFFFAOYSA-J
    # Small inorganic / proton carriers
    "XLYOFNOQVPJJNP": "water",      # XLYOFNOQVPJJNP-UHFFFAOYSA-N
    "GPRLSGONYQIRFK": "H+",         # GPRLSGONYQIRFK-UHFFFAOYSA-N
    "CURLTUGMZLYLDI": "CO2",        # CURLTUGMZLYLDI-UHFFFAOYSA-N
    "MYMOFIZGZYHOMD": "O2",         # MYMOFIZGZYHOMD-UHFFFAOYSA-N
    "UFHFLCQGNIYNRP": "H2",
    # Ubiquinones
    "QNTNKSLOFHEFPK": "ubiquinone",
    "RJZTYIYKBRSPGZ": "ubiquinol",
    # Conditional-currency cofactors. These are real metabolites that ARE
    # measured directly in many studies (so we want to keep them as features
    # when possible), but they also act as cosubstrates / cofactors in
    # hundreds of reactions — leaving them unflagged in BP causes massive
    # false-pooling at the cofactor hub. Listed here, excluded from
    # currency by default; pass include_borderline=True for BP / propagation
    # runs and False for ORA / single-feature scoring.
    "WHUUTDBJXJRKMK": "L-glutamate",          # transamination acceptor
    "KPGXRSRHYNQIFN": "2-oxoglutarate",       # α-KG, transamination donor + Fe(II)/α-KG dioxygenase cofactor
    "ZDXPYRJPNDTMRX": "L-glutamine",          # transamination + glutaminolysis hub
    "KDYFGRWQOYBRFD": "succinate",            # TCA + Fe(II)/α-KG dioxygenase product
    "MEFKEPWMEQBLKI": "S-adenosyl-L-methionine (SAM)",  # universal methyl donor
    "ZJUKTBDSGOFHSH": "S-adenosyl-L-homocysteine (SAH)",# methyl-transfer product
    "WSFSSNUMVMOOMR": "formaldehyde",         # HDM demethylation product
}

# Subset of InChIKeys that are borderline. Mirror BORDERLINE_CURRENCY_CHEBI
# — these are excluded from default canonical flagging but added when
# `include_borderline=True` is passed.
#
# IMPORTANT: for BP, prefer ``compute_conditional_currency_edges`` over
# the borderline-flag mechanism — the edge-level approach keeps these
# metabolites as features in their canonical pathways (TCA, glutamate
# metabolism, methionine cycle) while still removing the cofactor-hub
# artefact in α-KG dioxygenase / methyltransferase / transaminase contexts.
BORDERLINE_CURRENCY_INCHIKEY14: set[str] = {
    "WHUUTDBJXJRKMK",   # L-glutamate
    "KPGXRSRHYNQIFN",   # 2-oxoglutarate (α-KG)
    "ZDXPYRJPNDTMRX",   # L-glutamine
    "KDYFGRWQOYBRFD",   # succinate (HDM / TCA product)
    "MEFKEPWMEQBLKI",   # SAM
    "ZJUKTBDSGOFHSH",   # SAH
    "WSFSSNUMVMOOMR",   # formaldehyde
    "OIRDTQYFTABQOQ",   # adenosine
}


# ---------------------------------------------------------------------------
# Name-based currency matching
#
# Reactome compartmented metabolite nodes (e.g., "Pi [nucleoplasm]",
# "H2O [Golgi lumen]", "Mg2+ [cytosol]") frequently arrive without a
# resolved chebi_id or inchikey because the compartment-specific rendering
# is not cross-referenced. The structural-ID matching above misses these
# nodes entirely, leaving hundreds of cofactor instances unflagged.
#
# The set below captures canonical aliases for ubiquitous currency
# metabolites. _normalize_name() strips compartment suffixes and case;
# any node whose normalized name matches is flagged as currency.
# ---------------------------------------------------------------------------

CURRENCY_NAMES_NORMALIZED: set[str] = {
    # Phosphates
    "pi", "ppi", "phosphate", "orthophosphate", "diphosphate",
    "pyrophosphate", "hydrogenphosphate", "hydrogen phosphate",
    # Nucleotide phosphates (canonical + deoxy)
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "ctp", "cdp", "cmp",
    "utp", "udp", "ump", "ttp", "tdp", "tmp",
    "datp", "dadp", "damp", "dgtp", "dgdp", "dgmp", "dctp", "dcdp",
    "dcmp", "dttp", "dtdp", "dtmp", "dutp", "dudp", "dump",
    # Pyridine nucleotides
    "nad", "nad+", "nadh", "nadh2", "nadp", "nadp+", "nadph", "nadph2",
    # Flavin nucleotides
    "fad", "fadh", "fadh2", "fmn", "fmnh2", "riboflavin",
    # Coenzyme A
    "coa", "coa-sh", "coenzyme a", "hscoa",
    "acetyl-coa", "acetyl coa", "acetylcoa",
    # Small inorganics / proton carriers
    "h+", "h(+)", "proton", "h2o", "water", "h2o2", "hydrogen peroxide",
    "co2", "carbon dioxide", "co", "carbon monoxide",
    "o2", "oxygen", "dioxygen", "h2", "hydrogen",
    "no", "nitric oxide", "no2-", "nitrite", "no3-", "nitrate",
    "h2s", "hydrogen sulfide", "hs-", "sulfide",
    "hco3-", "hco3(-)", "bicarbonate", "carbonate", "co3(2-)", "co3 2-",
    "oh-", "hydroxide",
    # Ammonia / ammonium (nitrogen currency)
    "nh3", "ammonia", "nh4+", "nh4(+)", "ammonium",
    # Sulfate / sulfite
    "so4(2-)", "so4 2-", "sulfate", "so3(2-)", "sulfite",
    "ps", "paps", "3'-phosphoadenylyl sulfate",
    # Metal ions (currency cofactors)
    "mg2+", "mg(2+)", "mg++", "magnesium", "magnesium(2+)",
    "na+", "na(+)", "sodium", "sodium(1+)",
    "k+", "k(+)", "potassium", "potassium(1+)",
    "ca2+", "ca(2+)", "calcium", "calcium(2+)",
    "zn2+", "zn(2+)", "zinc", "zinc(2+)",
    "fe2+", "fe(2+)", "fe3+", "fe(3+)", "iron", "iron(2+)", "iron(3+)",
    "cu+", "cu(+)", "cu2+", "cu(2+)", "copper", "copper(1+)", "copper(2+)",
    "mn2+", "mn(2+)", "manganese", "manganese(2+)",
    "co2+", "co(2+)", "cobalt", "cobalt(2+)",
    "ni2+", "ni(2+)", "nickel",
    # Toxic/heavy metals (xenobiotic — never a focal metabolite)
    "cd2+", "cd(2+)", "cadmium", "cadmium(2+)",
    "hg2+", "hg(2+)", "mercury", "mercury(2+)", "hg+", "hg(+)",
    "pb2+", "pb(2+)", "lead", "lead(2+)",
    "as3+", "as(3+)", "arsenic", "arsenite", "arsenate",
    "al3+", "al(3+)", "aluminium", "aluminum",
    "cr3+", "cr(3+)", "cr6+", "cr(6+)", "chromium", "chromate",
    "se", "selenide",
    "tl+", "tl(+)", "thallium",
    "ba2+", "ba(2+)", "barium",
    "li+", "li(+)", "lithium",
    "cl-", "cl(-)", "chloride",
    "br-", "bromide", "i-", "iodide", "f-", "fluoride",
    # Ubiquinones / electron carriers
    "ubiquinone", "ubiquinol", "coq", "coqh2",
    "coenzyme q", "coenzyme q10", "ubiquinone-10", "ubiquinol-10",
    "cytochrome c", "ferricytochrome c", "ferrocytochrome c",
    # Glutathione (very high degree, often currency)
    "gsh", "gssg", "glutathione", "glutathione disulfide",
    # Tetrahydrofolate (one-carbon currency)
    "thf", "5-methyl-thf", "5-methyl thf", "tetrahydrofolate",
    "5-methyltetrahydrofolate", "5-formyl-thf",
    "5,10-methylene-thf", "5,10-methylenetetrahydrofolate",
    "10-formyl-thf", "10-formyltetrahydrofolate",
    "dihydrofolate", "dhf", "folate",
    # Selenium currency
    "selenocysteine",
    # Misc small-molecule cofactors
    "tpp", "thiamine pyrophosphate", "thiamine diphosphate",
    "biotin", "lipoamide", "lipoate", "lipoyl",
    "plp", "pyridoxal phosphate", "pyridoxal 5'-phosphate",
    "pmp", "pyridoxamine phosphate",
}


def _normalize_name(name: str) -> str:
    """Lowercase + strip ``[compartment]`` and trailing whitespace.

    Reactome rendering uses ``"Pi [nucleoplasm]"``, ``"H2O [Golgi lumen]"``,
    etc. — strip the bracketed compartment so the canonical name lookup
    works regardless of compartment.
    """
    if not name:
        return ""
    n = re.sub(r"\s*\[[^\]]+\]\s*$", "", name)
    return n.strip().lower()


# ---------------------------------------------------------------------------
# Reactome abstract / placeholder patterns
#
# Reactome models many "metabolites" as abstract placeholders rather than
# real chemical compounds:
#   - DNA elements: "Promotor region of beta-globin", "FoxO3a-binding Element"
#   - Immune abstracts: "Antigen", "Allergen", "Bacterial X surface pattern"
#   - Generic chemical classes: "NTP", "dNTP", "LCFA", "acyl-CoA" (no specific
#                                 acyl chain)
#   - Protein states: "Misfolded protein", "Unfolded protein",
#                      "K48polyUb-partially digested Ag"
#   - Peptide abstracts: "polypeptide", "oligopeptide fragment"
#   - Generic radicals: "ROS", "hydroperoxyl" (without specific lipid context)
#
# These are hub placeholders — they participate in many reactions but
# don't represent measurable chemistry. Treat as currency.
# ---------------------------------------------------------------------------

ABSTRACT_PATTERNS: list[re.Pattern] = [
    # DNA / response elements — pure regulatory placeholders, not chemistry.
    # These NEVER represent measurable biology, so always currency.
    re.compile(r"\b(promot[oe]r|response element|binding element|"
                r"polymerase\s+iii?\s+type\s+\d|"
                r"dna\s+containing|rdna\s+promoter|"
                r"mitochondrial\s+dna\s+promoter)\b", re.I),
    # Generic chemical-class placeholders (any-NTP, any-LCFA, etc.) —
    # represent a CLASS not a specific compound. Always currency.
    re.compile(r"^(?:ntp|dntp|ntp\(4-\)|dntp\(4-\)|"
                r"lcfa(?:s|\(-?\)|-coa)?|"
                r"short-chain\s+\w+|medium-chain\s+\w+|long-chain\s+\w+|"
                r"a\s+nucleotide\s+sugar|"
                r"a\s+(?:fatty|amino)\s+acid|"
                r"any\s+\w+|generic\s+\w+)$", re.I),
    # Generic "drug" placeholder
    re.compile(r"^(?:drug|exogenous\s+drug|substrate\s+\w+)$", re.I),
    # Note: antigens, allergens, peptide states, ubiquitinated forms,
    # ROS species, and surface-pattern recognition entities are LEFT
    # OUT — they're Reactome abstractions but they map to real
    # immunology / proteostasis / oxidative-stress biology that
    # cohorts can probe. Removing them strips real signal in
    # immune-disease cohorts (Crohn, COVID, RA).
]


def _matches_abstract_pattern(name: str) -> bool:
    """True if name matches a Reactome abstract / placeholder pattern."""
    if not name:
        return False
    n = re.sub(r"\s*\[[^\]]+\]\s*$", "", name).strip()
    return any(p.search(n) for p in ABSTRACT_PATTERNS)


# ---------------------------------------------------------------------------
# Conditional currency: edge-level skip rules
#
# A metabolite that is "currency" in one reaction class can be a feature
# in another. Rather than a global node flag, we compute the set of
# ``(metabolite_node, reaction_node)`` edges where the metabolite acts as
# a cofactor cosubstrate. BP skips just those edges.
# ---------------------------------------------------------------------------

# Reaction-name regex patterns → set of InChIKey14s to treat as currency
# in matching reactions. Patterns derived from Reactome reaction-name
# conventions ("X demethylates Y", "X transaminates Y to form Z", etc.).
EC_CLASS_PATTERNS: list[tuple[re.Pattern, set[str]]] = [
    # α-KG-dependent dioxygenases: HDMs, prolyl/asparaginyl hydroxylases,
    # JMJC family, TET enzymes, BBOX1, etc. Substrate (cofactor side):
    # α-KG + O2; products: succinate + CO2 + (formaldehyde for demethylases).
    (re.compile(r"\b(?:demethylates?|demethylation|"
                r"hydroxylates?(?!\s+(?:by|via)\s+CYP)|"
                r"prolyl[\s-]hydroxyl|asparaginyl[\s-]hydroxyl)\b", re.I),
     {"KPGXRSRHYNQIFN",   # α-KG
      "KDYFGRWQOYBRFD",   # succinate
      "WSFSSNUMVMOOMR",   # formaldehyde
      }),
    # Methyltransferases: KMTs (PKMTs), DNMTs, PEMT, GNMT, COMT, NNMT, etc.
    # Cofactors: SAM (donor), SAH (product).
    (re.compile(r"\b(?:methylates?|methylation|methyltransferase|"
                r"transmethylation)\b", re.I),
     {"MEFKEPWMEQBLKI",   # SAM
      "ZJUKTBDSGOFHSH",   # SAH
      }),
    # Transaminases: AST/GOT, ALT/GPT, BCAT, GPT2, OAT, etc.
    # Cofactors: α-KG ↔ L-glutamate (and sometimes L-glutamine).
    (re.compile(r"\b(?:transaminates?|transamination|"
                r"aminotransferase)\b", re.I),
     {"KPGXRSRHYNQIFN",   # α-KG
      "WHUUTDBJXJRKMK",   # L-glutamate
      "ZDXPYRJPNDTMRX",   # L-glutamine
      }),
    # Acetyltransferases (KATs, NATs, etc.): acetyl-CoA donor → CoA product.
    (re.compile(r"\b(?:acetylates?|acetyl[\s-]?transferase|"
                r"acetylation)\b", re.I),
     {"ZSLZBFCDCINBPY",   # acetyl-CoA
      "RGJOEKWQDUBAIZ",   # CoA
      }),
]


# Pathways where conditionally-currency metabolites are the focal FEATURE,
# not a cofactor. Inside these pathways, the EC-class skip rules above
# are reverted.
#
# Reactome IDs verified from the production graph; expand as needed.
PATHWAY_FEATURE_OVERRIDES: dict[str, set[str]] = {
    # Citric acid cycle (TCA / Krebs cycle)
    "R-HSA-71403": {
        "KPGXRSRHYNQIFN",   # 2-oxoglutarate (α-KG)
        "KDYFGRWQOYBRFD",   # succinate
        "VZCYOOQTPOCHFL",   # fumarate
        "BJEPYKJPYRNKOW",   # L-malate
        "KRKNYBCHXYNGOX",   # citrate
        "KHPXUQMNIQBQEV",   # oxaloacetate
        "GPRCLPUAAABXKL",   # cis-aconitate (for completeness)
        "GACDQMDRPRGCTN",   # isocitrate
    },
    # Glutamate and glutamine metabolism
    "R-HSA-8964539": {
        "WHUUTDBJXJRKMK",   # L-glutamate
        "ZDXPYRJPNDTMRX",   # L-glutamine
        "KPGXRSRHYNQIFN",   # α-KG (substrate of GDH etc.)
    },
    # Sulfur amino acid metabolism (methionine cycle)
    "R-HSA-1614635": {
        "MEFKEPWMEQBLKI",   # SAM
        "ZJUKTBDSGOFHSH",   # SAH
        "FFFHZYDWPBMWHY",   # L-homocysteine
        "FFEARJCKVFRZRR",   # L-methionine
    },
    # Krebs / Citric acid cycle (alt parent)
    "R-HSA-1428517": {
        "KPGXRSRHYNQIFN", "KDYFGRWQOYBRFD", "VZCYOOQTPOCHFL",
        "BJEPYKJPYRNKOW", "KRKNYBCHXYNGOX",
    },
}


def compute_conditional_currency_edges(mg: GizmoGraph) -> set[tuple[str, str]]:
    """
    Return the set of ``(metabolite_node_id, reaction_node_id)`` edges where
    the metabolite acts as a cofactor cosubstrate in this reaction's
    EC class, with pathway overrides for canonical-feature contexts.

    BP propagators should skip these edges — the metabolite node remains
    in the candidate set (it can still drive other reactions where it's
    a feature).

    Implementation: iterate reactions, match name → EC class patterns →
    set of currency InChIKey14s for this reaction, then subtract any
    overrides triggered by the reaction's pathway memberships.
    """
    g = mg.graph
    skip: set[tuple[str, str]] = set()

    for rxn_id, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "reaction":
            continue
        rxn_name = (attrs.get("name") or "").strip()
        if not rxn_name:
            continue

        # Step 1: which conditional-currency InChIKeys does this reaction's
        # EC class implicate?
        currency_iks: set[str] = set()
        for pattern, iks in EC_CLASS_PATTERNS:
            if pattern.search(rxn_name):
                currency_iks |= iks
        if not currency_iks:
            continue

        # Step 2: pathway override — if the reaction is in a "feature"
        # pathway, un-flag those metabolites for this reaction.
        rxn_pathways = set(attrs.get("pathways") or [])
        for path_id in rxn_pathways:
            override = PATHWAY_FEATURE_OVERRIDES.get(path_id)
            if override:
                currency_iks -= override
        if not currency_iks:
            continue

        # Step 3: emit (metabolite, reaction) edges for matching substrates/products
        for u, _, _ in g.in_edges(rxn_id, data=True):
            if g.nodes[u].get("node_type") != "metabolite":
                continue
            ik = (g.nodes[u].get("inchikey") or "")[:14]
            if ik in currency_iks:
                skip.add((u, rxn_id))
        for _, v, _ in g.out_edges(rxn_id, data=True):
            if g.nodes[v].get("node_type") != "metabolite":
                continue
            ik = (g.nodes[v].get("inchikey") or "")[:14]
            if ik in currency_iks:
                skip.add((v, rxn_id))

    return skip


def flag_currency_metabolites(
    mg: GizmoGraph,
    *,
    degree_threshold_k: Optional[float] = 3.0,
    include_borderline: bool = False,
) -> dict[str, list[str]]:
    """
    Flag currency metabolites in-place on mg.

    Parameters
    ----------
    mg : GizmoGraph
        Graph to annotate.
    degree_threshold_k : float | None
        Flag metabolites whose reaction-degree > mean + k*std as statistical
        currency. Set to None to skip heuristic detection.
    include_borderline : bool
        Whether to also flag borderline currency metabolites.

    Returns
    -------
    dict with keys:
      "canonical"   — node_ids flagged by known ChEBI list
      "statistical" — node_ids flagged by degree heuristic
      "total"       — all flagged node_ids
    """
    g = mg.graph
    canonical: list[str] = []
    name_matched: list[str] = []
    statistical: list[str] = []

    currency_chebi = set(KNOWN_CURRENCY_CHEBI.keys())
    currency_ik14 = set(KNOWN_CURRENCY_INCHIKEY14.keys())
    currency_names = set(CURRENCY_NAMES_NORMALIZED)
    if not include_borderline:
        currency_chebi -= BORDERLINE_CURRENCY_CHEBI
        currency_ik14 -= BORDERLINE_CURRENCY_INCHIKEY14

    # --- Pass 1: canonical structural matching (InChIKey + ChEBI) ---
    # Reliable whenever enrich_graph_inchikey populated nodes; ChEBI is a
    # fallback for graphs where the loader set chebi_id directly.
    #
    # Pass 1b: name-based matching for compartmented Reactome metabolites.
    # Many ubiquitous cofactors (Pi, H2O, Mg2+, ATP, etc.) appear in
    # Reactome as compartment-specific instances (``Pi [nucleoplasm]``)
    # without a resolved chebi_id/inchikey. Match by normalized name.
    for nid, data in g.nodes(data=True):
        if data.get("node_type") != "metabolite":
            continue
        ik = data.get("inchikey") or ""
        chebi = data.get("chebi_id")
        if (ik and ik[:14] in currency_ik14) or (chebi and chebi in currency_chebi):
            g.nodes[nid]["is_currency"] = True
            canonical.append(nid)
            continue
        # Name-based fallback (handles Reactome's compartmented forms)
        nm = data.get("name") or ""
        norm = _normalize_name(nm)
        if norm and norm in currency_names:
            g.nodes[nid]["is_currency"] = True
            g.nodes[nid].setdefault("currency_match", "name")
            name_matched.append(nid)
            continue
        # Abstract Reactome placeholder (DNA elements, antigens, generic
        # class names, protein states, ROS placeholders, etc.)
        if _matches_abstract_pattern(nm):
            g.nodes[nid]["is_currency"] = True
            g.nodes[nid].setdefault("currency_match", "abstract")
            name_matched.append(nid)

    # --- Pass 2: degree heuristic ---
    if degree_threshold_k is not None:
        met_nodes = mg.metabolite_nodes()
        # Degree = number of distinct reaction neighbours
        degrees = []
        for nid in met_nodes:
            rxn_neighbors = {
                n for n in list(g.predecessors(nid)) + list(g.successors(nid))
                if g.nodes[n].get("node_type") == "reaction"
            }
            degrees.append((nid, len(rxn_neighbors)))

        if len(degrees) > 2:
            vals = [d for _, d in degrees]
            mean = statistics.mean(vals)
            stdev = statistics.stdev(vals) if len(vals) > 1 else 0.0
            cutoff = mean + degree_threshold_k * stdev

            for nid, deg in degrees:
                if deg > cutoff and not g.nodes[nid].get("is_currency", False):
                    g.nodes[nid]["is_currency"] = True
                    statistical.append(nid)

    all_flagged = list(set(canonical + name_matched + statistical))
    return {
        "canonical": canonical,
        "name_matched": name_matched,
        "statistical": statistical,
        "total": all_flagged,
    }


def noncurrency_subgraph(mg: GizmoGraph) -> GizmoGraph:
    """
    Return a GizmoGraph with all currency-flagged metabolites removed.
    Reactions that become substrate-less or product-less are also removed.
    Use for metabolite co-occurrence / graph kernel analysis.
    """
    g = mg.graph
    keep_mets = {
        n for n in mg.metabolite_nodes()
        if not g.nodes[n].get("is_currency", False)
    }
    keep_rxns: set[str] = set()
    for rxn in mg.reaction_nodes():
        substrates = {
            n for n in g.predecessors(rxn)
            if g.nodes[n].get("node_type") == "metabolite"
        }
        products = {
            n for n in g.successors(rxn)
            if g.nodes[n].get("node_type") == "metabolite"
        }
        # Keep reaction only if it still has at least one non-currency participant on each side
        if substrates & keep_mets and products & keep_mets:
            keep_rxns.add(rxn)

    from gizmo.graph.network import GizmoGraph
    sub = GizmoGraph()
    sub._g = g.subgraph(keep_mets | keep_rxns).copy()
    return sub
