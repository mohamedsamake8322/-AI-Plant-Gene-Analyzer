"""Topic-based candidate protein search on the bioinfo_pipeline database.

The app reads the indexed ``v_gene_text`` materialized view created by
``scripts/create_trait_search_view.sql``. It never loads the legacy species
JSON files into memory.

Advanced layer: ten topic templates, free keywords, identifier search (Pfam, InterPro, GO,
MapMan, accession), multi-species comparison, evidence levels (convergent / functional /
positional), per-candidate explanations, charts, CSV and Word exports. The original
functions (search_candidates, score_candidates, generate_docx_report, ...) are unchanged.
"""

from __future__ import annotations

import html
import json
import os
import re
import time
import argparse
import unicodedata
import tempfile
from pathlib import Path

import requests
from i18n import translate

try:
    import streamlit as st
except ImportError:
    st = None


def _normalize_text(s: str) -> str:
    """Normalizes unicode (NFC) and case before any comparison of a string
    typed by a user.

    Without this: an accent typed via a terminal, a PowerShell heredoc, or
    pasted from certain sources can be encoded in decomposed form
    (NFD -- 'e' + separate combining accent) rather than composed (NFC --
    a single code point for 'é'). Visually identical on screen,
    but `"secheresse" in q` silently fails between the two
    forms, without any error at all. This is what caused
    match_topic("sécheresse") to fail when tested on the command line."""
    return unicodedata.normalize("NFC", s or "").strip().lower()

# ── Topic template library ──────────────────────────────────────────────────
# Each template defines keyword categories (used for the diversity
# score) and a tier A (strong signal) / tier B (weak signal) split
# for the specificity score. NOISY = keywords deliberately
# excluded because they're too generic (noise empirically observed on quinoa).

def resolve_species_filter(species_input: str | None) -> str | None:
    """Normalize a registry key supplied to the CLI.

    Display names and aliases are intentionally not maintained in the app;
    they belong in ``species_registry`` and are resolved by the UI query.
    """
    if not species_input:
        return None
    key = _normalize_text(species_input)
    return key.replace(" ", "_")



TOPIC_TEMPLATES: dict[str, dict] = {
    "verse": {
        "aliases": ["verse", "lodging", "tige", "rigidite"],
        "label": "Lodging / stem strength resistance",
        "pubmed_context": "lodging",
        "qtl_terms": ["lodging", "stem strength", "culm strength", "stalk strength", "breaking strength"],
        "keywords": {
            "lignification": [
                "lignin", "lignification", "cell wall", "secondary cell wall",
                "cellulose synthase", "laccase", "peroxidase",
                "cinnamyl alcohol", "phenylpropanoid", "4cl", "ccoaomt", "comt",
            ],
            "stem_rigidity": [
                "stem", "culm", "stalk", "lodging", "verse", "internode",
                "mechanical strength", "stem strength",
            ],
            "hormonal": [
                "gibberellin", "gibberellic", "della", "ga20ox", "ga3ox",
                "ga2ox", "brassinosteroid", "bri1",
            ],
        },
        "tier_a": {
            "lignin", "lignification", "laccase", "cellulose synthase",
            "secondary cell wall", "phenylpropanoid", "cinnamyl alcohol",
            "4cl", "comt", "ccoaomt", "della", "gibberellin", "gibberellic",
            "ga20ox", "ga3ox", "ga2ox", "brassinosteroid", "bri1",
            "verse", "lodging", "stalk", "mechanical strength",
            "stem strength",
        },
        "tier_b": {"peroxidase", "culm", "internode"},
        "noisy": {"stem", "cell wall"},
    },
    "secheresse": {
        "aliases": ["secheresse", "sécheresse", "drought", "hydrique", "eau"],
        "label": "Drought tolerance / water stress",
        "pubmed_context": "drought",
        "qtl_terms": ["drought", "water stress", "water deficit", "water use efficiency"],
        "keywords": {
            "aba_signaling": [
                "abscisic acid", "aba receptor", "pyr/pyl", "snrk2",
                "aba signaling",
            ],
            "osmotic_response": [
                "osmotic stress", "pyrroline-5-carboxylate synthase", "pyrroline-5-carboxylate synthetase", "p5cs", "proline accumulation", "proline biosynthesis", "proline dehydrogenase", "proline transporter", "dehydrin", "late embryogenesis",
                "lea protein", "osmotic adjustment",
            ],
            "stomatal_regulation": [
                "stomatal closure", "guard cell", "aquaporin", "water use efficiency",
            ],
            "transcription_factors": [
                "dreb", "nac transcription factor", "wrky", "myb drought",
            ],
            "direct_mention": ["drought", "water deficit", "dehydration"],
        },
        "tier_a": {
            "abscisic acid", "aba receptor", "pyr/pyl", "snrk2",
            "dehydrin", "lea protein", "dreb", "aquaporin",
            "stomatal closure", "drought",
        },
        "tier_b": {
            "aba signaling", "osmotic stress", "pyrroline-5-carboxylate synthase", "pyrroline-5-carboxylate synthetase", "p5cs", "proline accumulation", "proline biosynthesis", "proline dehydrogenase", "proline transporter",
            "late embryogenesis", "osmotic adjustment", "guard cell",
            "water use efficiency", "nac transcription factor", "wrky",
            "myb drought", "water deficit", "dehydration",
        },
        "noisy": set(),
    },
    "salinite": {
        "aliases": ["salinite", "salinité", "salinity", "salt", "sel", "salt stress", "tuz"],
        "label": "Salinity tolerance / salt stress",
        "pubmed_context": "salt stress",
        "qtl_terms": ["salt", "salinity", "sodium"],
        "keywords": {
            "ion_homeostasis": [
                "sodium", "na+/h+", "nhx", "hkt", "sos1", "salt overly sensitive",
                "potassium transporter", "cation/h+ exchanger",
            ],
            "osmoprotection": ["pyrroline-5-carboxylate synthase", "pyrroline-5-carboxylate synthetase", "p5cs", "proline accumulation", "proline biosynthesis", "proline dehydrogenase", "proline transporter", "betaine", "trehalose", "osmotic stress"],
            "ros_detox": ["superoxide dismutase", "catalase", "ascorbate peroxidase"],
            "direct_mention": ["salt stress", "salinity", "salt tolerance"],
        },
        "tier_a": {"nhx", "hkt", "sos1", "salt overly sensitive", "salt stress",
                   "salinity", "salt tolerance", "na+/h+"},
        "tier_b": {"sodium", "potassium transporter", "cation/h+ exchanger", "pyrroline-5-carboxylate synthase", "pyrroline-5-carboxylate synthetase", "p5cs", "proline accumulation", "proline biosynthesis", "proline dehydrogenase", "proline transporter",
                   "betaine", "trehalose", "superoxide dismutase", "catalase",
                   "ascorbate peroxidase", "osmotic stress"},
        "noisy": set(),
    },
    "froid": {
        "aliases": ["froid", "cold", "chilling", "freezing", "gel", "cold stress", "soguk"],
        "label": "Cold / chilling / freezing tolerance",
        "pubmed_context": "cold stress",
        "qtl_terms": ["cold", "chilling", "freezing"],
        "keywords": {
            "cbf_pathway": ["cbf", "dreb1", "ice1", "inducer of cbf", "cold-regulated", "cor15"],
            "membrane_protection": ["dehydrin", "fatty acid desaturase", "antifreeze"],
            "direct_mention": ["cold stress", "chilling", "freezing tolerance", "cold acclimation"],
        },
        "tier_a": {"cbf", "dreb1", "ice1", "inducer of cbf", "cold-regulated", "cor15",
                   "cold stress", "chilling", "freezing tolerance", "cold acclimation"},
        "tier_b": {"dehydrin", "fatty acid desaturase", "antifreeze"},
        "noisy": set(),
    },
    "chaleur": {
        "aliases": ["chaleur", "heat", "heat stress", "thermotolerance", "sicaklik"],
        "label": "Heat stress / thermotolerance",
        "pubmed_context": "heat stress",
        "qtl_terms": ["heat", "thermotolerance", "high temperature"],
        "keywords": {
            "chaperones": ["heat shock", "hsp70", "hsp90", "hsp101", "small heat shock", "chaperone"],
            "regulators": ["heat stress transcription factor", "hsf", "dreb2"],
            "direct_mention": ["heat stress", "thermotolerance", "high temperature"],
        },
        "tier_a": {"heat shock", "hsp70", "hsp90", "hsp101", "small heat shock",
                   "heat stress transcription factor", "hsf", "thermotolerance", "heat stress"},
        "tier_b": {"chaperone", "dreb2", "high temperature"},
        "noisy": set(),
    },
    "resistance_maladies": {
        "aliases": ["resistance", "résistance", "maladie", "disease", "disease resistance",
                    "pathogen", "blast", "mildiou", "hastalik"],
        "label": "Disease / pathogen resistance",
        "pubmed_context": "disease resistance",
        "qtl_terms": ["disease", "blast", "blight", "pathogen", "resistance to"],
        "keywords": {
            "r_genes": ["nb-arc", "nbs-lrr", "leucine-rich repeat", "resistance protein",
                        "disease resistance protein", "tir domain"],
            "defense": ["pathogenesis-related", "chitinase", "defensin", "thaumatin",
                        "plant defense", "hypersensitive response"],
            "signalling": ["receptor-like kinase", "salicylic acid", "jasmonate", "wrky"],
            "direct_mention": ["disease resistance", "pathogen", "blast"],
        },
        "tier_a": {"nb-arc", "nbs-lrr", "resistance protein", "disease resistance protein",
                   "pathogenesis-related", "hypersensitive response", "disease resistance"},
        "tier_b": {"leucine-rich repeat", "tir domain", "chitinase", "defensin", "thaumatin",
                   "plant defense", "receptor-like kinase", "salicylic acid", "jasmonate",
                   "wrky", "pathogen", "blast"},
        "noisy": set(),
    },
    "azote": {
        "aliases": ["azote", "nitrogen", "nitrate", "nue", "nitrogen use efficiency"],
        "label": "Nitrogen uptake and use efficiency",
        "pubmed_context": "nitrogen use efficiency",
        "qtl_terms": ["nitrogen", "nitrate"],
        "keywords": {
            "uptake": ["nitrate transporter", "ammonium transporter", "nrt1", "nrt2", "amt1"],
            "assimilation": ["nitrate reductase", "nitrite reductase", "glutamine synthetase",
                             "glutamate synthase", "asparagine synthetase"],
            "regulation": ["nlp", "nin-like", "nitrate regulat"],
            "direct_mention": ["nitrogen use efficiency", "nitrogen"],
        },
        "tier_a": {"nitrate transporter", "ammonium transporter", "nrt1", "nrt2", "amt1",
                   "nitrate reductase", "nitrite reductase", "glutamine synthetase",
                   "glutamate synthase", "nitrogen use efficiency"},
        "tier_b": {"asparagine synthetase", "nlp", "nin-like", "nitrate regulat", "nitrogen"},
        "noisy": set(),
    },
    "floraison": {
        "aliases": ["floraison", "flowering", "heading date", "flowering time", "cicek"],
        "label": "Flowering time / heading date",
        "pubmed_context": "flowering time",
        "qtl_terms": ["flowering", "heading date", "days to heading", "days to anthesis"],
        "keywords": {
            "photoperiod": ["constans", "gigantea", "phytochrome", "cryptochrome", "photoperiod"],
            "florigen": ["flowering locus t", "hd3a", "heading date", "florigen", "rft1"],
            "vernalization": ["vernalization", "flowering locus c", "vin3"],
            "direct_mention": ["flowering time", "flowering"],
        },
        "tier_a": {"constans", "gigantea", "flowering locus t", "hd3a", "heading date",
                   "florigen", "vernalization", "flowering locus c", "flowering time"},
        "tier_b": {"phytochrome", "cryptochrome", "photoperiod", "rft1", "vin3", "flowering"},
        "noisy": set(),
    },
    "rendement_grain": {
        "aliases": ["rendement", "yield", "grain", "grain size", "grain yield", "verim"],
        "label": "Grain size / yield components",
        "pubmed_context": "grain size",
        "qtl_terms": ["grain size", "grain weight", "grain yield", "grain number", "yield"],
        "keywords": {
            "grain_size": ["grain size", "grain weight", "gs3", "gw2", "gw5", "gif1", "tgw6"],
            "cytokinin": ["cytokinin oxidase", "ckx", "cytokinin"],
            "starch_filling": ["starch synthase", "sucrose transporter", "sucrose synthase",
                               "adp-glucose pyrophosphorylase"],
            "direct_mention": ["grain yield", "grain number"],
        },
        "tier_a": {"grain size", "grain weight", "gs3", "gw2", "gw5", "gif1", "tgw6",
                   "cytokinin oxidase", "ckx", "grain yield", "grain number"},
        "tier_b": {"cytokinin", "starch synthase", "sucrose transporter", "sucrose synthase",
                   "adp-glucose pyrophosphorylase"},
        "noisy": set(),
    },
    "nicotine_alcaloides": {
        "aliases": ["nicotine", "alcaloide", "alcaloides", "alkaloid", "nicotine biosynthesis",
                    "tobacco alkaloid", "nikotin"],
        "label": "Nicotine / alkaloid metabolism (tobacco)",
        "pubmed_context": "nicotine",
        "qtl_terms": ["nicotine", "alkaloid"],
        "keywords": {
            "pyrrolidine_ring": ["putrescine n-methyltransferase", "pmt", "ornithine decarboxylase",
                                 "n-methylputrescine oxidase", "mpo", "arginine decarboxylase"],
            "pyridine_ring": ["quinolinate phosphoribosyltransferase", "qpt", "aspartate oxidase",
                              "quinolinate synthase", "nicotinate"],
            "coupling_demethylation": ["berberine bridge", "a622", "nicotine demethylase",
                                       "cyp82e", "isoflavone reductase"],
            "regulation": ["erf189", "nic2", "jasmonate", "myc2"],
            "direct_mention": ["nicotine", "alkaloid", "tropane"],
        },
        "tier_a": {"putrescine n-methyltransferase", "pmt", "n-methylputrescine oxidase",
                   "quinolinate phosphoribosyltransferase", "qpt", "berberine bridge", "a622",
                   "nicotine demethylase", "cyp82e", "erf189", "nic2", "nicotine"},
        "tier_b": {"ornithine decarboxylase", "mpo", "arginine decarboxylase", "aspartate oxidase",
                   "quinolinate synthase", "nicotinate", "isoflavone reductase", "jasmonate",
                   "alkaloid", "tropane"},
        "noisy": set(),
    },
}



def match_topic(user_query: str) -> str | None:
    """Finds the topic template closest to the text typed by
    the user. Returns the template's key, or None if there is no match --
    in that case the app should indicate that the topic is not yet covered
    (future extension: keyword generation via LLM)."""
    q = _normalize_text(user_query)
    for key, tpl in TOPIC_TEMPLATES.items():
        if any(_normalize_text(alias) in q for alias in tpl["aliases"]):
            return key
    return None


# ── PostgreSQL-backed candidate search ─────────────────────────────────────

def _connect_bioinfo_pipeline():
    """Open the dedicated ESM-2/annotation database, never the main genes DB."""
    from core_engines.bioinfo_db import connect

    return connect()


def _view_setup_error(error: Exception) -> RuntimeError:
    return RuntimeError(
        "Trait search database is not initialized. Run "
        "scripts/create_trait_search_view.sql against bioinfo_pipeline, "
        "then refresh v_gene_text after data imports. "
        f"Database reported: {error}"
    )


def _decode_qtl_evidence(value: object) -> list[dict]:
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError as error:
            raise ValueError("Malformed QTL JSON returned by v_gene_text.") from error
    if value is None:
        return []
    if not isinstance(value, list) or any(not isinstance(item, dict) for item in value):
        raise ValueError("Expected qtl_evidence to be a JSON array of records.")
    return value


def _list_available_species() -> list[dict]:
    """Return only assemblies registered for the ESM-2 reference collection."""
    connection = _connect_bioinfo_pipeline()
    try:
        with connection.cursor() as cursor:
            cursor.execute(
                """
                SELECT species, display_name
                FROM species_registry
                WHERE enabled IS TRUE
                  AND EXISTS (
                      SELECT 1 FROM v_gene_text AS v WHERE v.species = species_registry.species
                  )
                ORDER BY display_name;
                """
            )
            return [
                {"species": species, "display_name": label}
                for species, label in cursor.fetchall()
            ]
    except Exception as error:
        raise _view_setup_error(error) from error
    finally:
        connection.close()


if st is not None:
    list_available_species = st.cache_data(ttl=300, show_spinner=False)(
        _list_available_species
    )
else:
    list_available_species = _list_available_species


def search_candidates(
    template: dict,
    species_filter: str,
    limit: int = 1000,
) -> list[dict]:
    """Search indexed annotations in SQL and score the bounded result set."""
    if not species_filter:
        raise ValueError("A registered species is required for trait search.")

    all_keywords = [
        (keyword, category)
        for category, keywords in template["keywords"].items()
        for keyword in keywords
        if keyword not in template["noisy"]
    ]
    if not all_keywords:
        return []

    fts_query = " OR ".join(
        f'"{keyword.replace(chr(34), "")}"' for keyword, _category in all_keywords
    )
    pool_size = max(100, min(limit, 5000))
    connection = _connect_bioinfo_pipeline()
    try:
        with connection.cursor() as cursor:
            cursor.execute(
                """
                SELECT species, assembly_accession, protein_accession, plaza_gene_id,
                       gene_symbol, description, go_text, pfam_text, interpro_text,
                       trait_text, qtl_evidence, search_text,
                       ts_rank_cd(search_vector, websearch_to_tsquery('simple', translate(%s, './', '  '))) AS text_rank,
                       COUNT(*) OVER () AS total_matches
                FROM v_gene_text
                WHERE species = %s
                  AND (
                      search_vector @@ websearch_to_tsquery('simple', translate(%s, './', '  '))
                      OR search_text ILIKE ANY(%s)
                  )
                ORDER BY text_rank DESC, protein_accession
                LIMIT %s;
                """,
                (
                    fts_query,
                    species_filter,
                    fts_query,
                    [f"%{keyword}%" for keyword, _category in all_keywords],
                    pool_size,
                ),
            )
            column_names = [column[0] for column in cursor.description]
            rows = [dict(zip(column_names, row)) for row in cursor.fetchall()]
    except Exception as error:
        raise _view_setup_error(error) from error
    finally:
        connection.close()

    results = []
    for row in rows:
        text = str(row.get("search_text") or "").lower()
        matches = {
            f"{category}:{keyword.strip()}"
            for keyword, category in all_keywords
            if keyword.lower() in text
        }
        if not matches:
            continue

        protein_accession = row["protein_accession"]
        gene = {
            "common_name": row.get("gene_symbol") or row.get("plaza_gene_id") or protein_accession,
            "gene_symbol": row.get("gene_symbol"),
            "description": row.get("description") or "",
            "organism": row["species"],
            "assembly_accession": row["assembly_accession"],
            "protein_accession": protein_accession,
            "annotation": {
                "go_text": row.get("go_text") or "",
                "pfam_text": row.get("pfam_text") or "",
                "interpro_text": row.get("interpro_text") or "",
            },
            "trait_text": row.get("trait_text") or "",
            "qtl_evidence": _decode_qtl_evidence(row.get("qtl_evidence")),
            "text_rank": float(row.get("text_rank") or 0.0),
            "total_matches": int(row.get("total_matches") or 0),
        }
        results.append({
            "gene_id": protein_accession,
            "candidate_key": (
                row["assembly_accession"],
                protein_accession,
            ),
            "gene": gene,
            "matches": matches,
            "text_rank": gene["text_rank"],
        })
    return results


def score_candidates(candidates: list[dict], template: dict) -> list[dict]:
    tier_a, tier_b = template["tier_a"], template["tier_b"]
    for c in candidates:
        kws = {m.split(":", 1)[1] for m in c["matches"]}
        n_tier_a = len(kws & tier_a)
        n_tier_b = len(kws & tier_b)
        n_categories = len({m.split(":")[0] for m in c["matches"]})
        score = n_tier_a * 3 + n_tier_b + max(0, n_categories - 1) * 2
        c["score"] = score
        c["n_categories"] = n_categories
        c["categories"] = sorted({m.split(":")[0] for m in c["matches"]})
    candidates.sort(
        key=lambda c: (c["score"], c.get("text_rank", 0.0)),
        reverse=True,
    )
    return candidates

# ═══════════════════════════════════════════════════════════════════════════
# ADVANCED SEARCH LAYER
# Evidence model, multi-species search, free keywords, identifiers, explanations,
# data frames for charts, CSV / Word exports. The functions above are unchanged.
# ═══════════════════════════════════════════════════════════════════════════

# Evidence fields, in display order. "qtl" is positional evidence (the gene lies inside a QTL
# interval); every other field is functional evidence attached to the protein itself.
EVIDENCE_LABELS = {
    "product": "Product / gene name",
    "description": "Curated description",
    "mapman": "MapMan bin",
    "interpro": "InterPro",
    "pfam": "Pfam domain",
    "go": "GO",
    "qtl": "QTL (position)",
}
FUNCTIONAL_FIELDS = ("product", "description", "mapman", "interpro", "pfam", "go")
# Curated, human-readable sources weigh more than identifier-only sources.
FIELD_WEIGHT = {"product": 2, "description": 2, "mapman": 2, "interpro": 1, "pfam": 1, "go": 1, "qtl": 0}

TIER_CONVERGENT = "convergent"
TIER_FUNCTIONAL = "functional"
TIER_POSITIONAL = "positional"
TIER_LABELS = {
    TIER_CONVERGENT: "★ Convergent (function + QTL)",
    TIER_FUNCTIONAL: "Functional evidence",
    TIER_POSITIONAL: "Positional only (QTL)",
}
TIER_RANK = {TIER_CONVERGENT: 3, TIER_FUNCTIONAL: 2, TIER_POSITIONAL: 1}

_KW_CACHE: dict[str, "re.Pattern[str]"] = {}


def _kw_regex(keyword: str) -> "re.Pattern[str]":
    """Prefix word-boundary match: 'nhx' matches 'nhx1' but not 'xnhx'; 'peroxidase' matches 'peroxidases'."""
    pattern = _KW_CACHE.get(keyword)
    if pattern is None:
        pattern = re.compile(r"(?<![a-z0-9])" + re.escape(keyword.lower()))
        _KW_CACHE[keyword] = pattern
    return pattern


def _t(key: str, default: str, **kwargs) -> str:
    """UI text through the app's i18n layer, falling back to the English default."""
    try:
        text = translate(key, default=default)
    except Exception:
        text = default
    if kwargs:
        try:
            text = text.format(**kwargs)
        except (KeyError, IndexError, ValueError):
            pass
    return text


def template_keywords(template: dict) -> list[tuple[str, str]]:
    noisy = template.get("noisy") or set()
    return [
        (keyword, category)
        for category, keywords in template["keywords"].items()
        for keyword in keywords
        if keyword not in noisy
    ]


def template_qtl_terms(template: dict) -> list[str]:
    terms = template.get("qtl_terms")
    if terms:
        return list(terms)
    return list(dict.fromkeys(list(template.get("aliases", [])) + [template.get("pubmed_context", "")]))


def build_custom_template(raw_terms: str, match_all: bool = False) -> dict | None:
    """Template from user-typed keywords (separated by comma, semicolon or new line)."""
    terms = []
    for part in re.split(r"[,;\n]+", raw_terms or ""):
        term = _normalize_text(part)
        if len(term) >= 2 and term not in terms:
            terms.append(term)
    if not terms:
        return None
    return {
        "aliases": [],
        "label": "Custom keywords: " + (" AND " if match_all else " OR ").join(terms),
        "pubmed_context": "",
        "qtl_terms": terms,
        "keywords": {"user_terms": terms},
        "tier_a": set(terms),
        "tier_b": set(),
        "noisy": set(),
        "match_all": match_all,
    }


_IDENTIFIER_PATTERNS = (
    ("Pfam", re.compile(r"^PF\d{5}(?:\.\d+)?$", re.IGNORECASE)),
    ("InterPro", re.compile(r"^IPR\d{6}$", re.IGNORECASE)),
    ("GO", re.compile(r"^GO:\d{7}$", re.IGNORECASE)),
    ("MapMan bin", re.compile(r"^\d{1,2}(?:\.\d+)+$")),
    ("Protein accession", re.compile(r"^(?:XP|NP|YP|WP)_\d+(?:\.\d+)?$", re.IGNORECASE)),
)


def classify_identifier(text: str) -> tuple[str, str]:
    """Return (kind, normalised identifier). Kind is 'Gene name / text' when nothing matches."""
    value = (text or "").strip()
    for kind, pattern in _IDENTIFIER_PATTERNS:
        if pattern.match(value):
            if kind == "Pfam":
                value = re.sub(r"\.\d+$", "", value.upper())
            elif kind in ("InterPro", "GO"):
                value = value.upper()
            return kind, value
    return "Gene name / text", value


def _decode_records(value: object) -> list[dict]:
    try:
        return _decode_qtl_evidence(value)
    except ValueError:
        return []


def _functional_text(row: dict) -> str:
    """search_text without its trailing QTL part (trait_text is always concatenated last in the view)."""
    text = str(row.get("search_text") or "")
    trait = str(row.get("trait_text") or "")
    if trait and text.endswith(trait):
        text = text[: len(text) - len(trait)].rstrip()
    return text


def _field_texts(row: dict) -> dict[str, str]:
    texts = {
        "description": str(row.get("description") or ""),
        "mapman": str(row.get("mapman_text") or ""),
        "interpro": str(row.get("interpro_text") or ""),
        "pfam": str(row.get("pfam_text") or ""),
        "go": str(row.get("go_text") or ""),
    }
    product = _functional_text(row)
    for text in texts.values():
        if text:
            product = product.replace(text, " ")
    texts["product"] = product
    return texts


def _qtl_text(record: dict) -> str:
    return " ".join(
        str(record.get(key) or "") for key in ("major_category", "trait_category", "trait_detail")
    ).lower()


def analyse_row(row: dict, template: dict, include_positional: bool) -> dict | None:
    """Turn one v_gene_text row into a candidate with field-level evidence, or None if nothing matches."""
    # Source descriptions contain HTML such as Na<sup>+</sup>; strip it so it is neither matched nor shown raw.
    row = dict(row)
    for _key in ("description", "go_text", "pfam_text", "interpro_text", "mapman_text", "trait_text"):
        if row.get(_key):
            row[_key] = html.unescape(re.sub(r"</?(?:sup|sub|i|b|em|strong|span)[^>]*>", "", str(row[_key])))
    keywords = template_keywords(template)
    texts = {field: text.lower() for field, text in _field_texts(row).items()}
    field_hits: dict[str, set[str]] = {}
    categories_by_keyword: dict[str, str] = {}
    for keyword, category in keywords:
        categories_by_keyword[keyword] = category
        regex = _kw_regex(keyword)
        for field in FUNCTIONAL_FIELDS:
            if regex.search(texts[field]):
                field_hits.setdefault(field, set()).add(keyword)

    qtl_records = _decode_records(row.get("qtl_evidence"))
    qtl_terms = [term.lower() for term in template_qtl_terms(template) if term]
    topic_qtls = [
        record for record in qtl_records
        if any(term in _qtl_text(record) for term in qtl_terms)
    ]
    trait_lower = str(row.get("trait_text") or "").lower()
    qtl_keywords = {keyword for keyword, _c in keywords if _kw_regex(keyword).search(trait_lower)}
    if qtl_keywords:
        field_hits["qtl"] = qtl_keywords

    functional_fields = [f for f in FUNCTIONAL_FIELDS if f in field_hits]
    has_qtl_signal = bool(topic_qtls) or bool(qtl_keywords)
    if functional_fields:
        tier = TIER_CONVERGENT if topic_qtls else TIER_FUNCTIONAL
    elif has_qtl_signal and include_positional:
        tier = TIER_POSITIONAL
    else:
        return None

    all_keywords_hit = set().union(*field_hits.values()) if field_hits else set()
    if template.get("match_all") and {k for k, _c in keywords} - all_keywords_hit:
        return None
    functional_keywords = set().union(*(field_hits[f] for f in functional_fields)) if functional_fields else set()
    matches = {f"{categories_by_keyword[k]}:{k}" for k in all_keywords_hit}
    categories = sorted({categories_by_keyword[k] for k in all_keywords_hit})

    gene_symbol = row.get("gene_symbol")
    protein_accession = row["protein_accession"]
    candidate = {
        "gene_id": protein_accession,
        "candidate_key": (row["assembly_accession"], protein_accession),
        "gene": {
            "common_name": gene_symbol or row.get("plaza_gene_id") or protein_accession,
            "gene_symbol": gene_symbol,
            "description": row.get("description") or "",
            "organism": row["species"],
            "assembly_accession": row["assembly_accession"],
            "protein_accession": protein_accession,
            "plaza_gene_id": row.get("plaza_gene_id"),
            "annotation": {
                "go_text": row.get("go_text") or "",
                "pfam_text": row.get("pfam_text") or "",
                "interpro_text": row.get("interpro_text") or "",
                "mapman_text": row.get("mapman_text") or "",
            },
            "trait_text": row.get("trait_text") or "",
            "qtl_evidence": qtl_records,
            "text_rank": float(row.get("text_rank") or 0.0),
            "total_matches": int(row.get("total_matches") or 0),
        },
        "matches": matches,
        "text_rank": float(row.get("text_rank") or 0.0),
        "tier": tier,
        "convergent": tier == TIER_CONVERGENT,
        "field_hits": {f: sorted(field_hits[f]) for f in EVIDENCE_LABELS if f in field_hits},
        "functional_fields": functional_fields,
        "functional_keywords": sorted(functional_keywords),
        "topic_qtls": topic_qtls,
        "categories": categories,
        "n_categories": len(categories),
        "references": [],
    }
    return candidate


def score_candidates_advanced(candidates: list[dict], template: dict) -> list[dict]:
    """Score and sort. Priority: convergent > functional > positional; then score; then full-text rank.

    Functional score = 3 x strong keywords + 1 x weak keywords + 2 x (extra categories)
                       + curated-source weights (product/description/MapMan 2, InterPro/Pfam/GO 1, each once)
                       + 4 if the protein also lies in a QTL of the same topic.
    Positional-only score = number of matching QTLs, capped at 3 (position is weak evidence).
    """
    tier_a, tier_b = template["tier_a"], template["tier_b"]
    for candidate in candidates:
        keywords = set(candidate["functional_keywords"])
        if candidate["tier"] == TIER_POSITIONAL:
            candidate["score"] = min(len(candidate["topic_qtls"]) or 1, 3)
        else:
            categories = {m.split(":", 1)[0] for m in candidate["matches"]
                          if m.split(":", 1)[1] in keywords}
            base = len(keywords & tier_a) * 3 + len(keywords & tier_b)
            base += max(0, len(categories) - 1) * 2
            base += sum(FIELD_WEIGHT[f] for f in candidate["functional_fields"])
            candidate["score"] = base + (4 if candidate["convergent"] else 0)
    candidates.sort(
        key=lambda c: (TIER_RANK[c["tier"]], c["score"], c.get("text_rank", 0.0)),
        reverse=True,
    )
    return candidates


def _prefix_tsquery(keywords: list[str], match_all: bool = False) -> str:
    """to_tsquery text with word-start prefix matching: 'dreb' also finds DREB1A, 'abscisic acid' is a phrase
    whose last word is a prefix. Mirrors _kw_regex (prefix word-boundary match)."""
    parts = []
    for keyword in keywords:
        words = re.findall(r"[a-z0-9]+", keyword.lower().replace("/", " ").replace(".", " "))
        if not words:
            continue
        words[-1] += ":*"
        parts.append("(" + " <-> ".join(words) + ")")
    return (" & " if match_all else " | ").join(parts)


def search_candidates_advanced(
    template: dict,
    species_filter: str,
    limit: int = 1000,
    include_positional: bool = False,
) -> dict:
    """Search one species. Functional hits are always ranked before QTL-only hits in the SQL pool,
    so a large QTL interval can never push functional candidates out of the pool.

    Returns {"candidates", "total_matches", "functional_matches", "pool_size"}."""
    if not species_filter:
        raise ValueError("A registered species is required for trait search.")
    keywords = template_keywords(template)
    if not keywords:
        return {"candidates": [], "total_matches": 0, "functional_matches": 0, "pool_size": 0}

    match_all = bool(template.get("match_all"))
    joiner = " " if match_all else " OR "
    quote = lambda term: '"' + term.replace('"', "") + '"'
    fq = joiner.join(quote(k) for k, _c in keywords)
    qq = " OR ".join(quote(t) for t in template_qtl_terms(template) if t)
    likes = [f"%{k}%" for k, _c in keywords]
    like_op = "ALL" if match_all else "ANY"
    pool_size = max(100, min(limit, 5000))

    legacy_sql = f"""
        WITH q AS (
            SELECT websearch_to_tsquery('simple', translate(%(fq)s, './', '  ')) AS fq,
                   websearch_to_tsquery('simple', translate(%(qq)s, './', '  ')) AS qq
        ),
        hits AS (
            SELECT v.species, v.assembly_accession, v.protein_accession, v.plaza_gene_id,
                   v.gene_symbol, v.description, v.go_text, v.pfam_text, v.interpro_text,
                   v.mapman_text, v.trait_text, v.qtl_evidence, v.search_text,
                   ts_rank_cd(v.search_vector, q.fq) AS text_rank,
                   (
                       to_tsvector('simple', translate(
                           CASE WHEN v.trait_text IS NULL THEN v.search_text
                                ELSE left(v.search_text,
                                          greatest(length(v.search_text) - length(v.trait_text) - 1, 0))
                           END, './', '  ')) @@ q.fq
                       OR (CASE WHEN v.trait_text IS NULL THEN v.search_text
                                ELSE left(v.search_text,
                                          greatest(length(v.search_text) - length(v.trait_text) - 1, 0))
                           END) ILIKE {like_op}(%(likes)s)
                   ) AS functional_hit
            FROM v_gene_text AS v, q
            WHERE v.species = %(species)s
              AND (
                  v.search_vector @@ q.fq
                  OR (%(use_qtl)s AND v.search_vector @@ q.qq)
                  OR v.search_text ILIKE {like_op}(%(likes)s)
              )
        )
        SELECT hits.*,
               COUNT(*) OVER () AS total_matches,
               COUNT(*) FILTER (WHERE functional_hit) OVER () AS functional_matches
        FROM hits
        ORDER BY functional_hit DESC, text_rank DESC, protein_accession
        LIMIT %(limit)s
    """
    legacy_params = {
        "fq": fq, "qq": qq or "", "likes": likes, "species": species_filter,
        "use_qtl": bool(include_positional and qq), "limit": pool_size,
    }

    # Fast path (needs the small materialised view v_gene_func, see scripts/create_gene_func_view.py):
    # prefix matching on word starts through the GIN index, same semantics as the keyword check in analyse_row.
    fast_sql = """
        WITH q AS (
            SELECT to_tsquery('simple', %(fq_ts)s) AS fq,
                   websearch_to_tsquery('simple', translate(%(qq)s, './', '  ')) AS qq
        ),
        func AS (
            SELECT f.species, f.assembly_accession, f.protein_accession,
                   ts_rank_cd(f.functional_vector, q.fq) AS text_rank
            FROM v_gene_func AS f, q
            WHERE f.species = %(species)s AND f.functional_vector @@ q.fq
        ),
        pos AS (
            SELECT v.species, v.assembly_accession, v.protein_accession
            FROM v_gene_text AS v, q
            WHERE v.species = %(species)s AND v.search_vector @@ (q.qq || q.fq)
              AND NOT EXISTS (SELECT 1 FROM func AS f2
                              WHERE f2.assembly_accession = v.assembly_accession
                                AND f2.protein_accession = v.protein_accession)
        ),
        counts AS (SELECT (SELECT count(*) FROM func) AS nf, (SELECT count(*) FROM pos) AS np),
        hits AS (
            SELECT species, assembly_accession, protein_accession, text_rank, TRUE AS functional_hit FROM func
            UNION ALL
            SELECT species, assembly_accession, protein_accession, 0::real, FALSE FROM pos WHERE %(use_qtl)s
        ),
        top AS (
            SELECT * FROM hits ORDER BY functional_hit DESC, text_rank DESC, protein_accession LIMIT %(limit)s
        )
        SELECT v.species, v.assembly_accession, v.protein_accession, v.plaza_gene_id,
               v.gene_symbol, v.description, v.go_text, v.pfam_text, v.interpro_text,
               v.mapman_text, v.trait_text, v.qtl_evidence, v.search_text,
               top.text_rank, top.functional_hit,
               counts.nf + counts.np AS total_matches, counts.nf AS functional_matches
        FROM top
        JOIN v_gene_text AS v ON v.species = top.species
                             AND v.assembly_accession = top.assembly_accession
                             AND v.protein_accession = top.protein_accession
        CROSS JOIN counts
        ORDER BY top.functional_hit DESC, top.text_rank DESC, top.protein_accession
    """
    fq_ts = _prefix_tsquery([k for k, _c in keywords], match_all)
    fast_params = {"fq_ts": fq_ts, "qq": qq or "", "species": species_filter,
                   "use_qtl": bool(include_positional and qq), "limit": pool_size}

    connection = _connect_bioinfo_pipeline()
    try:
        rows = None
        if fq_ts:
            try:
                with connection.cursor() as cursor:
                    cursor.execute(fast_sql, fast_params)
                    column_names = [column[0] for column in cursor.description]
                    rows = [dict(zip(column_names, row)) for row in cursor.fetchall()]
            except Exception as error:  # noqa: BLE001
                if "v_gene_func" not in str(error):
                    raise
                connection.rollback()  # fast view missing: use the slower legacy query
                rows = None
        if rows is None:
            with connection.cursor() as cursor:
                cursor.execute(legacy_sql, legacy_params)
                column_names = [column[0] for column in cursor.description]
                rows = [dict(zip(column_names, row)) for row in cursor.fetchall()]
    except Exception as error:
        raise _view_setup_error(error) from error
    finally:
        connection.close()

    candidates = []
    for row in rows:
        candidate = analyse_row(row, template, include_positional)
        if candidate is not None:
            candidates.append(candidate)
    total = max((int(r.get("total_matches") or 0) for r in rows), default=0)
    functional_total = max((int(r.get("functional_matches") or 0) for r in rows), default=0)
    return {
        "candidates": score_candidates_advanced(candidates, template),
        "total_matches": total,
        "functional_matches": functional_total,
        "pool_size": len(rows),
    }


def search_by_identifier(species_filter: str, query: str, limit: int = 1000) -> dict:
    """Exact-token search for a Pfam / InterPro / GO / MapMan identifier, a protein accession or a gene name."""
    kind, value = classify_identifier(query)
    if not value:
        return {"candidates": [], "total_matches": 0, "functional_matches": 0, "pool_size": 0, "kind": kind}
    if kind == "Protein accession":
        where, param = "v.protein_accession ILIKE %(p)s", value.rsplit(".", 1)[0] + "%"
    else:
        where, param = "v.search_text ~* %(p)s", r"(^|[^[:alnum:]])" + re.escape(value)
    sql = f"""
        SELECT v.species, v.assembly_accession, v.protein_accession, v.plaza_gene_id, v.gene_symbol,
               v.description, v.go_text, v.pfam_text, v.interpro_text, v.mapman_text, v.trait_text,
               v.qtl_evidence, v.search_text, 0.0 AS text_rank,
               COUNT(*) OVER () AS total_matches
        FROM v_gene_text AS v
        WHERE v.species = %(species)s AND {where}
        ORDER BY v.protein_accession
        LIMIT %(limit)s
    """
    connection = _connect_bioinfo_pipeline()
    try:
        with connection.cursor() as cursor:
            cursor.execute(sql, {"p": param, "species": species_filter, "limit": max(100, min(limit, 5000))})
            names = [c[0] for c in cursor.description]
            rows = [dict(zip(names, r)) for r in cursor.fetchall()]
    except Exception as error:
        raise _view_setup_error(error) from error
    finally:
        connection.close()

    template = {"keywords": {"identifier": [value.lower()]}, "tier_a": {value.lower()}, "tier_b": set(),
                "noisy": set(), "qtl_terms": []}
    candidates = []
    needle = value.lower()
    for row in rows:
        texts = {f: t.lower() for f, t in _field_texts(row).items()}
        fields = [f for f in FUNCTIONAL_FIELDS if _kw_regex(needle).search(texts[f])]
        if kind == "Protein accession":
            fields = fields or ["product"]
        if not fields:
            continue
        candidate = analyse_row(row, template, include_positional=False)
        if candidate is None:
            candidate = analyse_row({**row, "search_text": row.get("search_text") or ""}, template, False)
        if candidate is None:
            gene = {"common_name": row.get("gene_symbol") or row["protein_accession"]}
            candidate = {
                "gene_id": row["protein_accession"], "candidate_key": (row["assembly_accession"], row["protein_accession"]),
                "gene": {**gene, "gene_symbol": row.get("gene_symbol"), "description": row.get("description") or "",
                         "organism": row["species"], "assembly_accession": row["assembly_accession"],
                         "protein_accession": row["protein_accession"], "plaza_gene_id": row.get("plaza_gene_id"),
                         "annotation": {"go_text": row.get("go_text") or "", "pfam_text": row.get("pfam_text") or "",
                                        "interpro_text": row.get("interpro_text") or "",
                                        "mapman_text": row.get("mapman_text") or ""},
                         "trait_text": row.get("trait_text") or "",
                         "qtl_evidence": _decode_records(row.get("qtl_evidence")), "text_rank": 0.0, "total_matches": 0},
                "matches": {f"identifier:{needle}"}, "text_rank": 0.0, "tier": TIER_FUNCTIONAL, "convergent": False,
                "field_hits": {f: [needle] for f in fields}, "functional_fields": fields,
                "functional_keywords": [needle], "topic_qtls": [], "categories": ["identifier"], "n_categories": 1,
                "references": [],
            }
        candidate["score"] = 1
        candidates.append(candidate)
    total = max((int(r.get("total_matches") or 0) for r in rows), default=0)
    return {"candidates": candidates, "total_matches": total, "functional_matches": len(candidates),
            "pool_size": len(rows), "kind": kind, "identifier": value}


# ── Explanations ───────────────────────────────────────────────────────────

def _md_escape(text: str) -> str:
    return str(text).replace("$", r"\$")


def highlight_keywords(text: str, keywords: list[str]) -> str:
    """Markdown text with matched keywords in bold (case-insensitive, prefix word boundary)."""
    text = str(text or "")
    if not text or not keywords:
        return _md_escape(text)
    pattern = re.compile(
        "(?<![A-Za-z0-9])(" + "|".join(re.escape(k) for k in sorted(set(keywords), key=len, reverse=True)) + ")",
        re.IGNORECASE,
    )
    out, last = [], 0
    for match in pattern.finditer(text):
        out.append(_md_escape(text[last:match.start()]))
        out.append("**" + _md_escape(match.group(0)) + "**")
        last = match.end()
    out.append(_md_escape(text[last:]))
    return "".join(out)


def explain_candidate(candidate: dict) -> list[str]:
    """Plain-language explanation of why a protein was returned and how strong the evidence is."""
    gene = candidate["gene"]
    name = gene.get("common_name") or candidate["gene_id"]
    lines = []
    tier = candidate["tier"]
    fields = candidate.get("functional_fields", [])
    keywords = candidate.get("functional_keywords", [])
    if tier == TIER_POSITIONAL:
        lines.append(
            f"{name} has no functional annotation matching the topic. It was returned only because its gene lies "
            f"inside {len(candidate['topic_qtls'])} QTL interval(s) associated with the topic. Position alone does "
            f"not show that this gene acts on the trait."
        )
        return lines
    evidence = ", ".join(
        f"{EVIDENCE_LABELS[f]} ({', '.join(candidate['field_hits'][f][:3])})" for f in fields
    )
    lines.append(
        f"{name} matches {len(keywords)} topic keyword(s) across {len(fields)} evidence type(s): {evidence}."
    )
    if len(candidate.get("categories", [])) > 1:
        lines.append("Keywords come from several functional categories: " + ", ".join(candidate["categories"]) + ".")
    if candidate["convergent"]:
        lines.append(
            f"Convergent evidence: the gene also lies inside {len(candidate['topic_qtls'])} QTL interval(s) "
            f"reported for this topic. Function and position agree, which is the strongest case this search can make."
        )
    elif gene.get("qtl_evidence"):
        lines.append(
            f"The gene lies inside {len(gene['qtl_evidence'])} QTL interval(s), but none is labelled with this topic."
        )
    weak_only = set(fields) <= {"go"}
    if weak_only:
        lines.append("Only GO identifiers matched; GO terms in this database carry identifiers without descriptions, "
                     "so this is weaker evidence than a product, description, MapMan or InterPro match.")
    return lines


# ── GO names (optional) ────────────────────────────────────────────────────

_GO_NAMES_CACHE: dict[str, str] | None = None


def _go_names() -> dict[str, str]:
    """GO id -> name from a go-basic.obo file when available (env GO_OBO_PATH, or pipeline_data/reference)."""
    global _GO_NAMES_CACHE
    if _GO_NAMES_CACHE is not None:
        return _GO_NAMES_CACHE
    names: dict[str, str] = {}
    root = Path(__file__).resolve().parent
    candidates = [os.environ.get("GO_OBO_PATH", ""), str(root / "pipeline_data" / "reference" / "go-basic.obo"),
                  str(root / "reference" / "go-basic.obo")]
    for path in candidates:
        if path and Path(path).is_file():
            current = None
            try:
                with open(path, encoding="utf-8") as handle:
                    for line in handle:
                        line = line.strip()
                        if line.startswith("id: GO:"):
                            current = line[4:]
                        elif line.startswith("name: ") and current:
                            names[current] = line[6:]
                            current = None
            except OSError:
                names = {}
            break
    _GO_NAMES_CACHE = names
    return names


def go_labels(go_text: str) -> list[str]:
    ids = list(dict.fromkeys(re.findall(r"GO:\d{7}", go_text or "")))
    names = _go_names()
    return [f"{i} — {names[i]}" if i in names else i for i in ids]


# ── Data frames for tables and charts ──────────────────────────────────────

def _chr_sort_key(value: object) -> tuple:
    text = str(value)
    digits = re.findall(r"\d+", text)
    return (0, int(digits[0]), text) if digits else (1, 0, text)


def candidates_to_rows(candidates: list[dict]) -> list[dict]:
    rows = []
    for c in candidates:
        gene = c["gene"]
        annotation = gene.get("annotation") or {}
        pfam_ids = list(dict.fromkeys(re.findall(r"PF\d{5}", annotation.get("pfam_text", ""))))
        rows.append({
            "Priority": TIER_LABELS[c["tier"]],
            "Gene": gene.get("common_name") or c["gene_id"],
            "Species": gene.get("organism", ""),
            "Protein": gene.get("protein_accession", c["gene_id"]),
            "Score": c["score"],
            "Evidence": ", ".join(EVIDENCE_LABELS[f] for f in EVIDENCE_LABELS if f in c["field_hits"]),
            "Categories": ", ".join(c["categories"]),
            "Keywords matched": ", ".join(sorted({m.split(":", 1)[1] for m in c["matches"]})),
            "Pfam": ", ".join(pfam_ids[:4]),
            "Isoforms": c.get("isoforms", 1),
            "Topic QTLs": len(c["topic_qtls"]),
            "All QTLs": len(gene.get("qtl_evidence", [])),
            "Full-text relevance": round(c.get("text_rank", 0.0), 4),
            "PubMed refs": len(c.get("references", [])),
        })
    return rows


def species_summary_frame(per_species: dict[str, dict], labels: dict[str, str]):
    import pandas as pd
    rows = []
    for species, result in per_species.items():
        cands = result["candidates"]
        rows.append({
            "Species": labels.get(species, species),
            "Convergent": sum(c["tier"] == TIER_CONVERGENT for c in cands),
            "Functional only": sum(c["tier"] == TIER_FUNCTIONAL for c in cands),
            "Positional only (QTL)": sum(c["tier"] == TIER_POSITIONAL for c in cands),
        })
    return pd.DataFrame(rows).set_index("Species") if rows else pd.DataFrame()


def evidence_usage_frame(candidates: list[dict]):
    import pandas as pd
    counts = {EVIDENCE_LABELS[f]: sum(f in c["field_hits"] for c in candidates) for f in EVIDENCE_LABELS}
    return pd.DataFrame({"Candidates": counts})


def category_usage_frame(candidates: list[dict]):
    import pandas as pd
    counter: dict[str, int] = {}
    for c in candidates:
        for category in c["categories"]:
            counter[category] = counter.get(category, 0) + 1
    return pd.DataFrame({"Candidates": counter}).sort_values("Candidates", ascending=False) if counter else pd.DataFrame()


def heatmap_frame(candidates: list[dict], top_n: int = 25):
    import pandas as pd
    rows = []
    for c in candidates[:top_n]:
        label = f"{c['gene'].get('common_name') or c['gene_id']} · {c['gene_id']}"
        per_category: dict[str, int] = {}
        for match in c["matches"]:
            category = match.split(":", 1)[0]
            per_category[category] = per_category.get(category, 0) + 1
        for category in sorted({m.split(":", 1)[0] for cand in candidates[:top_n] for m in cand["matches"]}):
            rows.append({"Candidate": label, "Category": category, "Keywords": per_category.get(category, 0)})
    return pd.DataFrame(rows)


def qtl_frame(candidates: list[dict]):
    import pandas as pd
    rows = []
    for c in candidates:
        topic_ids = {id(r) for r in c["topic_qtls"]}
        for record in c["gene"].get("qtl_evidence", []):
            rows.append({
                "Candidate": c["gene"].get("common_name") or c["gene_id"],
                "Protein": c["gene_id"],
                "Chromosome": str(record.get("chr") or "?"),
                "Start": record.get("start"),
                "End": record.get("end"),
                "Trait": record.get("trait_detail") or record.get("trait_category") or "QTL",
                "Topic QTL": id(record) in topic_ids,
                "Explained variance (%)": record.get("explained_variance"),
                "Year": record.get("year"),
                "Reference": record.get("reference") or "",
            })
    return pd.DataFrame(rows)


def species_coverage_frame(rows: list[dict]):
    import pandas as pd
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows).set_index("Species")
    return frame


def _species_coverage() -> list[dict]:
    """Share of proteins carrying each evidence layer, per published species."""
    connection = _connect_bioinfo_pipeline()
    try:
        with connection.cursor() as cursor:
            cursor.execute(
                """
                SELECT r.display_name, r.assembly_accession, count(*),
                       count(v.description), count(v.mapman_text), count(v.interpro_text),
                       count(v.pfam_text), count(v.go_text), count(v.trait_text)
                FROM v_gene_text AS v
                JOIN species_registry AS r ON r.species = v.species
                GROUP BY r.display_name, r.assembly_accession
                ORDER BY r.display_name
                """
            )
            rows = []
            for name, assembly, total, desc, mapman, ipr, pfam, go, qtl in cursor.fetchall():
                pct = lambda n: round(100.0 * n / total, 1) if total else 0.0
                rows.append({
                    "Species": f"{name} ({assembly})", "Proteins": total,
                    "Description %": pct(desc), "MapMan %": pct(mapman), "InterPro %": pct(ipr),
                    "Pfam %": pct(pfam), "GO ids %": pct(go), "QTL %": pct(qtl),
                })
            return rows
    except Exception:
        return []
    finally:
        connection.close()


if st is not None:
    species_coverage = st.cache_data(ttl=600, show_spinner=False)(_species_coverage)
else:
    species_coverage = _species_coverage


# ── Orchestration (no Streamlit dependency, so it can be tested and scripted) ──

def run_advanced_search(
    mode: str,
    query: str,
    species_list: list[str],
    include_positional: bool = False,
    pool: int = 1000,
    match_all: bool = False,
) -> dict:
    """mode: 'topic' | 'keywords' | 'identifier'. Returns a result bundle for the UI / CLI."""
    notes: list[str] = []
    template: dict | None = None
    topic_key = None
    if mode == "topic":
        topic_key = match_topic(query) if query else None
        if not topic_key:
            return {"error": "topic_not_recognized", "notes": notes}
        template = TOPIC_TEMPLATES[topic_key]
        label = template["label"]
    elif mode == "keywords":
        template = build_custom_template(query, match_all)
        if template is None:
            return {"error": "no_keywords", "notes": notes}
        label = template["label"]
    else:
        kind, value = classify_identifier(query)
        if not value:
            return {"error": "no_identifier", "notes": notes}
        label = f"{kind}: {value}"

    per_species: dict[str, dict] = {}
    for species in species_list:
        if mode == "identifier":
            per_species[species] = search_by_identifier(species, query, pool)
        else:
            per_species[species] = search_candidates_advanced(template, species, pool, include_positional)
        if not per_species[species]["candidates"] and per_species[species].get("total_matches", 0) == 0:
            notes.append(f"No match in {species}.")
    return {
        "mode": mode, "query": query, "label": label, "template": template, "topic_key": topic_key,
        "species": per_species, "include_positional": include_positional, "match_all": match_all,
        "notes": notes, "error": None,
    }


def collapse_isoforms(candidates: list[dict]) -> list[dict]:
    """One row per gene: keep the best-ranked protein, record how many isoforms the gene has in the list."""
    groups: dict[tuple, list[dict]] = {}
    order: list[tuple] = []
    for c in candidates:
        gene = c["gene"]
        key = (gene.get("organism"), gene.get("gene_symbol") or gene.get("plaza_gene_id") or c["gene_id"])
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(c)
    out = []
    for key in order:
        members = groups[key]
        best = dict(members[0])
        best["isoforms"] = len(members)
        best["isoform_accessions"] = [m["gene"].get("protein_accession", m["gene_id"]) for m in members]
        out.append(best)
    return out


def merge_ranked(per_species: dict[str, dict]) -> list[dict]:
    merged = [c for result in per_species.values() for c in result["candidates"]]
    merged.sort(key=lambda c: (TIER_RANK[c["tier"]], c["score"], c.get("text_rank", 0.0)), reverse=True)
    return merged


# ── Exports ────────────────────────────────────────────────────────────────

def candidates_to_csv(candidates: list[dict], params: dict | None = None) -> bytes:
    import csv as _csv
    import io as _io
    buffer = _io.StringIO()
    writer = _csv.writer(buffer)
    if params:
        for key, value in params.items():
            writer.writerow([f"# {key}", value])
    columns = ["species", "protein_accession", "gene", "priority", "score", "evidence_types", "categories",
               "keywords", "description", "pfam", "interpro", "mapman", "go", "topic_qtls", "all_qtls",
               "pubmed_references"]
    writer.writerow(columns)
    for c in candidates:
        gene = c["gene"]
        annotation = gene.get("annotation") or {}
        writer.writerow([
            gene.get("organism", ""), gene.get("protein_accession", c["gene_id"]),
            gene.get("common_name") or c["gene_id"], TIER_LABELS[c["tier"]], c["score"],
            "; ".join(EVIDENCE_LABELS[f] for f in EVIDENCE_LABELS if f in c["field_hits"]),
            "; ".join(c["categories"]), "; ".join(sorted({m.split(":", 1)[1] for m in c["matches"]})),
            gene.get("description", ""), annotation.get("pfam_text", ""), annotation.get("interpro_text", ""),
            annotation.get("mapman_text", ""), "; ".join(go_labels(annotation.get("go_text", ""))),
            len(c["topic_qtls"]), len(gene.get("qtl_evidence", [])),
            " | ".join(f"{r['title']} ({r['year']}) {r['url']}" for r in c.get("references", [])),
        ])
    return ("﻿" + buffer.getvalue()).encode("utf-8")


def generate_advanced_docx(label: str, params: dict, candidates: list[dict], out_path: str, top_n: int = 30) -> str:
    from docx import Document

    doc = Document()
    doc.add_heading(f"Candidate genes — {label}", level=1)
    intro = doc.add_paragraph()
    intro.add_run(
        "Generated by Plant Gene Analyzer. Candidates are ranked by evidence level (convergent, functional, "
        "positional), then by score. The ranking is a prioritisation aid, not experimental validation. "
        "QTL membership reflects genomic position only."
    ).italic = True
    doc.add_heading("Search parameters", level=2)
    for key, value in params.items():
        doc.add_paragraph(f"{key}: {value}", style="List Bullet")

    shown = candidates[:top_n]
    doc.add_heading(f"Top {len(shown)} candidates", level=2)
    table = doc.add_table(rows=1, cols=7)
    table.style = "Light Grid Accent 1"
    for i, header in enumerate(["Gene", "Species", "Protein", "Priority", "Score", "Evidence", "Keywords"]):
        table.rows[0].cells[i].text = header
    for c in shown:
        gene = c["gene"]
        cells = table.add_row().cells
        cells[0].text = str(gene.get("common_name") or c["gene_id"])
        cells[1].text = str(gene.get("organism", ""))
        cells[2].text = str(gene.get("protein_accession", c["gene_id"]))
        cells[3].text = TIER_LABELS[c["tier"]]
        cells[4].text = str(c["score"])
        cells[5].text = ", ".join(EVIDENCE_LABELS[f] for f in EVIDENCE_LABELS if f in c["field_hits"])
        cells[6].text = ", ".join(sorted({m.split(":", 1)[1] for m in c["matches"]}))

    doc.add_heading("Why each candidate was selected", level=2)
    for c in shown:
        gene = c["gene"]
        doc.add_paragraph(f"{gene.get('common_name') or c['gene_id']} ({c['gene_id']})").runs[0].bold = True
        for line in explain_candidate(c):
            doc.add_paragraph(line)
        if gene.get("description"):
            doc.add_paragraph("Description: " + gene["description"])
        for ref in c.get("references", []):
            doc.add_paragraph(f"{ref['title']} ({ref['year']}) — {ref['url']}", style="List Bullet")
    doc.save(out_path)
    return out_path

# ── Automatic PubMed sourcing (NCBI E-utilities, free, no key required) ───

PUBMED_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_pubmed_references(gene_name: str, extra_terms: list[str] | None = None,
                             retmax: int = 3, sleep: float = 0.34) -> list[dict]:
    """Searches for up to `retmax` PubMed publications for a gene + context.
    `sleep` respects the ~3 req/s rate limit of the keyless NCBI API."""
    terms = [gene_name] + (extra_terms or [])
    query = " AND ".join(f'"{t}"[Title/Abstract]' for t in terms)

    try:
        search_resp = requests.get(
            f"{PUBMED_BASE}/esearch.fcgi",
            params={"db": "pubmed", "term": query, "retmax": retmax, "retmode": "json"},
            timeout=15,
        )
        pmids = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not pmids:
            return []
        time.sleep(sleep)

        summary_resp = requests.get(
            f"{PUBMED_BASE}/esummary.fcgi",
            params={"db": "pubmed", "id": ",".join(pmids), "retmode": "json"},
            timeout=15,
        )
        result = summary_resp.json().get("result", {})
        refs = []
        for pmid in pmids:
            doc = result.get(pmid, {})
            if not doc:
                continue
            refs.append({
                "pmid": pmid,
                "title": doc.get("title", ""),
                "journal": doc.get("fulljournalname", ""),
                "year": (doc.get("pubdate", "") or "")[:4],
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })
        time.sleep(sleep)
        return refs
    except Exception:
        return []  # a failing source should never break the whole run


# ── Word export (python-docx, to be integrated directly into the app) ─────

def generate_docx_report(topic_label: str, species: str, candidates: list[dict],
                          out_path: str, top_n: int = 30) -> str:
    from docx import Document
    from docx.shared import Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH

    doc = Document()
    doc.add_heading(f"Candidate genes — {topic_label} in {species}", level=1)
    p = doc.add_paragraph()
    p.add_run(
        f"Report automatically generated by Plant Gene Analyzer. "
        f"{len(candidates)} candidates found, top {min(top_n, len(candidates))} shown."
    ).italic = True

    table = doc.add_table(rows=1, cols=6)
    table.style = "Light Grid Accent 1"
    hdr = table.rows[0].cells
    for i, h in enumerate(
        ["Gene", "Assembly", "Accession", "Score", "Categories", "Reference(s)"]
    ):
        hdr[i].text = h

    for c in candidates[:top_n]:
        row = table.add_row().cells
        row[0].text = c["gene"].get("common_name", "") or c["gene_id"]
        row[1].text = c["gene"].get("assembly_accession", "")
        row[2].text = c["gene"].get("protein_accession", c["gene_id"])
        row[3].text = str(c["score"])
        row[4].text = ", ".join(c["categories"])
        refs = c.get("references", [])
        row[5].text = "\n".join(
            f"{r['title']} ({r['year']}) — {r['url']}" for r in refs
        ) or "Not found"

    doc.save(out_path)
    return out_path


# ── CLI entry point for testing the module outside Streamlit ──────────────



def _looks_like_systematic_id(name: str) -> bool:
    """Detects a systematic identifier (e.g. 'Zm00001e014008', a raw
    UniProt accession) rather than a real, readable gene name --
    these identifiers are useless as a PubMed search term
    (they almost never appear as-is in a title/abstract)."""
    if not name:
        return True
    normalized = name.strip()
    lowered = normalized.lower()
    if any(term in lowered for term in ("uncharacterized", "hypothetical")):
        return True
    if re.match(r"^(?:XP|NP|XM|NM)_\d+(?:\.\d+)?$", normalized, re.IGNORECASE):
        return True
    if re.match(r"^(?:Zm|Os|LOC_Os|Nitab|Nta|At)\w+$", normalized):
        return True
    if re.match(r"^[A-Za-z]{1,3}\d{4,}[A-Za-z0-9]*$", normalized):  # Zm00001e014008
        return True
    if re.match(r"^[A-Z][0-9][A-Z0-9]{3,8}$", normalized):  # raw UniProt accession (P41979...)
        return True
    return False


def _first_go_term(go_text: object) -> str | None:
    """Return one readable GO label instead of sending the whole GO blob to PubMed."""
    text = str(go_text or "").strip()
    if not text:
        return None
    first = re.split(r"\s*(?:;|\||\n)\s*", text, maxsplit=1)[0].strip()
    first = re.sub(r"^GO:\d+\s*[:\-]?\s*", "", first, flags=re.IGNORECASE)
    first = re.sub(r"\s*\[[^]]+\]\s*$", "", first).strip()
    return first or None


def pick_search_term(gene: dict, gene_id: str) -> str | None:
    """Prefer a real symbol; otherwise use a protein description or GO label."""
    for candidate in (
        gene.get("gene_symbol"),
        gene.get("common_name"),
        gene.get("description") if not _looks_like_systematic_id(
            str(gene.get("description") or "")
        ) else None,
        _first_go_term((gene.get("annotation") or {}).get("go_text")),
    ):
        term = str(candidate or "").strip()
        if term and not _looks_like_systematic_id(term):
            return term
    return None



# ── Streamlit interface (to be called from app.py, see integration guide) ──

def _render_trait_research_tab_classic() -> None:
    """Previous single-species interface, kept for reference and rollback."""
    """Render the SQL-backed topic search tab."""
    if st is None:
        raise RuntimeError("streamlit is not installed in this environment.")

    st.markdown("### 🌱 Candidate gene search by topic")
    st.markdown(
        "Search protein annotations and QTL evidence by biological topic. "
        "Published reference assemblies are loaded from species_registry."
    )

    try:
        species_records = list_available_species()
    except Exception as error:
        st.error(str(error))
        return
    if not species_records:
        st.warning("No enabled species are available in species_registry.")
        return
    display_names = {
        record["display_name"]: record["species"] for record in species_records
    }

    col1, col2 = st.columns([1, 2])
    with col1:
        species_label = st.selectbox("Species", options=list(display_names))
    with col2:
        topic_query = st.text_input(
            "Topic / problem studied",
            placeholder="e.g. lodging in maize, drought in rice...",
        )

    fetch_refs = st.checkbox(
        "Search PubMed references (slower)",
        value=True,
    )
    col_limit, col_display = st.columns(2)
    with col_limit:
        result_limit = st.slider("Candidate pool", 100, 2000, 1000, step=100)
    with col_display:
        top_n = st.slider("Number of candidates to display", 5, 50, 20)

    if st.button("🔍 Run search", type="primary", key="trait_search_run"):
        topic_key = match_topic(topic_query) if topic_query else None
        if not topic_key:
            st.session_state.pop("trait_search_result", None)
            st.warning(
                f"⚠ Topic not recognized. Currently available topics: "
                f"{', '.join(t['label'] for t in TOPIC_TEMPLATES.values())}. "
                f"To add a new topic, see TOPIC_TEMPLATES in trait_research.py."
            )
            return

        template = TOPIC_TEMPLATES[topic_key]
        species = display_names[species_label]
        try:
            with st.spinner(f"Searching indexed annotations for {species_label}..."):
                candidates = score_candidates(
                    search_candidates(template, species, result_limit),
                    template,
                )
        except Exception as error:
            st.session_state.pop("trait_search_result", None)
            st.error(str(error))
            return
        if fetch_refs:
            progress = st.progress(0, text="Searching PubMed...")
            for index, candidate in enumerate(candidates[:top_n]):
                term = pick_search_term(candidate["gene"], candidate["gene_id"])
                candidate["references"] = fetch_pubmed_references(
                    term, extra_terms=[template["pubmed_context"]]
                ) if term else []
                progress.progress((index + 1) / min(top_n, len(candidates)))
            progress.empty()
        else:
            for candidate in candidates[:top_n]:
                candidate["references"] = []

        st.session_state["trait_search_result"] = {
            "topic_key": topic_key,
            "species_label": species_label,
            "species": species,
            "candidates": candidates,
            "top_n": top_n,
        }

    search_result = st.session_state.get("trait_search_result")
    if not search_result:
        return

    topic_key = search_result["topic_key"]
    species_label = search_result["species_label"]
    species = search_result["species"]
    candidates = search_result["candidates"]
    top_n = search_result["top_n"]
    template = TOPIC_TEMPLATES[topic_key]
    if not candidates:
        st.info("No candidates found for this species/topic combination.")
        return

    total_matches = max(
        (candidate["gene"].get("total_matches", 0) for candidate in candidates),
        default=len(candidates),
    )
    st.success(
        f"Pool of {len(candidates)} candidates (total matching records: {total_matches}) — "
        f"top {min(top_n, len(candidates))} shown."
    )
    st.caption(
        "Ranking combines topic-keyword specificity, evidence-category diversity "
        "and PostgreSQL full-text relevance. It is a prioritization aid, not "
        "experimental validation."
    )

    table_rows = [{
        "Gene": c["gene"].get("common_name", "") or c["gene_id"],
        "Assembly": c["gene"].get("assembly_accession", ""),
        "Protein accession": c["gene"].get("protein_accession", c["gene_id"]),
        "Score": c["score"],
        "Categories": ", ".join(c["categories"]),
        "Full-text relevance": round(c.get("text_rank", 0.0), 4),
        "QTL records": len(c["gene"].get("qtl_evidence", [])),
        "References": len(c.get("references", [])),
    } for c in candidates[:top_n]]
    st.dataframe(table_rows, width="stretch", hide_index=True)

    with st.expander("View details of the PubMed references found"):
        for c in candidates[:top_n]:
            refs = c.get("references", [])
            if refs:
                name = c["gene"].get("common_name", "") or c["gene_id"]
                st.markdown(f"**{name}** ({c['gene_id']})")
                for r in refs:
                    st.markdown(f"- {r['title']} ({r['year']}) — [{r['pmid']}]({r['url']})")

    with st.expander("Functional and QTL evidence"):
        for candidate in candidates[:top_n]:
            gene = candidate["gene"]
            st.markdown(
                f"**{gene.get('common_name') or candidate['gene_id']}** "
                f"({candidate['gene_id']})"
            )
            st.caption(gene.get("description") or "No curated description in this database.")
            annotation = gene.get("annotation") or {}
            for label, key in (
                ("GO", "go_text"),
                ("Pfam", "pfam_text"),
                ("InterPro", "interpro_text"),
            ):
                if annotation.get(key):
                    st.markdown(f"**{label}:** {annotation[key]}")
            for qtl in gene.get("qtl_evidence", []):
                detail = qtl.get("trait_detail") or qtl.get("trait_category") or "QTL"
                reference = qtl.get("reference") or "Reference not supplied"
                evidence = []
                if qtl.get("chr") and qtl.get("start") is not None and qtl.get("end") is not None:
                    evidence.append(f"QTL interval chr {qtl['chr']}:{qtl['start']}-{qtl['end']}")
                if qtl.get("explained_variance") is not None:
                    evidence.append(f"explained variance {qtl['explained_variance']}%")
                if qtl.get("year"):
                    evidence.append(str(qtl["year"]))
                if qtl.get("resolution"):
                    evidence.append(f"{qtl['resolution']} resolution")
                summary = f" ({'; '.join(evidence)})" if evidence else ""
                st.markdown(f"- **{detail}**{summary} — {reference}")
            if not gene.get("qtl_evidence"):
                st.caption("No linked QTL record for this candidate in the current database.")

    candidate_by_key = {
        candidate["candidate_key"]: candidate for candidate in candidates[:top_n]
    }
    selected_candidate_key = st.selectbox(
        translate("ui.trait_esm2_candidate_label", default="Candidate for ESM-2 neighbor search"),
        options=list(candidate_by_key),
        format_func=lambda key: (
            f"{candidate_by_key[key]['gene'].get('common_name') or key[1]} "
            f"({key[0]} / {key[1]})"
        ),
        key="trait_esm2_candidate",
    )
    coverage_label = ", ".join(record["display_name"] for record in species_records)
    st.caption(
        "Uses the candidate's stored ESM-2 embedding; no 650M-parameter "
        f"model inference is needed. Similarity is separate from the topic "
        f"score. Published reference coverage: {coverage_label}."
    )
    if st.button(translate("ui.trait_esm2_button", default="Find ESM-2 neighbors"), key="trait_esm2_search"):
        candidate = candidate_by_key[selected_candidate_key]
        try:
            from core_engines.esm2_similarity import find_similar_by_stored_embedding

            with st.spinner("Searching the stored ESM-2 reference embeddings..."):
                hits = find_similar_by_stored_embedding(
                    candidate["gene"]["assembly_accession"],
                    candidate["gene"]["protein_accession"],
                    top_n=12,
                )
            st.session_state["trait_esm2_result"] = {
                "candidate_key": selected_candidate_key,
                "hits": hits,
            }
        except Exception as error:
            st.session_state.pop("trait_esm2_result", None)
            st.error(
                translate("ui.esm2_search_error", default="ESM-2 search is unavailable: {error}").format(error=error)
            )

    esm2_result = st.session_state.get("trait_esm2_result")
    if esm2_result and esm2_result.get("candidate_key") == selected_candidate_key:
        esm2_rows = []
        for hit in esm2_result["hits"]:
            annotations = []
            for source, terms in hit.get("annotation", {}).items():
                for term in terms:
                    label = term.get("term_desc") or term.get("term_id")
                    if label:
                        annotations.append(f"{source}: {label}")
            esm2_rows.append({
                translate("ui.esm2_accession_column", default="Accession"): hit["accession"],
                translate("ui.esm2_organism_column", default="Organism"): hit["organism"],
                translate("ui.esm2_similarity_column", default="ESM-2 cosine similarity"): hit["esm2_cosine_similarity"],
                translate("ui.esm2_annotation_column", default="Annotation"): " | ".join(annotations[:5]) or translate("ui.esm2_not_available", default="Not available"),
            })
        if esm2_rows:
            st.dataframe(esm2_rows, width="stretch", hide_index=True)
        else:
            st.info(translate("ui.trait_esm2_no_neighbors", default="No ESM-2 neighbors were returned for this sequence."))

    with tempfile.TemporaryDirectory() as temp_dir:
        out_path = Path(temp_dir) / f"report_{topic_key}_{species_label}.docx"
        generate_docx_report(template["label"], species_label, candidates, str(out_path), top_n)
        st.download_button(
            "📄 Download Word report",
            data=out_path.read_bytes(),
            file_name=f"candidates_{topic_key}_{species_label}.docx",
            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
            key="trait_search_download_report",
        )



# ── Streamlit interface: advanced search tab ───────────────────────────────

_TIER_HELP = (
    "**★ Convergent**: the protein has functional evidence for the topic AND lies inside a QTL labelled with that topic. "
    "**Functional**: annotation of the protein itself (product, description, MapMan, InterPro, Pfam, GO) matches. "
    "**Positional only**: the gene is merely inside a QTL interval; position does not show that the gene acts on the trait."
)


def _altair():
    try:
        import altair as alt
        return alt
    except ImportError:
        return None


def _candidate_label(candidate: dict) -> str:
    gene = candidate["gene"]
    return (f"{gene.get('common_name') or candidate['gene_id']} · {candidate['gene_id']} · "
            f"{gene.get('organism', '')} · {TIER_LABELS[candidate['tier']]}")


def _render_candidate_detail(candidate: dict) -> None:
    gene = candidate["gene"]
    annotation = gene.get("annotation") or {}
    keywords = sorted({m.split(":", 1)[1] for m in candidate["matches"]})
    st.markdown(f"#### {gene.get('common_name') or candidate['gene_id']}  `{candidate['gene_id']}`")
    st.caption(f"{gene.get('organism', '')} · assembly {gene.get('assembly_accession', '')} · "
               f"{TIER_LABELS[candidate['tier']]} · score {candidate['score']}")
    for line in explain_candidate(candidate):
        st.markdown("- " + _md_escape(line))

    if gene.get("description"):
        st.markdown("**Description**")
        st.markdown(highlight_keywords(gene["description"], keywords))
    if annotation.get("mapman_text"):
        st.markdown("**MapMan**")
        st.markdown(highlight_keywords(annotation["mapman_text"], keywords))
    if annotation.get("interpro_text"):
        st.markdown("**InterPro**")
        st.markdown(highlight_keywords(annotation["interpro_text"], keywords))
    if annotation.get("pfam_text"):
        st.markdown("**Pfam domains**")
        st.markdown(highlight_keywords(annotation["pfam_text"], keywords))
    labels = go_labels(annotation.get("go_text", ""))
    if labels:
        st.markdown(f"**GO** ({len(labels)} terms)")
        st.caption("Identifiers only unless a go-basic.obo file is available (set GO_OBO_PATH)." if not _go_names() else "")
        st.markdown(", ".join(_md_escape(x) for x in labels[:40]) + (" …" if len(labels) > 40 else ""))

    records = gene.get("qtl_evidence", [])
    if records:
        st.markdown(f"**QTL intervals containing this gene** ({len(records)}; "
                    f"{len(candidate['topic_qtls'])} labelled with the topic)")
        st.caption("Position only: being inside an interval is not evidence that this gene causes the trait.")
        topic_ids = {id(r) for r in candidate["topic_qtls"]}
        table = []
        for record in sorted(records, key=lambda r: (id(r) not in topic_ids, str(r.get("chr")))):
            table.append({
                "Topic QTL": "yes" if id(record) in topic_ids else "",
                "Trait": record.get("trait_detail") or record.get("trait_category") or "QTL",
                "Chromosome": record.get("chr"), "Start": record.get("start"), "End": record.get("end"),
                "Explained variance (%)": record.get("explained_variance"),
                "Year": record.get("year"), "Reference": record.get("reference") or "",
            })
        st.dataframe(table, width="stretch", hide_index=True)
    refs = candidate.get("references", [])
    if refs:
        st.markdown("**PubMed**")
        for ref in refs:
            st.markdown(f"- {ref['title']} ({ref['year']}) — [{ref['pmid']}]({ref['url']})")


def _render_esm2_tab(candidates: list[dict], selected_key) -> None:
    st.caption(
        "Uses the candidate's stored ESM-2 embedding; no model inference is needed. Similarity measures sequence "
        "context only and is independent of the topic score."
    )
    by_key = {c["candidate_key"]: c for c in candidates}
    if selected_key not in by_key:
        st.info("Select a candidate above.")
        return
    if st.button(_t("ui.trait_esm2_button", "Find ESM-2 neighbors"), key="trait_adv_esm2_btn"):
        candidate = by_key[selected_key]
        try:
            from core_engines.esm2_similarity import find_similar_by_stored_embedding
            with st.spinner("Searching the stored ESM-2 reference embeddings..."):
                hits = find_similar_by_stored_embedding(
                    candidate["gene"]["assembly_accession"], candidate["gene"]["protein_accession"], top_n=12)
            st.session_state["trait_adv_esm2"] = {"key": selected_key, "hits": hits}
        except Exception as error:
            st.session_state.pop("trait_adv_esm2", None)
            st.error(_t("ui.esm2_search_error", "ESM-2 search is unavailable: {error}", error=error))
    result = st.session_state.get("trait_adv_esm2")
    if result and result.get("key") == selected_key:
        rows = []
        for hit in result["hits"]:
            annotations = []
            for source, terms in hit.get("annotation", {}).items():
                for term in terms:
                    label = term.get("term_desc") or term.get("term_id")
                    if label:
                        annotations.append(f"{source}: {label}")
            rows.append({
                "Accession": hit["accession"], "Organism": hit["organism"],
                "ESM-2 cosine similarity": hit["esm2_cosine_similarity"],
                "Annotation": " | ".join(annotations[:5]) or "Not available",
            })
        if rows:
            st.dataframe(rows, width="stretch", hide_index=True, column_config={
                "ESM-2 cosine similarity": st.column_config.ProgressColumn(
                    "ESM-2 cosine similarity", min_value=0.0, max_value=1.0, format="%.3f"),
            })
        else:
            st.info(_t("ui.trait_esm2_no_neighbors", "No ESM-2 neighbors were returned for this sequence."))


def _render_charts_tab(bundle: dict, candidates: list[dict], labels: dict) -> None:
    alt = _altair()
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("**Candidates by evidence level, per species**")
        frame = species_summary_frame(bundle["species"], labels)
        if not frame.empty:
            st.bar_chart(frame)
        st.caption("Counts the whole candidate pool (after the filters of the Candidates tab).")
    with col_b:
        st.markdown("**Evidence types behind the candidates**")
        st.bar_chart(evidence_usage_frame(candidates))
        st.caption("A candidate can appear in several bars.")
    col_c, col_d = st.columns(2)
    with col_c:
        st.markdown("**Functional categories reached**")
        categories = category_usage_frame(candidates)
        if categories.empty:
            st.caption("No categories to show.")
        else:
            st.bar_chart(categories)
    with col_d:
        st.markdown("**Score distribution**")
        import pandas as pd
        scores = pd.Series([c["score"] for c in candidates], name="score")
        histogram = scores.value_counts().sort_index().rename("Candidates")
        st.bar_chart(histogram)
    st.markdown("**Which categories support each top candidate**")
    heat = heatmap_frame(candidates, top_n=25)
    if heat.empty:
        st.caption("Nothing to show.")
    elif alt is not None:
        chart = alt.Chart(heat).mark_rect().encode(
            x=alt.X("Category:N", title=None),
            y=alt.Y("Candidate:N", sort=None, title=None, axis=alt.Axis(labelLimit=260, labelOverlap=False)),
            color=alt.Color("Keywords:Q", scale=alt.Scale(scheme="greens"), title="Keywords hit"),
            tooltip=["Candidate", "Category", "Keywords"],
        ).properties(height=max(160, 22 * heat["Candidate"].nunique()))
        st.altair_chart(chart, width="stretch")
    else:
        st.dataframe(heat.pivot(index="Candidate", columns="Category", values="Keywords"), width="stretch")


def _render_qtl_tab(candidates: list[dict], bundle: dict) -> None:
    frame = qtl_frame(candidates)
    if frame.empty:
        st.info("None of the displayed candidates lies inside a QTL interval. QTL data exist only for species "
                "with a QTL layer (see the Data & method tab).")
        return
    alt = _altair()
    st.caption("Each bar is a QTL interval that contains at least one displayed candidate. Highlighted intervals are "
               "labelled with the searched topic. Being inside an interval is positional evidence only.")
    per_chr = frame.groupby(["Chromosome", "Topic QTL"]).size().unstack(fill_value=0)
    per_chr = per_chr.reindex(sorted(per_chr.index, key=_chr_sort_key))
    per_chr.columns = ["Topic QTL" if bool(c) else "Other QTL" for c in per_chr.columns]
    st.markdown("**QTL intervals per chromosome**")
    st.bar_chart(per_chr)
    valid = frame.dropna(subset=["Start", "End"]).drop_duplicates(["Chromosome", "Start", "End", "Trait"])
    if alt is not None and not valid.empty:
        order = sorted(valid["Chromosome"].unique(), key=_chr_sort_key)
        chart = alt.Chart(valid).mark_bar(size=9).encode(
            x=alt.X("Start:Q", title="Genomic position (bp)"), x2="End:Q",
            y=alt.Y("Chromosome:N", sort=order),
            color=alt.Color("Topic QTL:N", scale=alt.Scale(domain=[True, False], range=["#1b9e77", "#bdbdbd"]),
                            title="Labelled with topic"),
            tooltip=["Chromosome", "Start", "End", "Trait", "Candidate", "Reference"],
        ).properties(height=max(180, 26 * len(order)))
        st.altair_chart(chart, width="stretch")
    st.dataframe(frame, width="stretch", hide_index=True)


def _render_method_tab(species_records: list[dict]) -> None:
    coverage = species_coverage()
    st.markdown("**Evidence layers available per species**")
    if coverage:
        frame = species_coverage_frame(coverage)
        st.dataframe(frame, width="stretch")
        pct_columns = [c for c in frame.columns if c.endswith("%")]
        st.bar_chart(frame[pct_columns])
        st.caption("Share of each proteome that carries the layer. A species with QTL at 0 % cannot return "
                   "positional or convergent candidates.")
    else:
        st.caption("Coverage figures are unavailable.")
    st.markdown(
        """
**How candidates are found and ranked**
1. Topic keywords are looked up in the indexed text of every protein (product name, curated description, MapMan bins,
   InterPro, Pfam, GO identifiers and QTL labels). Separators `.` and `/` are treated as spaces so that
   “salt/drought” and “response.drought.stomatal” are searchable.
2. Each match is attributed to the field it comes from. Functional fields describe the protein itself;
   the QTL field only says the gene lies inside an interval.
3. Priority order: convergent, functional, positional. Inside a level: score = 3 × strong keywords + 1 × weak
   keywords + 2 × additional categories + source weights (product, description, MapMan 2; InterPro, Pfam, GO 1)
   + 4 when a topic QTL also contains the gene.
4. Functional hits always enter the candidate pool before QTL-only hits, so large QTL intervals cannot displace them.

**What the data can and cannot say**
- PLAZA-derived annotations (GO, InterPro, MapMan, QTL links) are used only for proteins whose PLAZA gene matches the
  RefSeq annotation by chromosome, strand and at least 50 % CDS overlap. Other proteins keep their RefSeq and
  Pfam annotations only.
- GO terms carry identifiers without descriptions in this database. They count as weaker evidence, and keyword search
  cannot match GO names unless a go-basic.obo file is provided.
- QTL data exist only for species with a QTL layer. For other species the search is purely functional.
- The ranking prioritises candidates for follow-up. It is not experimental validation.
        """
    )


def render_trait_research_tab() -> None:
    """Advanced candidate search: topics, free keywords, identifiers, multi-species, explained evidence."""
    if st is None:
        raise RuntimeError("streamlit is not installed in this environment.")

    st.markdown("### 🌱 " + _t("ui.trait_title", "Candidate gene search"))
    st.markdown(_t(
        "ui.trait_intro",
        "Search the published reference proteomes by biological topic, free keywords or an identifier. Candidates are "
        "ranked by evidence level, and every result explains why it was selected."))

    try:
        species_records = list_available_species()
    except Exception as error:
        st.error(str(error))
        return
    if not species_records:
        st.warning("No enabled species are available in species_registry.")
        return
    labels = {r["species"]: r["display_name"] for r in species_records}

    mode_labels = {
        "topic": "Topic",
        "keywords": "Free keywords",
        "identifier": "Identifier (Pfam, InterPro, GO, MapMan, accession)",
    }
    mode = st.radio("Search mode", list(mode_labels), format_func=mode_labels.get, horizontal=True,
                    key="trait_adv_mode")
    col1, col2 = st.columns([1, 2])
    with col1:
        chosen_species = st.multiselect(
            "Species", options=list(labels), default=list(labels)[:1], format_func=labels.get,
            key="trait_adv_species", help="Select several species to compare them in one search.")
    match_all = False
    with col2:
        if mode == "topic":
            query = st.text_input(
                "Topic / problem studied", key="trait_adv_topic",
                placeholder="e.g. drought, salinity, cold, lodging, nicotine, flowering time…")
            st.caption("Available topics: " + ", ".join(t["label"] for t in TOPIC_TEMPLATES.values()))
        elif mode == "keywords":
            query = st.text_area(
                "Keywords (comma, semicolon or one per line)", key="trait_adv_keywords", height=90,
                placeholder="dreb, aquaporin, nac transcription factor")
            match_all = st.radio("Combine keywords with", ["Any keyword (OR)", "All keywords (AND)"],
                                 horizontal=True, key="trait_adv_combine") == "All keywords (AND)"
        else:
            query = st.text_input(
                "Identifier or name", key="trait_adv_identifier",
                placeholder="PF00069, IPR000719, GO:0006629, 26.5.1.1, XP_015643239.1 …")
            if query.strip():
                kind, value = classify_identifier(query)
                st.caption(f"Recognised as: {kind} ({value})")

    with st.expander("Advanced options"):
        include_positional = st.checkbox(
            "Include positional (QTL-only) candidates", value=False, key="trait_adv_qtl",
            help="Adds genes that are merely located inside a QTL labelled with the topic. They are ranked last.")
        opt1, opt2 = st.columns(2)
        with opt1:
            pool = st.slider("Candidate pool per species", 100, 2000, 1000, step=100, key="trait_adv_pool")
        with opt2:
            top_n = st.slider("Candidates to display", 5, 100, 25, key="trait_adv_topn")
        fetch_refs = st.checkbox("Search PubMed references for the top candidates (slow)", value=False,
                                 key="trait_adv_refs")

    if st.button("🔍 Run search", type="primary", key="trait_adv_run"):
        if not chosen_species:
            st.warning("Select at least one species.")
            return
        if not (query or "").strip():
            st.warning("Enter a topic, keywords or an identifier.")
            return
        try:
            with st.spinner("Searching indexed annotations..."):
                bundle = run_advanced_search(mode, query, chosen_species, include_positional, pool, match_all)
        except Exception as error:
            st.session_state.pop("trait_adv_result", None)
            st.error(str(error))
            return
        if bundle.get("error") == "topic_not_recognized":
            st.session_state.pop("trait_adv_result", None)
            st.warning(
                "Topic not recognised. Available topics: "
                + ", ".join(t["label"] for t in TOPIC_TEMPLATES.values())
                + ". Switch to **Free keywords** to search any term.")
            return
        if bundle.get("error"):
            st.session_state.pop("trait_adv_result", None)
            st.warning("Nothing to search for.")
            return
        merged = merge_ranked(bundle["species"])
        if fetch_refs and merged:
            context = (bundle["template"] or {}).get("pubmed_context", "") if bundle.get("template") else ""
            progress = st.progress(0, text="Searching PubMed...")
            shown = merged[:top_n]
            for index, candidate in enumerate(shown):
                term = pick_search_term(candidate["gene"], candidate["gene_id"])
                candidate["references"] = fetch_pubmed_references(
                    term, extra_terms=[context] if context else []) if term else []
                progress.progress((index + 1) / len(shown))
            progress.empty()
        st.session_state["trait_adv_result"] = {"bundle": bundle, "top_n": top_n, "labels": labels}

    stored = st.session_state.get("trait_adv_result")
    if not stored:
        return
    bundle, top_n, labels = stored["bundle"], stored["top_n"], stored["labels"]
    merged = merge_ranked(bundle["species"])

    total = sum(r.get("total_matches", 0) for r in bundle["species"].values())
    functional_total = sum(r.get("functional_matches", 0) for r in bundle["species"].values())
    pool_total = sum(r.get("pool_size", 0) for r in bundle["species"].values())
    n_conv = sum(c["tier"] == TIER_CONVERGENT for c in merged)
    n_func = sum(c["tier"] == TIER_FUNCTIONAL for c in merged)
    n_pos = sum(c["tier"] == TIER_POSITIONAL for c in merged)

    st.markdown(f"#### Results — {bundle['label']}")
    metrics = st.columns(4)
    metrics[0].metric("Matching proteins in the database", f"{total:,}")
    validated = n_conv + n_func
    metrics[1].metric("With functional evidence", f"{validated:,}" if bundle["mode"] != "identifier" else f"{total:,}")
    if bundle["mode"] != "identifier" and functional_total > validated:
        st.caption(f"{functional_total:,} proteins matched the text query; {validated:,} were confirmed by an exact "
                   f"keyword in their own annotation (the others matched only loosely, e.g. as part of a longer word).")
    metrics[2].metric("★ Convergent in pool", n_conv)
    metrics[3].metric("Positional only in pool", n_pos if bundle["include_positional"] else "excluded")
    if pool_total and total > pool_total:
        st.caption(f"The pool holds the {pool_total:,} best-ranked proteins out of {total:,}; functional hits are "
                   "always ranked before QTL-only hits.")
    for note in bundle.get("notes", []):
        st.caption(note)
    if bundle["mode"] == "identifier":
        st.caption("Identifier search returns every protein carrying the identifier, ordered by accession.")
    if not merged:
        st.info("No candidates found for this search. Try fewer or broader keywords, another species, or enable "
                "positional candidates in the advanced options.")
        return
    st.caption(_TIER_HELP)
    if st.checkbox("One row per gene (group isoforms of the same gene)", value=True, key="trait_adv_collapse",
                   help="Several proteins (isoforms) often come from one gene; grouping avoids a list full of "
                        "near-identical rows. The Isoforms column gives the count."):
        n_before = len(merged)
        merged = collapse_isoforms(merged)
        st.caption(f"{n_before:,} proteins grouped into {len(merged):,} genes.")

    shown = merged[:top_n]  # replaced by the filtered list once the filters in the first tab are applied
    tabs = st.tabs(["Candidates", "Why this candidate", "Charts", "QTL map", "ESM-2 neighbours",
                    "Data & method", "Export"])

    with tabs[0]:
        f1, f2, f3 = st.columns([2, 1, 2])
        present = [t for t in (TIER_CONVERGENT, TIER_FUNCTIONAL, TIER_POSITIONAL) if any(c["tier"] == t for c in merged)]
        with f1:
            tier_filter = st.multiselect("Evidence level", present, default=present,
                                         format_func=TIER_LABELS.get, key="trait_adv_tier_filter")
        max_score = max(c["score"] for c in merged)
        with f2:
            min_score = st.slider("Minimum score", 0, max(max_score, 1), 0, key="trait_adv_min_score")
        with f3:
            text_filter = st.text_input("Filter by name, description or domain", key="trait_adv_text_filter")
        needle = text_filter.strip().lower()
        visible = [
            c for c in merged
            if c["tier"] in tier_filter and c["score"] >= min_score and (
                not needle or needle in " ".join([
                    str(c["gene"].get("common_name", "")), c["gene"].get("description", ""),
                    c["gene"].get("annotation", {}).get("pfam_text", ""), c["gene_id"]]).lower())
        ]
        shown = visible[:top_n]
        st.caption(f"{len(shown)} displayed of {len(visible)} matching the filters "
                   f"({n_conv} convergent, {n_func} functional, {n_pos} positional in the pool).")
        st.dataframe(candidates_to_rows(shown), width="stretch", hide_index=True, column_config={
            "Score": st.column_config.ProgressColumn("Score", min_value=0, max_value=max(max_score, 1), format="%d"),
        })

    options = [c["candidate_key"] for c in shown] or [merged[0]["candidate_key"]]
    by_key = {c["candidate_key"]: c for c in merged}
    st.session_state.setdefault("trait_adv_selected", options[0])
    if st.session_state["trait_adv_selected"] not in options:
        st.session_state["trait_adv_selected"] = options[0]

    with tabs[1]:
        selected_key = st.selectbox("Inspect a candidate", options, key="trait_adv_selected",
                                    format_func=lambda k: _candidate_label(by_key[k]))
        _render_candidate_detail(by_key[selected_key])

    with tabs[2]:
        _render_charts_tab(bundle, visible or merged[:top_n], labels)

    with tabs[3]:
        _render_qtl_tab(shown or merged[:top_n], bundle)

    with tabs[4]:
        _render_esm2_tab(shown or merged[:top_n], st.session_state.get("trait_adv_selected"))

    with tabs[5]:
        _render_method_tab(species_records)

    with tabs[6]:
        params = {
            "search": bundle["label"], "mode": bundle["mode"],
            "species": ", ".join(labels.get(s, s) for s in bundle["species"]),
            "positional candidates included": bundle["include_positional"],
            "matching proteins": total, "pool": pool_total,
        }
        st.download_button("⬇ Download CSV (all pool candidates)", data=candidates_to_csv(merged, params),
                           file_name="candidates.csv", mime="text/csv", key="trait_adv_csv")
        with tempfile.TemporaryDirectory() as temp_dir:
            out_path = Path(temp_dir) / "candidates.docx"
            generate_advanced_docx(bundle["label"], params, merged, str(out_path), top_n)
            st.download_button(
                "📄 Download Word report", data=out_path.read_bytes(), file_name="candidates.docx",
                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                key="trait_adv_docx")


def main():
    ap = argparse.ArgumentParser(description="Advanced candidate search on the bioinfo_pipeline database.")
    ap.add_argument("--species", help="Species key(s) from species_registry, comma separated (e.g. zea_mays,oryza_sativa)")
    ap.add_argument("--topic", help="Topic query, e.g. 'lodging' or 'drought'")
    ap.add_argument("--keywords", help="Free keywords, comma separated")
    ap.add_argument("--all-keywords", action="store_true", help="With --keywords: require every keyword (AND)")
    ap.add_argument("--identifier", help="Pfam / InterPro / GO / MapMan identifier, accession or gene name")
    ap.add_argument("--include-qtl", action="store_true", help="Also return QTL-only (positional) candidates")
    ap.add_argument("--list-species", action="store_true", help="List registered searchable species")
    ap.add_argument("--candidate-pool", type=int, default=1000)
    ap.add_argument("--fetch-refs", action="store_true", help="Query PubMed (slow, ~0.35s/gene)")
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--out", default="report_candidates.docx")
    ap.add_argument("--csv", help="Also write all pool candidates to this CSV file")
    args = ap.parse_args()

    if args.list_species:
        for record in list_available_species():
            print(f"{record['species']}\t{record['display_name']}")
        return
    if not args.species or not (args.topic or args.keywords or args.identifier):
        ap.error("--species and one of --topic / --keywords / --identifier are required unless --list-species is used")

    species_list = [resolve_species_filter(s) for s in args.species.split(",") if s.strip()]
    if args.topic:
        mode, query = "topic", args.topic
    elif args.keywords:
        mode, query = "keywords", args.keywords
    else:
        mode, query = "identifier", args.identifier
    bundle = run_advanced_search(mode, query, species_list, args.include_qtl, args.candidate_pool, args.all_keywords)
    if bundle.get("error"):
        print(f"⚠ Nothing to search: {bundle['error']}")
        if bundle["error"] == "topic_not_recognized":
            print(f"  Available topics: {list(TOPIC_TEMPLATES.keys())} — or use --keywords.")
        return

    merged = merge_ranked(bundle["species"])
    print(f"Search: {bundle['label']}")
    for species, result in bundle["species"].items():
        cands = result["candidates"]
        print(f"  {species}: {result.get('total_matches', 0)} matching proteins, pool {result.get('pool_size', 0)}, "
              f"{sum(c['tier'] == TIER_CONVERGENT for c in cands)} convergent, "
              f"{sum(c['tier'] == TIER_FUNCTIONAL for c in cands)} functional, "
              f"{sum(c['tier'] == TIER_POSITIONAL for c in cands)} positional")
    print()
    for c in merged[:args.top_n]:
        evidence = ", ".join(EVIDENCE_LABELS[f] for f in EVIDENCE_LABELS if f in c["field_hits"])
        print(f"  [{c['score']:2d}] {TIER_LABELS[c['tier']]:<32} {c['gene_id']} "
              f"{c['gene'].get('common_name', '')} ({c['gene'].get('organism', '')}) — {evidence}")

    if args.fetch_refs and bundle.get("template"):
        context = bundle["template"].get("pubmed_context", "")
        print("\nSearching PubMed...")
        for c in merged[:args.top_n]:
            term = pick_search_term(c["gene"], c["gene_id"])
            if not term:
                print(f"  {c['gene_id']}: no usable search term, skipped")
                continue
            c["references"] = fetch_pubmed_references(term, extra_terms=[context] if context else [])
            print(f"  {term}: {len(c['references'])} reference(s)")

    params = {"search": bundle["label"], "mode": mode, "species": ", ".join(species_list),
              "positional candidates included": args.include_qtl}
    out = generate_advanced_docx(bundle["label"], params, merged, args.out, args.top_n)
    print(f"\n✓ Report written: {out}")
    if args.csv:
        Path(args.csv).write_bytes(candidates_to_csv(merged, params))
        print(f"✓ CSV written: {args.csv}")


if __name__ == "__main__":
    main()
