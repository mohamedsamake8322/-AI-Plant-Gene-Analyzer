"""
trait_research.py — Topic-based candidate gene search engine,
generalized to any species/topic in the dataset (not just
lodging/quinoa). Designed to be called directly from app.py via a
new Streamlit "Topic search" tab.

Pipeline:
  1. The topic typed by the user ("lodging in quinoa", "drought
     in maize") is matched against a library of pre-written
     templates (TOPIC_TEMPLATES). If no template matches, a
     clear message indicates this (future extension: keyword
     generation via an LLM for uncovered topics).
  2. search_candidates() scans the <species>_all_sources.json file for
     the template's keywords, restricted to genes with origin=sequence_backed.
  3. score_candidates() sorts by specificity + category diversity.
  4. fetch_pubmed_references() queries the NCBI E-utilities API (free,
     no key required) for each retained gene.
  5. generate_docx_report() produces a downloadable .docx directly
     from the app.

Direct usage (outside Streamlit, for testing):
    python trait_research.py --species-file path.json --topic "lodging"
"""

from __future__ import annotations

import json
import re
import time
import argparse
import unicodedata
from pathlib import Path
from dataclasses import dataclass, field
from urllib.parse import quote

import requests


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

SPECIES_ALIASES: dict[str, str] = {
    "quinoa": "chenopodium quinoa",
    "riz": "oryza sativa", "rice": "oryza sativa",
    "mais": "zea mays", "maïs": "zea mays", "maize": "zea mays", "corn": "zea mays",
    "tomate": "solanum lycopersicum", "tomato": "solanum lycopersicum",
    "raisin": "vitis vinifera", "grape": "vitis vinifera", "vigne": "vitis vinifera",
    "tabac": "nicotiana tabacum", "tobacco": "nicotiana tabacum",
    "pomme de terre": "solanum tuberosum", "potato": "solanum tuberosum",
}


def resolve_species_filter(species_input: str | None) -> str | None:
    """Normalizes a species name typed by the user (common or
    scientific, in any common language) to the scientific name
    actually stored in the dataset's "organism" field.
    Without this, a common name ("maize") never matches "Zea mays" and
    silently filters out ALL genes -- this is exactly the bug
    encountered when testing the drought topic on maize (0 candidates)."""
    if not species_input:
        return None
    key = _normalize_text(species_input)
    return SPECIES_ALIASES.get(key, key)  # if not in the table, try as-is


TOPIC_TEMPLATES: dict[str, dict] = {
    "verse": {
        "aliases": ["verse", "lodging", "tige", "rigidite"],
        "label": "Lodging / stem strength resistance",
        "pubmed_context": "lodging",
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
        "keywords": {
            "aba_signaling": [
                "abscisic acid", "aba receptor", "pyr/pyl", "snrk2",
                "aba signaling",
            ],
            "osmotic_response": [
                "osmotic stress", "proline", "dehydrin", "late embryogenesis",
                "lea protein", "osmotic adjustment",
            ],
            "stomatal_regulation": [
                "stomatal closure", "guard cell", "aquaporin", "water use efficiency",
            ],
            "transcription_factors": [
                "dreb", "nac transcription factor", "wrky", "myb drought",
            ],
        },
        "tier_a": {
            "abscisic acid", "aba receptor", "pyr/pyl", "snrk2",
            "dehydrin", "lea protein", "dreb", "aquaporin",
            "stomatal closure",
        },
        "tier_b": {
            "aba signaling", "osmotic stress", "proline",
            "late embryogenesis", "osmotic adjustment", "guard cell",
            "water use efficiency", "nac transcription factor", "wrky",
            "myb drought",
        },
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


# ── Candidate search (generalized from search_lodging_candidates_v2.py) ────

def _collect_text(gene: dict) -> str:
    parts = []

    def add(v):
        if isinstance(v, str):
            parts.append(v.lower())
        elif isinstance(v, dict):
            for x in v.values():
                add(x)
        elif isinstance(v, list):
            for x in v:
                add(x)

    for field_name in ("description", "common_name", "traits"):
        add(gene.get(field_name))
    annotation = gene.get("annotation") or {}
    for field_name in ("go_terms", "kegg_pathways", "mapman", "tf_family"):
        add(annotation.get(field_name))
    add(gene.get("literature"))
    return " | ".join(parts)


def _has_real_sequence(gene: dict) -> bool:
    seq = gene.get("sequence") or {}
    return bool(seq.get("dna") or seq.get("rna") or seq.get("protein"))


def load_genes(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    if isinstance(raw, dict) and "genes" in raw:
        raw = raw["genes"]
    if isinstance(raw, dict) and all(isinstance(v, dict) for v in raw.values()):
        return raw
    if isinstance(raw, list):
        return {g.get("gene_id", str(i)): g for i, g in enumerate(raw)}
    raise ValueError("Unrecognized file format")


def search_candidates(genes: dict, template: dict, species_filter: str | None = None) -> list[dict]:
    all_kw = [(kw, cat) for cat, kws in template["keywords"].items() for kw in kws]
    noisy_tagged = {f"{cat}:{kw}" for cat, kws in template["keywords"].items()
                     for kw in kws if kw in template["noisy"]}

    results = []
    for gene_id, gene in genes.items():
        if not isinstance(gene, dict):
            continue
        if gene.get("origin") != "sequence_backed":
            continue
        if species_filter:
            organism = str(gene.get("organism", "")).lower()
            if species_filter.lower() not in organism:
                continue
        text = _collect_text(gene)
        matches = {f"{cat}:{kw.strip()}" for kw, cat in all_kw if kw in text} - noisy_tagged
        if not matches or not _has_real_sequence(gene):
            continue
        results.append({"gene_id": gene_id, "gene": gene, "matches": matches})
    return results


def score_candidates(candidates: list[dict], template: dict) -> list[dict]:
    tier_a, tier_b = template["tier_a"], template["tier_b"]
    for c in candidates:
        kws = {m.split(":", 1)[1] for m in c["matches"]}
        n_tier_a = len(kws & tier_a)
        n_tier_b = len(kws & tier_b)
        n_categories = len({m.split(":")[0] for m in c["matches"]})
        score = n_tier_a * 3 + n_tier_b * 1 + max(0, n_categories - 1) * 2
        c["score"] = score
        c["n_categories"] = n_categories
        c["categories"] = sorted({m.split(":")[0] for m in c["matches"]})
    candidates.sort(key=lambda c: c["score"], reverse=True)
    return candidates


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

    table = doc.add_table(rows=1, cols=5)
    table.style = "Light Grid Accent 1"
    hdr = table.rows[0].cells
    for i, h in enumerate(["Gene", "Accession", "Score", "Categories", "Reference(s)"]):
        hdr[i].text = h

    for c in candidates[:top_n]:
        row = table.add_row().cells
        row[0].text = c["gene"].get("common_name", "") or c["gene_id"]
        row[1].text = c["gene_id"]
        row[2].text = str(c["score"])
        row[3].text = ", ".join(c["categories"])
        refs = c.get("references", [])
        row[4].text = "\n".join(f"{r['title']} ({r['year']}) — {r['url']}" for r in refs) or "Not found"

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
    if re.match(r"^[A-Za-z]{1,3}\d{4,}[A-Za-z0-9]*$", name):  # Zm00001e014008
        return True
    if re.match(r"^[A-Z][0-9][A-Z0-9]{3,8}$", name):  # raw UniProt accession (P41979...)
        return True
    return False


def pick_search_term(gene: dict, gene_id: str) -> str | None:
    """Picks the best available PubMed search term for a
    gene: readable gene name in priority, otherwise a GO term
    (molecular_function/biological_process given priority, more
    informative than cellular_component), otherwise the last segment
    of a mapman description, otherwise give up."""
    common_name = (gene.get("common_name") or "").split(":")[0].strip()
    if common_name and not _looks_like_systematic_id(common_name):
        return common_name

    annotation = gene.get("annotation") or {}
    go_terms = annotation.get("go_terms") or []
    for pref_aspect in ("molecular_function", "biological_process"):
        for gt in go_terms:
            if isinstance(gt, dict) and gt.get("aspect") == pref_aspect and gt.get("term"):
                return gt["term"]
    if go_terms and isinstance(go_terms[0], dict) and go_terms[0].get("term"):
        return go_terms[0]["term"]

    mapman = annotation.get("mapman") or []
    if mapman and isinstance(mapman[0], dict) and mapman[0].get("description"):
        # e.g. "Cell wall organisation.lignin.monolignol conjugation..." -> last segment
        return mapman[0]["description"].split(".")[-1].strip()

    return None


# ── Streamlit interface (to be called from app.py, see integration guide) ──

try:
    import streamlit as st
except ImportError:
    st = None

SPECIES_FILES: dict[str, str] = {
    "quinoa": "chenopodium_quinoa_all_sources.json",
    "rice": "oryza_sativa_all_sources.json",
    "maize": "zea_mays_all_sources.json",
    "tomato": "solanum_lycopersicum_all_sources.json",
    "grape": "vitis_vinifera_all_sources.json",
    "tobacco": "nicotiana_tabacum_all_sources.json",
    "potato": "solanum_tuberosum_all_sources.json",
}


def render_trait_research_tab(species_dir: str) -> None:
    """Complete Streamlit "Topic search" section. To be called from
    app.py with the path to the folder containing the
    <species>_all_sources.json files (see integration guide).

    Does NOT depend on any sequence entered by the user -- works
    completely independently of the existing sequence analysis flow.
    """
    if st is None:
        raise RuntimeError("streamlit is not installed in this environment.")

    st.markdown("### 🌱 Candidate gene search by topic")
    st.markdown(
        "Choose a species and describe an agronomic problem "
        "(e.g. *lodging in quinoa*, *drought in maize*) to get "
        "a list of candidate genes sourced from the scientific literature."
    )

    col1, col2 = st.columns([1, 2])
    with col1:
        species_label = st.selectbox("Species", options=list(SPECIES_FILES.keys()))
    with col2:
        topic_query = st.text_input(
            "Topic / problem studied",
            placeholder="e.g. lodging, drought, cold resistance...",
        )

    fetch_refs = st.checkbox(
        "Search PubMed references (slower, ~0.7s per gene)",
        value=True,
    )
    top_n = st.slider("Number of candidates to display", 5, 50, 20)

    if not st.button("🔍 Run search", type="primary"):
        return

    topic_key = match_topic(topic_query) if topic_query else None
    if not topic_key:
        st.warning(
            f"⚠ Topic not recognized. Currently available topics: "
            f"{', '.join(t['label'] for t in TOPIC_TEMPLATES.values())}. "
            f"To add a new topic, see TOPIC_TEMPLATES in trait_research.py."
        )
        return

    template = TOPIC_TEMPLATES[topic_key]
    species_file = Path(species_dir) / SPECIES_FILES[species_label]

    with st.spinner(f"Loading {species_label} data..."):
        genes = _load_genes_cached(str(species_file))

    species_filter = resolve_species_filter(species_label)
    candidates = search_candidates(genes, template, species_filter)
    candidates = score_candidates(candidates, template)

    if not candidates:
        st.info("No candidates found for this species/topic combination.")
        return

    st.success(f"{len(candidates)} candidates found — top {min(top_n, len(candidates))} shown.")

    if fetch_refs:
        progress = st.progress(0, text="Searching PubMed...")
        for i, c in enumerate(candidates[:top_n]):
            term = pick_search_term(c["gene"], c["gene_id"])
            c["references"] = fetch_pubmed_references(term, extra_terms=[template["pubmed_context"]]) if term else []
            progress.progress((i + 1) / min(top_n, len(candidates)))
        progress.empty()
    else:
        for c in candidates[:top_n]:
            c["references"] = []

    table_rows = [{
        "Gene": c["gene"].get("common_name", "") or c["gene_id"],
        "Accession": c["gene_id"],
        "Score": c["score"],
        "Categories": ", ".join(c["categories"]),
        "References": len(c.get("references", [])),
    } for c in candidates[:top_n]]
    st.dataframe(table_rows, width="stretch")

    with st.expander("View details of the PubMed references found"):
        for c in candidates[:top_n]:
            refs = c.get("references", [])
            if refs:
                name = c["gene"].get("common_name", "") or c["gene_id"]
                st.markdown(f"**{name}** ({c['gene_id']})")
                for r in refs:
                    st.markdown(f"- {r['title']} ({r['year']}) — [{r['pmid']}]({r['url']})")

    out_path = f"/tmp/report_{topic_key}_{species_label}.docx"
    generate_docx_report(template["label"], species_label, candidates, out_path, top_n)
    with open(out_path, "rb") as f:
        st.download_button(
            "📄 Download Word report",
            data=f.read(),
            file_name=f"candidates_{topic_key}_{species_label}.docx",
            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
        )


def _load_genes_cached(species_file: str) -> dict:
    """Streamlit cache wrapper around load_genes -- avoids re-reading
    and re-parsing a file with tens of thousands of genes on
    every interaction with a page widget."""
    if st is not None:
        cached = st.cache_data(show_spinner=False)(load_genes)
        return cached(Path(species_file))
    return load_genes(Path(species_file))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-file", required=True)
    ap.add_argument("--topic", required=True)
    ap.add_argument("--species-name", default=None)
    ap.add_argument("--fetch-refs", action="store_true", help="Query PubMed (slow, ~0.35s/gene)")
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--out", default="report_candidates.docx")
    args = ap.parse_args()

    topic_key = match_topic(args.topic)
    if not topic_key:
        print(f"⚠ Topic '{args.topic}' not covered by the template library.")
        print(f"  Available topics: {list(TOPIC_TEMPLATES.keys())}")
        return

    template = TOPIC_TEMPLATES[topic_key]
    genes = load_genes(Path(args.species_file))
    species_filter = resolve_species_filter(args.species_name)
    candidates = search_candidates(genes, template, species_filter)
    candidates = score_candidates(candidates, template)

    print(f"{len(candidates)} candidates found for '{template['label']}'.\n")
    for c in candidates[:args.top_n]:
        print(f"  [{c['score']:2d}] {c['gene_id']} — {c['gene'].get('common_name', '')} ({', '.join(c['categories'])})")

    if args.fetch_refs:
        print("\nSearching PubMed...")
        for c in candidates[:args.top_n]:
            term = pick_search_term(c["gene"], c["gene_id"])
            if not term:
                c["references"] = []
                print(f"  {c['gene_id']}: no usable search term, skipped")
                continue
            c["references"] = fetch_pubmed_references(term, extra_terms=[template["pubmed_context"]])
            print(f"  {term}: {len(c['references'])} reference(s)")

    out = generate_docx_report(template["label"], args.species_name or "?", candidates, args.out, args.top_n)
    print(f"\n✓ Report written: {out}")


if __name__ == "__main__":
    main()