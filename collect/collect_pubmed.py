#!/usr/bin/env python3
"""
PubMed collector — fetch publications associated with plant genes/species.
Uses NCBI E-utilities (no key required, but EMAIL recommended for higher rate limits).
Set NCBI_EMAIL in your .env for better throughput.
"""

from __future__ import annotations

import os
import time
import xml.etree.ElementTree as ET
from typing import Optional
import requests
import request_utils as rq

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Optional: set your email for NCBI rate limits (3 req/s → 10 req/s)
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "researcher@example.com")
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")


def fetch_pubmed_for_species(
    species: str,
    retmax: int = 100,
    extra_terms: Optional[list[str]] = None,
) -> list[dict]:
    """
    Fetch PubMed publications for a plant species.

    Args:
        species: Scientific name (e.g. "Oryza sativa")
        retmax: Max publications to retrieve
        extra_terms: Additional MeSH/keyword terms (e.g. ["drought stress", "gene expression"])

    Returns:
        List of publication records
    """
    query_parts = [f'"{species}"[Organism]', "plant[filter]"]
    if extra_terms:
        query_parts += [f'"{t}"' for t in extra_terms]
    query = " AND ".join(query_parts)

    print(f"  [PubMed] Searching: {query}")
    pmids = _esearch(query, retmax)
    if not pmids:
        return []

    print(f"  [PubMed] Found {len(pmids)} publications, fetching details...")
    return _efetch_summaries(pmids, species)


def fetch_pubmed_for_gene(
    gene_symbol: str,
    species: str,
    retmax: int = 20,
) -> list[dict]:
    """
    Fetch publications for a specific gene in a species.

    Args:
        gene_symbol: Gene name/symbol (e.g. "OsGA20ox1")
        species: Scientific name
        retmax: Max publications

    Returns:
        List of publication records (to embed in gene record's publications field)
    """
    query = f'"{gene_symbol}"[Title/Abstract] AND "{species}"[Organism]'
    pmids = _esearch(query, retmax)
    if not pmids:
        return []
    return _efetch_summaries(pmids, species)


def fetch_pubmed_bulk_genes(
    gene_symbols: list[str],
    species: str,
    pubs_per_gene: int = 5,
) -> dict[str, list[dict]]:
    """
    Efficiently fetch publications for multiple genes in one species.
    Returns a dict mapping gene_symbol -> list of publication records.

    Args:
        gene_symbols: List of gene symbols to look up
        species: Scientific name
        pubs_per_gene: Max publications per gene

    Returns:
        Dict of gene_symbol -> [pub_records]
    """
    results: dict[str, list[dict]] = {}
    for symbol in gene_symbols:
        pubs = fetch_pubmed_for_gene(symbol, species, pubs_per_gene)
        results[symbol] = pubs
        time.sleep(0.15)  # NCBI rate limit
    return results


def _esearch(query: str, retmax: int) -> list[str]:
    """ESearch: get list of PMIDs for a query."""
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": min(retmax, 10000),
        "retmode": "json",
        "email": NCBI_EMAIL,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        resp = rq.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=20)
        data = resp.json()
        return data.get("esearchresult", {}).get("idlist", [])
    except (requests.RequestException, ValueError) as e:
        print(f"  [PubMed] ESearch error: {e}")
        return []


def _efetch_summaries(pmids: list[str], species: str) -> list[dict]:
    """EFetch: fetch full summaries for a list of PMIDs."""
    records = []
    batch_size = 100  # NCBI recommends ≤200 per request

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml",
            "email": NCBI_EMAIL,
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        try:
            resp = rq.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=30)
            batch_records = _parse_pubmed_xml(resp.text, species)
            records.extend(batch_records)
            time.sleep(0.1)

        except requests.RequestException as e:
            print(f"  [PubMed] EFetch error (batch {i}): {e}")

    return records


def _parse_pubmed_xml(xml_text: str, species: str) -> list[dict]:
    """Parse PubMed XML into normalized publication records."""
    records = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    for article in root.findall(".//PubmedArticle"):
        try:
            rec = _parse_article(article, species)
            if rec:
                records.append(rec)
        except Exception:
            continue

    return records


def _parse_article(article: ET.Element, species: str) -> dict | None:
    """Parse a single PubMed article XML element."""
    medline = article.find("MedlineCitation")
    if medline is None:
        return None

    pmid_el = medline.find("PMID")
    pmid = pmid_el.text if pmid_el is not None else ""
    if not pmid:
        return None

    art = medline.find("Article")
    if art is None:
        return None

    # Title
    title_el = art.find("ArticleTitle")
    title = title_el.text or "" if title_el is not None else ""

    # Abstract
    abstract_parts = []
    for ab_text in art.findall(".//AbstractText"):
        label = ab_text.get("Label", "")
        text = ab_text.text or ""
        if label:
            abstract_parts.append(f"{label}: {text}")
        else:
            abstract_parts.append(text)
    abstract = " ".join(abstract_parts)

    # Authors
    authors = []
    for author in art.findall(".//Author"):
        last = author.findtext("LastName", "")
        fore = author.findtext("ForeName", "")
        if last:
            authors.append(f"{last} {fore}".strip())

    # Journal
    journal_el = art.find(".//Journal")
    journal = ""
    pub_year = ""
    if journal_el is not None:
        journal = journal_el.findtext("Title", "") or journal_el.findtext("ISOAbbreviation", "")
        pub_year = journal_el.findtext(".//Year", "") or journal_el.findtext(".//MedlineDate", "")[:4]

    # DOI
    doi = ""
    for id_el in article.findall(".//ArticleId"):
        if id_el.get("IdType") == "doi":
            doi = id_el.text or ""
            break

    # MeSH terms
    mesh_terms = []
    for mesh in medline.findall(".//MeshHeading"):
        descriptor = mesh.findtext("DescriptorName", "")
        if descriptor:
            mesh_terms.append(descriptor)

    # Keywords
    keywords = [kw.text for kw in medline.findall(".//Keyword") if kw.text]

    return {
        "pmid": pmid,
        "title": title,
        "abstract": abstract[:1000],  # Truncate for DB storage
        "authors": authors[:10],
        "journal": journal,
        "year": pub_year,
        "doi": doi,
        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        "mesh_terms": mesh_terms[:20],
        "keywords": keywords[:15],
        "species": species,
        "source": "pubmed",
    }


def publications_to_gene_record(pubs: list[dict], species: str) -> dict:
    """
    Wrap publications into a gene record for storage.
    Use this when storing species-level publication summaries.
    """
    return {
        "gene_id": f"PUB_{species.replace(' ', '_').upper()}",
        "symbol": f"literature_{species.split()[0].lower()}",
        "organism": species,
        "sequence": "N",  # Placeholder — publications don't have sequences
        "sequence_type": "dna",
        "description": f"Publication metadata for {species}",
        "length": 1,
        "source": "pubmed",
        "publications": pubs,
        "annotations": {"publication_count": len(pubs)},
        "external_links": {},
        "traits": [],
        "expression_profiles": [],
        "pathways": [],
    }


if __name__ == "__main__":
    import json
    pubs = fetch_pubmed_for_species("Arabidopsis thaliana", retmax=5)
    print(json.dumps(pubs[:2], indent=2))
    print(f"Total pubs: {len(pubs)}")
