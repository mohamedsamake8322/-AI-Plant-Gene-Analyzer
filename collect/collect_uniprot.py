#!/usr/bin/env python3
"""
UniProt collector — reviewed (Swiss-Prot) + unreviewed (TrEMBL) plant proteomes.
Returns a list of gene records compatible with the pipeline's canonical schema.
"""

from __future__ import annotations

import time
import requests
import request_utils as rq
from typing import Optional

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"

# Map common plant names to UniProt taxonomy IDs
TAXON_MAP: dict[str, str] = {
    "arabidopsis thaliana": "3702",
    "oryza sativa": "39947",
    "zea mays": "4577",
    "triticum aestivum": "4565",
    "glycine max": "3847",
    "solanum lycopersicum": "4081",
    "solanum tuberosum": "4113",
    "vitis vinifera": "29760",
    "hordeum vulgare": "4513",
    "sorghum bicolor": "4558",
    "medicago sativa": "3879",
    "phaseolus vulgaris": "3885",
    "helianthus annuus": "4232",
    "daucus carota": "4039",
    "lactuca sativa": "4236",
    "allium cepa": "35883",
    "brassica oleracea": "3712",
    "cucumis sativus": "3659",
    "malus domestica": "3750",
    "prunus persica": "3760",
    "citrus sinensis": "2711",
    "fragaria ananassa": "3747",
    "olea europaea": "4146",
}


def fetch_uniprot(
    species: str,
    retmax: int = 300,
    reviewed_only: bool = False,
) -> list[dict]:
    """
    Fetch protein records from UniProt for a given plant species.

    Args:
        species: Scientific name (e.g. "Arabidopsis thaliana")
        retmax: Maximum number of records to retrieve
        reviewed_only: If True, only fetch Swiss-Prot (reviewed) entries

    Returns:
        List of normalized gene records
    """
    taxon_id = TAXON_MAP.get(species.lower())
    if taxon_id:
        query = f"taxonomy_id:{taxon_id}"
    else:
        query = f'organism_name:"{species}"'

    if reviewed_only:
        query += " AND reviewed:true"

    params = {
        "query": query,
        "format": "json",
        "size": min(retmax, 500),
        "fields": (
            "accession,gene_names,organism_name,protein_name,"
            "sequence,go_id,go,ft_domain,cc_function,"
            "xref_kegg,xref_ensembl,xref_refseq,reviewed,"
            "protein_existence,annotation_score,length,"
            "cc_subcellular_location,keyword,feature_count"
        ),
    }

    records: list[dict] = []
    next_link: Optional[str] = None
    fetched = 0

    while fetched < retmax:
        try:
            if next_link:
                resp = rq.get(next_link, timeout=30)
            else:
                resp = rq.get(UNIPROT_API, params=params, timeout=30)
            data = resp.json()
            results = data.get("results", [])

            for entry in results:
                if fetched >= retmax:
                    break
                rec = _parse_entry(entry, species)
                if rec:
                    records.append(rec)
                    fetched += 1

            # Pagination via Link header
            link_header = resp.headers.get("Link", "")
            next_link = _parse_next_link(link_header)
            if not next_link or not results:
                break

            time.sleep(0.3)  # Be polite to UniProt API

        except requests.RequestException as e:
            print(f"  [UniProt] Request error for {species}: {e}")
            break

    return records


def _parse_entry(entry: dict, species: str) -> dict | None:
    """Parse a single UniProt JSON entry into pipeline canonical format."""
    accession = entry.get("primaryAccession", "")
    if not accession:
        return None

    # Gene names
    gene_names = entry.get("genes", [])
    symbol = ""
    if gene_names:
        gn = gene_names[0]
        symbol = (
            gn.get("geneName", {}).get("value", "")
            or (gn.get("synonyms", [{}])[0].get("value", "") if gn.get("synonyms") else "")
            or accession
        )

    # Protein name
    pn = entry.get("proteinDescription", {})
    rec_name = pn.get("recommendedName", {})
    description = (
        rec_name.get("fullName", {}).get("value", "")
        or pn.get("submittedName", [{}])[0].get("fullName", {}).get("value", "")
        if pn.get("submittedName") else ""
    )

    # Sequence
    seq_obj = entry.get("sequence", {})
    seq = seq_obj.get("value", "")
    if not seq:
        return None

    # Organism
    organism = entry.get("organism", {}).get("scientificName", species)

    # GO terms
    go_terms = []
    uniProtKB_cross_refs = entry.get("uniProtKBCrossReferences", [])
    for ref in uniProtKB_cross_refs:
        if ref.get("database") == "GO":
            go_id = ref.get("id", "")
            props = {p["key"]: p["value"] for p in ref.get("properties", [])}
            go_terms.append({
                "id": go_id,
                "term": props.get("GoTerm", ""),
                "aspect": props.get("GoEvidenceType", ""),
            })

    # KEGG cross-refs
    kegg_ids = [
        ref.get("id") for ref in uniProtKB_cross_refs
        if ref.get("database") == "KEGG"
    ]

    # Ensembl cross-refs
    ensembl_ids = [
        ref.get("id") for ref in uniProtKB_cross_refs
        if ref.get("database") == "Ensembl"
    ]

    # Keywords
    keywords = [kw.get("name", "") for kw in entry.get("keywords", [])]

    # Subcellular location
    comments = entry.get("comments", [])
    subcell = []
    func_desc = ""
    for comment in comments:
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                loc_val = loc.get("location", {}).get("value", "")
                if loc_val:
                    subcell.append(loc_val)
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                func_desc = texts[0].get("value", "")

    reviewed = entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"

    return {
        "gene_id": accession,
        "symbol": symbol or accession,
        "organism": organism,
        "sequence": seq.upper(),
        "sequence_type": "protein",
        "description": description or func_desc,
        "length": len(seq),
        "source": "uniprot",
        "annotations": {
            "go_terms": go_terms,
            "keywords": keywords,
            "subcellular_location": subcell,
            "function": func_desc,
            "annotation_score": entry.get("annotationScore", 0),
            "protein_existence": entry.get("proteinExistence", ""),
            "reviewed": reviewed,
        },
        "external_links": {
            "uniprot": f"https://www.uniprot.org/uniprotkb/{accession}",
            "accession": accession,
            "kegg": kegg_ids[0] if kegg_ids else None,
            "ensembl": ensembl_ids[0] if ensembl_ids else None,
        },
        "traits": keywords[:10],  # top keywords as traits
        "expression_profiles": [],
        "pathways": [{"id": k, "source": "kegg"} for k in kegg_ids],
        "publications": [],
    }


def _parse_next_link(link_header: str) -> str | None:
    """Extract 'next' URL from Link header."""
    if not link_header:
        return None
    for part in link_header.split(","):
        part = part.strip()
        if 'rel="next"' in part:
            url = part.split(";")[0].strip().strip("<>")
            return url
    return None


if __name__ == "__main__":
    import json
    results = fetch_uniprot("Arabidopsis thaliana", retmax=5)
    print(json.dumps(results[:2], indent=2))
    print(f"Total fetched: {len(results)}")
