#!/usr/bin/env python3
"""
KEGG collector — metabolic pathways, orthologs (KO), and gene sequences for plant species.
Uses the KEGG REST API (free, no key required).
"""

from __future__ import annotations

import time
import requests
import request_utils as rq
from typing import Optional

KEGG_BASE = "https://rest.kegg.jp"

# KEGG organism codes for supported plant species
KEGG_ORG_MAP: dict[str, str] = {
    "arabidopsis thaliana": "ath",
    "oryza sativa": "osa",
    "zea mays": "zma",
    "triticum aestivum": "ta",
    "glycine max": "gmx",
    "solanum lycopersicum": "sly",
    "solanum tuberosum": "stu",
    "vitis vinifera": "vvi",
    "hordeum vulgare": "hvu",
    "sorghum bicolor": "sbi",
    "medicago sativa": "msa",
    "phaseolus vulgaris": "pvu",
    "helianthus annuus": "han",
    "daucus carota": "dcar",
    "brassica oleracea": "bol",
    "cucumis sativus": "csa",
    "malus domestica": "mdm",
    "prunus persica": "ppe",
    "citrus sinensis": "csi",
    "fragaria ananassa": "fan",
    "olea europaea": "oeu",
    "nicotiana tabacum": "nta",
    "chenopodium quinoa": "cqi",
}


def fetch_kegg(species: str, retmax: int = 300) -> list[dict]:
    """
    Fetch gene + pathway data from KEGG for a plant species.

    Args:
        species: Scientific name (e.g. "Oryza sativa")
        retmax: Maximum number of gene records to retrieve

    Returns:
        List of normalized gene records with pathway annotations
    """
    org_code = KEGG_ORG_MAP.get(species.lower())
    if not org_code:
        print(f"  [KEGG] No organism code for '{species}', skipping.")
        return []

    print(f"  [KEGG] Fetching gene list for {species} ({org_code})...")
    gene_ids = _get_gene_list(org_code, retmax)
    if not gene_ids:
        return []

    print(f"  [KEGG] Found {len(gene_ids)} genes, fetching details...")
    records = []
    for i, gene_id in enumerate(gene_ids[:retmax]):
        rec = _fetch_gene_entry(org_code, gene_id, species)
        if rec:
            records.append(rec)
        if (i + 1) % 50 == 0:
            print(f"    [KEGG] Processed {i+1}/{min(len(gene_ids), retmax)} genes...")
        time.sleep(0.2)  # KEGG rate limit: ~5 req/s

    return records


def _get_gene_list(org_code: str, limit: int) -> list[str]:
    """Retrieve list of gene IDs for an organism."""
    try:
        resp = rq.get(f"{KEGG_BASE}/list/{org_code}", timeout=30)
        # resp.raise_for_status() done in request_utils
        gene_ids = []
        for line in resp.text.strip().splitlines():
            parts = line.split("\t")
            if parts:
                gene_id = parts[0].split(":")[1] if ":" in parts[0] else parts[0]
                gene_ids.append(gene_id)
                if len(gene_ids) >= limit:
                    break
        return gene_ids
    except requests.RequestException as e:
        print(f"  [KEGG] Failed to get gene list: {e}")
        return []


def _fetch_gene_entry(org_code: str, gene_id: str, species: str) -> dict | None:
    """Fetch and parse a single KEGG gene entry."""
    try:
        resp = rq.get(f"{KEGG_BASE}/get/{org_code}:{gene_id}", timeout=30)
        return _parse_kegg_flat(resp.text, gene_id, org_code, species)
    except requests.RequestException as e:
        print(f"  [KEGG] Request failed for {org_code}:{gene_id}: {e}")
        return None
    except Exception as e:
        # A single malformed flat-file entry shouldn't abort the whole
        # species collection -- skip it and keep going.
        print(f"  [KEGG] Failed to parse {org_code}:{gene_id}: {e}")
        return None


def _parse_kegg_flat(text: str, gene_id: str, org_code: str, species: str) -> dict | None:
    """Parse KEGG flat-file format into pipeline schema."""
    fields: dict[str, list[str]] = {}
    current_key = ""

    for line in text.splitlines():
        if line.startswith("///"):
            break
        if line[:12].strip():
            current_key = line[:12].strip()
            fields.setdefault(current_key, [])
            fields[current_key].append(line[12:].strip())
        elif current_key:
            fields[current_key].append(line[12:].strip())

    # Extract fields
    name_lines = fields.get("NAME", [""])
    symbol = name_lines[0].split(",")[0].strip() if name_lines else gene_id
    description = " ".join(fields.get("DEFINITION", []))
    organism = " ".join(fields.get("ORGANISM", [species]))

    # Pathways
    pathways = []
    for line in fields.get("PATHWAY", []):
        parts = line.split(None, 1)
        if len(parts) == 2:
            pathways.append({"id": parts[0], "name": parts[1], "source": "kegg"})

    # Orthologs (KO)
    ko_ids = []
    for line in fields.get("ORTHOLOGY", []):
        parts = line.split(None, 1)
        if parts:
            ko_ids.append(parts[0])

    # Sequence (NTSEQ and/or AASEQ)
    # BUGFIX: KEGG gene entries very commonly carry BOTH the nucleotide
    # sequence (NTSEQ) and the protein sequence (AASEQ) in the same
    # flat-file record. The old code did
    #   fields.get("NTSEQ", fields.get("AASEQ", []))
    # which only ever reads ONE of the two blocks -- NTSEQ wins whenever
    # it's present, silently dropping the protein sequence for the vast
    # majority of entries (confirmed: 1000/1000 sampled records had
    # sequence.dna filled and sequence.protein empty). We now parse both
    # blocks independently.
    def _read_seq_block(lines: list[str]) -> str:
        # First line is the reported sequence length, not sequence data.
        return "".join(sl.replace(" ", "") for sl in lines[1:])

    nt_lines = fields.get("NTSEQ", [])
    aa_lines = fields.get("AASEQ", [])
    dna_sequence = _read_seq_block(nt_lines).upper() if nt_lines else ""
    protein_sequence = _read_seq_block(aa_lines).upper() if aa_lines else ""

    if not dna_sequence and not protein_sequence:
        # No NTSEQ/AASEQ in this KEGG entry -- still keep pathway/KO data
        # instead of discarding the whole record, same reasoning as the
        # PlantTFDB fix: metadata without a sequence is still useful.
        if not pathways and not ko_ids:
            return None

    # DBlinks
    external = {
        "kegg_gene": f"https://www.genome.jp/entry/{org_code}:{gene_id}",
        "accession": gene_id,
    }
    for line in fields.get("DBLINKS", []):
        parts = line.split(":", 1)
        if len(parts) == 2:
            db = parts[0].strip().lower().replace(" ", "_")
            val = parts[1].strip().split()[0]
            external[db] = val

    return {
        "gene_id": f"{org_code}:{gene_id}",
        "symbol": symbol,
        "organism": species,
        # BUGFIX: was a single ambiguous "sequence"/"sequence_type" pair.
        # Downstream (collect_all_sources.py KEGG block) now reads this
        # dict directly to fill sequence.dna AND sequence.protein.
        "sequences": {
            "dna": dna_sequence or None,
            "protein": protein_sequence or None,
        },
        "description": description,
        "length": len(dna_sequence) or len(protein_sequence),
        "source": "kegg",
        "pathways": pathways,
        "annotations": {
            "ko_ids": ko_ids,
            "kegg_org": org_code,
            "sequence_available": bool(dna_sequence or protein_sequence),
        },
        "external_links": external,
        # NOTE: this is a placeholder, not real phenotype data -- it just
        # echoes the first 5 pathway names so the field isn't empty. Real
        # trait data (verse/lodging, drought, etc.) needs the separate
        # manually-curated, PubMed-sourced trait table mentioned elsewhere
        # in this pipeline. Don't treat this as ground truth.
        "traits": [p["name"] for p in pathways[:5]],
        "expression_profiles": [],
        "publications": [],
    }


def fetch_pathway_genes(pathway_id: str, species: str = "") -> list[str]:
    """
    Fetch all gene IDs for a specific KEGG pathway.
    Useful for targeted pathway-based queries.
    Example: fetch_pathway_genes("map00010")  # Glycolysis
    """
    try:
        resp = rq.get(f"{KEGG_BASE}/link/genes/{pathway_id}", timeout=30)
        genes = []
        for line in resp.text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                genes.append(parts[1])
        return genes
    except requests.RequestException as e:
        print(f"  [KEGG] Failed to fetch pathway {pathway_id}: {e}")
        return []


if __name__ == "__main__":
    import json
    results = fetch_kegg("Arabidopsis thaliana", retmax=5)
    print(json.dumps(results[:2], indent=2))
    print(f"Total: {len(results)}")
