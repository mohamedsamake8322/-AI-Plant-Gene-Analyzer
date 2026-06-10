"""
similarityengine.py
-------------------
Database similarity search using global pairwise alignment (Needleman-Wunsch)
and percent identity over aligned non-gap columns.
"""

from __future__ import annotations

import json
import os
from typing import Optional

import alignment_engine as aln
import bioinformatics as bio
import sequence_loader


def _normalize_database(raw: object) -> dict:
    if isinstance(raw, dict):
        if "genes" in raw and isinstance(raw["genes"], list):
            entries = raw["genes"]
        elif all(isinstance(v, dict) for v in raw.values()):
            return raw
        else:
            entries = [raw] if raw.get("gene_id") or raw.get("symbol") else []
    elif isinstance(raw, list):
        entries = raw
    else:
        entries = []

    records: dict[str, dict] = {}
    for item in entries:
        if not isinstance(item, dict):
            continue
        key = item.get("gene_id") or item.get("symbol") or item.get("accession")
        if not key:
            continue
        records[key] = item
    return records


def load_gene_database(db_path: str = "genes_database.json") -> dict:
    """Load the gene database JSON file."""
    if isinstance(db_path, str) and (db_path == "postgres" or db_path.startswith("postgresql")):
        try:
            import sys
            from pathlib import Path

            script_root = Path(__file__).resolve().parent
            scripts_path = script_root / "scripts"
            if str(scripts_path) not in sys.path:
                sys.path.insert(0, str(scripts_path))
            from postgres_utils import load_gene_database_from_postgres

            return load_gene_database_from_postgres()
        except Exception:
            pass

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Gene database not found at: {db_path}")
    if db_path.lower().endswith((".fa", ".fasta")):
        return _load_database_from_fasta(db_path)

    with open(db_path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return _normalize_database(raw)


def _load_database_from_fasta(db_path: str) -> dict:
    records: dict[str, dict] = {}
    with open(db_path, "r", encoding="utf-8") as f:
        header_text = None
        header_metadata: dict[str, str] = {}
        sequence_parts: list[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header_text is not None:
                    seq = "".join(sequence_parts).upper()
                    entry_key = header_metadata.get("name", header_text)
                    records[entry_key] = {
                        "sequence": seq,
                        "trait": header_metadata.get("trait", "Unknown"),
                        "description": header_text,
                        "organism": header_metadata.get("organism", "Unknown"),
                        "accession": header_metadata.get("accession", "N/A"),
                        "symbol": header_metadata.get("name", entry_key),
                        "sequence_type": bio.detect_sequence_type(seq),
                    }
                header_text = line[1:].strip()
                header_metadata = sequence_loader._parse_header_metadata(header_text)
                sequence_parts = []
            else:
                sequence_parts.append(line)
        if header_text is not None:
            seq = "".join(sequence_parts).upper()
            entry_key = header_metadata.get("name", header_text)
            records[entry_key] = {
                "sequence": seq,
                "trait": header_metadata.get("trait", "Unknown"),
                "description": header_text,
                "organism": header_metadata.get("organism", "Unknown"),
                "accession": header_metadata.get("accession", "N/A"),
                "symbol": header_metadata.get("name", entry_key),
                "sequence_type": bio.detect_sequence_type(seq),
            }
    return records


def aligned_similarity(
    query: str,
    reference: str,
    query_type: str = "dna",
    reference_type: str = "dna",
) -> dict:
    """
    Compute global and local alignment identity between two sequences.
    For mixed DNA/protein pairs, tests all translation frames and keeps the best.
    """
    query = query.upper().replace(" ", "")
    reference = reference.upper().replace(" ", "")

    if query_type == "protein" and reference_type == "dna":
        best = {"similarity_score": 0.0, "alignment": {}, "frame": 0, "method": "protein_to_dna"}
        for frame in range(3):
            translated = bio.translate_dna(reference, frame)["protein"]
            global_aln = aln.needleman_wunsch(query, translated, seq_type="protein")
            if global_aln["identity_percent"] > best["similarity_score"]:
                best = {
                    "similarity_score": global_aln["identity_percent"],
                    "alignment": global_aln,
                    "frame": frame,
                    "method": "protein_to_dna_global",
                }
        return best

    if query_type == "dna" and reference_type == "protein":
        best = {"similarity_score": 0.0, "alignment": {}, "frame": 0, "method": "dna_to_protein"}
        for frame in range(3):
            translated = bio.translate_dna(query, frame)["protein"]
            global_aln = aln.needleman_wunsch(translated, reference, seq_type="protein")
            if global_aln["identity_percent"] > best["similarity_score"]:
                best = {
                    "similarity_score": global_aln["identity_percent"],
                    "alignment": global_aln,
                    "frame": frame,
                    "method": "dna_to_protein_global",
                }
        return best

    seq_type = "protein" if query_type == "protein" else "dna"
    global_aln = aln.needleman_wunsch(query, reference, seq_type=seq_type)
    local_aln = aln.smith_waterman(query, reference, seq_type=seq_type)
    return {
        "similarity_score": global_aln["identity_percent"],
        "global_identity": global_aln["identity_percent"],
        "local_identity": local_aln["identity_percent"],
        "alignment": global_aln,
        "method": "global",
        "local_alignment": local_aln,
    }


def compare_with_database(
    query: str,
    db_source: str | dict = "genes_database.json",
    top_n: int = 3,
) -> list[dict]:
    """Compare query against each gene using aligned percent identity."""
    database = db_source if isinstance(db_source, dict) else load_gene_database(db_source)
    query_type = bio.detect_sequence_type(query)
    results: list[dict] = []

    for gene_name, gene_info in database.items():
        ref_seq = gene_info["sequence"].upper().replace(" ", "")
        ref_type = gene_info.get("sequence_type") or bio.detect_sequence_type(ref_seq)
        match = aligned_similarity(query, ref_seq, query_type=query_type, reference_type=ref_type)
        alignment = match["alignment"]
        alignment_map = alignment.get("alignment_map", {}) if isinstance(alignment, dict) else {}

        results.append({
            "gene_name": gene_name,
            "trait": gene_info.get("trait", "Unknown"),
            "description": gene_info.get("description", ""),
            "organism": gene_info.get("organism", "Unknown"),
            "accession": gene_info.get("accession", "N/A"),
            "similarity_score": match["similarity_score"],
            "alignment_method": match.get("method", "global"),
            "alignment": {
                "alignment_map": alignment_map,
                "algorithm": alignment.get("algorithm") if isinstance(alignment, dict) else None,
                "alignment_score": alignment.get("alignment_score") if isinstance(alignment, dict) else None,
                "global_identity": match.get("global_identity"),
                "local_identity": match.get("local_identity"),
                "frame": match.get("frame"),
            },
            "reference_length": len(ref_seq),
            "query_length": len(query),
            "reference_type": ref_type,
            "query_type": query_type,
        })

    results.sort(key=lambda x: x["similarity_score"], reverse=True)
    return results[:top_n]


def get_best_match(query: str, db_path: str = "genes_database.json") -> Optional[dict]:
    matches = compare_with_database(query, db_path, top_n=1)
    return matches[0] if matches else None


def classify_similarity(score: float) -> dict[str, str]:
    if score >= 90:
        return {
            "level": "very_high",
            "label": "Very High Similarity",
            "color": "#00c853",
            "emoji": "🟢",
            "interpretation": "Near-identical aligned sequence — likely same gene or very close homolog.",
        }
    if score >= 75:
        return {
            "level": "high",
            "label": "High Similarity",
            "color": "#76ff03",
            "emoji": "🟩",
            "interpretation": "Strong homology after global alignment — likely functional equivalent.",
        }
    if score >= 55:
        return {
            "level": "moderate",
            "label": "Moderate Similarity",
            "color": "#ffd600",
            "emoji": "🟡",
            "interpretation": "Partial homology — may share conserved domains.",
        }
    if score >= 35:
        return {
            "level": "low",
            "label": "Low Similarity",
            "color": "#ff6d00",
            "emoji": "🟠",
            "interpretation": "Distant relationship — limited aligned identity.",
        }
    return {
        "level": "very_low",
        "label": "Very Low / No Significant Similarity",
        "color": "#d50000",
        "emoji": "🔴",
        "interpretation": "No significant aligned match in the reference database.",
    }
