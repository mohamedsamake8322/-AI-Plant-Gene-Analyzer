"""
similarityengine.py
-------------------
Similarity and alignment engine for the AI-Powered Plant Gene Analyzer.
Compares a query DNA sequence against the built-in gene database
using character-by-character alignment and a sliding-window approach.
"""

import json
import os
from typing import Optional

import bioinformatics as bio


# ─── Load gene database ────────────────────────────────────────────────────────

def load_gene_database(db_path: str = "genes_database.json") -> dict:
    """Load the gene database JSON file."""
    if not os.path.exists(db_path):
        raise FileNotFoundError(
            f"Gene database not found at: {db_path}"
        )
    if db_path.lower().endswith((".fa", ".fasta")):
        return _load_database_from_fasta(db_path)

    with open(db_path, "r", encoding="utf-8") as f:
        return json.load(f)


def _load_database_from_fasta(db_path: str) -> dict:
    records: dict[str, dict] = {}
    with open(db_path, "r", encoding="utf-8") as f:
        header = None
        sequence_parts: list[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(sequence_parts).upper()
                    records[header] = {
                        "sequence": seq,
                        "trait": "Unknown",
                        "description": "FASTA reference sequence",
                        "organism": "Unknown",
                        "accession": "N/A",
                        "sequence_type": bio.detect_sequence_type(seq),
                    }
                header = line[1:].split()[0] or f"sequence_{len(records)+1}"
                sequence_parts = []
            else:
                sequence_parts.append(line)
        if header is not None:
            seq = "".join(sequence_parts).upper()
            records[header] = {
                "sequence": seq,
                "trait": "Unknown",
                "description": "FASTA reference sequence",
                "organism": "Unknown",
                "accession": "N/A",
                "sequence_type": bio.detect_sequence_type(seq),
            }
    return records


# ─── Core alignment ────────────────────────────────────────────────────────────

def pairwise_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate the similarity percentage between two sequences
    by aligning them character by character over the shorter length.

    Returns a float in [0.0, 100.0].
    """
    if not seq1 or not seq2:
        return 0.0
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    return round(matches / min_len * 100, 2)


def sliding_window_similarity(
    query: str,
    reference: str,
    window_size: int = 30,
) -> dict:
    """
    Scan the query across the reference using a sliding window
    to find the best local alignment region.

    Args:
        query:       query DNA string
        reference:   reference DNA string
        window_size: length of the comparison window

    Returns:
        dict with best_score, best_position, best_window, alignment_map
    """
    if len(query) < window_size or len(reference) < window_size:
        score = pairwise_similarity(query, reference)
        return {
            "best_score": score,
            "best_position": 0,
            "best_window": reference[: len(query)],
            "alignment_map": _build_alignment_map(
                query[:window_size], reference[:window_size]
            ),
            "method": "direct",
        }

    best_score = 0.0
    best_pos = 0

    for i in range(len(reference) - window_size + 1):
        window = reference[i : i + window_size]
        score = pairwise_similarity(query[:window_size], window)
        if score > best_score:
            best_score = score
            best_pos = i

    best_window = reference[best_pos : best_pos + window_size]
    return {
        "best_score": best_score,
        "best_position": best_pos + 1,
        "best_window": best_window,
        "alignment_map": _build_alignment_map(query[:window_size], best_window),
        "method": "sliding_window",
    }


def _build_alignment_map(seq1: str, seq2: str) -> dict[str, str]:
    """
    Build a visual alignment representation with ASCII art.
    Returns a dict with 'query', 'match_line', 'reference', and ASCII alignment.
    """
    min_len = min(len(seq1), len(seq2))
    match_line = "".join(
        "|" if seq1[i] == seq2[i] else "X"
        for i in range(min_len)
    )
    
    # Build ASCII art alignment (formatted in blocks of 50)
    ascii_alignment = _build_ascii_alignment(seq1[:min_len], seq2[:min_len], match_line)
    
    return {
        "query": seq1[:min_len],
        "match_line": match_line,
        "reference": seq2[:min_len],
        "identity_count": match_line.count("|"),
        "mismatch_count": match_line.count("X"),
        "ascii_alignment": ascii_alignment,
    }


def _build_ascii_alignment(query: str, reference: str, match_line: str, block_size: int = 50) -> str:
    """
    Build a formatted ASCII alignment display, like BLAST output.

    Args:
        query: Query sequence
        reference: Reference sequence
        match_line: Line showing matches (|) and mismatches (X)
        block_size: Characters per line in the alignment

    Returns:
        Formatted ASCII alignment string
    """
    lines: list[str] = []
    for i in range(0, len(match_line), block_size):
        end = min(i + block_size, len(match_line))
        lines.append(f"\nQuery    {i + 1:6d}  {query[i:end]}  {end}")
        lines.append(f"         {'':6s}  {match_line[i:end]}")
        lines.append(f"Ref      {i + 1:6d}  {reference[i:end]}  {end}")
    return "\n".join(lines)


def _compare_protein_query_to_dna_reference(query_protein: str, reference_dna: str) -> dict:
    """Compare a protein query to a DNA reference by translation in all frames."""
    best_match = {
        "best_score": 0.0,
        "alignment": _build_alignment_map("", ""),
        "frame": 0,
    }
    for frame in range(3):
        translated = bio.translate_dna(reference_dna, frame)["protein"]
        score = pairwise_similarity(query_protein, translated)
        if score > best_match["best_score"]:
            best_match["best_score"] = score
            best_match["alignment"] = _build_alignment_map(query_protein, translated)
            best_match["frame"] = frame
    return best_match


def _compare_dna_query_to_protein_reference(query_dna: str, reference_protein: str) -> dict:
    """Compare a DNA query to a protein reference by translating the query."""
    best_match = {
        "best_score": 0.0,
        "alignment": _build_alignment_map("", ""),
        "frame": 0,
    }
    for frame in range(3):
        translated = bio.translate_dna(query_dna, frame)["protein"]
        score = pairwise_similarity(translated, reference_protein)
        if score > best_match["best_score"]:
            best_match["best_score"] = score
            best_match["alignment"] = _build_alignment_map(translated, reference_protein)
            best_match["frame"] = frame
    return best_match


# ─── Database comparison ───────────────────────────────────────────────────────

def compare_with_database(
    query: str,
    db_path: str = "genes_database.json",
    top_n: int = 3,
) -> list[dict]:
    """
    Compare the query sequence against every gene in the database.
    Returns the top_n best matches, sorted by similarity score.

    Each result dict contains:
        gene_name, trait, description, organism, accession,
        similarity_score, alignment, reference_length, query_length
    """
    database = load_gene_database(db_path)
    results: list[dict] = []

    query_type = bio.detect_sequence_type(query)

    for gene_name, gene_info in database.items():
        ref_seq = gene_info["sequence"].upper().replace(" ", "")
        ref_type = gene_info.get("sequence_type") or bio.detect_sequence_type(ref_seq)
        direct_score = 0.0
        best_score = 0.0
        alignment_payload: dict[str, object]

        if query_type == "protein" and ref_type == "dna":
            protein_result = _compare_protein_query_to_dna_reference(query, ref_seq)
            alignment_payload = {
                "alignment_map": protein_result["alignment"],
                "best_score": protein_result["best_score"],
                "method": "protein_to_dna",
                "frame": protein_result["frame"],
            }
            direct_score = protein_result["best_score"]
            best_score = protein_result["best_score"]
        elif query_type == "dna" and ref_type == "protein":
            protein_result = _compare_dna_query_to_protein_reference(query, ref_seq)
            alignment_payload = {
                "alignment_map": protein_result["alignment"],
                "best_score": protein_result["best_score"],
                "method": "dna_to_protein",
                "frame": protein_result["frame"],
            }
            direct_score = protein_result["best_score"]
            best_score = protein_result["best_score"]
        else:
            alignment_payload = sliding_window_similarity(query, ref_seq)
            direct_score = pairwise_similarity(query, ref_seq)
            best_score = max(alignment_payload["best_score"], direct_score)

        if isinstance(alignment_payload, dict) and "alignment_map" in alignment_payload:
            alignment_payload["ascii_alignment"] = alignment_payload["alignment_map"].get("ascii_alignment")

        results.append({
            "gene_name": gene_name,
            "trait": gene_info.get("trait", "Unknown"),
            "description": gene_info.get("description", ""),
            "organism": gene_info.get("organism", "Unknown"),
            "accession": gene_info.get("accession", "N/A"),
            "similarity_score": best_score,
            "direct_similarity": direct_score,
            "sliding_window_similarity": alignment_payload["best_score"],
            "alignment": alignment_payload,
            "reference_length": len(ref_seq),
            "query_length": len(query),
            "reference_type": ref_type,
            "query_type": query_type,
        })

    results.sort(key=lambda x: x["similarity_score"], reverse=True)
    return results[:top_n]


def get_best_match(
    query: str,
    db_path: str = "genes_database.json",
) -> Optional[dict]:
    """Return the single best matching gene from the database."""
    matches = compare_with_database(query, db_path, top_n=1)
    return matches[0] if matches else None


# ─── Similarity classification ─────────────────────────────────────────────────

def classify_similarity(score: float) -> dict[str, str]:
    """
    Map a similarity score to a biological classification label.

    Returns dict with 'level', 'label', 'color', and 'emoji'.
    """
    if score >= 90:
        return {
            "level": "very_high",
            "label": "Very High Similarity",
            "color": "#00c853",
            "emoji": "🟢",
            "interpretation": "Near-identical sequence — likely same gene or very close homolog.",
        }
    elif score >= 75:
        return {
            "level": "high",
            "label": "High Similarity",
            "color": "#76ff03",
            "emoji": "🟩",
            "interpretation": "Strong homology — likely functional equivalent across species.",
        }
    elif score >= 55:
        return {
            "level": "moderate",
            "label": "Moderate Similarity",
            "color": "#ffd600",
            "emoji": "🟡",
            "interpretation": "Partial homology — may share functional domains.",
        }
    elif score >= 35:
        return {
            "level": "low",
            "label": "Low Similarity",
            "color": "#ff6d00",
            "emoji": "🟠",
            "interpretation": "Distant relationship — possible convergent function.",
        }
    else:
        return {
            "level": "very_low",
            "label": "Very Low / No Significant Similarity",
            "color": "#d50000",
            "emoji": "🔴",
            "interpretation": "No significant match — likely novel or unrelated sequence.",
        }
