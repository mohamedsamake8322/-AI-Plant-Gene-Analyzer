"""
alignment_engine.py
-------------------
Sequence alignment: Needleman-Wunsch (global), Smith-Waterman (local),
star-based multiple sequence alignment, and alignment statistics.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

# Standard BLOSUM62 (Henikoff & Henikoff, 1992) — full 20×20 matrix
BLOSUM62_ALPHABET = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62_MATRIX: Dict[str, Dict[str, int]] = {
    "A": {"A": 4, "R": -1, "N": -2, "D": -2, "C": 0, "Q": -1, "E": -1, "G": 0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 0, "W": -3, "Y": -2, "V": 0},
    "R": {"A": -1, "R": 5, "N": 0, "D": -2, "C": -3, "Q": 1, "E": 0, "G": -2, "H": 0, "I": -3, "L": -2, "K": 2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
    "N": {"A": -2, "R": 0, "N": 6, "D": 1, "C": -3, "Q": 0, "E": 0, "G": 0, "H": 1, "I": -3, "L": -3, "K": 0, "M": -2, "F": -3, "P": -2, "S": 1, "T": 0, "W": -3, "Y": -2, "V": -3},
    "D": {"A": -2, "R": -2, "N": 1, "D": 6, "C": -3, "Q": 0, "E": 2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -3, "V": -3},
    "C": {"A": 0, "R": -3, "N": -3, "D": -3, "C": 9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
    "Q": {"A": -1, "R": 1, "N": 0, "D": 0, "C": -3, "Q": 5, "E": 2, "G": -2, "H": 0, "I": -3, "L": -2, "K": 1, "M": 0, "F": -3, "P": -1, "S": 0, "T": -1, "W": -2, "Y": -1, "V": -2},
    "E": {"A": -1, "R": 0, "N": 0, "D": 2, "C": -4, "Q": 2, "E": 5, "G": -2, "H": 0, "I": -3, "L": -3, "K": 1, "M": -2, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "G": {"A": 0, "R": -2, "N": 0, "D": -1, "C": -3, "Q": -2, "E": -2, "G": 6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S": 0, "T": -2, "W": -2, "Y": -3, "V": -3},
    "H": {"A": -2, "R": 0, "N": 1, "D": -1, "C": -3, "Q": 0, "E": 0, "G": -2, "H": 8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y": 2, "V": -3},
    "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I": 4, "L": 2, "K": -3, "M": 1, "F": 0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3},
    "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I": 2, "L": 4, "K": -2, "M": 2, "F": 0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1},
    "K": {"A": -1, "R": 2, "N": 0, "D": -1, "C": -3, "Q": 1, "E": 1, "G": -2, "H": -1, "I": -3, "L": -2, "K": 5, "M": -1, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q": 0, "E": -2, "G": -3, "H": -2, "I": 1, "L": 2, "K": -1, "M": 5, "F": 0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1},
    "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I": 0, "L": 0, "K": -3, "M": 0, "F": 6, "P": -4, "S": -2, "T": -2, "W": 1, "Y": 3, "V": -1},
    "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P": 7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"A": 1, "R": -1, "N": 1, "D": 0, "C": -1, "Q": 0, "E": 0, "G": 0, "H": -1, "I": -2, "L": -2, "K": 0, "M": -1, "F": -2, "P": -1, "S": 4, "T": 1, "W": -3, "Y": -2, "V": -2},
    "T": {"A": 0, "R": -1, "N": 0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 5, "W": -2, "Y": -2, "V": 0},
    "W": {"A": -3, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F": 1, "P": -4, "S": -3, "T": -2, "W": 11, "Y": 2, "V": -3},
    "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H": 2, "I": -1, "L": -1, "K": -2, "M": -1, "F": 3, "P": -3, "S": -2, "T": -2, "W": 2, "Y": 7, "V": -1},
    "V": {"A": 0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I": 3, "L": 1, "K": -2, "M": 1, "F": -1, "P": -2, "S": -2, "T": 0, "W": -3, "Y": -1, "V": 4},
}

DNA_MATRIX = {
    "A": {"A": 2, "T": -1, "G": -1, "C": -1, "N": 0},
    "T": {"A": -1, "T": 2, "G": -1, "C": -1, "N": 0},
    "G": {"A": -1, "T": -1, "G": 2, "C": -1, "N": 0},
    "C": {"A": -1, "T": -1, "G": -1, "C": 2, "N": 0},
    "N": {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0},
}


def get_score(char1: str, char2: str, seq_type: str = "dna") -> int:
    c1 = char1.upper()
    c2 = char2.upper()
    if seq_type == "dna":
        return DNA_MATRIX.get(c1, DNA_MATRIX["N"]).get(c2, 0)
    return BLOSUM62_MATRIX.get(c1, {}).get(c2, -4)


def alignment_statistics(aligned1: str, aligned2: str) -> Dict:
    """Compute identity and gap statistics from a gapped pairwise alignment."""
    matches = mismatches = gap_columns = 0
    alignment_length = len(aligned1)
    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            gap_columns += 1
            continue
        if a == b:
            matches += 1
        else:
            mismatches += 1
    # BLAST-style identity: identical columns / full alignment length (gaps reduce identity)
    identity = round(100.0 * matches / alignment_length, 2) if alignment_length else 0.0
    non_gap_cols = matches + mismatches
    non_gap_identity = round(100.0 * matches / non_gap_cols, 2) if non_gap_cols else 0.0
    return {
        "identity_percent": identity,
        "non_gap_identity_percent": non_gap_identity,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gap_columns,
        "aligned_length": alignment_length,
        "aligned_columns_without_gaps": non_gap_cols,
    }


def build_alignment_map(aligned1: str, aligned2: str) -> Dict[str, object]:
    stats = alignment_statistics(aligned1, aligned2)
    match_line = "".join(
        "|" if a == b and a != "-" else (" " if a == "-" or b == "-" else "X")
        for a, b in zip(aligned1, aligned2)
    )
    return {
        "query": aligned1,
        "reference": aligned2,
        "match_line": match_line,
        "identity_count": stats["matches"],
        "mismatch_count": stats["mismatches"],
        "gap_count": stats["gaps"],
        "identity_percent": stats["identity_percent"],
    }


def needleman_wunsch(seq1: str, seq2: str, gap_penalty: int = -2, seq_type: str = "dna") -> Dict:
    """Global alignment (Needleman-Wunsch)."""
    seq1 = seq1.upper().replace(" ", "")
    seq2 = seq2.upper().replace(" ", "")
    m, n = len(seq1), len(seq2)

    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)

    for i in range(1, m + 1):
        score_matrix[i, 0] = i * gap_penalty
        traceback_matrix[i, 0] = 1
    for j in range(1, n + 1):
        score_matrix[0, j] = j * gap_penalty
        traceback_matrix[0, j] = 2

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = get_score(seq1[i - 1], seq2[j - 1], seq_type)
            diagonal = score_matrix[i - 1, j - 1] + match_score
            up = score_matrix[i - 1, j] + gap_penalty
            left = score_matrix[i, j - 1] + gap_penalty
            score_matrix[i, j] = max(diagonal, up, left)
            if score_matrix[i, j] == diagonal:
                traceback_matrix[i, j] = 0
            elif score_matrix[i, j] == up:
                traceback_matrix[i, j] = 1
            else:
                traceback_matrix[i, j] = 2

    aligned_seq1, aligned_seq2 = _traceback(seq1, seq2, traceback_matrix, m, n, "nw")
    stats = alignment_statistics(aligned_seq1, aligned_seq2)
    return {
        "algorithm": "Needleman-Wunsch (Global)",
        "seq1_aligned": aligned_seq1,
        "seq2_aligned": aligned_seq2,
        "alignment_score": float(score_matrix[m, n]),
        "match_count": stats["matches"],
        "gap_count": stats["gaps"],
        "identity_percent": stats["identity_percent"],
        "alignment_map": build_alignment_map(aligned_seq1, aligned_seq2),
    }


def smith_waterman(seq1: str, seq2: str, gap_penalty: int = -2, seq_type: str = "dna") -> Dict:
    """Local alignment (Smith-Waterman)."""
    seq1 = seq1.upper().replace(" ", "")
    seq2 = seq2.upper().replace(" ", "")
    m, n = len(seq1), len(seq2)

    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = get_score(seq1[i - 1], seq2[j - 1], seq_type)
            diagonal = score_matrix[i - 1, j - 1] + match_score
            up = score_matrix[i - 1, j] + gap_penalty
            left = score_matrix[i, j - 1] + gap_penalty
            score_matrix[i, j] = max(0, diagonal, up, left)
            if score_matrix[i, j] == 0:
                traceback_matrix[i, j] = 3
            elif score_matrix[i, j] == diagonal:
                traceback_matrix[i, j] = 0
            elif score_matrix[i, j] == up:
                traceback_matrix[i, j] = 1
            else:
                traceback_matrix[i, j] = 2
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_pos = (i, j)

    aligned_seq1, aligned_seq2 = _traceback(seq1, seq2, traceback_matrix, max_pos[0], max_pos[1], "sw")
    stats = alignment_statistics(aligned_seq1, aligned_seq2)
    return {
        "algorithm": "Smith-Waterman (Local)",
        "seq1_aligned": aligned_seq1,
        "seq2_aligned": aligned_seq2,
        "alignment_score": float(max_score),
        "match_count": stats["matches"],
        "gap_count": stats["gaps"],
        "identity_percent": stats["identity_percent"],
        "alignment_map": build_alignment_map(aligned_seq1, aligned_seq2),
    }


def _traceback(
    seq1: str,
    seq2: str,
    traceback_matrix: np.ndarray,
    end_i: int,
    end_j: int,
    mode: str,
) -> Tuple[str, str]:
    aligned_seq1: List[str] = []
    aligned_seq2: List[str] = []
    i, j = end_i, end_j

    while i > 0 or j > 0:
        if mode == "sw" and traceback_matrix[i, j] == 3:
            break
        if i > 0 and j > 0 and traceback_matrix[i, j] == 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or traceback_matrix[i, j] == 1):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        elif j > 0:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        else:
            break

    return "".join(reversed(aligned_seq1)), "".join(reversed(aligned_seq2))


def star_alignment(sequences: List[str], seq_type: str = "dna") -> Dict:
    """
    Star alignment: align every sequence to the ungapped reference (first sequence),
    propagating gap columns across the full MSA.
    """
    if len(sequences) < 2:
        return {"error": "Need at least 2 sequences"}

    cleaned = [s.upper().replace(" ", "") for s in sequences]
    reference = cleaned[0]
    msa = [reference]

    for seq in cleaned[1:]:
        aln = needleman_wunsch(reference, seq, seq_type=seq_type)
        ref_aln = aln["seq1_aligned"]
        new_aln = aln["seq2_aligned"]

        expanded_rows: List[List[str]] = [[] for _ in msa]
        new_row: List[str] = []
        old_col = 0

        for rc, nc in zip(ref_aln, new_aln):
            if rc == "-":
                for row in expanded_rows:
                    row.append("-")
                new_row.append(nc if nc != "-" else "-")
            else:
                for r, row in enumerate(msa):
                    expanded_rows[r].append(row[old_col])
                new_row.append(nc)
                old_col += 1

        msa = ["".join(row) for row in expanded_rows] + ["".join(new_row)]

    width = len(msa[0])
    msa = [row.ljust(width, "-")[:width] for row in msa]
    return {
        "algorithm": "Star MSA (reference-guided)",
        "aligned_sequences": msa,
        "num_sequences": len(sequences),
        "alignment_length": width,
        "conservation_score": _calculate_conservation(msa),
        "reference_index": 0,
    }


def progressive_alignment(sequences: List[str], seq_type: str = "dna") -> Dict:
    """Backward-compatible alias for star alignment."""
    return star_alignment(sequences, seq_type=seq_type)


def _calculate_conservation(aligned_sequences: List[str]) -> float:
    if not aligned_sequences:
        return 0.0
    align_length = len(aligned_sequences[0])
    conserved = 0
    for pos in range(align_length):
        chars = {seq[pos] for seq in aligned_sequences if pos < len(seq) and seq[pos] != "-"}
        if len(chars) == 1:
            conserved += 1
    return round(conserved / align_length * 100, 2) if align_length else 0.0


def pairwise_align(
    seq1: str,
    seq2: str,
    mode: str = "global",
    seq_type: str = "dna",
) -> Dict:
    """Convenience wrapper returning alignment + identity statistics."""
    if mode == "local":
        return smith_waterman(seq1, seq2, seq_type=seq_type)
    return needleman_wunsch(seq1, seq2, seq_type=seq_type)
