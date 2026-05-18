"""
alignment_engine.py
-------------------
Advanced sequence alignment algorithms for bioinformatics analysis.
Implements Needleman-Wunsch (global) and Smith-Waterman (local) alignment.
"""

import numpy as np
from typing import Tuple, Dict, List


# ─── Scoring matrices ──────────────────────────────────────────────────────────

BLOSUM62_MATRIX = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0, '*': -4},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1, '*': -4},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -2, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1, '*': -4},
}

DNA_MATRIX = {
    'A': {'A': 1, 'T': 0, 'G': 0, 'C': 0, 'N': 0},
    'T': {'A': 0, 'T': 1, 'G': 0, 'C': 0, 'N': 0},
    'G': {'A': 0, 'T': 0, 'G': 1, 'C': 0, 'N': 0},
    'C': {'A': 0, 'T': 0, 'G': 0, 'C': 1, 'N': 0},
    'N': {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 1},
}


def get_score(char1: str, char2: str, matrix: Dict, seq_type: str = "dna") -> int:
    """Get alignment score for two characters using a scoring matrix."""
    matrix_to_use = DNA_MATRIX if seq_type == "dna" else BLOSUM62_MATRIX
    return matrix_to_use.get(char1, {}).get(char2, -1)


# ─── NEEDLEMAN-WUNSCH (GLOBAL ALIGNMENT) ──────────────────────────────────────

def needleman_wunsch(seq1: str, seq2: str, gap_penalty: int = -2, seq_type: str = "dna") -> Dict:
    """
    Global sequence alignment using Needleman-Wunsch algorithm.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        gap_penalty: Penalty for gaps (default -2)
        seq_type: "dna" or "protein"
    
    Returns:
        Dict with aligned sequences, score, and alignment details
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize matrices
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)
    
    # Initialize first row and column
    for i in range(m + 1):
        score_matrix[i, 0] = i * gap_penalty
        traceback_matrix[i, 0] = 1  # move up
    for j in range(n + 1):
        score_matrix[0, j] = j * gap_penalty
        traceback_matrix[0, j] = 2  # move left
    
    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = get_score(seq1[i-1], seq2[j-1], None, seq_type)
            diagonal = score_matrix[i-1, j-1] + match_score
            up = score_matrix[i-1, j] + gap_penalty
            left = score_matrix[i, j-1] + gap_penalty
            
            score_matrix[i, j] = max(diagonal, up, left)
            
            if score_matrix[i, j] == diagonal:
                traceback_matrix[i, j] = 0  # diagonal
            elif score_matrix[i, j] == up:
                traceback_matrix[i, j] = 1  # up
            else:
                traceback_matrix[i, j] = 2  # left
    
    # Traceback to get alignment
    aligned_seq1, aligned_seq2 = _traceback(seq1, seq2, traceback_matrix, m, n, "nw")
    
    return {
        "algorithm": "Needleman-Wunsch (Global)",
        "seq1_aligned": aligned_seq1,
        "seq2_aligned": aligned_seq2,
        "alignment_score": float(score_matrix[m, n]),
        "match_count": sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-'),
        "gap_count": sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == '-' or b == '-'),
    }


# ─── SMITH-WATERMAN (LOCAL ALIGNMENT) ──────────────────────────────────────────

def smith_waterman(seq1: str, seq2: str, gap_penalty: int = -2, seq_type: str = "dna") -> Dict:
    """
    Local sequence alignment using Smith-Waterman algorithm.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        gap_penalty: Penalty for gaps (default -2)
        seq_type: "dna" or "protein"
    
    Returns:
        Dict with aligned sequences, score, and alignment details
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize matrices
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)
    
    max_score = 0
    max_pos = (0, 0)
    
    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = get_score(seq1[i-1], seq2[j-1], None, seq_type)
            diagonal = score_matrix[i-1, j-1] + match_score
            up = score_matrix[i-1, j] + gap_penalty
            left = score_matrix[i, j-1] + gap_penalty
            
            score_matrix[i, j] = max(0, diagonal, up, left)
            
            if score_matrix[i, j] == 0:
                traceback_matrix[i, j] = 3  # stop
            elif score_matrix[i, j] == diagonal:
                traceback_matrix[i, j] = 0  # diagonal
            elif score_matrix[i, j] == up:
                traceback_matrix[i, j] = 1  # up
            else:
                traceback_matrix[i, j] = 2  # left
            
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_pos = (i, j)
    
    # Traceback from max position
    aligned_seq1, aligned_seq2 = _traceback(seq1, seq2, traceback_matrix, max_pos[0], max_pos[1], "sw")
    
    return {
        "algorithm": "Smith-Waterman (Local)",
        "seq1_aligned": aligned_seq1,
        "seq2_aligned": aligned_seq2,
        "alignment_score": float(max_score),
        "match_count": sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-'),
        "gap_count": sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == '-' or b == '-'),
    }


def _traceback(seq1: str, seq2: str, traceback_matrix: np.ndarray, end_i: int, end_j: int, mode: str) -> Tuple[str, str]:
    """Traceback to reconstruct alignment."""
    aligned_seq1, aligned_seq2 = [], []
    i, j = end_i, end_j
    
    while i > 0 or j > 0:
        if mode == "sw" and traceback_matrix[i, j] == 3:
            break
        
        if traceback_matrix[i, j] == 0:  # diagonal
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == 1:  # up
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback_matrix[i, j] == 2:  # left
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
        else:
            break
    
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


# ─── MULTIPLE SEQUENCE ALIGNMENT (MSA) ────────────────────────────────────────

def progressive_alignment(sequences: List[str], seq_type: str = "dna") -> Dict:
    """
    Simple progressive multiple sequence alignment.
    
    Args:
        sequences: List of sequences to align
        seq_type: "dna" or "protein"
    
    Returns:
        Dict with aligned sequences and alignment statistics
    """
    if len(sequences) < 2:
        return {"error": "Need at least 2 sequences"}
    
    # Start with first two sequences
    alignment_result = needleman_wunsch(sequences[0], sequences[1], seq_type=seq_type)
    aligned_seqs = [alignment_result["seq1_aligned"], alignment_result["seq2_aligned"]]
    
    # Progressively align remaining sequences
    for seq in sequences[2:]:
        profile = "".join(aligned_seqs[0])  # Use first as profile
        alignment_result = needleman_wunsch(profile, seq, seq_type=seq_type)
        aligned_seqs = [alignment_result["seq1_aligned"], alignment_result["seq2_aligned"]]
    
    # Collect all aligned sequences
    msa_seqs = [alignment_result["seq1_aligned"]] + [s for s in aligned_seqs[1:]]
    
    return {
        "algorithm": "Progressive MSA",
        "aligned_sequences": msa_seqs,
        "num_sequences": len(sequences),
        "alignment_length": len(msa_seqs[0]) if msa_seqs else 0,
        "conservation_score": _calculate_conservation(msa_seqs),
    }


def _calculate_conservation(aligned_sequences: List[str]) -> float:
    """Calculate sequence conservation score (0-100%)."""
    if not aligned_sequences or len(aligned_sequences[0]) == 0:
        return 0.0
    
    num_seqs = len(aligned_sequences)
    align_length = len(aligned_sequences[0])
    conserved = 0
    
    for pos in range(align_length):
        chars = [seq[pos] for seq in aligned_sequences if pos < len(seq)]
        if len(set(chars)) == 1:  # All same
            conserved += 1
    
    return round(conserved / align_length * 100, 2)
