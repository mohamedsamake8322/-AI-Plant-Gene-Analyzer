"""
distance_engine.py
------------------
Genetic distance calculation algorithms for phylogenetic analysis.
Implements Hamming, Jukes-Cantor, Kimura, and other evolutionary models.
"""

import numpy as np
from typing import Dict, List, Tuple
import math


# ─── HAMMING DISTANCE ─────────────────────────────────────────────────────────

def hamming_distance(seq1: str, seq2: str) -> Dict:
    """
    Calculate Hamming distance (number of differences).
    Works on aligned sequences of equal length.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
    
    Returns:
        Dict with distance metrics
    """
    if len(seq1) != len(seq2):
        return {"error": "Sequences must be same length"}
    
    differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')
    total = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    
    distance = differences / total if total > 0 else 0.0
    
    return {
        "algorithm": "Hamming Distance",
        "differences": differences,
        "total_positions": total,
        "distance": round(distance, 4),
        "similarity": round(1 - distance, 4),
        "percent_identity": round((1 - distance) * 100, 2),
    }


# ─── JUKES-CANTOR DISTANCE (DNA) ──────────────────────────────────────────────

def jukes_cantor_distance(seq1: str, seq2: str) -> Dict:
    """
    Jukes-Cantor distance for DNA sequences.
    Corrects for multiple substitutions assuming equal base frequencies.
    
    Formula: d = -3/4 * ln(1 - 4/3 * p)
    where p = proportion of differences
    
    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence
    
    Returns:
        Dict with distance metrics
    """
    if len(seq1) != len(seq2):
        return {"error": "Sequences must be same length"}
    
    differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')
    total = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    
    if total == 0:
        return {"error": "No valid positions"}
    
    p = differences / total  # proportion of differences
    
    # Check if p is valid for Jukes-Cantor
    if p >= 0.75:
        return {
            "algorithm": "Jukes-Cantor",
            "raw_distance": round(p, 4),
            "jukes_cantor_distance": None,
            "warning": "Sequences too divergent (p >= 0.75). Use raw distance.",
        }
    
    # Calculate JC distance
    jc_distance = -0.75 * math.log(1 - (4/3) * p) if p < 0.75 else None
    
    return {
        "algorithm": "Jukes-Cantor",
        "differences": differences,
        "total_positions": total,
        "proportion_differences": round(p, 4),
        "raw_distance": round(p, 4),
        "jukes_cantor_distance": round(jc_distance, 4) if jc_distance else None,
    }


# ─── KIMURA 2-PARAMETER (K2P) DISTANCE ────────────────────────────────────────

def kimura_distance(seq1: str, seq2: str) -> Dict:
    """
    Kimura 2-Parameter distance for DNA sequences.
    Distinguishes between transitions and transversions.
    
    Formula: d = -1/2 * ln((1 - 2P - Q) * sqrt(1 - 2Q))
    where P = transitions, Q = transversions
    
    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence
    
    Returns:
        Dict with distance metrics
    """
    if len(seq1) != len(seq2):
        return {"error": "Sequences must be same length"}
    
    transitions = 0
    transversions = 0
    total = 0
    
    for a, b in zip(seq1, seq2):
        if a == '-' or b == '-':
            continue
        if a != b:
            total += 1
            if _is_transition(a, b):
                transitions += 1
            else:
                transversions += 1
    
    if total == 0:
        return {"error": "No valid positions"}
    
    P = transitions / total
    Q = transversions / total
    
    # Check validity
    if (1 - 2*P - Q) <= 0 or (1 - 2*Q) <= 0:
        return {
            "algorithm": "Kimura 2-Parameter",
            "transitions": transitions,
            "transversions": transversions,
            "total_differences": total,
            "P": round(P, 4),
            "Q": round(Q, 4),
            "kimura_distance": None,
            "warning": "Sequences too divergent for K2P model",
        }
    
    k2p_dist = -0.5 * math.log((1 - 2*P - Q) * math.sqrt(1 - 2*Q))
    
    return {
        "algorithm": "Kimura 2-Parameter",
        "transitions": transitions,
        "transversions": transversions,
        "total_differences": total,
        "transition_ratio": round(transitions / total, 4) if total > 0 else 0,
        "transversion_ratio": round(transversions / total, 4) if total > 0 else 0,
        "P": round(P, 4),
        "Q": round(Q, 4),
        "kimura_distance": round(k2p_dist, 4),
    }


def _is_transition(base1: str, base2: str) -> bool:
    """Check if change is transition (purine↔purine or pyrimidine↔pyrimidine)."""
    transitions = {
        ('A', 'G'), ('G', 'A'),
        ('T', 'C'), ('C', 'T'),
    }
    return (base1, base2) in transitions


# ─── PAM DISTANCE (Protein) ────────────────────────────────────────────────────

def pam_distance(seq1: str, seq2: str) -> Dict:
    """
    Simple protein distance based on amino acid differences.
    
    Args:
        seq1: First protein sequence
        seq2: Second protein sequence
    
    Returns:
        Dict with distance metrics
    """
    if len(seq1) != len(seq2):
        return {"error": "Sequences must be same length"}
    
    differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')
    total = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    
    if total == 0:
        return {"error": "No valid positions"}
    
    distance = differences / total
    
    # Estimate PAM distance using Poisson correction
    if distance > 0:
        pam = -100 * math.log(1 - distance * 1.01) if distance < 0.99 else 200
    else:
        pam = 0
    
    return {
        "algorithm": "PAM Distance",
        "differences": differences,
        "total_positions": total,
        "proportion_differences": round(distance, 4),
        "pam_distance": round(pam, 2),
        "percent_identity": round((1 - distance) * 100, 2),
    }


# ─── PAIRWISE DISTANCE MATRIX ─────────────────────────────────────────────────

def distance_matrix(sequences: List[Dict], method: str = "kimura") -> Dict:
    """
    Calculate pairwise distance matrix for multiple sequences.
    
    Args:
        sequences: List of dicts with 'name' and 'sequence'
        method: "hamming", "jukes_cantor", "kimura", "pam"
    
    Returns:
        Dict with distance matrix, sequence names, and heatmap data
    """
    n = len(sequences)
    matrix = np.zeros((n, n))
    names = [seq["name"] for seq in sequences]
    
    distance_func = {
        "hamming": lambda s1, s2: hamming_distance(s1, s2)["distance"],
        "jukes_cantor": lambda s1, s2: jukes_cantor_distance(s1, s2).get("jukes_cantor_distance", 0),
        "kimura": lambda s1, s2: kimura_distance(s1, s2).get("kimura_distance", 0),
        "pam": lambda s1, s2: pam_distance(s1, s2).get("pam_distance", 0),
    }.get(method.lower(), hamming_distance)
    
    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i, j] = 0
            else:
                matrix[i, j] = distance_func(sequences[i]["sequence"], sequences[j]["sequence"])
    
    return {
        "method": method,
        "sequence_names": names,
        "distance_matrix": matrix.tolist(),
        "matrix_object": matrix,
    }
