"""
distance_engine.py
------------------
Evolutionary distance models computed on globally aligned sequence pairs.
"""

from __future__ import annotations

import math
from typing import Dict, List

import numpy as np

import alignment_engine as aln


def _pairwise_aligned_sequences(seq1: str, seq2: str, seq_type: str = "dna") -> tuple[str, str]:
    result = aln.needleman_wunsch(seq1, seq2, seq_type=seq_type)
    return result["seq1_aligned"], result["seq2_aligned"]


def hamming_distance(seq1: str, seq2: str, seq_type: str = "dna", prealigned: bool = False) -> Dict:
    if not prealigned:
        seq1, seq2 = _pairwise_aligned_sequences(seq1, seq2, seq_type=seq_type)

    differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != "-" and b != "-")
    total = sum(1 for a, b in zip(seq1, seq2) if a != "-" and b != "-")
    distance = differences / total if total > 0 else 0.0

    return {
        "algorithm": "Hamming Distance (aligned)",
        "differences": differences,
        "total_positions": total,
        "distance": round(distance, 6),
        "similarity": round(1 - distance, 6),
        "percent_identity": round((1 - distance) * 100, 2),
    }


def jukes_cantor_distance(seq1: str, seq2: str, seq_type: str = "dna", prealigned: bool = False) -> Dict:
    if not prealigned:
        seq1, seq2 = _pairwise_aligned_sequences(seq1, seq2, seq_type=seq_type)

    differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != "-" and b != "-")
    total = sum(1 for a, b in zip(seq1, seq2) if a != "-" and b != "-")
    if total == 0:
        return {"error": "No valid positions"}

    p = differences / total
    if p >= 0.75:
        return {
            "algorithm": "Jukes-Cantor (aligned)",
            "proportion_differences": round(p, 6),
            "jukes_cantor_distance": None,
            "warning": "Sequences too divergent (p >= 0.75).",
        }

    jc_distance = -0.75 * math.log(1 - (4 / 3) * p)
    return {
        "algorithm": "Jukes-Cantor (aligned)",
        "differences": differences,
        "total_positions": total,
        "proportion_differences": round(p, 6),
        "jukes_cantor_distance": round(jc_distance, 6),
    }


def kimura_distance(seq1: str, seq2: str, seq_type: str = "dna", prealigned: bool = False) -> Dict:
    if not prealigned:
        seq1, seq2 = _pairwise_aligned_sequences(seq1, seq2, seq_type=seq_type)

    transitions = transversions = 0
    non_gap = 0
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        non_gap += 1
        if a != b:
            if _is_transition(a, b):
                transitions += 1
            else:
                transversions += 1

    if non_gap == 0:
        return {"error": "No valid positions"}

    total = transitions + transversions
    p = transitions / non_gap
    q = transversions / non_gap
    if (1 - 2 * p - q) <= 0 or (1 - 2 * q) <= 0:
        return {
            "algorithm": "Kimura 2-Parameter (aligned)",
            "transitions": transitions,
            "transversions": transversions,
            "kimura_distance": None,
            "warning": "Sequences too divergent for K2P model.",
        }

    k2p_dist = -0.5 * math.log((1 - 2 * p - q) * math.sqrt(1 - 2 * q))
    return {
        "algorithm": "Kimura 2-Parameter (aligned)",
        "transitions": transitions,
        "transversions": transversions,
        "total_differences": total,
        "kimura_distance": round(k2p_dist, 6),
    }


def _is_transition(base1: str, base2: str) -> bool:
    return (base1, base2) in {("A", "G"), ("G", "A"), ("T", "C"), ("C", "T")}


def pam_distance(seq1: str, seq2: str, seq_type: str = "protein", prealigned: bool = False) -> Dict:
    if not prealigned:
        seq1, seq2 = _pairwise_aligned_sequences(seq1, seq2, seq_type="protein")

    differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != "-" and b != "-")
    total = sum(1 for a, b in zip(seq1, seq2) if a != "-" and b != "-")
    if total == 0:
        return {"error": "No valid positions"}

    proportion = differences / total
    if proportion > 0 and proportion < 0.99:
        pam = -100 * math.log(1 - proportion * 1.01)
    elif proportion == 0:
        pam = 0.0
    else:
        pam = 200.0

    return {
        "algorithm": "PAM Distance (aligned)",
        "differences": differences,
        "total_positions": total,
        "proportion_differences": round(proportion, 6),
        "pam_distance": round(pam, 2),
        "percent_identity": round((1 - proportion) * 100, 2),
    }


def distance_matrix(sequences: List[Dict], method: str = "kimura", seq_type: str = "dna") -> Dict:
    """
    Build a pairwise distance matrix.
    Sequences are star-aligned first so columns are comparable.
    """
    names = [seq["name"] for seq in sequences]
    raw_seqs = [seq["sequence"].upper().replace(" ", "") for seq in sequences]
    msa = aln.star_alignment(raw_seqs, seq_type=seq_type)
    aligned = msa.get("aligned_sequences", raw_seqs)

    n = len(aligned)
    matrix = np.zeros((n, n))

    key_by_method = {
        "hamming": ("hamming_distance", "distance"),
        "jukes_cantor": ("jukes_cantor_distance", "jukes_cantor_distance"),
        "kimura": ("kimura_distance", "kimura_distance"),
        "pam": ("pam_distance", "pam_distance"),
    }
    func_name, value_key = key_by_method.get(method.lower(), ("kimura_distance", "kimura_distance"))
    func = {
        "hamming_distance": hamming_distance,
        "jukes_cantor_distance": jukes_cantor_distance,
        "kimura_distance": kimura_distance,
        "pam_distance": pam_distance,
    }[func_name]
    stype = "protein" if method.lower() == "pam" else seq_type

    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i, j] = 0.0
            else:
                result = func(aligned[i], aligned[j], seq_type=stype, prealigned=True)
                value = result.get(value_key)
                matrix[i, j] = float(value) if value is not None else 0.0

    return {
        "method": method,
        "sequence_names": names,
        "distance_matrix": matrix.tolist(),
        "matrix_object": matrix,
        "alignment_method": msa.get("algorithm"),
        "aligned_sequences": aligned,
    }
