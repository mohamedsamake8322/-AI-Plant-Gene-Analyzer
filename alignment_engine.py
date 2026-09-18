"""
alignment_engine.py
-------------------
Sequence alignment: Needleman-Wunsch (global), Smith-Waterman (local),
star-based multiple sequence alignment, and alignment statistics.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, List, Tuple

import numpy as np

import config

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised only without optional acceleration
    njit = None

# "Negative infinity" sentinel for affine-gap DP. Must be far below any
# reachable real score (max ~15,000 * 11 for a full BLOSUM62 W-W run at the
# MAX_ALIGNMENT_SEQUENCE_LENGTH cap) but far enough from int32's limits
# (~2.1e9) that repeated gap_extend additions during traceback-adjacent
# arithmetic can't overflow.
NEG_INF = -10**7

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


def get_score(char1: str, char2: str, seq_type: str = "dna", matrix_name: str = "default") -> int:
    c1 = char1.upper()
    c2 = char2.upper()
    if seq_type == "dna":
        return DNA_MATRIX.get(c1, DNA_MATRIX["N"]).get(c2, 0)
    if matrix_name != "default":
        return _protein_matrix_score(c1, c2, matrix_name)
    return BLOSUM62_MATRIX.get(c1, {}).get(c2, -4)


def _protein_matrix_score(char1: str, char2: str, matrix_name: str) -> int:
    from Bio.Align import substitution_matrices

    matrix = substitution_matrices.load(matrix_name)
    try:
        return int(matrix[char1, char2])
    except (KeyError, IndexError):
        return -4


@lru_cache(maxsize=16)
def _score_lookup_table(seq_type: str, matrix_name: str = "default") -> tuple[dict[str, int], int, np.ndarray]:
    """Precompute a dense numeric substitution matrix + char->index map.

    Replaces get_score()'s per-cell dict-of-dict lookups (two hash lookups
    plus a Python function call, executed up to m*n times in the DP inner
    loop) with O(1) numpy array indexing. Built once per seq_type and
    cached — the matrices never change at runtime. The extra row/column
    (index = len(alphabet)) is the fallback for any character outside the
    known alphabet (e.g. an IUPAC ambiguity code); its score against
    everything is set to match get_score()'s existing fallback exactly:
    0 for DNA (same as falling back to "N", whose row is all zeros) and -4
    for protein (get_score's fixed fallback penalty).
    """
    if seq_type == "dna":
        alphabet, matrix, fallback = "ATGCN", DNA_MATRIX, 0
    else:
        if matrix_name == "default":
            alphabet, matrix, fallback = BLOSUM62_ALPHABET, BLOSUM62_MATRIX, -4
        else:
            from Bio.Align import substitution_matrices

            loaded = substitution_matrices.load(matrix_name)
            alphabet = "".join(str(letter) for letter in loaded.alphabet if str(letter).isalpha())
            matrix = {a: {b: int(loaded[a, b]) for b in alphabet} for a in alphabet}
            fallback = -4

    index = {ch: i for i, ch in enumerate(alphabet)}
    unknown_idx = len(alphabet)
    size = len(alphabet) + 1
    table = np.full((size, size), fallback, dtype=np.int32)
    for a in alphabet:
        for b in alphabet:
            table[index[a], index[b]] = matrix.get(a, {}).get(b, fallback)
    return index, unknown_idx, table


def _encode_sequence(sequence: str, index: dict[str, int], unknown_idx: int) -> np.ndarray:
    """Map a sequence string to an int array of lookup-table indices."""
    return np.array([index.get(ch, unknown_idx) for ch in sequence], dtype=np.int64)


if njit is not None:
    @njit(cache=True)
    def _needleman_wunsch_dp_numba(enc1, enc2, score_table, gap_open, gap_extend, neg_inf):
        m, n = len(enc1), len(enc2)
        M = np.full((m + 1, n + 1), neg_inf, dtype=np.int32)
        Ix = np.full((m + 1, n + 1), neg_inf, dtype=np.int32)
        Iy = np.full((m + 1, n + 1), neg_inf, dtype=np.int32)
        tb_M = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Ix = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Iy = np.zeros((m + 1, n + 1), dtype=np.int8)
        M[0, 0] = 0

        for i in range(1, m + 1):
            open_score = M[i - 1, 0] + gap_open
            extend_score = Ix[i - 1, 0] + gap_extend
            if open_score >= extend_score:
                Ix[i, 0], tb_Ix[i, 0] = open_score, 0
            else:
                Ix[i, 0], tb_Ix[i, 0] = extend_score, 1
        for j in range(1, n + 1):
            open_score = M[0, j - 1] + gap_open
            extend_score = Iy[0, j - 1] + gap_extend
            if open_score >= extend_score:
                Iy[0, j], tb_Iy[0, j] = open_score, 0
            else:
                Iy[0, j], tb_Iy[0, j] = extend_score, 1

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                best_prev = M[i - 1, j - 1]
                state = 0
                if Ix[i - 1, j - 1] > best_prev:
                    best_prev = Ix[i - 1, j - 1]
                    state = 1
                if Iy[i - 1, j - 1] > best_prev:
                    best_prev = Iy[i - 1, j - 1]
                    state = 2
                M[i, j] = score_table[enc1[i - 1], enc2[j - 1]] + best_prev
                tb_M[i, j] = state

                open_score = M[i - 1, j] + gap_open
                extend_score = Ix[i - 1, j] + gap_extend
                if open_score >= extend_score:
                    Ix[i, j], tb_Ix[i, j] = open_score, 0
                else:
                    Ix[i, j], tb_Ix[i, j] = extend_score, 1

                open_score = M[i, j - 1] + gap_open
                extend_score = Iy[i, j - 1] + gap_extend
                if open_score >= extend_score:
                    Iy[i, j], tb_Iy[i, j] = open_score, 0
                else:
                    Iy[i, j], tb_Iy[i, j] = extend_score, 1
        return M, Ix, Iy, tb_M, tb_Ix, tb_Iy


    @njit(cache=True)
    def _smith_waterman_dp_numba(enc1, enc2, score_table, gap_penalty):
        m, n = len(enc1), len(enc2)
        score_matrix = np.zeros((m + 1, n + 1), dtype=np.int32)
        traceback_matrix = np.zeros((m + 1, n + 1), dtype=np.int8)
        max_score = 0
        max_i, max_j = 0, 0
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                diagonal = score_matrix[i - 1, j - 1] + score_table[enc1[i - 1], enc2[j - 1]]
                up = score_matrix[i - 1, j] + gap_penalty
                left = score_matrix[i, j - 1] + gap_penalty
                cell = max(0, diagonal, up, left)
                score_matrix[i, j] = cell
                if cell == 0:
                    traceback_matrix[i, j] = 3
                elif cell == diagonal:
                    traceback_matrix[i, j] = 0
                elif cell == up:
                    traceback_matrix[i, j] = 1
                else:
                    traceback_matrix[i, j] = 2
                if cell > max_score:
                    max_score = cell
                    max_i, max_j = i, j
        return score_matrix, traceback_matrix, max_score, max_i, max_j


    @njit(cache=True)
    def _smith_waterman_affine_dp_numba(enc1, enc2, score_table, gap_open, gap_extend):
        m, n = len(enc1), len(enc2)
        M = np.zeros((m + 1, n + 1), dtype=np.int32)
        Ix = np.zeros((m + 1, n + 1), dtype=np.int32)
        Iy = np.zeros((m + 1, n + 1), dtype=np.int32)
        tb_M = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Ix = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Iy = np.zeros((m + 1, n + 1), dtype=np.int8)
        max_score = 0
        max_i, max_j, max_state = 0, 0, 0

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match_score = score_table[enc1[i - 1], enc2[j - 1]]
                diag_m = M[i - 1, j - 1] + match_score
                diag_ix = Ix[i - 1, j - 1] + match_score
                diag_iy = Iy[i - 1, j - 1] + match_score
                best_m = max(0, diag_m, diag_ix, diag_iy)
                M[i, j] = best_m
                if best_m == 0:
                    tb_M[i, j] = 3
                elif best_m == diag_m:
                    tb_M[i, j] = 0
                elif best_m == diag_ix:
                    tb_M[i, j] = 1
                else:
                    tb_M[i, j] = 2

                open_ix = M[i - 1, j] + gap_open
                extend_ix = Ix[i - 1, j] + gap_extend
                Ix[i, j] = max(0, open_ix, extend_ix)
                if Ix[i, j] == 0:
                    tb_Ix[i, j] = 3
                elif Ix[i, j] == open_ix:
                    tb_Ix[i, j] = 0
                else:
                    tb_Ix[i, j] = 1

                open_iy = M[i, j - 1] + gap_open
                extend_iy = Iy[i, j - 1] + gap_extend
                Iy[i, j] = max(0, open_iy, extend_iy)
                if Iy[i, j] == 0:
                    tb_Iy[i, j] = 3
                elif Iy[i, j] == open_iy:
                    tb_Iy[i, j] = 0
                else:
                    tb_Iy[i, j] = 1

                if best_m > max_score:
                    max_score, max_i, max_j, max_state = best_m, i, j, 0
                if Ix[i, j] > max_score:
                    max_score, max_i, max_j, max_state = Ix[i, j], i, j, 1
                if Iy[i, j] > max_score:
                    max_score, max_i, max_j, max_state = Iy[i, j], i, j, 2

        return M, Ix, Iy, tb_M, tb_Ix, tb_Iy, max_score, max_i, max_j, max_state


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


def needleman_wunsch(
    seq1: str,
    seq2: str,
    gap_penalty: int | None = None,
    seq_type: str = "dna",
    gap_open: int = -10,
    gap_extend: int = -1,
    matrix_name: str = "default",
) -> Dict:
    """Global alignment (Needleman-Wunsch) with affine gap penalties (Gotoh's algorithm).

    Affine gaps score opening a gap (gap_open) separately from extending an
    already-open one (gap_extend, applied to each additional consecutive
    gap position). This matches real indel biology far better than a flat
    per-position penalty: one long insertion/deletion is much more likely
    than many scattered single-base gaps, but a linear penalty scores both
    identically, which fragments alignments unrealistically.

    gap_penalty is kept for backward compatibility: if given explicitly, it
    is used as BOTH gap_open and gap_extend, reproducing this function's
    previous linear-gap behavior exactly (verified against it — see the
    test suite).

    Memory note: uses six (m+1)x(n+1) matrices (three DP score matrices +
    three traceback matrices) — O(n*m) space, not just O(n*m) time. dtype
    is pinned (int32 for scores, int8 for traceback codes) rather than
    relying on the platform-dependent numpy default. Callers should cap
    input length well before memory becomes an issue — see
    config.MAX_ALIGNMENT_SEQUENCE_LENGTH, enforced below.
    """
    if gap_penalty is not None:
        gap_open = gap_extend = gap_penalty

    seq1 = seq1.upper().replace(" ", "")
    seq2 = seq2.upper().replace(" ", "")
    m, n = len(seq1), len(seq2)

    if max(m, n) > config.MAX_ALIGNMENT_SEQUENCE_LENGTH:
        raise ValueError(
            f"Sequence too long for pairwise alignment ({max(m, n):,} > "
            f"{config.MAX_ALIGNMENT_SEQUENCE_LENGTH:,}). This guard applies "
            "regardless of caller, so any direct use of this function (not "
            "just the main analysis pipeline) is protected from excessive "
            "memory/time use."
        )

    char_index, unknown_idx, score_table = _score_lookup_table(seq_type, matrix_name)
    enc1 = _encode_sequence(seq1, char_index, unknown_idx)
    enc2 = _encode_sequence(seq2, char_index, unknown_idx)

    if njit is not None:
        M, Ix, Iy, tb_M, tb_Ix, tb_Iy = _needleman_wunsch_dp_numba(
            enc1, enc2, score_table, gap_open, gap_extend, NEG_INF
        )
    else:
        M = np.full((m + 1, n + 1), NEG_INF, dtype=np.int32)
        Ix = np.full((m + 1, n + 1), NEG_INF, dtype=np.int32)
        Iy = np.full((m + 1, n + 1), NEG_INF, dtype=np.int32)
        tb_M = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Ix = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Iy = np.zeros((m + 1, n + 1), dtype=np.int8)
        M[0, 0] = 0
        for i in range(1, m + 1):
            open_ = M[i - 1, 0] + gap_open
            extend_ = Ix[i - 1, 0] + gap_extend
            if open_ >= extend_:
                Ix[i, 0], tb_Ix[i, 0] = open_, 0
            else:
                Ix[i, 0], tb_Ix[i, 0] = extend_, 1
        for j in range(1, n + 1):
            open_ = M[0, j - 1] + gap_open
            extend_ = Iy[0, j - 1] + gap_extend
            if open_ >= extend_:
                Iy[0, j], tb_Iy[0, j] = open_, 0
            else:
                Iy[0, j], tb_Iy[0, j] = extend_, 1
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                s = int(score_table[enc1[i - 1], enc2[j - 1]])
                prevs = (M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1])
                best_prev = max(prevs)
                M[i, j] = s + best_prev
                tb_M[i, j] = prevs.index(best_prev)
                open_ = M[i - 1, j] + gap_open
                extend_ = Ix[i - 1, j] + gap_extend
                if open_ >= extend_:
                    Ix[i, j], tb_Ix[i, j] = open_, 0
                else:
                    Ix[i, j], tb_Ix[i, j] = extend_, 1
                open_ = M[i, j - 1] + gap_open
                extend_ = Iy[i, j - 1] + gap_extend
                if open_ >= extend_:
                    Iy[i, j], tb_Iy[i, j] = open_, 0
                else:
                    Iy[i, j], tb_Iy[i, j] = extend_, 1

    finals = {"M": M[m, n], "Ix": Ix[m, n], "Iy": Iy[m, n]}
    end_state = max(finals, key=finals.get)
    final_score = int(finals[end_state])

    aligned_seq1, aligned_seq2 = _affine_traceback(seq1, seq2, tb_M, tb_Ix, tb_Iy, m, n, end_state)
    stats = alignment_statistics(aligned_seq1, aligned_seq2)
    return {
        "algorithm": "Needleman-Wunsch (Global, affine gap)",
        "seq1_aligned": aligned_seq1,
        "seq2_aligned": aligned_seq2,
        "alignment_score": float(final_score),
        "match_count": stats["matches"],
        "gap_count": stats["gaps"],
        "identity_percent": stats["identity_percent"],
        "alignment_map": build_alignment_map(aligned_seq1, aligned_seq2),
    }


def _affine_traceback(
    seq1: str,
    seq2: str,
    tb_M: np.ndarray,
    tb_Ix: np.ndarray,
    tb_Iy: np.ndarray,
    end_i: int,
    end_j: int,
    end_state: str,
) -> Tuple[str, str]:
    aligned1: List[str] = []
    aligned2: List[str] = []
    i, j, state = end_i, end_j, end_state

    while i > 0 or j > 0:
        if state == "M":
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            state = ("M", "Ix", "Iy")[tb_M[i, j]]
            i -= 1
            j -= 1
        elif state == "Ix":
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            state = "M" if tb_Ix[i, j] == 0 else "Ix"
            i -= 1
        else:  # Iy
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            state = "M" if tb_Iy[i, j] == 0 else "Iy"
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2))


def smith_waterman(
    seq1: str,
    seq2: str,
    gap_penalty: int | None = None,
    seq_type: str = "dna",
    gap_open: int = -10,
    gap_extend: int = -1,
    matrix_name: str = "default",
) -> Dict:
    """Local alignment (Smith-Waterman) with affine gap penalties by default.

    `gap_penalty` is retained for backward compatibility: when provided, it is
    applied as both the opening and extension penalty, matching the historical
    linear-gap implementation. New call sites should prefer the explicit
    `gap_open` / `gap_extend` pair, which is both more biologically realistic and
    consistent with the affine-gap model already used by `needleman_wunsch()`.
    """
    if gap_penalty is not None:
        gap_open = gap_extend = gap_penalty

    seq1 = seq1.upper().replace(" ", "")
    seq2 = seq2.upper().replace(" ", "")
    m, n = len(seq1), len(seq2)

    if max(m, n) > config.MAX_ALIGNMENT_SEQUENCE_LENGTH:
        raise ValueError(
            f"Sequence too long for pairwise alignment ({max(m, n):,} > "
            f"{config.MAX_ALIGNMENT_SEQUENCE_LENGTH:,}). This guard applies "
            "regardless of caller, so any direct use of this function (not "
            "just the main analysis pipeline) is protected from excessive "
            "memory/time use."
        )

    char_index, unknown_idx, score_table = _score_lookup_table(seq_type, matrix_name)
    enc1 = _encode_sequence(seq1, char_index, unknown_idx)
    enc2 = _encode_sequence(seq2, char_index, unknown_idx)

    if njit is not None:
        if gap_penalty is not None:
            score_matrix, traceback_matrix, max_score, max_i, max_j = _smith_waterman_dp_numba(
                enc1, enc2, score_table, gap_penalty
            )
            max_pos = (max_i, max_j)
        else:
            score_matrix, Ix, Iy, tb_M, tb_Ix, tb_Iy, max_score, max_i, max_j, max_state = _smith_waterman_affine_dp_numba(
                enc1, enc2, score_table, gap_open, gap_extend
            )
            max_pos = (max_i, max_j)
    else:
        score_matrix = np.zeros((m + 1, n + 1), dtype=np.int32)
        Ix = np.zeros((m + 1, n + 1), dtype=np.int32)
        Iy = np.zeros((m + 1, n + 1), dtype=np.int32)
        tb_M = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Ix = np.zeros((m + 1, n + 1), dtype=np.int8)
        tb_Iy = np.zeros((m + 1, n + 1), dtype=np.int8)
        max_score = 0
        max_pos = (0, 0)
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match_score = int(score_table[enc1[i - 1], enc2[j - 1]])
                diag_M = score_matrix[i - 1, j - 1] + match_score
                diag_Ix = Ix[i - 1, j - 1] + match_score
                diag_Iy = Iy[i - 1, j - 1] + match_score
                score_matrix[i, j] = max(0, diag_M, diag_Ix, diag_Iy)
                if score_matrix[i, j] == 0:
                    tb_M[i, j] = 3
                elif score_matrix[i, j] == diag_M:
                    tb_M[i, j] = 0
                elif score_matrix[i, j] == diag_Ix:
                    tb_M[i, j] = 1
                else:
                    tb_M[i, j] = 2

                from_m = score_matrix[i - 1, j] + gap_open
                from_i = Ix[i - 1, j] + gap_extend
                Ix[i, j] = max(0, from_m, from_i)
                tb_Ix[i, j] = 0 if Ix[i, j] == from_m else 1 if Ix[i, j] > 0 else 3

                from_m = score_matrix[i, j - 1] + gap_open
                from_i = Iy[i, j - 1] + gap_extend
                Iy[i, j] = max(0, from_m, from_i)
                tb_Iy[i, j] = 0 if Iy[i, j] == from_m else 1 if Iy[i, j] > 0 else 3

                if score_matrix[i, j] > max_score:
                    max_score = score_matrix[i, j]
                    max_pos = (i, j)

        traceback_matrix = tb_M

    if gap_penalty is not None:
        aligned_seq1, aligned_seq2 = _traceback(seq1, seq2, traceback_matrix, max_pos[0], max_pos[1], "sw")
    else:
        i, j = max_pos[0], max_pos[1]
        aligned_seq1: List[str] = []
        aligned_seq2: List[str] = []
        state = "M"
        while i > 0 and j > 0 and score_matrix[i, j] > 0:
            if state == "M":
                if tb_M[i, j] == 0:
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append(seq2[j - 1])
                    i -= 1
                    j -= 1
                elif tb_M[i, j] == 1:
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append(seq2[j - 1])
                    i -= 1
                    j -= 1
                    state = "Ix"
                else:
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append(seq2[j - 1])
                    i -= 1
                    j -= 1
                    state = "Iy"
            elif state == "Ix":
                if tb_Ix[i, j] == 0:
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append("-")
                    i -= 1
                    state = "M"
                else:
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append("-")
                    i -= 1
            else:
                if tb_Iy[i, j] == 0:
                    aligned_seq1.append("-")
                    aligned_seq2.append(seq2[j - 1])
                    j -= 1
                    state = "M"
                else:
                    aligned_seq1.append("-")
                    aligned_seq2.append(seq2[j - 1])
                    j -= 1

        aligned_seq1 = "".join(reversed(aligned_seq1))
        aligned_seq2 = "".join(reversed(aligned_seq2))

    stats = alignment_statistics(aligned_seq1, aligned_seq2)
    return {
        "algorithm": "Smith-Waterman (Local, affine gap)" if gap_penalty is None else "Smith-Waterman (Local)",
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


def star_alignment(
    sequences: List[str],
    seq_type: str = "dna",
    gap_open: int = -10,
    gap_extend: int = -1,
    matrix_name: str = "default",
) -> Dict:
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
        aln = needleman_wunsch(
            reference,
            seq,
            seq_type=seq_type,
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=matrix_name,
        )
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
        column = [seq[pos] if pos < len(seq) else "-" for seq in aligned_sequences]
        # A column is conserved only when every displayed sequence contains
        # the same residue. Gaps must not be ignored because the chart claims
        # conservation across all sequences shown.
        if column[0] != "-" and len(set(column)) == 1:
            conserved += 1
    return round(conserved / align_length * 100, 2) if align_length else 0.0


def consensus_profile(aligned_sequences: List[str]) -> Dict:
    """Return consensus, conservation and variable-column data for an MSA."""
    if not aligned_sequences:
        return {"consensus": "", "columns": [], "variable_columns": [], "conservation_score": 0.0}
    width = max(len(sequence) for sequence in aligned_sequences)
    columns = []
    consensus = []
    variable_columns = []
    for index in range(width):
        column = [sequence[index] if index < len(sequence) else "-" for sequence in aligned_sequences]
        counts = {char: column.count(char) for char in set(column) if char != "-"}
        if not counts:
            consensus_char, conservation = "-", 0.0
        else:
            consensus_char, count = max(counts.items(), key=lambda item: (item[1], item[0]))
            conservation = round(count / len(column) * 100, 2)
        consensus.append(consensus_char)
        row = {
            "position": index + 1,
            "consensus": consensus_char,
            "conservation_percent": conservation,
            "counts": counts,
            "variable": len(counts) > 1,
        }
        columns.append(row)
        if row["variable"]:
            variable_columns.append(index + 1)
    return {
        "consensus": "".join(consensus),
        "columns": columns,
        "variable_columns": variable_columns,
        "conservation_score": round(
            sum(row["conservation_percent"] == 100.0 for row in columns) / width * 100,
            2,
        ) if width else 0.0,
    }


def pairwise_align(
    seq1: str,
    seq2: str,
    mode: str = "global",
    seq_type: str = "dna",
    gap_open: int = -10,
    gap_extend: int = -1,
    matrix_name: str = "default",
) -> Dict:
    """Convenience wrapper returning alignment + identity statistics."""
    if mode == "local":
        return smith_waterman(
            seq1, seq2, seq_type=seq_type, gap_open=gap_open,
            gap_extend=gap_extend, matrix_name=matrix_name,
        )
    return needleman_wunsch(
        seq1, seq2, seq_type=seq_type, gap_open=gap_open,
        gap_extend=gap_extend, matrix_name=matrix_name,
    )