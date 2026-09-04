"""Shared sequence-quality rules used by ingestion and analysis."""

from __future__ import annotations

VALID_CHARS = {
    "dna": set("ACGTURYSWKMBDHVN"),
    "rna": set("ACGTURYSWKMBDHVN"),
    "protein": set("ACDEFGHIKLMNPQRSTVWYXBZJUO*"),
}
MIN_SEQUENCE_LENGTH = 50
MAX_N_RATIO = 0.05


def validate_sequence_quality(
    sequence: str | None,
    sequence_type: str | None = "dna",
    min_length: int = MIN_SEQUENCE_LENGTH,
    max_n_ratio: float = MAX_N_RATIO,
) -> tuple[bool, str]:
    """Apply the same quality gate used before PostgreSQL ingestion."""
    if not sequence:
        return False, "empty"
    seq = sequence.upper().strip()
    seq_type = (sequence_type or "dna").lower()
    if seq_type == "protein":
        return True, "not_applicable_protein"
    if len(seq) < min_length:
        return False, f"too_short(<{min_length})"
    n_ratio = seq.count("N") / len(seq)
    if n_ratio > max_n_ratio:
        return False, f"too_many_n({n_ratio:.1%})"
    allowed = VALID_CHARS.get(seq_type, VALID_CHARS["dna"])
    invalid = set(seq) - allowed
    if invalid:
        return False, f"invalid_chars({''.join(sorted(invalid))})"
    return True, ""


def quality_report(
    sequence: str | None,
    sequence_type: str | None = "dna",
    min_length: int = MIN_SEQUENCE_LENGTH,
    max_n_ratio: float = MAX_N_RATIO,
) -> dict[str, object]:
    """Return UI-friendly quality status and the N percentage."""
    seq = (sequence or "").upper().strip()
    seq_type = (sequence_type or "dna").lower()
    n_pct = (seq.count("N") / len(seq) * 100) if seq else 0.0
    valid, reason = validate_sequence_quality(seq, seq_type, min_length, max_n_ratio)
    return {
        "valid": valid,
        "reason": None if not reason or reason == "not_applicable_protein" else reason,
        "n_pct": round(n_pct, 2),
        "threshold_pct": round(max_n_ratio * 100, 2),
        "min_length": min_length,
        "applicable": seq_type != "protein",
    }
