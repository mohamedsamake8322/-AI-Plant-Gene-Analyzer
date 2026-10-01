"""Compatibility exports for shared sequence-quality rules."""

from scripts.quality_rules import (
    MAX_N_RATIO,
    MIN_SEQUENCE_LENGTH,
    quality_report,
    validate_sequence_quality,
)

__all__ = [
    "MAX_N_RATIO",
    "MIN_SEQUENCE_LENGTH",
    "quality_report",
    "validate_sequence_quality",
]