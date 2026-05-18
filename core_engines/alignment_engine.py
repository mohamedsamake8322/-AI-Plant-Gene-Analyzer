"""
Copied alignment_engine into core_engines package for refactor.
"""

# For simplicity duplicate of root module to allow import from core_engines
import alignment_engine as _orig

# Re-export functions
needleman_wunsch = _orig.needleman_wunsch
smith_waterman = _orig.smith_waterman
progressive_alignment = _orig.progressive_alignment

__all__ = [
    'needleman_wunsch',
    'smith_waterman',
    'progressive_alignment',
]
