"""
Wrapper for distance_engine inside core_engines package.
"""
import distance_engine as _orig

hamming_distance = _orig.hamming_distance
jukes_cantor_distance = _orig.jukes_cantor_distance
kimura_distance = _orig.kimura_distance
pam_distance = _orig.pam_distance
distance_matrix = _orig.distance_matrix

__all__ = [
    'hamming_distance',
    'jukes_cantor_distance',
    'kimura_distance',
    'pam_distance',
    'distance_matrix',
]
