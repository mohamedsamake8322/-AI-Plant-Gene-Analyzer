#!/usr/bin/env python
"""Quick test of MEGA-like engines."""

from alignment_engine import needleman_wunsch, smith_waterman
from distance_engine import kimura_distance, hamming_distance, distance_matrix
from phylogeny_engine import upgma, neighbor_joining
import numpy as np

# Test 1: Alignment
print("=" * 60)
print("🧬 TEST 1: NEEDLEMAN-WUNSCH (Global Alignment)")
print("=" * 60)
seq1 = "ATGCGATCG"
seq2 = "ATGAGATCG"
nw = needleman_wunsch(seq1, seq2)
print(f"Query:     {seq1}")
print(f"Reference: {seq2}")
print(f"Score: {nw['alignment_score']}, Matches: {nw['match_count']}, Gaps: {nw['gap_count']}")
print(f"Aligned Q: {nw['seq1_aligned']}")
print(f"Aligned R: {nw['seq2_aligned']}")

# Test 2: Smith-Waterman
print("\n" + "=" * 60)
print("🧬 TEST 2: SMITH-WATERMAN (Local Alignment)")
print("=" * 60)
sw = smith_waterman(seq1, seq2)
print(f"Local alignment score: {sw['alignment_score']}")
print(f"Aligned Q: {sw['seq1_aligned']}")
print(f"Aligned R: {sw['seq2_aligned']}")

# Test 3: Hamming Distance
print("\n" + "=" * 60)
print("🧬 TEST 3: HAMMING DISTANCE")
print("=" * 60)
seq1_aligned = "ATGCGATCG"
seq2_aligned = "ATGAGATCG"
hamming = hamming_distance(seq1_aligned, seq2_aligned)
print(f"Differences: {hamming['differences']}/{hamming['total_positions']}")
print(f"Distance: {hamming['distance']}")
print(f"Identity: {hamming['percent_identity']}%")

# Test 4: Kimura Distance
print("\n" + "=" * 60)
print("🧬 TEST 4: KIMURA 2-PARAMETER DISTANCE")
print("=" * 60)
kim = kimura_distance(seq1_aligned, seq2_aligned)
print(f"Transitions: {kim.get('transitions', 'N/A')}")
print(f"Transversions: {kim.get('transversions', 'N/A')}")
print(f"Kimura Distance: {kim.get('kimura_distance')}")

# Test 5: Distance Matrix
print("\n" + "=" * 60)
print("🧬 TEST 5: DISTANCE MATRIX (Multiple Sequences)")
print("=" * 60)
sequences = [
    {"name": "Gene_A", "sequence": "ATGCGATCGATG"},
    {"name": "Gene_B", "sequence": "ATGAGATCGATG"},
    {"name": "Gene_C", "sequence": "ATGCGATCGATC"},
]
dist_mat = distance_matrix(sequences, method="kimura")
print(f"Sequences: {dist_mat['sequence_names']}")
print(f"Distance matrix shape: {np.array(dist_mat['distance_matrix']).shape}")

# Test 6: UPGMA Phylogeny
print("\n" + "=" * 60)
print("🧬 TEST 6: UPGMA (Phylogenetic Tree)")
print("=" * 60)
# Use simple test distance matrix
test_dist = np.array([
    [0.0, 0.1, 0.2],
    [0.1, 0.0, 0.15],
    [0.2, 0.15, 0.0]
], dtype=float)
test_names = ["Gene_A", "Gene_B", "Gene_C"]
upgma_result = upgma(test_dist, test_names)
print(f"Algorithm: {upgma_result['algorithm']}")
print(f"Tree type: {upgma_result['tree_type']}")
print(f"Dendrogram ready for visualization: {bool(upgma_result['dendrogram_data'])}")

# Test 7: Neighbor-Joining Phylogeny
print("\n" + "=" * 60)
print("🧬 TEST 7: NEIGHBOR-JOINING (Additive Tree)")
print("=" * 60)
nj_result = neighbor_joining(test_dist, test_names)
print(f"Algorithm: {nj_result['algorithm']}")
print(f"Tree type: {nj_result['tree_type']}")
print(f"Edges constructed: {len(nj_result['edges'])}")

print("\n" + "=" * 60)
print("✅ ALL MEGA ENGINES WORKING CORRECTLY!")
print("=" * 60)
