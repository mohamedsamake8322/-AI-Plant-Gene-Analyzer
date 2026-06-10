"""Tests for alignment-based bioinformatics engines (no simulated metrics)."""

import alignment_engine as aln
import bioinformatics as bio
import distance_engine as dist
import phylogeny_engine as phylo
import similarityengine as sim
import numpy as np


def test_needleman_wunsch_identity():
    s = "ATGCGATCGATCG"
    result = aln.needleman_wunsch(s, s)
    assert result["identity_percent"] == 100.0
    assert result["gap_count"] == 0


def test_needleman_wunsch_with_indel():
    result = aln.needleman_wunsch("ATGC", "ATGGC")
    assert "-" in result["seq1_aligned"] or "-" in result["seq2_aligned"]
    assert 0 < result["identity_percent"] < 100


def test_star_alignment_preserves_all_sequences():
    seqs = ["ATGCGATCG", "ATGCGATCA", "ATGCGATCGA"]
    msa = aln.star_alignment(seqs)
    assert msa["num_sequences"] == 3
    assert len(msa["aligned_sequences"]) == 3
    assert len({len(s) for s in msa["aligned_sequences"]}) == 1


def test_detect_mutations_uses_alignment():
    query = "ATGCGATCA"
    ref = "ATGCGATCG"
    report = bio.detect_mutations(query, ref)
    assert report["identity_percent"] < 100
    assert report["total_mutations"] >= 1
    assert "alignment" in report
    assert report["alignment"]["query_aligned"]


def test_translate_reverse_frames():
    seq = "ATGAAATAA"
    frames = bio.translate_all_frames(seq, include_reverse=True)
    assert "Frame +1" in frames
    assert "Frame -1" in frames


def test_distance_matrix_after_alignment():
    sequences = [
        {"name": "A", "sequence": "ATGCGATCG"},
        {"name": "B", "sequence": "ATGCGATCA"},
        {"name": "C", "sequence": "ATGCGATCC"},
    ]
    dm = dist.distance_matrix(sequences, method="hamming")
    assert dm["distance_matrix"][0][0] == 0.0
    assert dm["distance_matrix"][0][1] > 0.0
    assert len(dm["aligned_sequences"]) == 3


def test_similarity_uses_global_alignment():
    db = {
        "GENE1": {
            "sequence": "ATGCGATCGATCG",
            "trait": "Test",
            "description": "Test gene",
            "organism": "Test",
            "accession": "T1",
        }
    }
    matches = sim.compare_with_database("ATGCGATCGATCA", db, top_n=1)
    assert matches
    assert matches[0]["similarity_score"] < 100
    assert matches[0]["alignment"]["algorithm"]


def test_upgma_produces_newick():
    matrix = np.array([
        [0.0, 0.1, 0.2],
        [0.1, 0.0, 0.15],
        [0.2, 0.15, 0.0],
    ])
    tree = phylo.upgma(matrix, ["A", "B", "C"])
    assert tree["newick"].endswith(";")
    assert "A" in tree["newick"] and "B" in tree["newick"]
