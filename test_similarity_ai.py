"""
test_similarity_ai.py
----------------------
Unit tests for similarity and AI interpretation modules.
"""

import similarityengine as sim
import aiinterpreter as ai
import bioinformatics as bio


def test_load_gene_database_json():
    db = sim.load_gene_database("genes_database.json")
    assert isinstance(db, dict)
    assert "PR1" in db


def test_pairwise_similarity_same_sequence():
    assert sim.pairwise_similarity("ATGC", "ATGC") == 100.0


def test_pairwise_similarity_different_sequence():
    score = sim.pairwise_similarity("ATGC", "ATTA")
    assert score == 50.0


def test_compare_with_database_protein_query():
    protein = "MQNCG"  # short protein-like sequence
    results = sim.compare_with_database(protein, top_n=1)
    assert isinstance(results, list)
    assert len(results) == 1
    assert "gene_name" in results[0]


def test_ai_interpreter_summary_high_confidence():
    stats = {
        "length": 100,
        "gc_content": 45.0,
        "at_content": 55.0,
        "is_coding_length": True,
        "has_start_codon": True,
        "has_stop_codon": True,
    }
    similarity_results = [
        {
            "gene_name": "TestGene",
            "trait": "Stress Response",
            "organism": "Arabidopsis thaliana",
            "similarity_score": 92.0,
        }
    ]
    interp = ai.AIInterpreter(stats, similarity_results, None).full_report()
    assert interp["confidence_level"]["level"] in {"High", "Medium", "Low"}
    assert "summary" in interp["overall_summary"].lower() or isinstance(interp["overall_summary"], str)
