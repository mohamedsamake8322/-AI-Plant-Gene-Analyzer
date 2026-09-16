import config
import pipeline
import similarityengine


def test_similarity_skip_reason_for_overlong_sequence():
    result = pipeline.analyze_sequence_record(
        {"header": "long-sequence", "sequence": "A" * (config.MAX_ALIGNMENT_SEQUENCE_LENGTH + 1)},
        "dna",
        reading_frame=0,
        db={},
    )

    assert result["similarity_skipped_reason"] == "sequence_too_long"
    assert result["similarity_results"] == []
    assert any("alignment threshold" in warning for warning in result["metadata_warnings"])


def test_empty_similarity_results_have_no_skip_reason():
    result = pipeline.analyze_sequence_record(
        {"header": "short-sequence", "sequence": "ATGAAATAGCCC"},
        "dna",
        reading_frame=0,
        db={},
    )

    assert result["similarity_skipped_reason"] is None
    assert result["similarity_results"] == []


def test_distant_best_match_does_not_create_mutation_report():
    result = pipeline.analyze_sequence_record(
        {"header": "distant", "sequence": "A" * 12},
        "DNA",
        reading_frame=0,
        db={"reference": {"sequence": "C" * 12, "sequence_type": "dna"}},
        enable_length_prefilter=False,
    )

    assert result["mutation_report"] is None
    assert result["variant_report"] is None
    assert any("No close reference found" in warning for warning in result["metadata_warnings"])


def test_explicit_reference_bypasses_similarity_threshold():
    result = pipeline.analyze_sequence_record(
        {"header": "explicit-reference", "sequence": "ATGAAATAGCCC"},
        "DNA",
        reading_frame=0,
        db={},
        reference_sequence="ATGAGA TAGCCC".replace(" ", ""),
    )

    assert result["mutation_reference_source"] == "explicit_reference"
    assert result["mutation_report"]["total_mutations"] == 1
    assert result["variant_report"]["substitutions"][0]["consequence"] == "silent"


def test_similarity_skip_reason_for_candidate_cost(monkeypatch):
    monkeypatch.setattr(config, "MAX_ALIGNMENT_CELL_BUDGET", 1)
    result = pipeline.analyze_sequence_record(
        {"header": "candidate-cost", "sequence": "ATGAAATAGCCC"},
        "dna",
        reading_frame=0,
        db={"candidate": {"sequence": "ATGAAATAGCCC"}},
    )

    assert result["similarity_skipped_reason"] == "alignment_cost_too_high"
    assert result["similarity_results"] == []


def test_long_query_candidate_pool_respects_alignment_budget(monkeypatch):
    monkeypatch.setattr(config, "MAX_ALIGNMENT_CELL_BUDGET", 300)

    assert similarityengine._budgeted_candidate_pool_size(10, 45) == 3
    assert similarityengine._budgeted_candidate_pool_size(100, 45) == 1


def test_compare_with_database_uses_traits_list_of_dicts():
    db = {
        "gene_trait_test": {
            "sequence": "ATGAAATAGCCC",
            "sequence_type": "dna",
            "organism": "Test species",
            "traits": [
                {"trait": "ATP-binding", "source": "planttfdb"},
                {"trait": "Chromatin regulator", "source": "planttfdb"},
                {"trait": "Coiled coil", "source": "planttfdb"},
            ],
        }
    }

    result = similarityengine.compare_with_database("ATGAAATAGCCC", db, top_n=1, enable_length_prefilter=False)

    assert result
    assert result[0]["trait"] == "ATP-binding, Chromatin regulator, Coiled coil"


def test_similarity_candidates_track_requested_vs_actual_pool_size():
    candidates = similarityengine.SimilarityCandidates(
        {"g1": {"sequence": "ATG"}},
        source="kmer_prefilter",
        requested_pool_size=75,
    )

    assert candidates.requested_pool_size == 75
    assert candidates.actual_pool_size == 1
    assert candidates.pool_reduced is True


def test_pipeline_keeps_candidate_pool_metadata_without_crashing():
    sequence = "ATGAAATAGCCC"
    db = similarityengine.SimilarityCandidates(
        {"gene": {"sequence": sequence, "sequence_type": "dna"}},
        source="kmer_prefilter",
        requested_pool_size=75,
    )

    result = pipeline.analyze_sequence_record(
        {"header": "candidate-pool-check", "sequence": sequence},
        "dna",
        reading_frame=0,
        db=db,
    )

    assert result["similarity_candidate_pool_requested"] == 75
    assert result["similarity_candidate_pool_reduced"] is True


def test_compare_with_database_extracts_trait_labels_from_traits_list():
    db = {
        "A0A803N8F1": {
            "sequence": "ATGC" * 25,
            "sequence_type": "dna",
            "organism": "Chenopodium quinoa",
            "traits": [
                {"trait": "ATP-binding", "source": "planttfdb"},
                {"trait": "Chromatin regulator", "source": "planttfdb"},
                {"trait": "Coiled coil", "source": "planttfdb"},
            ],
        }
    }

    results = similarityengine.compare_with_database("ATGC" * 25, db, top_n=1, logger=None)

    assert results
    assert results[0]["trait"] == "ATP-binding, Chromatin regulator, Coiled coil"


def test_similarity_exposes_global_primary_and_local_complementary_scores():
    match = similarityengine.aligned_similarity("ATGCATGC", "ATGCATGC", compute_local=True)
    assert match["method"] == "global"
    assert match["similarity_score"] == 100.0
    assert match["global_identity"] == 100.0
    assert match["local_identity"] == 100.0
    assert match["local_coverage_percent"] == 100.0

    global_match = similarityengine.aligned_similarity("ATGCATGC", "ATGCATGC", compute_local=False)
    assert global_match["method"] == "global"
    assert global_match["similarity_score"] == 100.0


def test_compare_with_database_uses_global_scoring_by_default():
    db = {
        "gene_1": {"sequence": "ATGCAATGCGTA", "sequence_type": "dna"},
        "gene_2": {"sequence": "ATGCCATGCGTA", "sequence_type": "dna"},
    }
    results = similarityengine.compare_with_database("ATGCAATGCGTA", db, top_n=2, enable_length_prefilter=False)
    assert results
    assert results[0]["alignment_method"] == "global"
    assert results[0]["similarity_score"] >= 90.0


def test_global_score_remains_primary_when_local_finds_perfect_subsegment():
    reference = "ATGC" * 25
    query = list(reference)
    for index in (0, 4, 8, 12, 16, 20, 24, 28):
        query[index] = {"A": "C", "C": "G", "G": "T", "T": "A"}[query[index]]

    match = similarityengine.aligned_similarity("".join(query), reference, compute_local=True)

    assert match["similarity_score"] == match["global_identity"]
    assert match["similarity_score"] < match["local_identity"]
    assert 90.0 < match["local_coverage_percent"] < 100.0


def test_local_alignment_estimates_small_evalue_for_near_perfect_match():
    db = {
        "gene_1": {"sequence": "ATGCATGCATGC", "sequence_type": "dna"},
    }
    results = similarityengine.compare_with_database("ATGCATGCATGC", db, top_n=1, enable_length_prefilter=False)
    assert results
    assert "evalue" not in results[0]
    assert results[0]["approximate_local_significance"] > 0.0
    assert results[0]["approximate_local_significance"] < 1.0
