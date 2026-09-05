import config
import pipeline


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
