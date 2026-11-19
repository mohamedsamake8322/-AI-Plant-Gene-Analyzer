from analysis_state import build_analysis_signature, should_reuse_cached_results


def test_analysis_signature_changes_when_input_changes():
    sig1 = build_analysis_signature("ATGC", "", "DNA", 1, 10, False)
    sig2 = build_analysis_signature("ATGC", "", "DNA", 1, 10, False)
    sig3 = build_analysis_signature("GGGG", "", "DNA", 1, 10, False)
    sig4 = build_analysis_signature("GGGG", "", "DNA", 2, 10, False)

    assert sig1 == sig2
    assert sig1 != sig3
    assert sig3 != sig4


def test_cached_results_only_reused_for_matching_signature():
    current = build_analysis_signature("ATGC", "", "DNA", 1, 10, False)
    other = build_analysis_signature("GGGG", "", "DNA", 1, 10, False)

    assert should_reuse_cached_results(current, current, "ATGC") is True
    assert should_reuse_cached_results(other, current, "GGGG") is False
    assert should_reuse_cached_results(current, other, "ATGC") is False
