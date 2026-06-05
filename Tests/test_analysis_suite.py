def test_analysis_suite_imports():
    import scripts.run_analysis_suite as suite

    report = {
        "generated_at": "2026-01-01T00:00:00Z",
        "records_analyzed": 0,
        "type_counts": {},
        "top_alignments": [],
        "sequence_analysis": [],
        "analysis_notes": [],
    }

    markdown = suite.format_markdown_report(report)
    assert "Bioinformatics Analysis Suite Report" in markdown
