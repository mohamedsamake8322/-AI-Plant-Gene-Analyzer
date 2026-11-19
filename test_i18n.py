import importlib


def test_translation_fallback_and_glossary():
    i18n = importlib.import_module("i18n")

    assert i18n.translate("ui.app_title", lang="fr") == "Analyseur de gènes végétaux"
    assert i18n.translate("ui.app_title", lang="en") == "Plant Gene Analyzer"
    assert i18n.translate("ui.app_title", lang="tr") == "Bitki Gen Analizörü"

    assert i18n.translate("ui.missing_key", lang="tr") == "ui.missing_key"
    assert i18n.get_glossary("en", "gene") == "gene"
    assert i18n.get_glossary("fr", "gene") == "gène"


def test_esm2_ui_translations_cover_all_supported_languages():
    i18n = importlib.import_module("i18n")
    keys = (
        "ui.esm2_similarity_metric",
        "ui.esm2_distance_note",
        "ui.esm2_tree_warning",
        "ui.esm2_accession_column",
        "ui.esm2_organism_column",
        "ui.esm2_similarity_column",
        "ui.esm2_annotation_column",
        "ui.esm2_scope_note",
        "ui.esm2_pairwise_error",
        "ui.esm2_search_error",
        "ui.distance_matrix_compute_error",
        "ui.trait_esm2_candidate_label",
        "ui.trait_esm2_scope_note",
        "ui.trait_esm2_button",
    )

    for language in i18n.SUPPORTED_LANGS:
        for key in keys:
            assert i18n.translate(key, lang=language) != key
        assert "database unavailable" in i18n.translate(
            "ui.esm2_search_error", lang=language
        ).format(error="database unavailable").lower()


def test_ai_interpretation_localizes_output_for_french():
    i18n = importlib.import_module("i18n")
    ai = importlib.import_module("aiinterpreter")

    stats = {
        "length": 1200,
        "gc_content": 55,
        "has_start_codon": True,
        "has_stop_codon": True,
        "sequence_type": "dna",
    }
    similarity_results = [{
        "gene_name": "P5CS",
        "trait": "drought",
        "organism": "Arabidopsis thaliana",
        "similarity_score": 92,
    }]

    report = ai.AIInterpreter(stats, similarity_results, None).full_report(lang="fr")

    assert isinstance(report["overall_summary"], str)
    assert "drought" not in report["overall_summary"].lower()
    assert "séquence" in report["overall_summary"].lower() or "analyse" in report["overall_summary"].lower()
    assert "Aucune" in i18n.translate("ai.no_similarity_data", lang="fr") or "Aucune" in i18n.translate("ai.no_mutation_data", lang="fr")
