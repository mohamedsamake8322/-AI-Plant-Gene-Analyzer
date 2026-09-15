import importlib


def test_translation_fallback_and_glossary():
    i18n = importlib.import_module("i18n")

    assert i18n.translate("ui.app_title", lang="fr") == "Analyseur de gènes végétaux"
    assert i18n.translate("ui.app_title", lang="en") == "Plant Gene Analyzer"
    assert i18n.translate("ui.app_title", lang="tr") == "Bitki Gen Analizörü"

    assert i18n.translate("ui.missing_key", lang="tr") == "ui.missing_key"
    assert i18n.get_glossary("en", "gene") == "gene"
    assert i18n.get_glossary("fr", "gene") == "gène"


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
