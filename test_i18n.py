import importlib


def test_translation_fallback_and_glossary():
    i18n = importlib.import_module("i18n")

    assert i18n.translate("ui.app_title", lang="fr") == "Analyseur de gènes végétaux"
    assert i18n.translate("ui.app_title", lang="en") == "Plant Gene Analyzer"
    assert i18n.translate("ui.app_title", lang="tr") == "Bitki Gen Analizörü"

    assert i18n.translate("ui.missing_key", lang="tr") == "ui.missing_key"
    assert i18n.get_glossary("en", "gene") == "gene"
    assert i18n.get_glossary("fr", "gene") == "gène"
