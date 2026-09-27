import json

import pytest

from collect import collect_all_sources
from collect.collect_all_sources import (
    _build_accession_index,
    _build_uniprot_index,
    _accept_ncbi_record,
    _preflight_quality_gate,
    _source_quality_report,
    restructure_to_schema,
)
from scripts.collect_ncbi import make_record_from_fasta
from collect.collect_kegg import _parse_kegg_flat


def test_restructure_preserves_uniprot_functional_annotations():
    flat_record = {
        "gene_id": "P12345",
        "symbol": "GENE1",
        "organism": "Oryza sativa",
        "description": "Example transcription factor",
        "source": "uniprot",
        "external_links": {"accession": "P12345"},
        "_raw_sequences": {"protein": "MPEPTIDE"},
        "annotations": {
            "function": "Regulates grain development.",
            "keywords": ["Transcription", "DNA-binding"],
            "keywords_detailed": [{"name": "Transcription", "category": "Biological process"}],
            "subcellular_location": ["Nucleus"],
            "annotation_score": 5,
            "protein_existence": "Evidence at protein level",
            "reviewed": True,
            "go_terms": [{"id": "GO:0000001", "evidence_code": "IDA"}],
        },
    }

    record = restructure_to_schema("P12345", flat_record)

    assert record["common_name"] == "Example transcription factor"
    assert record["external_links"]["uniprot_accession"] == "P12345"
    assert record["annotation"]["protein_name"] == "Example transcription factor"
    assert record["annotation"]["function"] == "Regulates grain development."
    assert record["annotation"]["keywords"] == ["Transcription", "DNA-binding"]
    assert record["annotation"]["keywords_detailed"] == [
        {"name": "Transcription", "category": "Biological process"}
    ]
    assert record["annotation"]["subcellular_location"] == ["Nucleus"]
    assert record["annotation"]["annotation_score"] == 5
    assert record["annotation"]["protein_existence"] == "Evidence at protein level"
    assert record["annotation"]["reviewed"] is True
    assert record["annotation"]["go_terms"][0]["id"] == "GO:0000001"


def test_ncbi_records_keep_sequence_origin_and_evidence():
    genomic = make_record_from_fasta(
        "NC_000001.1 Arabidopsis thaliana locus",
        "ATGC",
        resolved_gene_id="1",
        organism="Arabidopsis thaliana",
    )
    transcript = make_record_from_fasta(
        "NM_000001.1 Arabidopsis thaliana transcript",
        "ATGC",
        organism="Arabidopsis thaliana",
    )

    assert genomic["sequence_origin"] == "genomic_locus"
    assert genomic["evidence_code"] == "NCBI_GENE_LINK"
    assert transcript["sequence_origin"] == "transcript_only"


def test_source_quality_report_explains_zero_and_failure():
    report = _source_quality_report(
        ["ncbi", "kegg", "ensembl"],
        {"ncbi": 3, "kegg": 0},
        ["kegg: request failed"],
    )

    assert report["ncbi"]["status"] == "success"
    assert report["kegg"]["status"] == "request_failed"
    assert report["ensembl"]["status"] == "unsupported"
    assert report["ncbi"]["records_accepted"] == 3


def test_preflight_keeps_sequence_types_separate():
    gate = _preflight_quality_gate({
        "g1": {"_raw_sequences": {"dna": "ATGC", "protein": "MPEP"}},
    })

    assert gate["passed"] is True
    assert gate["sequence_records_by_type"] == {"dna": 1, "rna": 0, "protein": 1}


def test_ncbi_accession_only_dna_is_not_treated_as_gene_locus():
    record = make_record_from_fasta("AB123456.1 plant sequence", "ATGC")
    assert not _accept_ncbi_record(record, "dna")


def test_accession_index_matches_versioned_nucleotide_cross_references():
    records = {
        "GeneID:1": {
            "gene_id": "GeneID:1",
            "external_links": {"accession": "NM_001234.2"},
            "_source_accessions": {"rna": ["NM_001234.2"], "dna": ["NC_000001.1"]},
        },
    }

    index = _build_accession_index(records)

    assert index["NM_001234"] == "GeneID:1"
    assert index["NC_000001"] == "GeneID:1"


def test_uniprot_index_retains_alias_after_merging_onto_ncbi_gene():
    records = {
        "GeneID:1": {
            "gene_id": "GeneID:1",
            "external_links": {
                "accession": "NC_000001.1",
                "uniprot_accession": "A0A123",
            },
        },
    }

    index = _build_uniprot_index(records)

    assert index["A0A123"] == "GeneID:1"
    assert index["NC_000001.1"] == "GeneID:1"


def test_kegg_pathways_are_not_mislabeled_as_traits():
    entry = """ENTRY       osa:LOC_Os01g00010  CDS\nNAME        test_gene\nDEFINITION  example enzyme\nPATHWAY     osa00010  Glycolysis / Gluconeogenesis\nNTSEQ       6\n            ATGGCT\n///\n"""

    record = _parse_kegg_flat(entry, "LOC_Os01g00010", "osa", "Oryza sativa")

    assert record is not None
    assert record["pathways"] == [
        {"id": "osa00010", "name": "Glycolysis / Gluconeogenesis", "source": "kegg"}
    ]
    assert record.get("traits", []) == []

    normalized = restructure_to_schema(
        "osa:LOC_Os01g00010",
        {
            **record,
            "traits": ["Glycolysis / Gluconeogenesis", "TF:MYB", "transcription_factor"],
        },
    )
    assert [trait["trait"] for trait in normalized["traits"]] == ["TF:MYB", "transcription_factor"]
    assert "kegg" in normalized["sources_summary"]
    assert "uniprot" not in normalized["sources_summary"]


def test_species_output_writer_replaces_atomically(tmp_path, monkeypatch):
    output_path = tmp_path / "species.json"
    original = '{"metadata":{"version":"old"},"genes":[]}'
    output_path.write_text(original, encoding="utf-8")

    def fail_restructure(*args):
        raise RuntimeError("synthetic serialization failure")

    monkeypatch.setattr(collect_all_sources, "restructure_to_schema", fail_restructure)
    with pytest.raises(RuntimeError, match="synthetic serialization failure"):
        collect_all_sources._write_species_json_atomic(
            output_path,
            {"version": "new"},
            {"g1": {"gene_id": "g1"}},
        )

    assert output_path.read_text(encoding="utf-8") == original
    assert json.loads(output_path.read_text(encoding="utf-8"))["metadata"]["version"] == "old"
    assert list(tmp_path.glob(".species.json.*.tmp")) == []
