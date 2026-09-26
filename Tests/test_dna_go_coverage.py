import json

from scripts.analyze_dna_go_coverage import analyze_file


def test_dna_go_coverage_finds_unmerged_accession_candidates(tmp_path):
    data = {
        "metadata": {"plant": "Test plant"},
        "genes": [
            {
                "gene_id": "GeneID:1",
                "organism": "Test plant",
                "sequence": {"dna": "ATGC", "rna": "AUGC"},
                "external_links": {
                    "accession": "NM_MATCH.1",
                    "ncbi_dna_accessions": ["NM_MATCH.1"],
                    "ncbi_rna_accessions": ["NM_MATCH.1"],
                },
                "sources_summary": ["ncbi", "uniprot"],
                "annotation": {"go_terms": [{"id": "GO:DNA", "evidence_code": "IEA"}]},
            },
            {
                "gene_id": "P1",
                "organism": "Test plant",
                "sequence": {"protein": "MPEP"},
                "external_links": {"refseq_nucleotide": "NM_MATCH.1"},
                "sources_summary": ["uniprot"],
                "annotation": {"go_terms": [{"id": "GO:P1", "evidence_code": "EXP"}]},
            },
            {
                "gene_id": "P2",
                "organism": "Test plant",
                "sequence": {"protein": "MKPE"},
                "external_links": {"refseq_nucleotide": "NM_NOT_COLLECTED.1"},
                "sources_summary": ["uniprot"],
                "annotation": {"go_terms": [{"id": "GO:P2", "evidence_code": "IEA"}]},
            },
            {
                "gene_id": "G3",
                "organism": "Test plant",
                "sequence": {},
                "sources_summary": ["plaza"],
                "annotation": {"go_terms": [{"id": "GO:NONE", "evidence_code": "IEA"}]},
            },
        ],
    }
    input_path = tmp_path / "test_species_all_sources.json"
    input_path.write_text(json.dumps(data), encoding="utf-8")

    result = analyze_file(input_path)

    assert result["records"]["records"] == 4
    assert result["records"]["go_with_dna"] == 1
    assert result["records"]["go_without_dna_with_protein"] == 2
    assert result["records"]["go_without_any_sequence"] == 1
    assert result["cross_reference_recovery"]["GO_protein_records_without_DNA_with_nucleotide_xref"] == 2
    assert result["cross_reference_recovery"]["those_with_xref_matching_a_collected_DNA_accession"] == 1
    assert result["cross_reference_recovery"]["those_with_xref_matching_a_collected_RNA_accession"] == 1
    assert result["cross_reference_recovery"]["candidate_not_currently_merged_to_nucleotide"] == 1
