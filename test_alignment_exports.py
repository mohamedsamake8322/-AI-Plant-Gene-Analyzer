import alignment_engine as aln
import alignment_exports as exports


def test_protein_alignment_supports_standard_matrices():
    result = aln.needleman_wunsch("MKTLL", "MKTLV", seq_type="protein", matrix_name="BLOSUM80")
    assert result["identity_percent"] == 80.0
    assert result["algorithm"].startswith("Needleman-Wunsch")


def test_msa_consensus_reports_variable_columns():
    result = aln.star_alignment(["ATGC", "ATGT", "ATGA"], seq_type="dna")
    profile = aln.consensus_profile(result["aligned_sequences"])
    assert profile["consensus"] == "ATGT"
    assert profile["variable_columns"] == [4]


def test_alignment_exports_embed_reproducibility_metadata():
    sequences = ["ATGC", "ATGT"]
    labels = ["ref", "variant"]
    metadata = exports.reproducibility_metadata(
        sequences, labels, "Star MSA (reference-guided)", "dna", "DNA", -10, -1
    )
    fasta = exports.aligned_fasta(sequences, labels, metadata)
    clustal = exports.clustal(sequences, labels, metadata)
    assert metadata["input_sha256"]
    assert "gap_open" in fasta
    assert "input_sha256" in clustal
    assert ">ref" in fasta
