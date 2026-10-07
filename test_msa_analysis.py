import alignment_engine as aln
import alignment_exports as exports
from slider_utils import _alignment_window_slider_bounds
from msa_analysis import (
    alignment_column_to_residue,
    analyze_alignment,
    iupac_consensus_char,
    qc_unaligned,
    residue_to_alignment_column,
    trim_alignment,
)


def test_iupac_and_variable_definition_ignore_gaps():
    profile = analyze_alignment(["ATGC", "ATGT", "ATGA"])
    assert profile["iupac_consensus"] == "ATGH"
    assert profile["consensus"] == "ATGH"
    assert profile["variable_columns"] == [4]
    assert profile["site_counts"]["conserved"] == 3
    assert profile["site_counts"]["variable"] == 1
    assert profile["site_counts"]["singletons"] == 1
    assert profile["site_counts"]["parsimony_informative"] == 0

    gapped = analyze_alignment(["AT-C", "ATAC", "ATAC"])
    assert 3 not in gapped["variable_columns"]
    assert gapped["columns"][2]["class"] == "conserved"
    assert 3 in gapped["gapped_columns"]


def test_parsimony_informative_requires_two_states_twice():
    profile = analyze_alignment(["AA", "AT", "AT", "AA"])
    assert profile["columns"][1]["class"] == "parsimony_informative"
    assert profile["variable_columns"] == [2]


def test_terminal_gaps_are_missing_not_indels():
    profile = analyze_alignment(["--ATGC", "CGATGC", "CGATGC"])
    first = profile["columns"][0]
    assert first["class"] == "conserved"
    assert first["missing_fraction"] > 0
    assert first["gap_fraction"] == 0
    assert profile["sequence_qc"][0]["internal_gaps"] == 0


def test_coordinates_round_trip_with_real_headers():
    sequences = ["ATG-C", "ATGCC"]
    labels = ["A_110688118", "B_110684354"]
    assert alignment_column_to_residue(sequences[0], 4) is None
    assert alignment_column_to_residue(sequences[0], 5) == 4
    assert residue_to_alignment_column(sequences[0], 4) == 5
    metadata = exports.reproducibility_metadata(
        ["ATGC", "ATGCC"],
        labels,
        "Star MSA (fast / approximate)",
        "dna",
        "DNA",
        -10,
        -1,
        aligned_sequences=sequences,
    )
    assert metadata["input_sha256"]
    mega = exports.mega(sequences, labels, metadata)
    nexus = exports.nexus(sequences, labels, metadata)
    assert mega.startswith("#MEGA")
    assert "CHARSET codon1" in nexus
    assert "#A_110688118" in mega


def test_qc_detects_duplicates_and_mixed_types():
    issues = qc_unaligned(["ATGCATGCATGCATGCATGCATGCATGCAT", "ATGCATGCATGCATGCATGCATGCATGCAT"], ["a", "b"], "dna")
    assert any(issue["code"] == "duplicate" for issue in issues)
    mixed = qc_unaligned(["ATGATGATGATGATGATGATGATGATGATG", "MKTLLILAV"], ["dna", "prot"], "dna")
    assert any(issue["code"] == "mixed_types" for issue in mixed)


def test_codon_trimming_keeps_multiples_of_three():
    aligned = ["ATG---CCC", "ATGAAACCC"]
    trimmed = trim_alignment(aligned, mode="gt", gt=0.9, codon=True)
    assert all(len(sequence) % 3 == 0 for sequence in trimmed["aligned_sequences"])
    assert trimmed["removed_columns"] == [4, 5, 6]


def test_alignment_window_slider_bounds_never_equal():
    assert _alignment_window_slider_bounds(20) == (1, 20, 20)
    assert _alignment_window_slider_bounds(100) == (20, 100, 60)


def test_star_alignment_uses_same_site_classes():
    result = aln.star_alignment(["ATGC", "ATGT", "ATGA"], seq_type="dna")
    profile = aln.consensus_profile(result["aligned_sequences"])
    assert profile["consensus"] == "ATGH"
    assert iupac_consensus_char(["T", "G", "A"]) == "D"
