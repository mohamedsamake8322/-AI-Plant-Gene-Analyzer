"""
test_bioinformatics.py
----------------------
Unit tests for bioinformatics module.
Run with: pytest test_bioinformatics.py -v
"""

import pytest
import bioinformatics as bio
import sequence_loader as loader


class TestSequenceCleaning:
    """Test sequence cleaning and validation."""
    
    def test_clean_sequence_basic(self):
        """Test basic sequence cleaning."""
        raw = "ATGC\nGATA"
        cleaned = bio.clean_sequence(raw)
        assert cleaned == "ATGCGATA"
    
    def test_clean_sequence_with_fasta_header(self):
        """Test removal of FASTA headers."""
        fasta = ">chr1\nATGCGATA\n>chr2\nTTAAGC"
        cleaned = bio.clean_sequence(fasta)
        assert cleaned == "ATGCGATATTAAGC"
    
    def test_clean_sequence_lowercase(self):
        """Test conversion to uppercase."""
        raw = "atgc"
        cleaned = bio.clean_sequence(raw)
        assert cleaned == "ATGC"
    
    def test_clean_sequence_removes_invalid_chars(self):
        """Test removal of invalid characters."""
        raw = "ATG@C#GA$TA"
        cleaned = bio.clean_sequence(raw)
        assert cleaned == "ATGCGATA"
    
    def test_clean_sequence_preserves_n(self):
        """Test that N (unknown nucleotide) is preserved."""
        raw = "ATGCNGATA"
        cleaned = bio.clean_sequence(raw)
        assert "N" in cleaned
    
    def test_validate_sequence_valid(self):
        """Test validation of valid sequence."""
        is_valid, msg = bio.validate_sequence("ATGCGATATATGC")
        assert is_valid is True
    
    def test_validate_sequence_too_short(self):
        """Test validation of too short sequence."""
        is_valid, msg = bio.validate_sequence("AT")
        assert is_valid is False
        assert "too short" in msg.lower()
    
    def test_validate_sequence_empty(self):
        """Test validation of empty sequence."""
        is_valid, msg = bio.validate_sequence("")
        assert is_valid is False
        assert "empty" in msg.lower()
    
    def test_validate_sequence_invalid_chars(self):
        """Test validation with invalid characters."""
        is_valid, msg = bio.validate_sequence("ATGCXYZ123")
        assert is_valid is False
        assert "invalid" in msg.lower()


class TestGCContent:
    """Test GC content calculation."""
    
    def test_calculate_gc_content_basic(self):
        """Test basic GC content calculation."""
        seq = "ATGC"
        gc = bio.calculate_gc_content(seq)
        assert gc == 50.0
    
    def test_calculate_gc_content_all_gc(self):
        """Test sequence with only G and C."""
        seq = "GGCCGGCC"
        gc = bio.calculate_gc_content(seq)
        assert gc == 100.0
    
    def test_calculate_gc_content_no_gc(self):
        """Test sequence with no G or C."""
        seq = "AAATTTT"
        gc = bio.calculate_gc_content(seq)
        assert gc == 0.0
    
    def test_calculate_gc_content_empty(self):
        """Test GC content of empty sequence."""
        gc = bio.calculate_gc_content("")
        assert gc == 0.0


class TestNucleotideDistribution:
    """Test nucleotide counting and distribution."""
    
    def test_nucleotide_distribution_equal(self):
        """Test distribution with equal nucleotides."""
        seq = "ATGC"
        dist = bio.nucleotide_distribution(seq)
        assert dist["counts"]["A"] == 1
        assert dist["counts"]["T"] == 1
        assert dist["counts"]["G"] == 1
        assert dist["counts"]["C"] == 1
        assert all(pct == 25.0 for pct in dist["percentages"].values() if pct > 0)
    
    def test_nucleotide_distribution_with_n(self):
        """Test distribution with unknown nucleotides."""
        seq = "ATGCN"
        dist = bio.nucleotide_distribution(seq)
        assert dist["counts"]["N"] == 1
        assert dist["percentages"]["N"] == 20.0


class TestProteinTranslation:
    """Test DNA to protein translation."""
    
    def test_translate_dna_basic(self):
        """Test basic translation."""
        # ATG = M, TTA = L, AAA = K, TAA = stop
        seq = "ATGTTAAAA"
        result = bio.translate_dna(seq)
        assert result["protein"] == "MLK"
    
    def test_translate_dna_with_stop_codon(self):
        """Test translation stopping at stop codon."""
        seq = "ATGAAATAATTAG"  # M, K, stop
        result = bio.translate_dna(seq)
        assert result["status"] == "complete"
        assert result["stop_position"] is not None
    
    def test_translate_dna_no_stop_codon(self):
        """Test translation without stop codon."""
        seq = "ATGAAACCC"  # M, K, P
        result = bio.translate_dna(seq)
        assert result["status"] == "no_stop_codon"
    
    def test_translate_all_frames(self):
        """Test translation of all three frames."""
        seq = "ATGATGATG"
        result = bio.translate_all_frames(seq)
        assert len(result) == 3
        assert "Frame +1" in result
        assert "Frame +2" in result
        assert "Frame +3" in result


class TestMutationDetection:
    """Test mutation detection."""
    
    def test_detect_mutations_identical(self):
        """Test mutation detection with identical sequences."""
        query = "ATGCATGC"
        reference = "ATGCATGC"
        result = bio.detect_mutations(query, reference)
        assert result["total_mutations"] == 0
        assert result["identity_percent"] == 100.0
    
    def test_detect_mutations_single_substitution(self):
        """Test detection of single mutation."""
        query = "ATGCATGC"
        reference = "ATGAATGC"  # C -> A at position 3
        result = bio.detect_mutations(query, reference)
        assert result["total_mutations"] == 1
        assert result["mutations"][0]["position"] == 4
    
    def test_detect_mutations_transition(self):
        """Test classification of transition mutation."""
        # A -> G is a purine-to-purine transition
        result = bio._classify_mutation("A", "G")
        assert result == "transition"
    
    def test_detect_mutations_transversion(self):
        """Test classification of transversion mutation."""
        # A -> T is a purine-to-pyrimidine transversion
        result = bio._classify_mutation("A", "T")
        assert result == "transversion"


class TestComplementarySequences:
    """Test complementary and reverse complement sequences."""
    
    def test_complement_basic(self):
        """Test complementary sequence generation."""
        seq = "ATGC"
        comp = bio.complement(seq)
        assert comp == "TACG"
    
    def test_reverse_complement(self):
        """Test reverse complement generation."""
        seq = "ATGC"
        rev_comp = bio.reverse_complement(seq)
        assert rev_comp == "GCAT"


class TestSequenceStatistics:
    """Test aggregate sequence statistics."""
    
    def test_sequence_statistics_basic(self):
        """Test basic sequence statistics."""
        seq = "ATGCATGC"
        stats = bio.sequence_statistics(seq)
        assert stats["length"] == 8
        assert "gc_content" in stats
        assert "at_content" in stats
        assert "has_start_codon" in stats

    def test_sequence_statistics_with_start_codon(self):
        """Test statistics for sequence with ATG start."""
        seq = "ATGCATGC"
        stats = bio.sequence_statistics(seq)
        assert stats["has_start_codon"] is True


class TestSequenceLoader:
    """Test FASTA parsing and header metadata extraction."""

    def test_parse_fasta_header_metadata(self):
        fasta = ">geneX | GC=50% | trait=drought\nATGCGC"
        records = loader.parse_fasta(fasta)
        assert len(records) == 1
        assert records[0]["header"] == "geneX | GC=50% | trait=drought"
        assert records[0]["sequence"] == "ATGCGC"
        assert records[0]["metadata"]["name"] == "geneX"
        assert records[0]["metadata"]["gc"] == "50%"
        assert records[0]["metadata"]["trait"] == "drought"

    def test_parse_fasta_without_header(self):
        records = loader.parse_fasta("ATGCATGC")
        assert len(records) == 1
        assert records[0]["header"] == "Sequence 1"
        assert records[0]["sequence"] == "ATGCATGC"
        assert records[0]["metadata"] == {}


class TestMotifSearch:
    """Test motif detection."""
    
    def test_find_motifs_tata_box(self):
        """Test TATA-box detection."""
        seq = "GCTAGCTATAAATAGCTAG"  # Contains TATAAA
        motifs = bio.find_motifs(seq)
        assert "TATA-box" in [m["name"] for m in motifs]
    
    def test_find_motifs_empty_sequence(self):
        """Test motif search on empty sequence."""
        motifs = bio.find_motifs("")
        assert motifs == []


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
