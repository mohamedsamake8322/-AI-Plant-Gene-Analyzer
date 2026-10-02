"""
test_bioinformatics.py
----------------------
Unit tests for bioinformatics module.
Run with: pytest test_bioinformatics.py -v
"""

import pytest
import bioinformatics as bio
import sequence_loader as loader
import visualization as viz
import external_tools
import numpy as np
import phylogeny_engine as phylo
from Bio import Phylo
from io import StringIO
from unittest.mock import MagicMock, patch
import phylo_view


class TestVisualizationRenderingRegression:
    """Ensure Plotly figures keep valid coordinate arrays after the legend-only fix."""

    def test_upgma_plot_uses_cluster_leaf_order_and_evolutionary_heights(self):
        dendro = {
            "icoord": [[5, 5, 15, 15], [25, 25, 35, 35], [10, 10, 30, 30]],
            "dcoord": [[0, 0.0246, 0.0246, 0], [0, 0.0246, 0.0246, 0], [0.0246, 0.0745, 0.0745, 0.0246]],
            "leaves": [0, 2, 1, 3],
            "color_list": ["#00d9a3", "#00d9a3", "#00d9a3"],
        }

        fig = viz.plot_dendrogram(dendro, labels=["A", "B", "C", "D"], method="K2P")

        label_trace = next(trace for trace in fig.data if getattr(trace, "mode", "") == "text")
        assert list(label_trace.text) == ["A", "C", "B", "D"]
        assert list(label_trace.customdata) == ["A", "C", "B", "D"]
        assert max(value for trace in fig.data for value in (trace.y or []) if value is not None) == pytest.approx(0.03725)
        assert "K2P" in fig.layout.yaxis.title.text

    def test_plot_tree_figures_do_not_emit_none_coordinates(self):
        dendro = {
            "icoord": [[0, 5, 10, 15], [0, 5, 10, 15]],
            "dcoord": [[0, 1, 1, 2], [0, 1, 1, 2]],
            "color_list": ["#00d9a3", "#00d9a3"],
        }
        dendro_fig = viz.plot_dendrogram(dendro, labels=["A", "B"])
        nj_fig = viz.plot_neighbor_joining(
            [{"parent": "Root", "child_1": "A", "child_2": "B", "branch_1": 1.0, "branch_2": 1.0}],
            ["A", "B"],
        )
        newick_fig = viz.plot_newick_tree("((A:0.1,B:0.2)95:0.3,(C:0.15,D:0.25)88:0.35);")

        assert any("A" in (trace.text or []) for trace in newick_fig.data if hasattr(trace, "text"))
        assert len(newick_fig.data) == 2

        for fig in (dendro_fig, nj_fig, newick_fig):
            for trace in fig.data:
                if getattr(trace, "visible", None) == "legendonly":
                    assert list(getattr(trace, "x", [])) == [], "Legend-only trace should be empty"
                    assert list(getattr(trace, "y", [])) == [], "Legend-only trace should be empty"
                else:
                    if hasattr(trace, "x"):
                        assert any(v is not None for v in trace.x), "Data trace has no valid x coordinates"
                    if hasattr(trace, "y"):
                        assert any(v is not None for v in trace.y), "Data trace has no valid y coordinates"

    def test_iqtree_midpoint_plot_and_restored_labels(self):
        original_header = "cqi:110688118 | organism=Chenopodium quinoa"
        original_headers = {
            "S0001": original_header,
            "S0002": "cqi:110684354 | organism=Chenopodium quinoa",
            "S0003": "test copy (A+15) | organism=Chenopodium quinoa",
            "S0004": "test copy (B+15) | organism=Chenopodium quinoa",
        }
        iqtree_newick = "(S0001:0.0000025118,(S0002:0.0000020379,S0004:0.0246483335)100:0.0501366272,S0003:0.0246269447);"

        restored = external_tools.restore_newick_labels(
            iqtree_newick,
            original_headers,
        )
        restored_report = external_tools.restore_iqtree_report_labels(
            "S0001\tS0002\nS00010 should not match S0001",
            {"S0001": original_header, "S0002": "B_110684354"},
        )
        fig = viz.plot_newick_tree(restored)

        terminal_trace = fig.data[1]
        assert set(terminal_trace.customdata) == set(original_headers.values())
        assert set(terminal_trace.text) == {
            "cqi:110688118",
            "cqi:110684354",
            "test copy (A+15)",
            "test copy (B+15)",
        }
        support_annotation = next(annotation for annotation in fig.layout.annotations if annotation.text == "UFBoot 100")
        assert support_annotation.bgcolor
        assert "cqi:110688118" in restored
        assert "test copy (A+15) | organism=Chenopodium quinoa" in restored
        assert "test copy (B+15) | organism=Chenopodium quinoa" in restored
        assert not restored.endswith(":0.0000000000;")
        assert original_header in restored_report
        assert "S0002" not in restored_report
        assert "S00010" in restored_report

    def test_neighbor_joining_keeps_precision_and_expected_terminal_branches(self):
        distances = np.array([
            [0.0, 0.050281, 0.024641, 0.074886],
            [0.050281, 0.0, 0.076584, 0.024635],
            [0.024641, 0.076584, 0.0, 0.096626],
            [0.074886, 0.024635, 0.096626, 0.0],
        ])
        names = ["A", "B", "C", "D"]

        tree = phylo.neighbor_joining(distances, names)
        terminal_branches = {
            child: edge[f"branch_{child_index}"]
            for edge in tree["edges"]
            for child_index, child in ((1, edge["child_1"]), (2, edge["child_2"]))
            if child in names
        }
        root_edge = tree["edges"][-1]
        fig = viz.plot_neighbor_joining(tree["edges"], names)
        tip_trace = next(trace for trace in fig.data if trace.mode == "markers+text")
        tree_trace = next(trace for trace in fig.data if trace.mode == "lines")
        horizontal_segments = [
            (x1, x2, y1, y2)
            for x1, x2, y1, y2 in zip(tree_trace.x, tree_trace.x[1:], tree_trace.y, tree_trace.y[1:])
            if x1 is not None and x2 is not None and y1 is not None and y2 is not None
        ]

        assert terminal_branches["A"] == pytest.approx(0.00030975)
        assert terminal_branches["C"] == pytest.approx(0.02433125)
        assert terminal_branches["B"] == pytest.approx(0.00115575)
        assert terminal_branches["D"] == pytest.approx(0.02347925)
        assert root_edge["branch_1"] + root_edge["branch_2"] == pytest.approx(0.04995625)
        assert list(tip_trace.customdata) == ["A", "C", "B", "D"]
        assert "midpoint-rooted for display" in fig.layout.title.text
        assert any(x1 != x2 and y1 == y2 for x1, x2, y1, y2 in horizontal_segments)

        midpoint_tree = Phylo.read(StringIO(tree["newick"]), "newick")
        midpoint_tree.root_at_midpoint()
        assert sorted(child.branch_length for child in midpoint_tree.root.clades) == pytest.approx(
            [0.024552, 0.025404]
        )

    def test_neighbor_joining_newick_preserves_real_fasta_headers(self):
        names = [
            "cqi:110688118 | organism=Chenopodium quinoa",
            "cqi:110684354 | organism=Chenopodium quinoa",
            "test copy (A+15) | organism=Chenopodium quinoa",
            "test copy (B+15) | organism=Chenopodium quinoa",
        ]
        distances = np.array([
            [0.0, 0.050281, 0.024641, 0.074886],
            [0.050281, 0.0, 0.076584, 0.024635],
            [0.024641, 0.076584, 0.0, 0.096626],
            [0.074886, 0.024635, 0.096626, 0.0],
        ])

        tree = phylo.neighbor_joining(distances, names)
        parsed = Phylo.read(StringIO(tree["newick"]), "newick")

        assert {leaf.name for leaf in parsed.get_terminals()} == set(names)


class TestProfessionalPhylogenyView:
    def test_newick_parser_preserves_quoted_headers_and_ufboot(self):
        newick = "('cqi:1 | organism=Chenopodium quinoa':0.1,(B:0.2,C:0.3)95.2/100:0.4);"

        root = phylo_view.parse_newick(newick)

        assert phylo_view._leaves(root)[0].name == "cqi:1 | organism=Chenopodium quinoa"
        assert root.children[1].support == 100

    def test_neighbor_joining_midpoint_roots_bifurcating_newick(self):
        newick = "((A:0.000310,C:0.024331):0.024978,(B:0.001156,D:0.023479):0.024978);"

        root, was_rerooted = phylo_view.prepare_tree(newick, midpoint=True)
        root_groups = [set(leaf.name for leaf in phylo_view._leaves(child)) for child in root.children]

        assert was_rerooted
        assert sorted(child.length for child in root.children) == pytest.approx([0.024552, 0.025404])
        assert {frozenset(group) for group in root_groups} == {frozenset({"A", "C"}), frozenset({"B", "D"})}

    def test_upgma_figure_uses_node_height_axis_without_rerooting(self):
        newick = "((B:0.012318,D:0.012318):0.024980,(A:0.012320,C:0.012320):0.024977);"

        fig = phylo_view.make_tree_figure(newick, axis_mode="height")

        assert fig.layout.xaxis.visible
        assert fig.layout.xaxis.title.text == "Node height (substitutions / site)"
        assert fig.layout.annotations == ()

    def test_iqtree_midpoint_split_shows_one_support_marker(self):
        newick = "(A:0.0000025118,(B:0.0000020379,D:0.0246483335)100:0.0501366272,C:0.0246269447);"

        fig = phylo_view.make_tree_figure(newick, support_label="UFBoot", midpoint=True)
        support_trace = next(
            trace for trace in fig.data
            if trace.mode == "markers+text" and list(trace.text) == ["100"]
        )

        assert list(support_trace.text) == ["100"]
        assert list(support_trace.marker.color) == ["#2E9E6B"]

    def test_four_taxon_split_warning_mentions_single_split(self):
        ins = phylo_view.tree_insights("((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);")

        warning = next(w for w in ins["warnings"] if "sequences" in w.lower())
        assert "single internal split" in warning.lower()
        assert "3 possible groupings" in warning.lower()
        assert "only one support value" not in warning.lower()

    def test_upgma_main_split_follows_display_order(self):
        newick = "((B:0.012318,D:0.012318):0.024980,(A:0.012320,C:0.012320):0.024977);"

        ins = phylo_view.tree_insights(newick, method="upgma")
        rendered_order = phylo_view.leaf_order(newick)
        expected_top_group = [name for name in rendered_order if name in {"A", "C"}]

        assert ins["split"][0][0] == expected_top_group == ["A", "C"]
        assert ins["split"][0][1] == ["B", "D"]

    def test_distance_heatmap_uses_tree_leaf_order(self):
        names = ["A", "B", "C", "D"]
        matrix = [[float(i == j) for j in range(4)] for i in range(4)]
        order = phylo_view.leaf_order("((B:0.1,D:0.1):0.2,(A:0.1,C:0.1):0.2);")

        fig = phylo_view.make_distance_heatmap(names, matrix, order=order)

        assert list(fig.data[0].x) == order
        assert list(fig.data[0].y) == order

    def test_streamlit_phylogeny_panel_renders_without_errors(self):
        newick = "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);"
        columns = lambda count: [MagicMock() for _ in range(count if isinstance(count, int) else len(count))]

        with (
            patch("streamlit.columns", side_effect=columns),
            patch("streamlit.plotly_chart"),
            patch("streamlit.caption"),
            patch("streamlit.markdown"),
            patch("streamlit.warning"),
            patch("streamlit.code"),
            patch("streamlit.expander"),
        ):
            phylo_view.render_phylo_result(
                newick,
                method="nj",
                meta={"distance_method": "K2P", "n_sites": 619},
                dist_names=["A", "B", "C", "D"],
                dist_matrix=[[0.0] * 4 for _ in range(4)],
            )


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
        """Test translation of forward and reverse frames."""
        seq = "ATGATGATG"
        result = bio.translate_all_frames(seq, include_reverse=True)
        assert len(result) == 6
        assert "Frame +1" in result
        assert "Frame -1" in result

    def test_translation_metadata_and_codon_rows(self):
        """Expose coordinates and stop information for the Translation tab."""
        result = bio.translate_dna("ATGAAATAG", frame=0)
        assert result["protein_with_stop"] == "MK*"
        assert result["stop_position_nt"] == 9
        assert result["remainder_nucleotides"] == 0

        rows = bio.translation_codon_rows("ATGAAATAG", frame=0)
        assert rows[0]["codon"] == "ATG"
        assert rows[0]["start"] == 1
        assert rows[-1]["amino_acid"] == "*"
        assert rows[-1]["is_stop"] is True


class TestProteinProperties:
    """Regression tests for reference protein-property calculations."""

    def test_instability_index_uses_guruprasad_diwv(self):
        sequence = "MKWVTFISLLFLFSSAYS"

        reference_index = bio.calculate_instability_index(sequence)
        proxy_index = bio.calculate_instability_proxy(sequence)

        assert reference_index == pytest.approx(17.5666666667, abs=1e-10)
        assert proxy_index == pytest.approx(27.4705882353, abs=1e-10)
        assert reference_index != pytest.approx(proxy_index)
        assert bio.generate_protein_statistics(sequence)["instability_index"] == pytest.approx(17.57)

    def test_profiles_and_motifs_use_explicit_positions(self):
        sequence = "MNNSTACNPSYTC"

        hydro = bio.hydrophobicity_profile(sequence, window=5)
        charge = bio.charge_profile(sequence, step=1.0)
        motifs = bio.detect_protein_motifs(sequence)

        assert hydro["points"][0]["start"] == 1
        assert hydro["points"][0]["end"] == 5
        assert hydro["points"][0]["position"] == 3.0
        assert charge["points"][0]["ph"] == 0.0
        assert charge["points"][-1]["ph"] == 14.0
        assert [item["position"] for item in motifs["n_glycosylation_candidates"]] == [2, 3]
        assert [item["position"] for item in motifs["phosphorylation_candidates"]] == [4, 5, 10, 11, 12]


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
        assert result["mutations"][0]["position_reference"] == 4
    
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


class TestReadingFrameSummary:
    """Test explicit complete/truncated ORF reporting."""

    def test_separates_complete_and_truncated_orfs(self):
        rows = bio.all_frames_summary("ATGAAATAGATG")
        frame = next(row for row in rows if row["label"] == "+1")

        assert frame["orfs_complete"] == 1
        assert frame["orfs_truncated"] == 1
        assert frame["orf_count"] == 2
        assert frame["longest_orf_length"] == 9


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
