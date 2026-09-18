"""
views/independent_tools.py
---------------------------
Independent analysis tools — usable without running the main sequence
pipeline on the home page. Registered as its own page via st.navigation
in app.py, so it's always reachable regardless of whether a main
analysis is in progress.
"""

import os
import streamlit as st

import bioinformatics as bio
import alignment_engine as aln
import alignment_exports as align_export
import visualization as viz
import trait_research as tr
import config
from i18n import translate, language_selector


# ─── Load custom CSS (same stylesheet as the main page) ────────────────────────
def load_css(css_file: str = "style.css") -> None:
    try:
        if os.path.exists(css_file):
            with open(css_file, "r", encoding="utf-8") as f:
                st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    except Exception:
        pass


load_css()

with st.sidebar:
    st.markdown(f"## 🧬 {translate('ui.app_title')}")
    st.markdown("---")
    st.markdown(translate('ui.standalone_tools_message'))

# Keep the language control beside the Streamlit header actions, rather than
# making it part of the sidebar's scrollable content.
top_spacer, top_language = st.columns([5, 1])
with top_language:
    language_selector(key="top_language_selector")

st.markdown(
    f"""
    <div class="hero-panel">
        <h1>🧪 {translate('ui.independent_tools_title')}</h1>
        <p class="hero-subtitle">{translate('ui.independent_tools_subtitle')}</p>
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("---")


tool_tabs = st.tabs([
    translate('ui.tab_alignments'),
    translate('ui.tab_distance_matrix'),
    translate('ui.tab_phylogeny'),
    translate('ui.tab_protein_analysis'),
    translate('ui.tab_trait_search'),
])

with tool_tabs[0]:
    st.markdown(f"#### {translate('ui.tab_alignments')}")
    st.info(
        "Global alignment compares sequences end to end. Local alignment finds the best shared region. "
        "MSA aligns several sequences in a common reference-guided layout."
    )
    if st.button("Load alignment example", key="independent_msa_example"):
        st.session_state["independent_msa_input"] = (
            ">Reference\nATGCGATCGATC\n"
            ">Variant_A\nATGCGATCAATC\n"
            ">Variant_B\nATGCGATCGGTC"
        )
    msa_input = st.text_area(translate('ui.msa_input_hint', default="Paste multiple FASTA sequences or one per line:"), height=160, key="independent_msa_input")
    msa_type = st.selectbox("MSA sequence type", ["Auto", "DNA", "Protein"], key="independent_msa_type")
    msa_reference = st.number_input("Star-MSA reference sequence", min_value=1, value=1, step=1, key="independent_msa_reference")
    msa_matrix = st.selectbox(
        "Protein substitution matrix",
        ["default", "BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM30", "PAM70", "PAM250"],
        index=2,
        key="independent_msa_matrix",
        help="Used for protein alignments. DNA alignments use the nucleotide matrix.",
    )
    msa_gap_open = st.number_input("MSA gap-open penalty", value=-10, step=1, key="independent_msa_gap_open")
    msa_gap_extend = st.number_input("MSA gap-extension penalty", value=-1, step=1, key="independent_msa_gap_extend")
    if st.button(translate('ui.run_msa', default="Run MSA"), key="independent_msa_run") and msa_input:
        from core_engines.alignment_engine import star_alignment
        from sequence_loader import parse_fasta
        records = parse_fasta(msa_input)
        sequences = [record["sequence"] for record in records]
        if len(sequences) < 2:
            st.warning(translate('ui.pairwise_missing', default="Provide at least 2 sequences for MSA."))
        else:
            detected_types = [bio.detect_sequence_type(sequence) for sequence in sequences]
            seq_type = detected_types[0] if msa_type == "Auto" else msa_type.lower()
            if any(detected != seq_type for detected in detected_types):
                st.error("All MSA sequences must have the same type (DNA or protein).")
            elif msa_reference > len(sequences):
                st.error("The selected MSA reference sequence does not exist.")
            else:
                ordered_sequences = sequences[msa_reference - 1:] + sequences[:msa_reference - 1]
                estimated_cells = sum(
                    len(ordered_sequences[0]) * len(sequence)
                    for sequence in ordered_sequences[1:]
                )
                if max(map(len, ordered_sequences)) > config.MAX_ALIGNMENT_SEQUENCE_LENGTH:
                    st.error(
                        f"MSA sequence length exceeds the alignment limit of "
                        f"{config.MAX_ALIGNMENT_SEQUENCE_LENGTH:,} characters."
                    )
                elif estimated_cells > config.MAX_ALIGNMENT_CELL_BUDGET:
                    st.error("This MSA exceeds the configured alignment budget. Use fewer or shorter sequences.")
                else:
                    result = star_alignment(
                        ordered_sequences,
                        seq_type=seq_type,
                        gap_open=msa_gap_open,
                        gap_extend=msa_gap_extend,
                        matrix_name=msa_matrix if seq_type == "protein" else "default",
                    )
                    st.warning(
                        f"Star MSA uses sequence {msa_reference} as its reference. "
                        "Changing the reference can change gap placement and conservation."
                    )
                    st.success(
                        f"MSA complete — {result.get('num_sequences')} sequences, "
                        f"{result.get('alignment_length')} aligned columns, "
                        f"{result.get('conservation_score', 0):.1f}% conservation"
                    )
                    labels = [record.get("header", f"Seq{i + 1}") for i, record in enumerate(records)]
                    labels = labels[msa_reference - 1:] + labels[:msa_reference - 1]
                    complexity_warnings = [
                        index + 1
                        for index, sequence in enumerate(ordered_sequences)
                        if bio.detect_low_complexity_regions(sequence).get("regions")
                    ]
                    if complexity_warnings:
                        st.warning(
                            "Low-complexity or repetitive regions were detected in sequence(s): "
                            + ", ".join(map(str, complexity_warnings))
                            + ". Conservation and gap placement may be less reliable there."
                        )
                    aligned_sequences = result.get("aligned_sequences", [])
                    st.plotly_chart(viz.plot_msa_table(aligned_sequences, labels=labels), width="stretch")
                    profile = aln.consensus_profile(aligned_sequences)
                    st.markdown("**Consensus sequence**")
                    st.code(profile["consensus"], language=None)
                    variable_positions = profile["variable_columns"]
                    st.caption(
                        f"Variable columns: {len(variable_positions)} | "
                        f"Fully conserved columns: {profile['conservation_score']:.1f}%"
                    )
                    window_start = st.number_input(
                        "MSA display window start",
                        min_value=1,
                        max_value=max(1, result.get("alignment_length", 1)),
                        value=1,
                        step=1,
                        key="independent_msa_window_start",
                    )
                    window_width = st.slider("MSA display window width", 20, 200, 60, key="independent_msa_window_width")
                    window_end = window_start + window_width - 1
                    st.dataframe(
                        [
                            {"Sequence": label, "Aligned sequence": sequence[window_start - 1:window_end]}
                            for label, sequence in zip(labels, aligned_sequences)
                        ],
                        hide_index=True,
                        width="stretch",
                    )
                    metadata = align_export.reproducibility_metadata(
                        ordered_sequences,
                        labels,
                        result.get("algorithm", "Star MSA"),
                        seq_type,
                        msa_matrix if seq_type == "protein" else "DNA",
                        msa_gap_open,
                        msa_gap_extend,
                    )
                    export_cols = st.columns(5)
                    export_cols[0].download_button("Aligned FASTA", align_export.aligned_fasta(aligned_sequences, labels, metadata), "alignment.fasta", "text/plain", key="msa_export_fasta")
                    export_cols[1].download_button("Clustal ALN", align_export.clustal(aligned_sequences, labels, metadata), "alignment.aln", "text/plain", key="msa_export_clustal")
                    export_cols[2].download_button("Metrics CSV", align_export.alignment_metrics_csv(aligned_sequences, labels, metadata), "alignment_metrics.csv", "text/csv", key="msa_export_csv")
                    export_cols[3].download_button("PHYLIP", align_export.phylip(aligned_sequences, labels, metadata), "alignment.phy", "text/plain", key="msa_export_phylip")
                    export_cols[4].download_button("NEXUS", align_export.nexus(aligned_sequences, labels, metadata), "alignment.nex", "text/plain", key="msa_export_nexus")

    pairwise_left, pairwise_right = st.columns(2)
    with pairwise_left:
        pairwise_seq1 = st.text_area(translate('ui.sequence_1', default="Sequence 1"), height=80, key="independent_pw1")
    with pairwise_right:
        pairwise_seq2 = st.text_area(translate('ui.sequence_2', default="Sequence 2"), height=80, key="independent_pw2")
    pairwise_type = st.selectbox("Pairwise sequence type", ["Auto", "DNA", "Protein"], key="independent_pairwise_type")
    pairwise_matrix = st.selectbox("Pairwise protein matrix", ["default", "BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM30", "PAM70", "PAM250"], index=2, key="independent_pairwise_matrix")
    pairwise_gap_open = st.number_input("Pairwise gap-open penalty", value=-10, step=1, key="independent_pairwise_gap_open")
    pairwise_gap_extend = st.number_input("Pairwise gap-extension penalty", value=-1, step=1, key="independent_pairwise_gap_extend")
    if st.button(translate('ui.align_pairwise', default="Align pairwise"), key="independent_pw_align"):
        if not pairwise_seq1 or not pairwise_seq2:
            st.warning(translate('ui.pairwise_missing', default="Provide two sequences for pairwise alignment."))
        else:
            from core_engines.alignment_engine import needleman_wunsch, smith_waterman
            seq1 = bio.clean_sequence(pairwise_seq1, sequence_type="protein" if pairwise_type == "Protein" else "dna")
            seq2 = bio.clean_sequence(pairwise_seq2, sequence_type="protein" if pairwise_type == "Protein" else "dna")
            seq_type = bio.detect_sequence_type(seq1) if pairwise_type == "Auto" else pairwise_type.lower()
            if pairwise_type == "Auto" and bio.detect_sequence_type(seq2) != seq_type:
                st.error("The two sequences must have the same detected type.")
                st.stop()
            if max(len(seq1), len(seq2)) > config.MAX_ALIGNMENT_SEQUENCE_LENGTH:
                st.error("Pairwise alignment exceeds the configured sequence-length limit.")
                st.stop()
            if len(seq1) * len(seq2) > config.MAX_ALIGNMENT_CELL_BUDGET:
                st.error("Pairwise alignment exceeds the configured alignment budget.")
                st.stop()
            matrix_name = pairwise_matrix if seq_type == "protein" else "default"
            global_result = needleman_wunsch(seq1, seq2, seq_type=seq_type, gap_open=pairwise_gap_open, gap_extend=pairwise_gap_extend, matrix_name=matrix_name)
            local_result = smith_waterman(seq1, seq2, seq_type=seq_type, gap_open=pairwise_gap_open, gap_extend=pairwise_gap_extend, matrix_name=matrix_name)
            st.markdown("**Needleman-Wunsch (global)**")
            st.code(global_result["seq1_aligned"] + "\n" + global_result["seq2_aligned"])
            global_stats = aln.alignment_statistics(global_result["seq1_aligned"], global_result["seq2_aligned"])
            metrics = st.columns(5)
            metrics[0].metric("Identity", f"{global_stats['identity_percent']:.2f}%")
            metrics[1].metric("Matches", global_stats["matches"])
            metrics[2].metric("Mismatches", global_stats["mismatches"])
            metrics[3].metric("Gaps", global_stats["gaps"])
            metrics[4].metric("Score", global_result["alignment_score"])
            st.caption(
                f"Coverage sequence 1: {global_stats['aligned_columns_without_gaps'] / len(seq1) * 100:.1f}% | "
                f"Coverage sequence 2: {global_stats['aligned_columns_without_gaps'] / len(seq2) * 100:.1f}%"
            )
            st.markdown("**Smith-Waterman (local)**")
            st.code(local_result["seq1_aligned"] + "\n" + local_result["seq2_aligned"])
            local_stats = aln.alignment_statistics(local_result["seq1_aligned"], local_result["seq2_aligned"])
            st.caption(
                f"Local identity: {local_stats['identity_percent']:.2f}% | "
                f"Local aligned length: {local_stats['aligned_length']} | "
                f"Local score: {local_result['alignment_score']}"
            )

with tool_tabs[1]:
    st.markdown(f"#### {translate('ui.distance_matrix_title', default='Compute Pairwise Distance Matrix')}")
    distance_input = st.text_area(translate('ui.distance_input_hint', default="Paste FASTA or one sequence per line:"), height=160, key="independent_distance_input")
    distance_method = st.selectbox(translate('ui.method', default="Method"), ["hamming", "jukes_cantor", "kimura", "pam"], index=2, key="independent_distance_method")
    if st.button(translate('ui.compute_distance_matrix', default="Compute Distance Matrix"), key="independent_dm_compute"):
        from sequence_loader import parse_fasta
        from core_engines.distance_engine import distance_matrix
        import pandas as pd
        records = parse_fasta(distance_input)
        sequences = [{"name": record.get("header", f"Seq{i + 1}"), "sequence": record["sequence"]} for i, record in enumerate(records)]
        if len(sequences) < 2:
            st.warning(translate('ui.pairwise_missing', default="Provide at least 2 sequences to build distance matrix."))
        else:
            result = distance_matrix(sequences, method=distance_method)
            names = result["sequence_names"]
            frame = pd.DataFrame(result["distance_matrix"], index=names, columns=names)
            st.dataframe(frame)
            st.download_button("Download CSV", frame.to_csv().encode("utf-8"), file_name="distance_matrix.csv", key="independent_dm_download")

with tool_tabs[2]:
    st.markdown(f"#### {translate('ui.phylogeny_title', default='Build Phylogenetic Tree')}")
    phylogeny_input = st.text_area(translate('ui.phylogeny_input_hint', default="Paste sequences for phylogeny (FASTA or lines):"), height=160, key="independent_phylogeny_input")
    phylogeny_method = st.selectbox(translate('ui.tree_algorithm', default="Tree algorithm"), ["upgma", "neighbor_joining"], key="independent_phylogeny_method")
    if st.button(translate('ui.build_tree', default="Build Tree"), key="independent_build_tree"):
        from sequence_loader import parse_fasta
        from core_engines.distance_engine import distance_matrix
        from core_engines.phylogeny_engine import upgma, neighbor_joining
        import numpy as np
        records = parse_fasta(phylogeny_input)
        sequences = [{"name": record.get("header", f"Seq{i + 1}"), "sequence": record["sequence"]} for i, record in enumerate(records)]
        if len(sequences) < 2:
            st.warning(translate('ui.pairwise_missing', default="Provide at least 2 sequences for a simple tree."))
        else:
            distances = distance_matrix(sequences, method="kimura")
            builder = upgma if phylogeny_method == "upgma" else neighbor_joining
            tree = builder(np.array(distances["distance_matrix"]), distances["sequence_names"])
            st.write("**Tree metadata**", {"algorithm": tree.get("algorithm")})
            if tree.get("newick"):
                st.code(tree["newick"])
                st.download_button("Download Newick", tree["newick"], file_name="phylogeny_tree.nwk", mime="text/plain", key="independent_newick_download")

with tool_tabs[3]:
    st.markdown(f"#### {translate('ui.protein_biochemical_analysis', default='Protein biochemical analysis')}")
    protein_input = st.text_area(translate('ui.paste_protein_sequence', default="Paste protein sequence:"), height=120, key="independent_protein_input")
    if st.button(translate('ui.analyze_protein', default="Analyze protein"), key="independent_protein_analyze"):
        cleaned = bio.clean_sequence(protein_input.strip(), sequence_type="protein")
        valid, message = bio.validate_sequence(cleaned, sequence_type="protein")
        if not valid:
            st.error(message)
        else:
            result = bio.generate_protein_statistics(cleaned)
            st.write({"length_aa": result["length"], "molecular_weight": result["molecular_weight"], "isoelectric_point": result["isoelectric_point"], "hydrophobicity": result["hydrophobicity"]})
            st.plotly_chart(viz.plot_amino_acid_bar(result["amino_acid_distribution"]), width="stretch")

with tool_tabs[4]:
    tr.render_trait_research_tab("Data/clean/species")