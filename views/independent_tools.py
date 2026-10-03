"""
views/independent_tools.py
---------------------------
Independent analysis tools — usable without running the main sequence
pipeline on the home page. Registered as its own page via st.navigation
in app.py, so it's always reachable regardless of whether a main
analysis is in progress.
"""

import json
import os
import subprocess
import streamlit as st

import bioinformatics as bio
import alignment_engine as aln
import alignment_exports as align_export
import external_tools
import msa_analysis
import visualization as viz
import trait_research as tr
import config
from i18n import translate, language_selector

STAR_MSA_ENGINE = "Star MSA (fast / approximate)"
MSA_ENGINES = ["MAFFT", "MUSCLE", "ClustalW", STAR_MSA_ENGINE]


def _hydrate_uploaded_fasta(uploaded_file, target_key: str, label: str = "Uploaded FASTA") -> None:
    if uploaded_file is None:
        return
    try:
        text = uploaded_file.read().decode("utf-8", errors="ignore").strip()
    except Exception:
        text = ""
    if text:
        st.session_state[target_key] = text
        st.success(f"{label}: {uploaded_file.name}")


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
    st.info(translate('ui.alignment_instructions', default="Use the MSA section for 2 or more homologous sequences. Use the pairwise section for exactly 2 sequences. Global alignment compares sequences end-to-end; local alignment finds the best shared region. These are sequence-similarity measurements, not proof of gene identity or biological function."))
    if st.button(translate('ui.load_alignment_example', default="Load alignment example"), key="independent_msa_example"):
        st.session_state["independent_msa_input"] = (
            ">Reference\nATGCGATCGATC\n"
            ">Variant_A\nATGCGATCAATC\n"
            ">Variant_B\nATGCGATCGGTC"
        )
    uploaded_msa = st.file_uploader(
        translate('ui.upload_fasta_file', default="Upload FASTA file"),
        type=["fasta", "fa", "txt"],
        key="independent_msa_upload",
    )
    if uploaded_msa is not None:
        _hydrate_uploaded_fasta(uploaded_msa, "independent_msa_input", translate('ui.uploaded_fasta', default="Uploaded FASTA"))
    msa_input = st.text_area(translate('ui.msa_input_hint', default="Paste multiple FASTA sequences or one per line:"), height=160, key="independent_msa_input")
    msa_type = st.selectbox(
        translate('ui.msa_sequence_type', default="MSA sequence type"),
        ["Auto", "DNA", "Protein"],
        key="independent_msa_type",
        format_func=lambda value: {
            "Auto": translate('ui.auto_detect', default="Auto"),
            "DNA": translate('ui.dna', default="DNA"),
            "Protein": translate('ui.protein', default="Protein"),
        }[value],
    )
    msa_engine = st.selectbox(
        translate('ui.msa_engine', default="MSA engine"),
        MSA_ENGINES,
        key="independent_msa_engine",
        help=translate('ui.msa_engine_help', default="MAFFT, MUSCLE and ClustalW run in WSL2. Star MSA is a fast approximate fallback, not the default publication aligner."),
        format_func=lambda value: {
            STAR_MSA_ENGINE: translate('ui.internal_star_msa', default="Star MSA (fast / approximate)"),
            "MAFFT": "MAFFT",
            "MUSCLE": "MUSCLE",
            "ClustalW": "ClustalW",
        }[value],
    )
    msa_reference = 1
    msa_matrix = "BLOSUM62"
    msa_gap_open = -10
    msa_gap_extend = -1
    mafft_strategy = "auto"
    mafft_op = 1.53
    mafft_ep = 0.123
    mafft_thread = 0
    mafft_adjust = False
    muscle_mode = "align"
    if msa_engine == "MAFFT":
        mafft_strategy = st.selectbox(
            translate("ui.mafft_strategy", default="MAFFT strategy"),
            list(external_tools.MAFFT_STRATEGIES),
            key="independent_mafft_strategy",
        )
        mafft_cols = st.columns(3)
        mafft_op = mafft_cols[0].number_input(translate("ui.mafft_op", default="MAFFT --op"), value=1.53, step=0.1, key="independent_mafft_op")
        mafft_ep = mafft_cols[1].number_input(translate("ui.mafft_ep", default="MAFFT --ep"), value=0.123, step=0.01, key="independent_mafft_ep")
        mafft_thread = int(mafft_cols[2].number_input(translate("ui.mafft_thread", default="MAFFT --thread"), value=0, step=1, key="independent_mafft_thread"))
        mafft_adjust = st.checkbox(translate("ui.mafft_adjustdirection", default="MAFFT --adjustdirection"), key="independent_mafft_adjust")
    elif msa_engine == "MUSCLE":
        muscle_mode = st.selectbox(
            translate("ui.muscle_mode", default="MUSCLE 5 mode"),
            ["align", "super5"],
            key="independent_muscle_mode",
        )
    elif msa_engine == STAR_MSA_ENGINE:
        st.caption(translate("ui.star_msa_fallback_help", default="Fast / approximate internal engine. Use MAFFT or MUSCLE for a tree or a paper."))
        msa_reference = st.number_input(translate('ui.star_msa_reference', default="Star-MSA reference sequence"), min_value=1, value=1, step=1, key="independent_msa_reference")
        msa_gap_open = st.number_input(translate('ui.msa_gap_open_penalty', default="MSA gap-open penalty"), value=-10, step=1, key="independent_msa_gap_open")
        msa_gap_extend = st.number_input(translate('ui.msa_gap_extension_penalty', default="MSA gap-extension penalty"), value=-1, step=1, key="independent_msa_gap_extend")
    if msa_type == "Protein":
        msa_matrix = st.selectbox(
            translate('ui.protein_substitution_matrix', default="Protein substitution matrix"),
            ["default", "BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM30", "PAM70", "PAM250"],
            index=2,
            key="independent_msa_matrix",
            help=translate('ui.protein_substitution_matrix_help', default="Used only for protein alignments."),
        )
    if st.button(translate('ui.run_msa', default="Run MSA"), key="independent_msa_run") and msa_input:
        from core_engines.alignment_engine import star_alignment
        from sequence_loader import parse_fasta
        records = parse_fasta(msa_input)
        sequences = [record["sequence"] for record in records]
        labels = [record.get("header", f"Seq{i + 1}") for i, record in enumerate(records)]
        if len(sequences) < 2:
            st.warning(translate('ui.pairwise_missing', default="Provide at least 2 sequences for MSA."))
        else:
            detected_types = [bio.detect_sequence_type(sequence) for sequence in sequences]
            seq_type = detected_types[0] if msa_type == "Auto" else msa_type.lower()
            if any(detected != seq_type for detected in detected_types):
                st.error(translate('ui.msa_type_mismatch', default="All MSA sequences must have the same type (DNA or protein)."))
            elif msa_engine == STAR_MSA_ENGINE and msa_reference > len(sequences):
                st.error(translate('ui.msa_reference_missing', default="The selected MSA reference sequence does not exist."))
            else:
                if msa_engine == STAR_MSA_ENGINE:
                    ordered_sequences = sequences[msa_reference - 1:] + sequences[:msa_reference - 1]
                    ordered_labels = labels[msa_reference - 1:] + labels[:msa_reference - 1]
                else:
                    ordered_sequences = sequences
                    ordered_labels = labels
                issues = msa_analysis.qc_unaligned(ordered_sequences, ordered_labels, seq_type)
                blocking = [issue for issue in issues if issue["level"] == "error"]
                if blocking:
                    for issue in blocking:
                        st.error(issue["message"])
                    st.stop()
                for issue in issues:
                    st.warning(issue["message"])
                estimated_cells = sum(
                    len(ordered_sequences[0]) * len(sequence)
                    for sequence in ordered_sequences[1:]
                )
                if max(map(len, ordered_sequences)) > config.MAX_ALIGNMENT_SEQUENCE_LENGTH:
                    st.error(
                        translate('ui.msa_length_limit', default="MSA sequence length exceeds the alignment limit of {limit:,} characters.").format(limit=config.MAX_ALIGNMENT_SEQUENCE_LENGTH)
                    )
                elif estimated_cells > config.MAX_ALIGNMENT_CELL_BUDGET:
                    st.error(translate('ui.msa_budget_limit', default="This MSA exceeds the configured alignment budget. Use fewer or shorter sequences."))
                else:
                    msa_options = {
                        "strategy": mafft_strategy,
                        "op": mafft_op,
                        "ep": mafft_ep,
                        "thread": mafft_thread or None,
                        "adjustdirection": mafft_adjust,
                        "muscle_mode": muscle_mode,
                    }
                    if msa_engine == STAR_MSA_ENGINE:
                        result = star_alignment(
                            ordered_sequences,
                            seq_type=seq_type,
                            gap_open=msa_gap_open,
                            gap_extend=msa_gap_extend,
                            matrix_name=msa_matrix if seq_type == "protein" else "default",
                        )
                        result["command"] = "star_alignment"
                        result["engine_version"] = "internal"
                        result["parameters"] = {
                            "gap_open": msa_gap_open,
                            "gap_extend": msa_gap_extend,
                            "reference_index": msa_reference,
                        }
                    else:
                        try:
                            if not external_tools.tool_status().get(msa_engine, False):
                                raise RuntimeError(f"{msa_engine} is not available in WSL2.")
                            result = external_tools.run_external_msa(
                                ordered_sequences,
                                ordered_labels,
                                msa_engine,
                                options=msa_options if msa_engine != "ClustalW" else {},
                            )
                        except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
                            st.error(f"{msa_engine} failed: {error}")
                            st.stop()
                    aligned_sequences = result.get("aligned_sequences", [])
                    profile = aln.consensus_profile(aligned_sequences, seq_type=seq_type)
                    result["num_sequences"] = len(aligned_sequences)
                    result["alignment_length"] = len(aligned_sequences[0]) if aligned_sequences else 0
                    result["conservation_score"] = profile["conservation_score"]
                    metadata = align_export.reproducibility_metadata(
                        ordered_sequences,
                        ordered_labels,
                        result.get("algorithm", STAR_MSA_ENGINE),
                        seq_type,
                        msa_matrix if seq_type == "protein" else "DNA",
                        msa_gap_open,
                        msa_gap_extend,
                        command=str(result.get("command", "")),
                        engine_version=str(result.get("engine_version", "")),
                        extra_parameters=dict(result.get("parameters") or msa_options),
                        aligned_sequences=aligned_sequences,
                    )
                    st.session_state["independent_msa_artifact"] = {
                        "unaligned": ordered_sequences,
                        "labels": ordered_labels,
                        "aligned": aligned_sequences,
                        "seq_type": seq_type,
                        "engine": msa_engine,
                        "result": result,
                        "metadata": metadata,
                        "pre_qc": issues,
                    }
                    st.rerun()

    artifact = st.session_state.get("independent_msa_artifact")
    if artifact:
        aligned_sequences = artifact["aligned"]
        labels = artifact["labels"]
        seq_type = artifact["seq_type"]
        metadata = artifact["metadata"]
        profile = aln.consensus_profile(aligned_sequences, seq_type=seq_type)
        alignment_length = len(aligned_sequences[0]) if aligned_sequences else 1
        if artifact["engine"] == STAR_MSA_ENGINE:
            st.warning(
                translate('ui.msa_reference_warning', default="Star MSA is a fast approximate fallback. Changing the reference can change gap placement.")
            )
        st.success(
            translate('ui.msa_complete_summary', default="MSA complete — {seqs} sequences, {cols} aligned columns, {score:.1f}% conservation").format(
                seqs=len(aligned_sequences),
                cols=alignment_length,
                score=profile['conservation_score'],
            )
        )
        counts = profile["site_counts"]
        st.caption(
            translate(
                "ui.msa_site_classes",
                default="Conserved: {conserved} | Variable: {variable} | Parsimony-informative: {pi} | Singletons: {singletons} | Gapped: {gapped}",
            ).format(
                conserved=counts.get("conserved", 0),
                variable=counts.get("variable", 0),
                pi=counts.get("parsimony_informative", 0),
                singletons=counts.get("singletons", 0),
                gapped=counts.get("gapped", 0),
            )
        )
        window_start = st.number_input(
            translate('ui.msa_window_start', default="MEGA-like alignment view: start column"),
            min_value=1,
            max_value=max(1, alignment_length),
            value=1,
            step=1,
            key="independent_msa_window_start",
        )
        window_width = st.slider(
            translate('ui.msa_window_visible', default="Visible alignment columns"),
            20,
            min(200, max(20, alignment_length)),
            min(60, max(20, alignment_length)),
            key="independent_msa_window_width",
        )
        window_end = min(alignment_length, window_start + window_width - 1)
        visible_sequences = [sequence[window_start - 1:window_end] for sequence in aligned_sequences]
        st.caption(
            translate('ui.msa_window_caption', default="Showing aligned columns {start}–{end} of {total}. The complete alignment remains available in the exports below.").format(start=window_start, end=window_end, total=alignment_length)
        )
        st.plotly_chart(
            viz.plot_msa_table(visible_sequences, labels=labels),
            width="stretch",
            key="independent_msa_window_plot",
        )
        st.markdown(f"**{translate('ui.consensus_iupac', default='IUPAC consensus')}**")
        st.code(profile.get("iupac_consensus") or profile["consensus"], language=None)
        coord_cols = st.columns(3)
        ref_seq = int(coord_cols[0].number_input(translate("ui.msa_ref_sequence", default="Reference sequence"), min_value=1, max_value=len(labels), value=1, key="independent_msa_ref_seq"))
        jump_col = int(coord_cols[1].number_input(translate("ui.msa_jump_column", default="Alignment column"), min_value=1, max_value=max(1, alignment_length), value=window_start, key="independent_msa_jump_col"))
        residue_at = msa_analysis.alignment_column_to_residue(aligned_sequences[ref_seq - 1], jump_col)
        coord_cols[2].metric(
            translate("ui.msa_residue_coordinate", default="Residue in reference"),
            residue_at if residue_at is not None else "—",
        )
        jump_res = st.number_input(translate("ui.msa_jump_residue", default="Jump to residue in reference"), min_value=1, value=1, key="independent_msa_jump_res")
        mapped_col = msa_analysis.residue_to_alignment_column(aligned_sequences[ref_seq - 1], int(jump_res))
        st.caption(
            translate("ui.msa_coordinate_caption", default="{label}: alignment column {column} ↔ residue {residue}.").format(
                label=labels[ref_seq - 1],
                column=mapped_col or "—",
                residue=jump_res,
            )
        )
        st.dataframe(profile["sequence_qc"], hide_index=True, width="stretch")
        trim_mode = st.selectbox(
            translate("ui.msa_trim_mode", default="Trimming"),
            ["none", "gt", "cons", "automated1", "smart-gap", "kpic-smart-gap"],
            key="independent_msa_trim_mode",
        )
        trim_gt = st.slider(translate("ui.msa_trim_gt", default="Minimum residue occupancy (gt)"), 0.1, 1.0, 0.8, 0.05, key="independent_msa_trim_gt")
        codon_trim = st.checkbox(translate("ui.msa_trim_codon", default="Trim as codon triplets"), key="independent_msa_trim_codon")
        display_alignment = aligned_sequences
        trim_result = None
        if trim_mode != "none":
            trim_result = msa_analysis.trim_alignment(aligned_sequences, mode=trim_mode, gt=trim_gt, codon=codon_trim)
            display_alignment = trim_result["aligned_sequences"]
            before = profile["site_counts"]
            after = aln.consensus_profile(display_alignment, seq_type=seq_type)["site_counts"]
            st.info(
                translate(
                    "ui.msa_trim_summary",
                    default="Trim {mode}: {before} → {after} columns (variable {vb} → {va}, PI {pb} → {pa}). Removed columns: {removed}",
                ).format(
                    mode=trim_mode,
                    before=before.get("columns"),
                    after=after.get("columns"),
                    vb=before.get("variable"),
                    va=after.get("variable"),
                    pb=before.get("parsimony_informative"),
                    pa=after.get("parsimony_informative"),
                    removed=", ".join(map(str, trim_result["removed_columns"][:40])) + ("…" if len(trim_result["removed_columns"]) > 40 else ""),
                )
            )
        export_alignment = display_alignment
        export_meta = dict(metadata)
        if trim_result:
            export_meta["trim_steps"] = [{
                "tool": "internal trim",
                "mode": trim_result["mode"],
                "removed_columns": trim_result["removed_columns"],
            }]
        st.markdown(f"**{translate('ui.msa_methods', default='Methods')}**")
        st.write(msa_analysis.methods_paragraph(export_meta))
        if st.button(translate("ui.send_msa_to_phylogeny", default="Send trimmed MSA to Phylogeny"), key="independent_msa_to_phylo"):
            st.session_state["independent_validated_msa"] = {
                "labels": labels,
                "aligned": export_alignment,
                "seq_type": seq_type,
                "metadata": export_meta,
            }
            st.success(translate("ui.msa_linked_to_phylogeny", default="Validated alignment linked to Phylogeny. The tree will not re-align pasted sequences."))
        export_cols = st.columns(4)
        export_cols[0].download_button(translate('ui.export_aligned_fasta', default="Aligned FASTA"), align_export.aligned_fasta(export_alignment, labels, export_meta), "alignment.fasta", "text/plain", key="msa_export_fasta")
        export_cols[1].download_button(translate('ui.export_clustal', default="Clustal ALN"), align_export.clustal(export_alignment, labels, export_meta), "alignment.aln", "text/plain", key="msa_export_clustal")
        export_cols[2].download_button(translate('ui.export_metrics_csv', default="Metrics CSV"), align_export.alignment_metrics_csv(export_alignment, labels, export_meta, seq_type=seq_type), "alignment_metrics.csv", "text/csv", key="msa_export_csv")
        export_cols[3].download_button(translate('ui.export_phylip', default="PHYLIP"), align_export.phylip(export_alignment, labels, export_meta), "alignment.phy", "text/plain", key="msa_export_phylip")
        export_cols2 = st.columns(4)
        export_cols2[0].download_button(translate('ui.export_nexus', default="NEXUS"), align_export.nexus(export_alignment, labels, export_meta, seq_type=seq_type), "alignment.nex", "text/plain", key="msa_export_nexus")
        export_cols2[1].download_button(translate('ui.export_mega', default="MEGA"), align_export.mega(export_alignment, labels, export_meta, seq_type=seq_type), "alignment.meg", "text/plain", key="msa_export_mega")
        export_cols2[2].download_button(translate('ui.export_stockholm', default="Stockholm"), align_export.stockholm(export_alignment, labels, export_meta), "alignment.sto", "text/plain", key="msa_export_sto")
        export_cols2[3].download_button(translate('ui.export_provenance', default="Provenance JSON"), align_export.provenance_json(export_meta), "alignment_provenance.json", "application/json", key="msa_export_prov")
        st.download_button(translate("ui.export_phylip_map", default="PHYLIP name map"), align_export.phylip_name_table(labels), "phylip_names.csv", "text/csv", key="msa_export_phy_map")

    pairwise_left, pairwise_right = st.columns(2)
    st.markdown(f"#### {translate('ui.pairwise_comparison', default='Pairwise comparison of two sequences')}")
    st.caption(translate('ui.pairwise_comparison_help', default="Paste one sequence in each field. Use DNA with DNA or protein with protein. If both cleaned inputs are identical, the 100% result is an intentional self-comparison, not biological evidence."))
    with pairwise_left:
        pairwise_seq1 = st.text_area(translate('ui.sequence_1', default="Sequence 1"), height=80, key="independent_pw1")
    with pairwise_right:
        pairwise_seq2 = st.text_area(translate('ui.sequence_2', default="Sequence 2"), height=80, key="independent_pw2")
    pairwise_type = st.selectbox(
        translate('ui.pairwise_sequence_type', default="Pairwise sequence type"),
        ["Auto", "DNA", "Protein"],
        key="independent_pairwise_type",
        format_func=lambda value: {
            "Auto": translate('ui.auto_detect', default="Auto"),
            "DNA": translate('ui.dna', default="DNA"),
            "Protein": translate('ui.protein', default="Protein"),
        }[value],
    )
    pairwise_matrix = st.selectbox(
        translate('ui.pairwise_protein_matrix', default="Pairwise protein matrix"),
        ["default", "BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM30", "PAM70", "PAM250"],
        index=2,
        key="independent_pairwise_matrix",
    )
    pairwise_gap_open = st.number_input(translate('ui.pairwise_gap_open_penalty', default="Pairwise gap-open penalty"), value=-10, step=1, key="independent_pairwise_gap_open")
    pairwise_gap_extend = st.number_input(translate('ui.pairwise_gap_extension_penalty', default="Pairwise gap-extension penalty"), value=-1, step=1, key="independent_pairwise_gap_extend")
    if st.button(translate('ui.align_pairwise', default="Align pairwise"), key="independent_pw_align"):
        if not pairwise_seq1 or not pairwise_seq2:
            st.warning(translate('ui.pairwise_missing', default="Provide two sequences for pairwise alignment."))
        else:
            from core_engines.alignment_engine import needleman_wunsch, smith_waterman
            seq1 = bio.clean_sequence(pairwise_seq1, sequence_type="protein" if pairwise_type == "Protein" else "dna")
            seq2 = bio.clean_sequence(pairwise_seq2, sequence_type="protein" if pairwise_type == "Protein" else "dna")
            if not seq1 or not seq2:
                st.error(translate('ui.invalid_pairwise_characters', default="Both sequences must contain valid characters for the selected type."))
                st.stop()
            if seq1 == seq2:
                st.warning(
                    translate('ui.identical_pairwise_input', default="The two cleaned inputs are identical ({count:,} characters). This is a self-comparison control, so 100% identity is expected.").format(count=len(seq1))
                )
            seq_type = bio.detect_sequence_type(seq1) if pairwise_type == "Auto" else pairwise_type.lower()
            if pairwise_type == "Auto" and bio.detect_sequence_type(seq2) != seq_type:
                st.error(translate('ui.pairwise_type_mismatch', default="The two sequences must have the same detected type."))
                st.stop()
            if max(len(seq1), len(seq2)) > config.MAX_ALIGNMENT_SEQUENCE_LENGTH:
                st.error(translate('ui.pairwise_length_limit', default="Pairwise alignment exceeds the configured sequence-length limit."))
                st.stop()
            if len(seq1) * len(seq2) > config.MAX_ALIGNMENT_CELL_BUDGET:
                st.error(translate('ui.pairwise_budget_limit', default="Pairwise alignment exceeds the configured alignment budget."))
                st.stop()
            matrix_name = pairwise_matrix if seq_type == "protein" else "default"
            global_result = needleman_wunsch(seq1, seq2, seq_type=seq_type, gap_open=pairwise_gap_open, gap_extend=pairwise_gap_extend, matrix_name=matrix_name)
            local_result = smith_waterman(seq1, seq2, seq_type=seq_type, gap_open=pairwise_gap_open, gap_extend=pairwise_gap_extend, matrix_name=matrix_name)
            st.markdown(f"**{translate('ui.global_alignment', default='Needleman-Wunsch (global)')}**")
            st.code(global_result["seq1_aligned"] + "\n" + global_result["seq2_aligned"])
            global_stats = aln.alignment_statistics(global_result["seq1_aligned"], global_result["seq2_aligned"])
            metrics = st.columns(5)
            metrics[0].metric(translate('ui.identity', default='Identity'), f"{global_stats['identity_percent']:.2f}%")
            metrics[1].metric(translate('ui.matches', default='Matches'), global_stats["matches"])
            metrics[2].metric(translate('ui.mismatches', default='Mismatches'), global_stats["mismatches"])
            metrics[3].metric(translate('ui.gaps', default='Gaps'), global_stats["gaps"])
            metrics[4].metric(translate('ui.score', default='Score'), global_result["alignment_score"])
            st.caption(
                translate('ui.pairwise_coverage_summary', default="Coverage sequence 1: {first:.1f}% | Coverage sequence 2: {second:.1f}% | Identity without gaps: {identity:.2f}%").format(
                    first=global_stats['aligned_columns_without_gaps'] / len(seq1) * 100,
                    second=global_stats['aligned_columns_without_gaps'] / len(seq2) * 100,
                    identity=global_stats['non_gap_identity_percent'],
                )
            )
            st.caption(translate('ui.pairwise_interpretation', default="Interpretation: identity, matches, mismatches and gaps describe this alignment. They do not by themselves confirm homology, annotation or biological function."))
            st.markdown(f"**{translate('ui.local_alignment', default='Smith-Waterman (local)')}**")
            st.code(local_result["seq1_aligned"] + "\n" + local_result["seq2_aligned"])
            local_stats = aln.alignment_statistics(local_result["seq1_aligned"], local_result["seq2_aligned"])
            st.caption(
                translate('ui.local_alignment_summary', default="Local identity: {identity:.2f}% | Local aligned length: {length} | Local score: {score}").format(
                    identity=local_stats['identity_percent'],
                    length=local_stats['aligned_length'],
                    score=local_result['alignment_score'],
                )
            )

with tool_tabs[1]:
    st.markdown(f"#### {translate('ui.distance_matrix_title', default='Compute Pairwise Distance Matrix')}")
    st.caption(translate('ui.distance_matrix_help', default="The matrix is the shared input for the phylogenetic tree. Compute it here, then send it to Phylogeny."))
    uploaded_distance = st.file_uploader(
        translate('ui.upload_distance_matrix_input', default="Upload FASTA or sequence file"),
        type=["fasta", "fa", "txt"],
        key="independent_distance_upload",
    )
    if uploaded_distance is not None:
        _hydrate_uploaded_fasta(uploaded_distance, "independent_distance_input", translate('ui.uploaded_fasta', default="Uploaded FASTA"))
    distance_input = st.text_area(translate('ui.distance_input_hint', default="Paste FASTA or one sequence per line:"), height=160, key="independent_distance_input")
    distance_method = st.selectbox(
        translate('ui.method', default="Method"),
        ["hamming", "jukes_cantor", "kimura", "pam"],
        index=2,
        key="independent_distance_method",
        format_func=lambda value: {
            "hamming": translate('ui.distance_method_hamming', default="Hamming"),
            "jukes_cantor": translate('ui.distance_method_jukes_cantor', default="Jukes-Cantor"),
            "kimura": translate('ui.distance_method_kimura', default="Kimura"),
            "pam": translate('ui.distance_method_pam', default="PAM"),
        }[value],
    )
    st.info({
        "hamming": translate('ui.distance_method_hamming_help', default="Raw fraction of observed differences."),
        "jukes_cantor": translate('ui.distance_method_jukes_cantor_help', default="Corrects nucleotide distances for repeated substitutions."),
        "kimura": translate('ui.distance_method_kimura_help', default="Separates transitions and transversions in a nucleotide model."),
        "pam": translate('ui.distance_method_pam_help', default="Protein distance; use amino-acid sequences, not DNA."),
    }[distance_method])
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
            st.session_state["independent_distance_result"] = result
            st.session_state["independent_distance_names"] = names
            st.session_state["independent_distance_method_used"] = distance_method
            st.session_state["independent_distance_sequences"] = sequences
            st.plotly_chart(viz.plot_distance_heatmap(result["distance_matrix"], names, distance_method), width="stretch")
            st.caption(translate('ui.distance_legend', default="Legend: green means a smaller model distance, yellow an intermediate distance, and red a larger model distance. Hover a cell for the exact value."))
            st.dataframe(frame)
            pairs = [(float(frame.iloc[i, j]), names[i], names[j]) for i in range(len(names)) for j in range(i + 1, len(names))]
            if pairs:
                closest = min(pairs)
                farthest = max(pairs)
                st.success(translate('ui.closest_pair', default="Closest pair: {first} and {second} (distance {distance:.6f}).").format(first=closest[1], second=closest[2], distance=closest[0]))
                st.warning(translate('ui.farthest_pair', default="Most distant pair: {first} and {second} (distance {distance:.6f}).").format(first=farthest[1], second=farthest[2], distance=farthest[0]))
            st.caption(translate('ui.distance_reference_note', default="Distances use one shared star alignment. The first sequence is the alignment reference."))
            if st.button(translate('ui.use_matrix_for_tree', default="Use this matrix to build the tree"), key="independent_use_matrix_for_tree"):
                st.session_state["independent_phylogeny_matrix_ready"] = True
                st.success(translate('ui.matrix_linked_to_phylogeny', default="Matrix linked to Phylogeny. Open the Phylogeny tab and build the tree."))
            st.download_button(translate('ui.download_csv', default="Download CSV"), frame.to_csv().encode("utf-8"), file_name="distance_matrix.csv", key="independent_dm_download")

with tool_tabs[2]:
    st.markdown(f"#### {translate('ui.phylogeny_title', default='Build Phylogenetic Tree')}")
    linked_result = st.session_state.get("independent_distance_result")
    linked_method = st.session_state.get("independent_distance_method_used")
    validated_msa = st.session_state.get("independent_validated_msa")
    if validated_msa:
        st.success(
            translate(
                "ui.using_validated_msa",
                default="Using the validated/trimmed MSA from Alignments ({n} sequences, {cols} columns). IQ-TREE will not re-align pasted sequences.",
            ).format(n=len(validated_msa["aligned"]), cols=len(validated_msa["aligned"][0]) if validated_msa["aligned"] else 0)
        )
    uploaded_phylogeny = st.file_uploader(
        translate('ui.upload_phylogeny_file', default="Upload FASTA or sequence file"),
        type=["fasta", "fa", "txt"],
        key="independent_phylogeny_upload",
    )
    if uploaded_phylogeny is not None:
        _hydrate_uploaded_fasta(uploaded_phylogeny, "independent_phylogeny_input", translate('ui.uploaded_fasta', default="Uploaded FASTA"))
    phylogeny_input = st.text_area(translate('ui.phylogeny_input_hint', default="Paste sequences for phylogeny (FASTA or lines):"), height=160, key="independent_phylogeny_input")
    if linked_result and st.session_state.get("independent_phylogeny_matrix_ready"):
        st.success(translate('ui.using_linked_matrix', default="Using the linked {method} distance matrix from Distance Matrix.").format(method=linked_method))
    elif not phylogeny_input.strip() and not validated_msa:
        st.info(translate('ui.phylogeny_input_missing', default="Paste FASTA sequences above or link a distance matrix from the Distance Matrix tab."))
    phylogeny_engine = st.selectbox(
        translate('ui.phylogeny_engine', default="Phylogeny engine"),
        ["Internal distance tree", "IQ-TREE (ModelFinder + bootstrap)"],
        key="independent_phylogeny_engine",
        format_func=lambda value: {
            "Internal distance tree": translate('ui.internal_distance_tree', default="Internal distance tree"),
            "IQ-TREE (ModelFinder + bootstrap)": translate('ui.iqtree_model_bootstrap', default="IQ-TREE (ModelFinder + bootstrap)"),
        }[value],
    )
    phylogeny_method = "upgma"
    if phylogeny_engine == "Internal distance tree":
        phylogeny_method = st.selectbox(
            translate('ui.tree_algorithm', default="Tree algorithm"),
            ["upgma", "neighbor_joining"],
            key="independent_phylogeny_method",
            format_func=lambda value: {
                "upgma": translate('ui.upgma_method', default="UPGMA"),
                "neighbor_joining": translate('ui.neighbor_joining_method', default="Neighbor-Joining"),
            }[value],
        )
    bootstrap_replicates = 1000
    if phylogeny_engine == "IQ-TREE (ModelFinder + bootstrap)":
        bootstrap_replicates = st.select_slider(
            translate('ui.bootstrap_replicates', default="Bootstrap replicates"),
            options=[1000, 2000, 5000],
            value=1000,
            key="independent_bootstrap_replicates",
        )
    tree_theme = st.radio(
        translate('ui.tree_appearance', default="Tree appearance"),
        ["Dark", "Light"],
        horizontal=True,
        key="independent_tree_theme",
        format_func=lambda value: {
            "Dark": translate('ui.dark_theme', default="Dark"),
            "Light": translate('ui.light_theme', default="Light"),
        }[value],
    ).lower()
    if st.button(translate('ui.build_tree', default="Build Tree"), key="independent_build_tree"):
        from sequence_loader import parse_fasta
        from core_engines.distance_engine import distance_matrix
        from core_engines.phylogeny_engine import upgma, neighbor_joining
        import numpy as np

        if phylogeny_engine == "IQ-TREE (ModelFinder + bootstrap)":
            validated_msa = st.session_state.get("independent_validated_msa")
            if validated_msa:
                iq_sequences = [
                    {"name": label, "sequence": sequence, "aligned": True}
                    for label, sequence in zip(validated_msa["labels"], validated_msa["aligned"])
                ]
            elif linked_result and st.session_state.get("independent_phylogeny_matrix_ready"):
                iq_sequences = st.session_state.get("independent_distance_sequences", [])
            else:
                iq_records = parse_fasta(phylogeny_input)
                iq_sequences = [
                    {"name": record.get("header", f"Seq{i + 1}"), "sequence": record["sequence"]}
                    for i, record in enumerate(iq_records)
                ]
            if len(iq_sequences) < 3:
                st.warning(translate('ui.iqtree_min_sequences', default="IQ-TREE requires at least three sequences."))
            else:
                iq_ids = [f"S{index + 1:04d}" for index in range(len(iq_sequences))]
                labels_by_id = {
                    sequence_id: item["name"]
                    for sequence_id, item in zip(iq_ids, iq_sequences)
                }
                try:
                    if validated_msa:
                        spinner = translate("ui.iqtree_on_msa_spinner", default="IQ-TREE model selection on the validated MSA...")
                    else:
                        spinner = translate('ui.iqtree_spinner', default="MAFFT alignment and IQ-TREE model selection in progress...")
                    with st.spinner(spinner):
                        if validated_msa:
                            iq_msa = {
                                "aligned_sequences": [item["sequence"] for item in iq_sequences],
                                "labels": iq_ids,
                            }
                            aligner_label = str((validated_msa.get("metadata") or {}).get("engine") or (validated_msa.get("metadata") or {}).get("algorithm") or "validated MSA")
                        else:
                            iq_msa = external_tools.run_external_msa(
                                [item["sequence"] for item in iq_sequences],
                                iq_ids,
                                "MAFFT",
                            )
                            aligner_label = "MAFFT (--auto)"
                        iq_result = external_tools.run_iqtree(
                            iq_msa["aligned_sequences"],
                            iq_msa["labels"],
                            bootstrap=bootstrap_replicates,
                        )
                    exported_newick = external_tools.restore_newick_labels(
                        iq_result["newick"],
                        labels_by_id,
                    )
                    exported_report = external_tools.restore_iqtree_report_labels(
                        iq_result["report"],
                        labels_by_id,
                    )
                    st.success(
                        translate('ui.iqtree_complete_summary', default="IQ-TREE complete — model: {model} — UFBoot: {bootstrap} replicates").format(
                            model=iq_result['model'],
                            bootstrap=iq_result['bootstrap'],
                        )
                    )
                    from phylo_view import render_phylo_result

                    render_phylo_result(
                        exported_newick,
                        method="iqtree",
                        meta={
                            "model": iq_result["model"],
                            "bootstrap": iq_result["bootstrap"],
                            "n_sites": len(iq_msa["aligned_sequences"][0]),
                            "aligner": aligner_label,
                            "engine": "IQ-TREE 2",
                        },
                        theme=tree_theme,
                        key="independent_iqtree",
                    )
                    st.code(exported_newick, language="text")
                    st.download_button(
                        translate('ui.download_iqtree_newick', default="Download IQ-TREE Newick"),
                        exported_newick,
                        file_name="iqtree_bootstrap.treefile",
                        mime="text/plain",
                        key="independent_iqtree_newick_download",
                    )
                    st.download_button(
                        translate('ui.download_iqtree_report', default="Download IQ-TREE report"),
                        exported_report,
                        file_name="iqtree_model_report.txt",
                        mime="text/plain",
                        key="independent_iqtree_report_download",
                    )
                except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
                    st.error(translate('ui.iqtree_failed', default="IQ-TREE failed: {error}").format(error=error))
            st.stop()

        if linked_result and st.session_state.get("independent_phylogeny_matrix_ready"):
            distances = linked_result
            names = st.session_state["independent_distance_names"]
        else:
            records = parse_fasta(phylogeny_input)
            sequences = [{"name": record.get("header", f"Seq{i + 1}"), "sequence": record["sequence"]} for i, record in enumerate(records)]
            if len(sequences) < 2:
                st.warning(translate("ui.phylogeny_missing", default="Provide at least 2 sequences to build a phylogenetic tree."))
                distances = None
            else:
                distances = distance_matrix(sequences, method="kimura")
            names = distances["sequence_names"] if distances else []
        if distances:
            builder = upgma if phylogeny_method == "upgma" else neighbor_joining
            tree = builder(np.array(distances["distance_matrix"]), names)
            from phylo_view import render_phylo_result

            aligned_sequences = distances.get("aligned_sequences", [])
            render_phylo_result(
                tree["newick"],
                method="upgma" if phylogeny_method == "upgma" else "nj",
                meta={
                    "distance_method": distances.get("method", "kimura"),
                    "n_sites": len(aligned_sequences[0]) if aligned_sequences else None,
                    "aligner": distances.get("alignment_method"),
                    "engine": tree.get("algorithm"),
                },
                theme=tree_theme,
                dist_names=names,
                dist_matrix=distances["distance_matrix"],
                key=f"independent_{phylogeny_method}",
            )
            if tree.get("newick"):
                st.code(tree["newick"])
                st.download_button(translate('ui.download_newick', default="Download Newick"), tree["newick"], file_name="phylogeny_tree.nwk", mime="text/plain", key="independent_newick_download")

with tool_tabs[3]:
    st.markdown(f"#### {translate('ui.protein_biochemical_analysis', default='Protein biochemical analysis')}")
    demo_protein = (
        "MAKQVTLSLVVLASLLALSSAGPVSAQNPAEALKAAGCPSAVWKC"
        "AAAKAGCEAAGLKCLADPKCGAGCKAVCDKAGCDKGSCRKFCS"
    )
    import_col, demo_col = st.columns(2)
    with import_col:
        uploaded_fasta = st.file_uploader(
            translate("ui.import_fasta", default="Import a FASTA file"),
            type=["fasta", "fa", "txt"],
            key="independent_protein_fasta_upload",
        )
    with demo_col:
        if st.button(translate("ui.load_demo_protein", default="Load a demo protein"), key="independent_protein_demo"):
            st.session_state["independent_protein_input"] = demo_protein

    if uploaded_fasta is not None:
        raw_upload = uploaded_fasta.read().decode("utf-8", errors="ignore")
        header, sequence_from_file = bio.parse_fasta_input(raw_upload)
        st.session_state["independent_protein_input"] = sequence_from_file
        if header:
            st.caption(f"FASTA header: {header}")

    protein_input = st.text_area(
        translate("ui.paste_protein_sequence", default="Paste protein sequence:"),
        height=120,
        key="independent_protein_input",
    )
    if protein_input:
        detected_header, sequence_preview = bio.parse_fasta_input(protein_input)
        preview_length = len(sequence_preview.replace("\n", "").replace(" ", ""))
        st.caption(f"{preview_length} {translate('ui.characters', default='characters')}" + (f" - header: {detected_header}" if detected_header else ""))

    if st.button(translate('ui.analyze_protein', default="Analyze protein"), key="independent_protein_analyze"):
        header, sequence_only = bio.parse_fasta_input(protein_input.strip())
        cleaned = bio.clean_sequence(sequence_only, sequence_type="protein")
        valid, message = bio.validate_sequence(cleaned, sequence_type="protein")
        if not valid:
            st.error(message)
        else:
            result = bio.generate_protein_statistics(cleaned)
            categories = result["biochemical_categories"]
            cysteines = result["cysteine_analysis"]
            signal_candidate = result["n_terminal_signal_candidate"]
            st.markdown(f"##### {translate('ui.summary', default='Summary')}")
            metrics = st.columns(4)
            metrics[0].metric(translate("ui.length", default="Length"), f"{result['length']} aa")
            metrics[1].metric(translate("ui.molecular_weight", default="Molecular weight"), f"{result['molecular_weight'] / 1000:.1f} kDa")
            metrics[2].metric(translate("ui.isoelectric_point", default="Estimated pI"), f"{result['isoelectric_point']:.2f}")
            metrics[3].metric(translate("ui.hydrophobicity", default="Avg. hydrophobicity"), f"{result['hydrophobicity']:.2f}")
            metrics = st.columns(4)
            metrics[0].metric(translate("ui.unique_residues", default="Unique residues"), result["unique_residues"])
            metrics[1].metric(translate("ui.cysteines", default="Cysteines"), cysteines["count"])
            metrics[2].metric(translate("ui.charged_residues", default="Charged residues"), result["charged_residues_count"])
            metrics[3].metric(
                translate("ui.instability_index", default="Guruprasad instability index"),
                f"{result['instability_index']:.1f}",
            )
            st.caption(translate(
                "ui.instability_index_help",
                default="Reference DIWV index (Guruprasad et al., 1990; ExPASy ProtParam). Values above 40 suggest a tendency toward instability.",
            ))

            st.markdown(f"##### {translate('ui.detailed_composition', default='Detailed composition')}")
            composition_rows = []
            for amino_acid, count in result["amino_acid_distribution"]["counts"].items():
                if not count:
                    continue
                amino_acid_categories = [
                    category for category, residues in bio.BIOCHEMICAL_CATEGORIES.items()
                    if amino_acid in residues
                ]
                composition_rows.append({
                    "Residue": amino_acid,
                    "Count": count,
                    "Percentage": result["amino_acid_distribution"]["percentages"][amino_acid],
                    "Categories": ", ".join(amino_acid_categories),
                })
            st.dataframe(composition_rows, width="stretch", hide_index=True)

            st.markdown(f"##### {translate('ui.protein_profiles', default='Sequence-level profiles')}")
            profile_col1, profile_col2 = st.columns(2)
            with profile_col1:
                st.plotly_chart(
                    viz.plot_hydrophobicity_profile(result["hydrophobicity_profile"]),
                    width="stretch",
                )
                st.caption(result["hydrophobicity_profile"]["note"])
            with profile_col2:
                st.plotly_chart(
                    viz.plot_charge_profile(result["charge_profile"]),
                    width="stretch",
                )
                st.caption(result["charge_profile"]["note"])

            st.markdown(f"##### {translate('ui.cysteine_motifs', default='Cysteine positions and sequence motifs')}")
            cysteine_rows = [
                {
                    "Cysteine": index + 1,
                    "Position": position,
                    "Distance to next C": cysteines["distances_between_consecutive"][index] if index < len(cysteines["distances_between_consecutive"]) else "-",
                }
                for index, position in enumerate(cysteines["positions"])
            ]
            if cysteine_rows:
                st.dataframe(cysteine_rows, width="stretch", hide_index=True)
            st.caption(cysteines["note"])

            motif_col1, motif_col2 = st.columns(2)
            with motif_col1:
                st.markdown(f"**{translate('ui.glycosylation_candidates', default='N-glycosylation candidates')}**")
                glycosylation_rows = result["protein_motifs"]["n_glycosylation_candidates"]
                st.dataframe(glycosylation_rows or [{"position": translate('ui.none_detected', default='None detected')}], width="stretch", hide_index=True)
            with motif_col2:
                st.markdown(f"**{translate('ui.phosphorylation_candidates', default='Phosphorylation candidates')}**")
                phosphorylation_rows = result["protein_motifs"]["phosphorylation_candidates"]
                st.dataframe(phosphorylation_rows or [{"position": translate('ui.none_detected', default='None detected')}], width="stretch", hide_index=True)
            st.caption(result["protein_motifs"]["note"])

            st.markdown(f"##### {translate('ui.interpretation', default='Interpretation')}")
            interpretation = []
            if result["isoelectric_point"] > 7.5:
                interpretation.append(translate("ui.interp_basic", default="High estimated pI suggests an overall basic protein."))
            elif result["isoelectric_point"] < 6.0:
                interpretation.append(translate("ui.interp_acidic", default="Low estimated pI suggests an overall acidic protein."))
            if result["hydrophobicity"] > 0.5:
                interpretation.append(translate("ui.interp_hydrophobic", default="Elevated average hydrophobicity suggests possible hydrophobic regions."))
            if result["aliphatic_index"] > 80:
                interpretation.append(translate("ui.interp_aliphatic", default="High aliphatic index suggests an aliphatically enriched protein."))
            if result["instability_index"] > 40:
                interpretation.append(translate("ui.interp_unstable", default="The reference instability index is above 40, suggesting a tendency toward instability."))
            else:
                interpretation.append(translate("ui.interp_stable", default="The reference instability index is at or below 40; this does not establish experimental stability."))
            if cysteines["count"] >= 2:
                interpretation.append(translate(
                    "ui.interp_cysteines",
                    default="{count} cysteine(s) detected; at most {pairs} disulfide pair(s) are possible, but bonds are not predicted.",
                    count=cysteines["count"],
                    pairs=cysteines["max_possible_disulfide_pairs"],
                ))
            if signal_candidate["candidate"]:
                interpretation.append(translate(
                    "ui.interp_signal",
                    default="The N-terminal region is strongly hydrophobic and is compatible with a signal peptide; this heuristic does not predict a cleavage site.",
                ))
            glyco_count = len(result["protein_motifs"]["n_glycosylation_candidates"])
            phospho_count = len(result["protein_motifs"]["phosphorylation_candidates"])
            if glyco_count:
                interpretation.append(translate("ui.interp_glycosylation", default="{count} N-glycosylation candidate motif(s) detected; occupancy is not predicted.", count=glyco_count))
            if phospho_count:
                interpretation.append(translate("ui.interp_phosphorylation", default="{count} Ser/Thr/Tyr phosphorylation candidate residue(s) detected; modification is not predicted.", count=phospho_count))
            if not interpretation:
                interpretation.append(translate("ui.interp_none", default="No strongly distinctive biochemical feature was detected from these estimates alone."))
            for line in interpretation:
                st.markdown(f"- {line}")
            st.caption(translate("ui.interp_disclaimer", default="These are sequence-based estimates, not experimental measurements."))

            st.plotly_chart(viz.plot_amino_acid_bar(result["amino_acid_distribution"]), width="stretch")
            import json
            import pandas as pd
            export_result = dict(result)
            export_result["guruprasad_instability_index"] = export_result.pop("instability_index")
            export_result.pop("instability_proxy", None)
            csv_export_rows = [
                {
                    "Residue": "Guruprasad Instability Index (DIWV)",
                    "Count": result["instability_index"],
                    "Percentage": "",
                    "Categories": "Reference index; >40 suggests instability",
                },
                *composition_rows,
            ]
            export_col1, export_col2 = st.columns(2)
            with export_col1:
                st.download_button(
                    translate("ui.export_json", default="Export JSON"),
                    json.dumps(export_result, ensure_ascii=False, indent=2),
                    file_name="protein_analysis.json",
                    mime="application/json",
                    key="independent_protein_export_json",
                )
            with export_col2:
                st.download_button(
                    translate("ui.export_csv", default="Export CSV"),
                    pd.DataFrame(csv_export_rows).to_csv(index=False),
                    file_name="protein_composition.csv",
                    mime="text/csv",
                    key="independent_protein_export_csv",
                )

with tool_tabs[4]:
    tr.render_trait_research_tab("Data/clean/species")