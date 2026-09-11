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
    language_selector(key="independent_lang_selector")
    st.markdown(f"## 🧬 {translate('ui.app_title')}")
    st.markdown("---")
    st.markdown(translate('ui.standalone_tools_message'))

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
    msa_input = st.text_area(translate('ui.msa_input_hint', default="Paste multiple FASTA sequences or one per line:"), height=160, key="independent_msa_input")
    if st.button(translate('ui.run_msa', default="Run MSA"), key="independent_msa_run") and msa_input:
        from core_engines.alignment_engine import star_alignment
        from sequence_loader import parse_fasta
        records = parse_fasta(msa_input)
        sequences = [record["sequence"] for record in records]
        if len(sequences) < 2:
            st.warning(translate('ui.pairwise_missing', default="Provide at least 2 sequences for MSA."))
        else:
            result = star_alignment(sequences, seq_type="dna")
            st.success(f"MSA complete — {result.get('num_sequences')} sequences")
            labels = [record.get("header", f"Seq{i + 1}") for i, record in enumerate(records)]
            st.plotly_chart(viz.plot_msa_table(result.get("aligned_sequences", []), labels=labels), width="stretch")

    pairwise_left, pairwise_right = st.columns(2)
    with pairwise_left:
        pairwise_seq1 = st.text_area(translate('ui.sequence_1', default="Sequence 1"), height=80, key="independent_pw1")
    with pairwise_right:
        pairwise_seq2 = st.text_area(translate('ui.sequence_2', default="Sequence 2"), height=80, key="independent_pw2")
    if st.button(translate('ui.align_pairwise', default="Align pairwise"), key="independent_pw_align"):
        if not pairwise_seq1 or not pairwise_seq2:
            st.warning(translate('ui.pairwise_missing', default="Provide two sequences for pairwise alignment."))
        else:
            from core_engines.alignment_engine import needleman_wunsch, smith_waterman
            global_result = needleman_wunsch(pairwise_seq1.strip(), pairwise_seq2.strip())
            local_result = smith_waterman(pairwise_seq1.strip(), pairwise_seq2.strip())
            st.markdown("**Needleman-Wunsch (global)**")
            st.code(global_result["seq1_aligned"] + "\n" + global_result["seq2_aligned"])
            st.write(f"Score: {global_result['alignment_score']} — Matches: {global_result['match_count']} — Gaps: {global_result['gap_count']}")
            st.markdown("**Smith-Waterman (local)**")
            st.code(local_result["seq1_aligned"] + "\n" + local_result["seq2_aligned"])

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