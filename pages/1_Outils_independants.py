"""
pages/1_Outils_independants.py
-------------------------------
Independent analysis tools — usable without running the main sequence
pipeline on the home page. Always reachable from the sidebar navigation,
regardless of whether a main analysis is in progress.
"""

import os
import streamlit as st

import bioinformatics as bio
import visualization as viz
import trait_research as tr
import config

# ─── Page configuration (each Streamlit page sets this independently) ──────────
st.set_page_config(
    page_title=f"{config.PAGE_TITLE} · Independent Tools",
    page_icon=config.PAGE_ICON,
    layout=config.DEFAULT_LAYOUT,
    initial_sidebar_state=config.DEFAULT_SIDEBAR_STATE,
)


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
    st.markdown("## 🧬 AI Plant Gene Analyzer")
    st.markdown("---")
    st.markdown(
        "Standalone tools that work on any sequences you paste here — "
        "no need to run the main analysis first."
    )

st.markdown(
    """
    <div class="hero-panel">
        <h1>🧪 Independent Analysis Tools</h1>
        <p class="hero-subtitle">Alignments · Distance Matrix · Phylogeny · Protein Analysis</p>
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("---")


tool_tabs = st.tabs(["Alignments", "Distance Matrix", "Phylogeny", "Protein Analysis", "Recherche par thème"])

with tool_tabs[0]:
    st.markdown("#### Multiple and pairwise alignment")
    msa_input = st.text_area("Paste multiple FASTA sequences or one per line:", height=160, key="independent_msa_input")
    if st.button("Run MSA", key="independent_msa_run") and msa_input:
        from core_engines.alignment_engine import star_alignment
        from sequence_loader import parse_fasta
        records = parse_fasta(msa_input)
        sequences = [record["sequence"] for record in records]
        if len(sequences) < 2:
            st.warning("Provide at least 2 sequences for MSA.")
        else:
            result = star_alignment(sequences, seq_type="dna")
            st.success(f"MSA complete — {result.get('num_sequences')} sequences")
            labels = [record.get("header", f"Seq{i + 1}") for i, record in enumerate(records)]
            st.plotly_chart(viz.plot_msa_table(result.get("aligned_sequences", []), labels=labels), width="stretch")

    pairwise_left, pairwise_right = st.columns(2)
    with pairwise_left:
        pairwise_seq1 = st.text_area("Sequence 1", height=80, key="independent_pw1")
    with pairwise_right:
        pairwise_seq2 = st.text_area("Sequence 2", height=80, key="independent_pw2")
    if st.button("Align pairwise", key="independent_pw_align"):
        if not pairwise_seq1 or not pairwise_seq2:
            st.warning("Provide two sequences for pairwise alignment.")
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
    st.markdown("#### Compute Pairwise Distance Matrix")
    distance_input = st.text_area("Paste FASTA or one sequence per line:", height=160, key="independent_distance_input")
    distance_method = st.selectbox("Method", ["hamming", "jukes_cantor", "kimura", "pam"], index=2, key="independent_distance_method")
    if st.button("Compute Distance Matrix", key="independent_dm_compute"):
        from sequence_loader import parse_fasta
        from core_engines.distance_engine import distance_matrix
        import pandas as pd
        records = parse_fasta(distance_input)
        sequences = [{"name": record.get("header", f"Seq{i + 1}"), "sequence": record["sequence"]} for i, record in enumerate(records)]
        if len(sequences) < 2:
            st.warning("Provide at least 2 sequences to build distance matrix.")
        else:
            result = distance_matrix(sequences, method=distance_method)
            names = result["sequence_names"]
            frame = pd.DataFrame(result["distance_matrix"], index=names, columns=names)
            st.dataframe(frame)
            st.download_button("Download CSV", frame.to_csv().encode("utf-8"), file_name="distance_matrix.csv", key="independent_dm_download")

with tool_tabs[2]:
    st.markdown("#### Build Phylogenetic Tree")
    phylogeny_input = st.text_area("Paste sequences for phylogeny (FASTA or lines):", height=160, key="independent_phylogeny_input")
    phylogeny_method = st.selectbox("Tree algorithm", ["upgma", "neighbor_joining"], key="independent_phylogeny_method")
    if st.button("Build Tree", key="independent_build_tree"):
        from sequence_loader import parse_fasta
        from core_engines.distance_engine import distance_matrix
        from core_engines.phylogeny_engine import upgma, neighbor_joining
        import numpy as np
        records = parse_fasta(phylogeny_input)
        sequences = [{"name": record.get("header", f"Seq{i + 1}"), "sequence": record["sequence"]} for i, record in enumerate(records)]
        if len(sequences) < 2:
            st.warning("Provide at least 2 sequences for a simple tree.")
        else:
            distances = distance_matrix(sequences, method="kimura")
            builder = upgma if phylogeny_method == "upgma" else neighbor_joining
            tree = builder(np.array(distances["distance_matrix"]), distances["sequence_names"])
            st.write("**Tree metadata**", {"algorithm": tree.get("algorithm")})
            if tree.get("newick"):
                st.code(tree["newick"])
                st.download_button("Download Newick", tree["newick"], file_name="phylogeny_tree.nwk", mime="text/plain", key="independent_newick_download")

with tool_tabs[3]:
    st.markdown("#### Protein biochemical analysis")
    protein_input = st.text_area("Paste protein sequence:", height=120, key="independent_protein_input")
    if st.button("Analyze protein", key="independent_protein_analyze"):
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
