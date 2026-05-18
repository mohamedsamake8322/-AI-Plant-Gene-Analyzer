"""
app.py
------
AI-Powered Plant Gene Analyzer — Streamlit frontend.

Run with:
    streamlit run app.py
"""

import streamlit as st
import json
import os
import io
import logging

# ── Local modules ──────────────────────────────────────────────────────────────
import bioinformatics as bio
import similarityengine as sim
import aiinterpreter as ai_interp
import visualization as viz
import export_utils as export_util
import sequence_loader as loader
import config

# ─── Configure logging ─────────────────────────────────────────────────────────
logger = config.get_logger(__name__)


# ─── Page configuration ────────────────────────────────────────────────────────
st.set_page_config(
    page_title=config.PAGE_TITLE,
    page_icon=config.PAGE_ICON,
    layout=config.DEFAULT_LAYOUT,
    initial_sidebar_state=config.DEFAULT_SIDEBAR_STATE,
)


# ─── Load custom CSS ────────────────────────────────────────────────────────────
def load_css(css_file: str = "style.css") -> None:
    """Load custom CSS with error handling."""
    try:
        if os.path.exists(css_file):
            with open(css_file, "r", encoding="utf-8") as f:
                st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
            logger.info(f"CSS loaded successfully from {css_file}")
        else:
            logger.warning(f"CSS file not found: {css_file}")
    except Exception as e:
        logger.error(f"Error loading CSS: {e}")
        st.warning("⚠️ Could not load custom styling (CSS file error)")


def analyze_sequence_record(record: dict[str, str], input_type: str, reading_frame: int) -> dict:
    """Analyze a single sequence record and return a structured result."""
    seq_type = input_type
    if seq_type == "Auto detect":
        seq_type = loader.detect_sequence_type(record["sequence"])
        if seq_type == "unknown":
            seq_type = "dna"

    sequence = bio.clean_sequence(record["sequence"], sequence_type="protein" if seq_type == "protein" else "dna")
    is_valid, validation_msg = bio.validate_sequence(sequence, sequence_type=seq_type)
    if not is_valid:
        raise ValueError(validation_msg)

    if seq_type == "protein":
        stats = bio.generate_protein_statistics(sequence)
        dist = bio.amino_acid_distribution(sequence)
        translation = None
        protein_props = bio.protein_properties(sequence)
        orfs = []
    else:
        stats = bio.sequence_statistics(sequence)
        dist = bio.nucleotide_distribution(sequence)
        translation = bio.translate_dna(sequence, frame=reading_frame)
        protein_props = None
        orfs = bio.find_orfs(sequence)

    motifs = bio.find_motifs(sequence)
    similarity_results = []
    best_match = None
    mutation_report = None

    try:
        similarity_results = sim.compare_with_database(sequence, top_n=top_n_matches)
        best_match = similarity_results[0] if similarity_results else None
    except Exception as e:
        logger.warning(f"Database comparison failed: {e}")
        similarity_results = []

    if best_match and db and seq_type == "dna":
        try:
            ref_seq = db[best_match["gene_name"]]["sequence"].upper().replace(" ", "")
            mutation_report = bio.detect_mutations(sequence, ref_seq)
        except Exception as e:
            logger.warning(f"Mutation detection failed: {e}")
            mutation_report = None

    interpretation = {}
    try:
        interpretation = ai_interp.interpret(stats, similarity_results, mutation_report)
    except Exception as e:
        logger.warning(f"AI interpretation failed: {e}")

    return {
        "header": record.get("header", "Sequence"),
        "sequence": sequence,
        "stats": stats,
        "protein_stats": protein_props,
        "dist": dist,
        "translation": translation,
        "motifs": motifs,
        "similarity_results": similarity_results,
        "best_match": best_match,
        "mutation_report": mutation_report,
        "interpretation": interpretation,
        "sequence_type": seq_type,
        "orfs": orfs,
    }


def get_alignment_map(match: dict) -> dict[str, str] | None:
    if isinstance(match.get("alignment"), dict):
        return match["alignment"].get("alignment_map")
    return None


# ─── Load gene database with caching ────────────────────────────────────────────────────────────
@st.cache_data
def load_gene_database_cached(db_path: str = "genes_database.json") -> dict:
    """
    Load gene database with Streamlit caching to improve performance.
    Caches the result to avoid reloading on every rerun.
    """
    try:
        if not os.path.exists(db_path):
            logger.warning(f"Database not found at {db_path}")
            return {}
        
        with open(db_path, "r", encoding="utf-8") as f:
            db = json.load(f)
        
        logger.info(f"Loaded {len(db)} genes from database")
        return db
    
    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing error: {e}")
        st.error("❌ Error parsing genes_database.json")
        return {}
    except Exception as e:
        logger.error(f"Error loading database: {e}")
        st.error(f"❌ Error loading database: {e}")
        return {}


load_css()


# ─── Demo sequences ────────────────────────────────────────────────────────────
DEMO_SEQUENCES: dict[str, dict] = config.DEMO_SEQUENCES


# ─── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🧬 AI Plant Gene Analyzer")
    st.markdown("---")

    st.markdown("### About")
    st.markdown(
        "Analyze plant DNA sequences for:\n"
        "- GC content & nucleotide stats\n"
        "- Gene database similarity\n"
        "- Mutation detection\n"
        "- AI biological interpretation\n"
        "- Agricultural recommendations"
    )
    st.markdown("---")

    st.markdown("### Settings")
    top_n_matches = st.slider(
        "Top database matches to show",
        min_value=1,
        max_value=8,
        value=3,
        help="Number of best-matching genes to display",
    )
    window_size = st.slider(
        "Sliding window (GC profile)",
        min_value=5,
        max_value=60,
        value=20,
        step=5,
        help="Window size (bp) for the GC content profile chart",
    )
    reading_frame = st.selectbox(
        "Reading frame for translation",
        options=[0, 1, 2],
        format_func=lambda x: f"+{x + 1}",
    )
    sequence_input_type = st.selectbox(
        "Input type",
        options=config.SUPPORTED_INPUT_TYPES,
        help="Choose the sequence type or let the app detect it automatically.",
    )
    st.markdown("---")

    st.markdown("### Database")
    db = load_gene_database_cached("genes_database.json")
    if db:
        st.success(f"✅ {len(db)} genes loaded")
        with st.expander("View gene names"):
            for gene_name, info in db.items():
                st.markdown(f"- **{gene_name}** — {info.get('trait', '')}")
    else:
        st.error("❌ No genes available in database")


# ─── Main header ───────────────────────────────────────────────────────────────
st.markdown(
    "<h1>🧬 AI-Powered Plant Gene Analyzer</h1>",
    unsafe_allow_html=True,
)
st.markdown(
    "<p style='text-align:center;color:#a5d6a7;font-family:Courier New;'>"
    "Bioinformatics · AI Interpretation · Agricultural Insights"
    "</p>",
    unsafe_allow_html=True,
)
st.markdown("---")


# ─── Input section ─────────────────────────────────────────────────────────────
col_input, col_demo = st.columns([2, 1])

with col_input:
    st.markdown("### Sequence Input")

    if sequence_input_type != "DNA":
        st.markdown(
            f"This app will accept **{sequence_input_type}** input and adjust analysis accordingly."
        )

    uploaded_file = st.file_uploader(
        "Upload a sequence file (.fasta / .fa / .txt)",
        type=["fasta", "fa", "txt"],
        help="FASTA files starting with '>' header lines are supported.",
    )

    raw_sequence = st.text_area(
        "Or paste your sequence here:",
        height=140,
        placeholder=(
            "Paste raw DNA, protein sequence, or FASTA format here.\n"
            "Example: ATGCGTAGCTAGCGATCGATCGAATTCG..."
        ),
    )

    records: list[dict[str, str]] = []
    analyze_all = False
    selected_index = 0

    if uploaded_file is not None:
        content = uploaded_file.read().decode("utf-8", errors="ignore")
        records = loader.parse_fasta(content)
        st.info(f"File loaded: {uploaded_file.name}")
    elif raw_sequence:
        records = loader.parse_fasta(raw_sequence)
    
    if records:
        if len(records) > 1:
            record_options = [f"{idx + 1}. {r['header']}" for idx, r in enumerate(records)]
            selected_index = st.selectbox(
                "Choose a sequence to analyze",
                options=list(range(len(records))),
                format_func=lambda i: record_options[i],
            )
            analyze_all = st.checkbox(
                "Analyze all sequences in this FASTA input",
                value=False,
                help="If checked, all parsed FASTA records will be analyzed in batch.",
            )
            if analyze_all:
                st.success(f"{len(records)} sequences will be analyzed as a batch.")
            else:
                raw_sequence = records[selected_index]["sequence"]
                st.info(f"Selected: {record_options[selected_index]}")
        elif len(records) == 1:
            raw_sequence = records[0]["sequence"]
    
with col_demo:
    st.markdown("### Quick Demo")
    selected_demo = st.selectbox(
        "Load a demo sequence",
        options=list(DEMO_SEQUENCES.keys()),
        label_visibility="collapsed",
    )
    if selected_demo != "Select a demo…":
        demo = DEMO_SEQUENCES[selected_demo]
        st.markdown(f"*{demo['desc']}*")
        if st.button("Load Demo Sequence", use_container_width=True):
            raw_sequence = demo["seq"]
            st.session_state["loaded_demo"] = demo["seq"]

    if "loaded_demo" in st.session_state and not raw_sequence:
        raw_sequence = st.session_state["loaded_demo"]

analyze_btn = st.button("🔬 Analyze Sequence", use_container_width=True, type="primary")

st.markdown("---")


# ─── Analysis pipeline ─────────────────────────────────────────────────────────
if analyze_btn or (raw_sequence and "last_result" in st.session_state):

    if analyze_btn and not raw_sequence:
        st.warning("⚠️ Please enter or paste a DNA sequence before analyzing.")
        st.stop()

    if analyze_btn and raw_sequence:
        try:
            analysis_targets: list[dict[str, str]] = []
            if records and len(records) > 1:
                if analyze_all:
                    analysis_targets = records
                else:
                    analysis_targets = [records[selected_index]]
            else:
                analysis_targets = [{"header": "Sequence 1", "sequence": raw_sequence}]

            analyzed_results: list[dict] = []
            with st.spinner("🧬 Running bioinformatics analysis…"):
                for idx, record in enumerate(analysis_targets):
                    logger.info(
                        f"Starting analysis for record {idx + 1}/{len(analysis_targets)}: {record.get('header', 'Sequence')}"
                    )
                    analyzed_results.append(
                        analyze_sequence_record(record, sequence_input_type, reading_frame)
                    )

            st.session_state["last_results"] = analyzed_results
            st.session_state["last_result"] = analyzed_results[0] if analyzed_results else None
            logger.info("Analysis session state saved")

        except Exception as e:
            logger.error(f"Unexpected error during analysis: {e}")
            st.error(f"❌ Analysis failed: {e}")
            st.stop()

    last_results = st.session_state.get("last_results")
    if not last_results:
        st.stop()

    result = last_results[0]
    batch_mode = len(last_results) > 1
    sequence = result["sequence"]
    stats = result["stats"]
    protein_stats = result.get("protein_stats")
    dist = result["dist"]
    translation = result["translation"]
    motifs = result["motifs"]
    similarity_results = result["similarity_results"]
    best_match = result["best_match"]
    mutation_report = result["mutation_report"]
    interpretation = result["interpretation"]
    sequence_type = result.get("sequence_type", "dna")

    if batch_mode:
        average_gc = round(sum(r["stats"].get("gc_content", 0) for r in last_results if r["sequence_type"] == "dna") / max(1, sum(1 for r in last_results if r["sequence_type"] == "dna")), 2)
        st.success(f"✅ Batch analysis complete — {len(last_results)} sequences processed.")
        st.markdown(
            f"**Batch summary:** {len(last_results)} sequences, "
            f"average DNA GC content {average_gc}% (protein sequences excluded from GC average)."
        )
        if len(last_results) > 1:
            summary_rows = []
            for idx, item in enumerate(last_results, start=1):
                similarity_value = "—"
                if item["best_match"]:
                    similarity_score = item["best_match"].get("similarity_score", 0.0)
                    similarity_value = f"{similarity_score:.1f}"
                summary_rows.append({
                    "Sequence": item["header"],
                    "Type": item["sequence_type"].upper(),
                    "Length": item["stats"]["length"],
                    "Best match": item["best_match"]["gene_name"] if item["best_match"] else "—",
                    "Similarity (%)": similarity_value,
                })
            st.table(summary_rows)
            st.markdown("#### Batch analysis details")
            for idx, item in enumerate(last_results, start=1):
                with st.expander(f"{idx}. {item['header']} — {item['sequence_type'].upper()} ({item['stats']['length']} { 'aa' if item['sequence_type'] == 'protein' else 'bp'})"):
                    st.markdown(f"- **Best match:** {item['best_match']['gene_name'] if item['best_match'] else '—'}")
                    st.markdown(f"- **Similarity:** {item['best_match']['similarity_score']:.1f}%" if item['best_match'] else "- **Similarity:** —")
                    if item['sequence_type'] == 'dna':
                        st.markdown(f"- **ORFs found:** {len(item['orfs'])}")
                        if item['orfs']:
                            st.markdown(f"- **Longest ORF:** {item['orfs'][0]['length']} bp in frame {item['orfs'][0]['frame']}")
                    else:
                        st.markdown(f"- **Protein weight:** {item['stats'].get('molecular_weight', 'N/A')} Da")
                        st.markdown(f"- **Estimated pI:** {item['stats'].get('isoelectric_point', 'N/A')}")
                        st.markdown(f"- **Hydrophobicity:** {item['stats'].get('hydrophobicity', 'N/A')}")
    else:
        length_unit = "aa" if sequence_type == "protein" else "bp"
        st.success(f"✅ Analysis complete — {stats['length']:,} {length_unit} sequence processed.")
    
    # ── Export Options ─────────────────────────────────────────────────────────
    st.markdown("---")
    st.markdown("### 📥 Export Results")
    export_col1, export_col2, export_col3, export_col4 = st.columns(4)
    
    with export_col1:
        if st.button("📄 Download JSON", use_container_width=True):
            try:
                json_path = export_util.export_results_json(result)
                with open(json_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 JSON Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.json",
                        mime="application/json",
                        use_container_width=True,
                    )
                logger.info(f"JSON export created: {json_path}")
                st.success("✅ JSON exported successfully")
            except Exception as e:
                logger.error(f"JSON export failed: {e}")
                st.error(f"Export failed: {e}")
    
    with export_col2:
        if st.button("📊 Download CSV", use_container_width=True):
            try:
                csv_path = export_util.export_results_csv(result)
                with open(csv_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 CSV Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.csv",
                        mime="text/csv",
                        use_container_width=True,
                    )
                logger.info(f"CSV export created: {csv_path}")
                st.success("✅ CSV exported successfully")
            except Exception as e:
                logger.error(f"CSV export failed: {e}")
                st.error(f"Export failed: {e}")
    
    with export_col3:
        if st.button("🌐 Download HTML", use_container_width=True):
            try:
                html_path = export_util.export_results_html(result)
                with open(html_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 HTML Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.html",
                        mime="text/html",
                        use_container_width=True,
                    )
                logger.info(f"HTML export created: {html_path}")
                st.success("✅ HTML exported successfully")
            except Exception as e:
                logger.error(f"HTML export failed: {e}")
                st.error(f"Export failed: {e}")
    with export_col4:
        if st.button("📑 Download XLSX", use_container_width=True):
            try:
                xlsx_path = export_util.export_results_xlsx(result)
                with open(xlsx_path, "rb") as f:
                    st.download_button(
                        "📥 XLSX Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                    )
                logger.info(f"XLSX export created: {xlsx_path}")
                st.success("✅ XLSX exported successfully")
            except Exception as e:
                logger.error(f"XLSX export failed: {e}")
                st.error(f"Export failed: {e}")

    st.markdown("---")

    # ── KPI Metrics ────────────────────────────────────────────────────────────
    # ── KPI Metrics ───────────────────────────────────────────────────────────
    st.markdown("### Sequence Overview")
    if sequence_type == "protein":
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Length (aa)", f"{stats['length']:,}")
        m2.metric("Unique Residues", f"{stats.get('unique_residues', 'N/A')}")
        m3.metric("Most Abundant", f"{max(dist['counts'], key=dist['counts'].get)}")
        m4.metric(
            "Best Match",
            best_match["gene_name"] if best_match else "—",
            f"{best_match['similarity_score']:.1f}%" if best_match else None,
        )
        m5.metric(
            "Mutations",
            mutation_report["total_mutations"] if mutation_report else "—",
        )
    else:
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Length (bp)", f"{stats['length']:,}")
        m2.metric("GC Content", f"{stats['gc_content']}%")
        m3.metric("AT Content", f"{stats['at_content']}%")
        m4.metric(
            "Best Match",
            best_match["gene_name"] if best_match else "—",
            f"{best_match['similarity_score']:.1f}%" if best_match else None,
        )
        m5.metric(
            "Mutations",
            mutation_report["total_mutations"] if mutation_report else "—",
        )

    st.markdown("---")

    # ── Tabs ───────────────────────────────────────────────────────────────────
    tabs = st.tabs([
        "📊 Statistics",
        "🔍 Similarity",
        "🧪 Mutations",
        "🔬 Translation",
        "🤖 AI Interpretation",
        "📋 Raw Sequence",
        "🧩 Alignments",
        "📐 Distance Matrix",
        "🌳 Phylogeny",
        "🔬 Protein Analysis",
    ])

    # ── Tab 1: Statistics ──────────────────────────────────────────────────────
    with tabs[0]:
        if sequence_type == "protein":
            st.markdown("#### Amino Acid Composition")

            col1, col2 = st.columns([1, 1])
            with col1:
                st.plotly_chart(viz.plot_amino_acid_bar(dist), use_container_width=True)
            with col2:
                st.markdown("#### Protein Statistics")
                st.markdown(f"**Sequence Length:** {stats['length']} aa")
                st.markdown(f"**Unique residues:** {stats.get('unique_residues', 'N/A')}")
                st.markdown(f"**Residue diversity:** {len([v for v in dist['counts'].values() if v > 0])} / {len(dist['counts'])}")
                st.markdown(f"**Most abundant residue:** {max(dist['counts'], key=dist['counts'].get)}")

            if motifs:
                st.markdown("#### Motifs Found")
                for motif in motifs:
                    st.markdown(
                        f"- **{motif['name']}** (`{motif['motif']}`) — "
                        f"Position {motif['start']}–{motif['end']}  "
                        f"Match: `{motif['match']}`"
                    )
            else:
                st.info("No known protein motifs or regulatory elements detected.")

        else:
            st.markdown("#### Nucleotide Composition & GC Profile")

            col1, col2, col3 = st.columns([1, 1, 1])
            with col1:
                st.plotly_chart(viz.plot_nucleotide_pie(dist), use_container_width=True)
            with col2:
                st.plotly_chart(viz.plot_nucleotide_bar(dist), use_container_width=True)
            with col3:
                st.plotly_chart(viz.plot_gc_gauge(stats["gc_content"]), use_container_width=True)

            st.plotly_chart(
                viz.plot_gc_sliding_window(sequence, window=window_size),
                use_container_width=True,
            )

            st.markdown("#### Detailed Statistics")
            stat_col1, stat_col2 = st.columns(2)
            with stat_col1:
                st.markdown(f"""
| Property | Value |
|---|---|
| Sequence Length | `{stats['length']} bp` |
| GC Content | `{stats['gc_content']}%` |
| AT Content | `{stats['at_content']}%` |
| GC/AT Ratio | `{stats.get('gc_ratio', 'N/A')}` |
                """)
            with stat_col2:
                st.markdown(f"""
| Property | Value |
|---|---|
| Is coding length (×3) | `{'Yes' if stats['is_coding_length'] else 'No'}` |
| Has start codon (ATG) | `{'Yes' if stats['has_start_codon'] else 'No'}` |
| Has in-frame stop codon | `{'Yes' if stats['has_stop_codon'] else 'No'}` |
| A count | `{dist['counts']['A']}` |
| T count | `{dist['counts']['T']}` |
| G count | `{dist['counts']['G']}` |
| C count | `{dist['counts']['C']}` |
                """)

            if motifs:
                st.markdown("#### Regulatory Motifs Found")
                for motif in motifs:
                    st.markdown(
                        f"- **{motif['name']}** (`{motif['motif']}`) — "
                        f"Position {motif['start']}–{motif['end']}  "
                        f"Match: `{motif['match']}`"
                    )
            else:
                st.info("No known regulatory motifs detected in this sequence.")

    # ── Tab 2: Similarity ──────────────────────────────────────────────────────
    with tabs[1]:
        st.markdown("#### Database Similarity Search")

        if not similarity_results:
            st.warning("No similarity results available.")
        else:
            st.plotly_chart(
                viz.plot_similarity_scores(similarity_results),
                use_container_width=True,
            )

            for i, match in enumerate(similarity_results):
                classification = sim.classify_similarity(match["similarity_score"])
                with st.expander(
                    f"{classification['emoji']}  "
                    f"{match['gene_name']} — {match['similarity_score']:.1f}% similarity"
                ):
                    c1, c2 = st.columns(2)
                    with c1:
                        st.markdown(f"**Gene:** {match['gene_name']}")
                        st.markdown(f"**Trait:** {match['trait']}")
                        st.markdown(f"**Organism:** {match['organism']}")
                        st.markdown(f"**Accession:** {match['accession']}")
                    with c2:
                        st.markdown(f"**Similarity:** {match['similarity_score']:.1f}%")
                        st.markdown(f"**Level:** {classification['label']}")
                        st.markdown(f"**Interpretation:** {classification['interpretation']}")
                        st.markdown(f"**Description:** {match['description']}")

                    if match.get("alignment"):
                        st.markdown("**Alignment Map:**")
                        st.plotly_chart(
                            viz.plot_alignment(match["alignment"]["alignment_map"]),
                            use_container_width=True,
                            key=f"alignment_{i}",
                        )

    # ── Tab 3: Mutations ───────────────────────────────────────────────────────
    with tabs[2]:
        st.markdown("#### Mutation Analysis")

        if not mutation_report:
            st.info("No mutation report — run analysis with a database match first.")
        else:
            mc1, mc2, mc3, mc4 = st.columns(4)
            mc1.metric("Total Mutations", mutation_report["total_mutations"])
            mc2.metric("Mutation Rate", f"{mutation_report['mutation_rate_percent']}%")
            mc3.metric("Identity", f"{mutation_report['identity_percent']}%")
            mc4.metric("Compared Length", f"{mutation_report['compared_length']} bp")

            st.plotly_chart(
                viz.plot_mutation_map(mutation_report, mutation_report["compared_length"]),
                use_container_width=True,
            )

            mutations = mutation_report.get("mutations", [])
            if mutations:
                st.markdown("#### Point Mutations Table")
                st.markdown(
                    "| Position | Reference | Query | Type |\n"
                    "|---|---|---|---|\n" +
                    "\n".join(
                        f"| {m['position']} | `{m['reference']}` | `{m['query']}` | {m['type'].capitalize()} |"
                        for m in mutations[:50]
                    )
                )
                if len(mutations) > 50:
                    st.info(f"Showing first 50 of {len(mutations)} mutations.")
            else:
                st.success("No point mutations detected — sequence matches reference perfectly.")

    # ── Tab 4: Translation ─────────────────────────────────────────────────────
    with tabs[3]:
        if sequence_type == "protein":
            st.markdown("#### Protein input detected — translation not applicable")
            st.info(
                "The uploaded sequence is interpreted as a protein sequence. "
                "DNA translation and nucleotide complement calculations are skipped."
            )
            st.markdown("#### Protein Properties")
            st.markdown(f"**Sequence length:** {stats['length']} aa")
            st.markdown(f"**Unique residues:** {stats.get('unique_residues', 'N/A')}")
            st.markdown(f"**Most abundant residue:** {max(dist['counts'], key=dist['counts'].get)}")
        else:
            st.markdown(f"#### Protein Translation (Frame +{reading_frame + 1})")

            tl = translation
            st.markdown(f"**Protein length:** {tl['length']} amino acids")
            st.markdown(
                f"**Status:** {'Complete ORF (stop codon found)' if tl['status'] == 'complete' else 'No stop codon in frame'}"
            )

            if tl["protein"]:
                st.code(tl["protein"], language=None)
            else:
                st.warning("No protein sequence translated — check reading frame or sequence length.")

            st.markdown("#### All Reading Frames")
            all_frames = bio.translate_all_frames(sequence)
            for frame_name, frame_result in all_frames.items():
                with st.expander(f"{frame_name} — {frame_result['length']} aa"):
                    st.code(frame_result["protein"] or "(empty)", language=None)
                    st.caption(f"Status: {frame_result['status']}")

            st.markdown("#### Complementary Sequences")
            comp_col1, comp_col2 = st.columns(2)
            with comp_col1:
                st.markdown("**5'→3' Complement:**")
                st.code(bio.complement(sequence[:80]) + ("…" if len(sequence) > 80 else ""), language=None)
            with comp_col2:
                st.markdown("**Reverse Complement:**")
                st.code(bio.reverse_complement(sequence[:80]) + ("…" if len(sequence) > 80 else ""), language=None)

    # ── Tab 5: AI Interpretation ───────────────────────────────────────────────
    with tabs[4]:
        st.markdown("#### AI Biological Interpretation")

        interp = interpretation

        # Overall summary
        st.info(f"**Summary:** {interp['overall_summary']}")

        # Confidence badge
        conf = interp["confidence_level"]
        conf_colors = {"High": "🟢", "Medium": "🟡", "Low": "🔴"}
        st.markdown(
            f"**Confidence:** {conf_colors.get(conf['level'], '⚪')} {conf['level']} — {conf['note']}"
        )

        st.markdown("---")

        # Two-column layout
        left, right = st.columns(2)

        with left:
            st.markdown("##### Sequence Profile")
            profile = interp["sequence_profile"]
            for note in profile["notes"]:
                st.markdown(f"- {note}")
            st.markdown(f"*Coding potential: **{profile['coding_potential'].upper()}***")

            st.markdown("##### GC Content Analysis")
            gc_interp = interp["gc_interpretation"]
            for line in gc_interp["interpretation"]:
                st.markdown(f"- {line}")
            st.markdown(f"*{gc_interp['stress_implication']}*")

            st.markdown("##### Functional Prediction")
            func = interp["functional_prediction"]
            for p in func["predictions"]:
                st.markdown(f"- {p}")

        with right:
            st.markdown("##### Similarity Interpretation")
            sim_interp = interp["similarity_interpretation"]
            for line in sim_interp.get("interpretation", ["—"]):
                st.markdown(f"- {line}")

            st.markdown("##### Mutation Interpretation")
            mut_interp = interp["mutation_interpretation"]
            for line in mut_interp.get("interpretation", ["—"]):
                st.markdown(f"- {line}")

            st.markdown("##### Stress Resistance Assessment")
            stress = interp["stress_resistance"]
            detected = stress.get("detected_resistance", {})
            if detected:
                for stress_type, detail in detected.items():
                    st.markdown(f"- **{stress_type.upper()}:** {detail}")
            else:
                st.markdown("- No specific stress resistance profile detected.")

        st.markdown("---")
        st.markdown("#### Agricultural Recommendations")

        recs = interp["agricultural_recommendations"]
        priority_colors = {"HIGH": "🔴", "MEDIUM": "🟡", "LOW": "🔵"}
        for rec in recs:
            priority_icon = priority_colors.get(rec["priority"], "⚪")
            with st.expander(f"{priority_icon} [{rec['priority']}] {rec['category']}"):
                st.markdown(rec["recommendation"])

    # ── Tab 6: Raw Sequence ────────────────────────────────────────────────────
    with tabs[5]:
        st.markdown("#### Cleaned Sequence")
        if sequence_type == "protein":
            st.markdown(
                f"**Length:** {len(sequence)} aa  |  "
                f"**Valid amino acids only** (standard residues + X/B/Z/*)"
            )
        else:
            st.markdown(
                f"**Length:** {len(sequence)} bp  |  "
                f"**GC:** {stats['gc_content']}%  |  "
                f"**Valid nucleotides only** (ATGCN)"
            )
        st.code(sequence, language=None)

        st.markdown("#### Download")
        if sequence_type == "protein":
            fasta_content = f">Query_sequence | length={len(sequence)}aa\n{sequence}\n"
        else:
            fasta_content = f">Query_sequence | length={len(sequence)}bp | GC={stats['gc_content']}%\n{sequence}\n"
        st.download_button(
            label="Download as FASTA",
            data=fasta_content,
            file_name="query_sequence.fasta",
            mime="text/plain",
        )

        report_lines = [
            "AI-Powered Plant Gene Analyzer — Analysis Report",
            "=" * 50,
            f"Sequence Length: {stats['length']} {'aa' if sequence_type == 'protein' else 'bp'}",
        ]
        if sequence_type != "protein":
            report_lines += [
                f"GC Content: {stats['gc_content']}%",
                f"AT Content: {stats['at_content']}%",
            ]
        report_lines += [
            "",
            "DATABASE MATCHES",
            "-" * 30,
        ]
        for match in similarity_results:
            report_lines.append(
                f"{match['gene_name']}: {match['similarity_score']:.1f}% ({match['trait']})"
            )
        report_lines += [
            "",
            "AI SUMMARY",
            "-" * 30,
            interpretation.get("overall_summary", ""),
            "",
            "AGRICULTURAL RECOMMENDATIONS",
            "-" * 30,
        ]
        for rec in interpretation.get("agricultural_recommendations", []):
            report_lines.append(f"[{rec['priority']}] {rec['category']}: {rec['recommendation']}")

        report_text = "\n".join(report_lines)
        st.download_button(
            label="Download Analysis Report (.txt)",
            data=report_text,
            file_name="gene_analysis_report.txt",
            mime="text/plain",
        )
        st.markdown("#### Biological Annotation")
        if st.button("Run annotation", key="run_annotation"):
            from core_engines.annotation_engine import annotate_sequence
            try:
                anns = annotate_sequence(sequence, db=db)
                st.success("Annotation complete")
                st.json(anns)
            except Exception as e:
                logger.error(f"Annotation failed: {e}")
                st.error(f"Annotation failed: {e}")

    # ── Tab 7: Alignments (MSA + pairwise) ─────────────────────────────────
    with tabs[6]:
        st.markdown("#### Multiple Sequence Alignment")
        seqs_input = st.text_area("Paste multiple FASTA sequences or one per line:", height=160)
        msa_btn = st.button("Run MSA", key="msa_run")
        if msa_btn and seqs_input:
            from core_engines.alignment_engine import progressive_alignment, needleman_wunsch

            # parse simple input (one sequence per line or FASTA)
            from sequence_loader import parse_fasta
            records_msa = parse_fasta(seqs_input)
            sequences = [r['sequence'] for r in records_msa]
            if len(sequences) < 2:
                st.warning("Provide at least 2 sequences for MSA.")
            else:
                with st.spinner("Running progressive MSA..."):
                    msa_result = progressive_alignment(sequences, seq_type="dna")
                st.success(f"MSA complete — {msa_result.get('num_sequences')} sequences aligned")
                aligned = msa_result.get('aligned_sequences', [])
                labels = [r.get('header', f"Seq{i+1}") for i, r in enumerate(records_msa)]
                try:
                    fig_msa = viz.plot_msa_table(aligned, labels=labels)
                    st.plotly_chart(fig_msa, use_container_width=True)
                except Exception:
                    for aseq in aligned:
                        st.code(aseq, language=None)

        st.markdown("#### Pairwise Alignment")
        col_a, col_b = st.columns(2)
        with col_a:
            p_seq1 = st.text_area("Sequence 1", height=80, key="pw1")
        with col_b:
            p_seq2 = st.text_area("Sequence 2", height=80, key="pw2")
        if st.button("Align pairwise", key="pw_align"):
            from core_engines.alignment_engine import needleman_wunsch, smith_waterman
            if not p_seq1 or not p_seq2:
                st.warning("Provide two sequences for pairwise alignment.")
            else:
                nw = needleman_wunsch(p_seq1.strip(), p_seq2.strip())
                sw = smith_waterman(p_seq1.strip(), p_seq2.strip())
                st.markdown("**Needleman-Wunsch (global)**")
                # Show aligned text and interactive alignment map
                st.code(nw['seq1_aligned'] + "\n" + nw['seq2_aligned'])
                # Build match line
                match_line_nw = ''.join(['|' if a==b and a!='-' else ' ' for a,b in zip(nw['seq1_aligned'], nw['seq2_aligned'])])
                try:
                    st.plotly_chart(viz.plot_alignment({'query': nw['seq1_aligned'], 'reference': nw['seq2_aligned'], 'match_line': match_line_nw}), use_container_width=True)
                except Exception:
                    pass
                st.markdown(f"Score: {nw['alignment_score']} — Matches: {nw['match_count']} — Gaps: {nw['gap_count']}")
                st.markdown("**Smith-Waterman (local)**")
                st.code(sw['seq1_aligned'] + "\n" + sw['seq2_aligned'])
                match_line_sw = ''.join(['|' if a==b and a!='-' else ' ' for a,b in zip(sw['seq1_aligned'], sw['seq2_aligned'])])
                try:
                    st.plotly_chart(viz.plot_alignment({'query': sw['seq1_aligned'], 'reference': sw['seq2_aligned'], 'match_line': match_line_sw}), use_container_width=True)
                except Exception:
                    pass

    # ── Tab 8: Distance Matrix ───────────────────────────────────────────────
    with tabs[7]:
        st.markdown("#### Compute Pairwise Distance Matrix")
        dm_input = st.text_area("Paste FASTA or one sequence per line:", height=160)
        dm_method = st.selectbox("Method", options=["hamming", "jukes_cantor", "kimura", "pam"], index=2)
        if st.button("Compute Distance Matrix", key="dm_compute"):
            from sequence_loader import parse_fasta
            from core_engines.distance_engine import distance_matrix
            records_dm = parse_fasta(dm_input)
            sequences = [{"name": r.get('header', f"Seq{i+1}"), "sequence": r['sequence']} for i, r in enumerate(records_dm)]
            if len(sequences) < 2:
                st.warning("Provide at least 2 sequences to build distance matrix.")
            else:
                with st.spinner("Calculating distances..."):
                    dm_res = distance_matrix(sequences, method=dm_method)
                st.write("**Sequence names**", dm_res['sequence_names'])
                import pandas as pd
                df = pd.DataFrame(dm_res['distance_matrix'], index=dm_res['sequence_names'], columns=dm_res['sequence_names'])
                st.dataframe(df)
                st.download_button("Download CSV", df.to_csv().encode('utf-8'), file_name="distance_matrix.csv")

    # ── Tab 9: Phylogeny ─────────────────────────────────────────────────────
    with tabs[8]:
        st.markdown("#### Build Phylogenetic Tree")
        ph_input = st.text_area("Paste sequences for phylogeny (FASTA or lines):", height=160)
        ph_method = st.selectbox("Tree algorithm", options=["upgma", "neighbor_joining"], index=0)
        if st.button("Build Tree", key="build_tree"):
            from sequence_loader import parse_fasta
            from core_engines.distance_engine import distance_matrix
            from core_engines.phylogeny_engine import upgma, neighbor_joining, phylo_to_newick, newick_to_plotly_tree
            records_ph = parse_fasta(ph_input)
            seqs = [{"name": r.get('header', f"Seq{i+1}"), "sequence": r['sequence']} for i, r in enumerate(records_ph)]
            if len(seqs) < 2:
                st.warning("Provide at least 2 sequences for a simple tree; 3+ sequences are recommended for more meaningful phylogeny.")
            else:
                with st.spinner("Computing distance matrix..."):
                    dm = distance_matrix(seqs, method="kimura")
                mat = dm['distance_matrix']
                import numpy as np
                mat_np = np.array(mat)
                with st.spinner("Building tree..."):
                    if ph_method == "upgma":
                        tree = upgma(mat_np, dm['sequence_names'])
                    else:
                        tree = neighbor_joining(mat_np, dm['sequence_names'])
                st.write("**Tree metadata**", {"algorithm": tree.get('algorithm'), "tree_type": tree.get('tree_type')})
                # If dendrogram data available, plot interactive dendrogram
                try:
                    if tree.get('dendrogram_data'):
                        fig = viz.plot_dendrogram(tree['dendrogram_data'], labels=tree.get('sequence_names'))
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.info('Dendrogram data not available for this method; showing edge list instead.')
                        if tree.get('edges'):
                            st.table(tree['edges'])
                except Exception as e:
                    logger.warning(f"Failed to render dendrogram: {e}")

                # Simple Newick fallback (leaf list)
                try:
                    newick = "(" + ",".join(tree.get('sequence_names', [])) + ");"
                except Exception:
                    newick = None
                st.markdown("#### Newick format")
                st.code(newick or "Newick not available")

    # ── Tab 10: Protein Analysis ─────────────────────────────────────────────
    with tabs[9]:
        st.markdown("#### Protein properties (simple)")
        prot_seq = st.text_area("Paste protein sequence:", height=120)
        if st.button("Analyze protein", key="prot_analyze"):
            from bioinformatics import generate_protein_statistics
            if not prot_seq:
                st.warning("Paste a protein sequence to analyze.")
            else:
                stats = generate_protein_statistics(prot_seq.strip())
                st.json(stats)

else:
    # ── Welcome screen ──────────────────────────────────────────────────────────
    st.markdown(
        """
        <div style='text-align:center;padding:3rem 1rem;'>
            <p style='font-size:3rem;'>🧬</p>
            <h3 style='color:#69f0ae;font-family:Courier New;'>Ready to Analyze</h3>
            <p style='color:#a5d6a7;font-family:Courier New;'>
                Paste a plant DNA sequence above or load a demo,<br>
                then click <b>Analyze Sequence</b>.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.expander("What does this app analyze?"):
        st.markdown(
            """
            | Feature | Description |
            |---|---|
            | GC Content | % of G and C nucleotides |
            | Nucleotide Distribution | Count and % of A, T, G, C, N |
            | Sliding Window GC | GC content profile along the sequence |
            | Database Similarity | Comparison against 8 known plant genes |
            | Alignment | Visual match/mismatch display |
            | Mutation Detection | Point mutations vs. best reference |
            | Protein Translation | All 3 reading frames, codon table |
            | Motif Search | Known plant regulatory elements |
            | AI Interpretation | Rule-based biological explanation |
            | Recommendations | Agronomic insights and research guidance |
            """
        )

    with st.expander("Supported input formats"):
        st.markdown(
            """
            - **Raw DNA**: paste directly (e.g., `ATGCGTAGCTAG...`)
            - **FASTA**: with or without header lines starting with `>`
            - **Upload**: `.fasta`, `.fa`, or `.txt` files
            - **Valid nucleotides**: A, T, G, C, N (case-insensitive)
            """
        )
