"""
app.py
------
Plant Gene Analyzer — Streamlit frontend.

Run with:
    streamlit run app.py
"""

import streamlit as st
import pandas as pd
import json
import os
import io
import base64
import logging
import sys
import time
from pathlib import Path

# ── Local modules ──────────────────────────────────────────────────────────────
import bioinformatics as bio
from organism_reference import get_organism_reference
import similarityengine as sim
import visualization as viz
import export_utils as export_util
import re
import sequence_loader as loader
import config
import pipeline
import trait_research as tr

SCRIPT_ROOT = Path(__file__).resolve().parent.parent  # project root (this file now lives in views/)
sys.path.insert(0, str(SCRIPT_ROOT / "scripts"))

try:
    from scripts.postgres_utils import (
        load_gene_database_from_postgres,
        load_gene_database_metadata_from_postgres,
        get_gene_count,
        search_gene_metadata,
        count_gene_metadata_matches,
        get_gc_content_stats_for_organism,
        get_codon_usage_for_organism,
        get_codon_reference_for_organism,
        get_length_stats_for_organism,
    )
except ImportError:
    load_gene_database_from_postgres = None
    load_gene_database_metadata_from_postgres = None
    get_gene_count = None
    search_gene_metadata = None
    count_gene_metadata_matches = None
    get_gc_content_stats_for_organism = None
    get_codon_usage_for_organism = None
    get_codon_reference_for_organism = None
    get_length_stats_for_organism = None

# ─── Configure logging ─────────────────────────────────────────────────────────
logger = config.get_logger(__name__)


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


# ─── Load background video ──────────────────────────────────────────────────────
def load_video_background(video_path: str = "assets/images.mp4", max_mb: float = 30.0) -> None:
    """
    Inject a fixed, full-screen, looping, muted video as the app background.

    The video is base64-embedded directly into the page as a data URI, so no
    Streamlit static-file-serving configuration or specific folder name
    (e.g. `static/`) is required — it works from wherever `video_path` points,
    such as the existing `assets/` folder.
    """
    path = Path(video_path)
    try:
        if not path.exists():
            logger.warning(f"Background video not found: {path.resolve()}")
            return

        size_mb = path.stat().st_size / (1024 * 1024)
        if size_mb > max_mb:
            logger.warning(
                f"Background video is {size_mb:.1f} MB — base64 embedding will "
                f"slow down page load. Consider compressing it (e.g. with "
                f"HandBrake, target < 10 MB, 1080p, no audio track)."
            )

        video_b64 = base64.b64encode(path.read_bytes()).decode("utf-8")
        st.markdown(
            f"""
            <video autoplay loop muted playsinline class="bg-video">
                <source src="data:video/mp4;base64,{video_b64}" type="video/mp4">
            </video>
            <div class="bg-video-overlay"></div>
            """,
            unsafe_allow_html=True,
        )
        logger.info(f"Background video loaded from {path} ({size_mb:.1f} MB)")
    except Exception as e:
        logger.error(f"Error loading background video: {e}")


# analyze_sequence_record() and get_alignment_map() now live in pipeline.py
# (see import below) — this keeps the analysis orchestration testable and
# reusable independently of the Streamlit UI.


@st.cache_data(show_spinner=False)
def _cached_analyze(
    record_json: str,
    input_type: str,
    reading_frame: int,
    top_n_matches: int,
    similarity_deep_search: bool,
    _db: dict,
) -> dict:
    """Streamlit-cached wrapper around pipeline.analyze_sequence_record.

    Avoids recomputing GC%, ORFs, alignments, mutation detection, etc. when
    Streamlit re-runs the script for an unrelated widget interaction (e.g.
    toggling a chart option) with the exact same sequence and settings.
    The record is passed as a JSON string (not a dict) because st.cache_data
    needs a hashable argument. `_db` is prefixed with an underscore, a
    Streamlit convention meaning "don't hash this for the cache key" — for
    a fixed record + settings, sim.find_similar_genes() should return the
    same small candidate set deterministically, so skipping the hash of a
    (potentially large-ish) candidate dict on every call is safe and saves
    work.
    """
    record = json.loads(record_json)
    return pipeline.analyze_sequence_record(
        record,
        input_type,
        reading_frame,
        db=_db,
        top_n_matches=top_n_matches,
        enable_length_prefilter=not similarity_deep_search,
        logger=logger,
    )


# ─── Load gene database with caching ────────────────────────────────────────────────────────────
@st.cache_data
def load_gene_database_cached(db_path: str = str(config.DATABASE_PATH)) -> dict:
    """
    Load gene database with Streamlit caching to improve performance.
    Prefer PostgreSQL if the helper is available and configured.
    """
    try:
        # Try PostgreSQL first if available
        if load_gene_database_from_postgres is not None:
            try:
                db = load_gene_database_from_postgres()
                if db:
                    logger.info(f"Loaded {len(db)} genes from PostgreSQL")
                    return db
                # PostgreSQL returned no records - log warning and fall through to JSON
                logger.warning("PostgreSQL database returned no records. Falling back to JSON.")
            except Exception as e:
                logger.warning(f"PostgreSQL load failed: {e}. Falling back to JSON.")
                # Continue to JSON fallback instead of returning empty
                pass

        # JSON fallback. The large generated dataset is intentionally ignored
        # by Git, so deployed environments may only have the small tracked
        # fallback database available when PostgreSQL is unavailable.
        fallback_path = Path(db_path)
        if not fallback_path.exists():
            tracked_fallback = SCRIPT_ROOT / "genes_database.json"
            if tracked_fallback.exists():
                logger.warning(
                    f"Configured database not found at {db_path}; using tracked fallback {tracked_fallback}"
                )
                fallback_path = tracked_fallback
            else:
                logger.warning(f"Database not found at {db_path}")
                return {}

        db = sim.load_gene_database(str(fallback_path))
        logger.info(f"Loaded {len(db)} genes from database")
        try:
            # Build a k-mer index once per cached load to accelerate
            # similarity prefilters. This mutates `db` in-place so the
            # cached object contains the precomputed `_kmers` sets.
            if isinstance(db, dict):
                try:
                    sim._ensure_kmer_index(db)
                    logger.info("K-mer index built for database (cached)")
                except Exception as e:
                    logger.warning(f"K-mer index build failed: {e}")
        except Exception:
            # Keep original behavior if anything goes wrong here.
            pass
        return db

    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing error: {e}")
        st.error(f"❌ Error parsing gene database: {db_path}")
        return {}
    except Exception as e:
        logger.error(f"Error loading database: {e}")
        st.error(f"❌ Error loading database: {e}")
        return {}


@st.cache_data(ttl=300, show_spinner=False)
def get_gene_count_cached() -> int:
    """Cheap total-row count for the sidebar header. Cached for 5 minutes
    so it isn't re-queried on every widget interaction/rerun."""
    return get_gene_count()


@st.cache_data(ttl=60, show_spinner=False)
def search_gene_metadata_cached(query: str, limit: int = 20, offset: int = 0) -> list[dict]:
    """Cached, server-side search — only `limit` rows ever get pulled from
    Postgres and only `limit` rows ever get built into Python dicts,
    regardless of how large the `genes` table is. A short TTL (rather than
    the default indefinite cache) keeps results from going stale if the
    table is being actively ingested into, while still absorbing the
    repeated calls a Streamlit rerun triggers for an unchanged query."""
    return search_gene_metadata(query or None, limit=limit, offset=offset)


@st.cache_data(ttl=60, show_spinner=False)
def count_gene_metadata_matches_cached(query: str) -> int:
    return count_gene_metadata_matches(query)


@st.cache_data(ttl=3600, show_spinner=False)
def get_gc_content_stats_cached(organism: str) -> dict:
    if get_gc_content_stats_for_organism is None:
        return {"mean_gc": 0.0, "stdev_gc": 0.0, "n_sequences": 0}
    return get_gc_content_stats_for_organism(organism)


@st.cache_data(ttl=3600, show_spinner=False)
def get_codon_usage_cached(organism: str) -> dict[str, float]:
    if get_codon_usage_for_organism is None:
        return {}
    return get_codon_usage_for_organism(organism)


@st.cache_data(ttl=3600, show_spinner=False)
def get_codon_reference_cached(organism: str) -> dict:
    if get_codon_reference_for_organism is None:
        return {"value": {}, "n": 0}
    return get_codon_reference_for_organism(organism)


@st.cache_data(ttl=3600, show_spinner=False)
def get_length_stats_cached(organism: str) -> dict:
    if get_length_stats_for_organism is None:
        return {"mean_length": 0.0, "stdev_length": 0.0, "n_sequences": 0}
    return get_length_stats_for_organism(organism)


def build_methods_paragraph(result: dict, references: dict[str, dict]) -> str:
    """Format available computed metrics as a publication-ready sentence."""
    stats = result.get("stats", {})
    length = stats.get("length")
    unit = "aa" if result.get("sequence_type") == "protein" else "bp"
    parts = [f"The {length:,} {unit} sequence"] if length else ["The sequence"]
    if result.get("sequence_type") == "protein":
        props = result.get("protein_stats") or {}
        if props.get("isoelectric_point") is not None:
            parts.append(f"had an estimated pI of {props['isoelectric_point']:.2f}")
        if props.get("gravy") is not None:
            parts.append(f"and a GRAVY score of {props['gravy']:.2f}")
    else:
        if stats.get("gc_content") is not None:
            sentence = f"exhibited a GC content of {stats['gc_content']:.2f}%"
            gc_ref = references.get("gc", {})
            if gc_ref.get("available"):
                sentence += f" (species average: {gc_ref['value']:.2f}%, n={gc_ref['n']})"
            parts.append(sentence)
        length_ref = references.get("length", {})
        if length and length_ref.get("available"):
            parts.append(
                f"with a length {length - float(length_ref['value']):+.0f} bp from the species mean"
                f" ({float(length_ref['value']):.0f} bp, n={length_ref['n']})"
            )
        methylation = result.get("methylation_context") or {}
        if methylation.get("cg") and methylation.get("chg") and methylation.get("chh"):
            parts.append(
                "and methylation-context proportions of "
                f"CG {methylation['cg']['pct']:.2f}%, "
                f"CHG {methylation['chg']['pct']:.2f}%, and "
                f"CHH {methylation['chh']['pct']:.2f}%"
            )
        if stats.get("has_complete_orf"):
            parts.append("consistent with a complete open reading frame (start-to-stop, same frame)")
        else:
            parts.append("without a complete start-to-stop open reading frame")
    return " ".join(parts) + "."


load_css()
if os.getenv("ENABLE_VIDEO_BACKGROUND", "false").lower() in {"1", "true", "yes"}:
    load_video_background()


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
    similarity_deep_search = st.checkbox(
        "Enable deep similarity search",
        value=False,
        help=(
            "Disable the alignment length prefilter and evaluate more candidates. "
            "This is slower, but increases sensitivity for short or divergent queries."
        ),
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

    db = None
    metadata = None
    metadata_available = False

    if get_gene_count is not None and search_gene_metadata is not None:
        try:
            # A cheap COUNT(*) rather than materializing all ~56k rows just
            # to call len() on them. Cached for 5 minutes so a widget
            # interaction elsewhere on the page doesn't re-issue it.
            total_genes = get_gene_count_cached()
            metadata_available = total_genes > 0
            st.success(f"✅ {total_genes} gene metadata records available")
            st.markdown(
                "The app loads lightweight gene metadata first for search and filtering. "
                "Full sequence data is loaded only when an analysis is run."
            )

            gene_search = st.text_input(
                "Search gene ID, symbol, or trait",
                value="",
                help="Filter the loaded gene database by gene_id, symbol, or trait.",
            )
            if gene_search:
                query = gene_search.strip()
                # Server-side ILIKE search (see postgres_utils.search_gene_metadata)
                # -- only the ~20 rows actually shown ever leave Postgres,
                # instead of pulling all ~56k rows into Python on every
                # keystroke and filtering them in a list comprehension.
                match_count = count_gene_metadata_matches_cached(query)
                filtered = search_gene_metadata_cached(query, limit=20)
                st.write(f"Showing {len(filtered)} of {match_count} matching gene metadata records")
            else:
                filtered = search_gene_metadata_cached("", limit=20)
                st.info("Showing a sample of 20 gene metadata records. Use search to filter specific genes.")

            with st.expander("Preview gene metadata"):
                for gene in filtered:
                    symbol = gene.get("symbol", "Unknown")
                    gene_id = gene.get("gene_id", "n/a")
                    trait = ", ".join(gene.get("traits", [])[:3]) or "No trait specified"
                    description = gene.get("description", "No description")
                    st.markdown(f"- **{symbol}** (`{gene_id}`) — {trait} — {description}")

            st.info("Full gene database with sequences will be loaded when you start an analysis.")

        except Exception as e:
            logger.warning(f"Lightweight gene metadata load failed: {e}")
            st.warning("Could not load gene metadata preview. Falling back to full database load.")
            db = load_gene_database_cached(str(config.DATABASE_PATH))
    else:
        db = load_gene_database_cached(str(config.DATABASE_PATH))

    if db is not None:
        if not db:
            st.error(
                "❌ No genes available: PostgreSQL could not be loaded and the tracked fallback database is missing."
            )
        elif isinstance(db, dict) and db:
            fallback_note = " (local fallback; PostgreSQL unavailable)" if not metadata_available else ""
            st.success(f"✅ {len(db)} genes loaded{fallback_note}")
    elif metadata_available:
        st.info("✅ Lightweight gene metadata available. Full database will load on analysis.")
    else:
        st.error("❌ No genes available: configure the PostgreSQL connection or add the fallback database.")


# ─── Main header ───────────────────────────────────────────────────────────────
# ─── Main header ───────────────────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-panel">
        <h1>🧬 Plant Gene Analyzer</h1>
        <p class="hero-subtitle">Bioinformatics · AI Interpretation · Agricultural Insights</p>
    </div>
    """,
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
        if uploaded_file.size > config.MAX_UPLOAD_SIZE_BYTES:
            st.error(
                f"❌ File too large ({uploaded_file.size / 1024 / 1024:.1f} MB). "
                f"Maximum allowed is {config.MAX_UPLOAD_SIZE_BYTES / 1024 / 1024:.0f} MB."
            )
            content = None
        else:
            raw_bytes = uploaded_file.read()
            try:
                content = raw_bytes.decode("utf-8")
            except UnicodeDecodeError:
                try:
                    content = raw_bytes.decode("latin-1")
                    st.warning(
                        "⚠️ File is not valid UTF-8; decoded as Latin-1 instead. "
                        "Double-check the sequence for unexpected characters."
                    )
                except UnicodeDecodeError:
                    st.error("❌ Could not decode file — unsupported text encoding.")
                    content = None
        if content is not None:
            records = loader.parse_fasta(content)
            st.info(f"File loaded: {uploaded_file.name}")
    elif raw_sequence:
        records = loader.parse_fasta(raw_sequence)
    
    if records:
        if len(records) > 1:
            record_options = [
                f"{idx + 1}. {r.get('metadata', {}).get('name', r['header'])}"
                for idx, r in enumerate(records)
            ]
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
        if st.button("Load Demo Sequence"):
            raw_sequence = demo["seq"]
            st.session_state["loaded_demo"] = demo["seq"]

    if "loaded_demo" in st.session_state and not raw_sequence:
        raw_sequence = st.session_state["loaded_demo"]

if raw_sequence and sequence_input_type != "Protein":
    preview_dna = bio.clean_sequence(raw_sequence, sequence_type="dna")
    if preview_dna:
        with st.expander("Reading-frame preview"):
            st.caption("Six-frame summary to guide the reading-frame choice before analysis.")
            st.dataframe(
                pd.DataFrame(bio.all_frames_summary(preview_dna))
                .drop(columns=["frame"])
                .rename(columns={
                    "label": "Frame", "strand": "Strand", "has_start_codon": "Start",
                    "has_stop_codon": "Stop",
                    "longest_orf_length": "Longest ORF (bp, stop codon included)",
                    "orfs_complete": "Complete ORFs",
                    "orfs_truncated": "Truncated ORFs",
                }),
                hide_index=True,
                width="stretch",
            )

analyze_btn = st.button("🔬 Analyze Sequence", type="primary")

st.markdown("---")

st.info(
    "🧪 Looking for alignments, distance matrices, phylogeny, or standalone "
    "protein analysis? They now live on their own page — see **Outils "
    "indépendants** in the sidebar navigation, always available."
)


# ─── Analysis pipeline ─────────────────────────────────────────────────────────
if analyze_btn or (raw_sequence and "last_result" in st.session_state):

    if analyze_btn and not raw_sequence and not records:
        st.warning("⚠️ Please enter or paste a DNA sequence before analyzing.")
        st.stop()

    if analyze_btn and (raw_sequence or records):

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
                    similarity_started_at = time.perf_counter()

                    if db is not None:
                        # JSON-file deployment (no Postgres configured) —
                        # the whole database already lives in memory, same
                        # as before.
                        target_db = db
                    elif metadata_available:
                        # Postgres-backed deployment: never load the full
                        # ~56k-gene database. Depending on search mode:
                        # - Balanced (default): sim.find_similar_genes() uses
                        #   compact pg_trgm with length prefilter
                        # - Deep search: sim.find_similar_genes_deep() scans
                        #   ALL genes by trigram (no length filter) then aligns
                        #   top ~500 candidates (more thorough, takes 30-60s)
                        if similarity_deep_search:
                            if logger:
                                logger.info("Deep Search mode: exhaustive trigram scan + precision alignment...")
                            with st.spinner("Deep Search in progress... (30-60s) Scanning all genes and aligning top candidates"):
                                target_db = sim.find_similar_genes_deep(
                                    record.get("sequence", ""),
                                    top_n=top_n_matches,
                                    alignment_limit=500,
                                    logger=logger,
                                )
                        else:
                            # Balanced (default) mode
                            target_db = sim.find_similar_genes(
                                record.get("sequence", ""),
                                top_n=top_n_matches,
                                logger=logger,
                            )
                        
                        if not target_db:
                            logger.warning(
                                f"No candidate genes found for record {idx + 1} "
                                "(database search returned no results)"
                            )
                    else:
                        target_db = {}

                    analyzed_result = _cached_analyze(
                            json.dumps(record, sort_keys=True),
                            sequence_input_type,
                            reading_frame,
                            top_n_matches,
                            similarity_deep_search,
                            _db=target_db,
                        )
                    analyzed_result["similarity_elapsed_seconds"] = round(
                        time.perf_counter() - similarity_started_at, 3
                    )
                    analyzed_result["similarity_candidate_pool_count"] = len(target_db or {})
                    analyzed_results.append(analyzed_result)

            st.session_state["last_results"] = analyzed_results
            st.session_state["last_result"] = analyzed_results[0] if analyzed_results else None
            logger.info("Analysis session state saved")

        except Exception as e:
            error_msg = str(e)
            logger.error(f"Unexpected error during analysis: {e}")
            
            # Provide user-friendly messages for common database errors
            if "SSL connection" in error_msg and "closed" in error_msg:
                st.error(
                    "⚠️ **Database connection interrupted** — The server temporarily lost connection to the gene database. "
                    "This sometimes happens with high volume. Please try again in a few seconds.\n\n"
                    "_Technical: SSL connection to PostgreSQL pooler was closed unexpectedly._"
                )
            elif "consuming input" in error_msg:
                st.error(
                    "⚠️ **Database query timeout** — The similarity search took too long. "
                    "Try using fewer top matches or disable deep search mode, then retry.\n\n"
                    f"_Error: {error_msg[:100]}_"
                )
            elif "connection" in error_msg.lower():
                st.error(
                    "⚠️ **Could not connect to the gene database** — Check your internet connection and try again.\n\n"
                    f"_Technical: {error_msg[:150]}_"
                )
            else:
                st.error(f"❌ Analysis failed: {e}")
            st.stop()

    last_results = st.session_state.get("last_results")
    if not last_results:
        st.stop()

    batch_mode = len(last_results) > 1
    selected_batch_index = 0
    if batch_mode:
        record_options = [
            f"{idx + 1}. {item.get('header_metadata', {}).get('name', item['header'])}"
            for idx, item in enumerate(last_results)
        ]
        selected_batch_index = st.selectbox(
            "Select a sequence to inspect in this batch:",
            options=list(range(len(last_results))),
            format_func=lambda i: record_options[i],
            help="Choose a sequence to display its detailed statistics and charts.",
        )

    result = last_results[selected_batch_index]
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
    organism = result.get("organism") or result.get("header_metadata", {}).get("organism")
    gc_reference = get_organism_reference(
        organism, "gc", fetcher=get_gc_content_stats_cached if organism else None
    )
    codon_reference = get_organism_reference(
        organism, "codon usage", fetcher=get_codon_reference_cached if organism else None
    )
    length_reference = get_organism_reference(
        organism, "length", fetcher=get_length_stats_cached if organism else None
    )
    organism_codon_usage = codon_reference.get("value") or {}
    low_complexity = result.get("low_complexity", {"regions": [], "coverage_pct": 0.0})

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
                    similarity_score = item["best_match"]["similarity_score"]
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
    export_col1, export_col2, export_col3, export_col4, export_col5, export_col6, export_col7 = st.columns(7)
    
    with export_col1:
        if st.button("📄 Download JSON"):
            try:
                json_path = export_util.export_results_json(result)
                with open(json_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 JSON Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.json",
                        mime="application/json",
                    )
                logger.info(f"JSON export created: {json_path}")
                st.success("✅ JSON exported successfully")
            except Exception as e:
                logger.error(f"JSON export failed: {e}")
                st.error(f"Export failed: {e}")
    
    with export_col2:
        if st.button("📊 Download CSV"):
            try:
                csv_path = export_util.export_results_csv(result)
                with open(csv_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 CSV Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.csv",
                        mime="text/csv",
                    )
                logger.info(f"CSV export created: {csv_path}")
                st.success("✅ CSV exported successfully")
            except Exception as e:
                logger.error(f"CSV export failed: {e}")
                st.error(f"Export failed: {e}")
    
    with export_col3:
        if st.button("🌐 Download HTML"):
            try:
                html_path = export_util.export_results_html(result)
                with open(html_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 HTML Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.html",
                        mime="text/html",
                    )
                logger.info(f"HTML export created: {html_path}")
                st.success("✅ HTML exported successfully")
            except Exception as e:
                logger.error(f"HTML export failed: {e}")
                st.error(f"Export failed: {e}")
    with export_col4:
        if st.button("📑 Download XLSX"):
            try:
                xlsx_path = export_util.export_results_xlsx(result)
                with open(xlsx_path, "rb") as f:
                    st.download_button(
                        "📥 XLSX Report",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    )
                logger.info(f"XLSX export created: {xlsx_path}")
                st.success("✅ XLSX exported successfully")
            except Exception as e:
                logger.error(f"XLSX export failed: {e}")
                st.error(f"Export failed: {e}")

    with export_col5:
        if st.button("🧬 Download FASTA"):
            try:
                fasta_path = export_util.export_results_fasta(result)
                with open(fasta_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 FASTA sequence",
                        f.read(),
                        file_name=f"analysis_{stats['length']}{'aa' if sequence_type == 'protein' else 'bp'}.fasta",
                        mime="text/plain",
                    )
                logger.info(f"FASTA export created: {fasta_path}")
                st.success("✅ FASTA exported successfully")
            except Exception as e:
                logger.error(f"FASTA export failed: {e}")
                st.error(f"Export failed: {e}")

    with export_col6:
        if sequence_type == "protein":
            st.caption("GFF3 is available for DNA features only.")
        elif st.button("🧭 Download GFF3"):
            try:
                gff3_path = export_util.export_results_gff3(result)
                with open(gff3_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        "📥 GFF3 annotations",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.gff3",
                        mime="text/plain",
                    )
                logger.info(f"GFF3 export created: {gff3_path}")
                st.success("✅ GFF3 exported successfully")
            except Exception as e:
                logger.error(f"GFF3 export failed: {e}")
                st.error(f"Export failed: {e}")

    with export_col7:
        methods_paragraph = build_methods_paragraph(
            result, {"gc": gc_reference, "codon": codon_reference, "length": length_reference}
        )
        st.download_button(
            "Copy methods paragraph",
            methods_paragraph,
            file_name="methods_paragraph.txt",
            mime="text/plain",
            help="Copy the generated methods sentence from the downloaded text.",
        )

    st.markdown("---")

    # ── KPI Metrics ────────────────────────────────────────────────────────────
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

    if result.get("header_metadata"):
        header_meta = result["header_metadata"]
        header_notes = []
        if header_meta.get("gc"):
            header_notes.append(f"Header GC: {header_meta['gc']}")
        if header_meta.get("trait"):
            header_notes.append(f"Header trait: {header_meta['trait']}")
        if header_notes:
            st.info("**Header annotations:** " + ", ".join(header_notes))

    if result.get("metadata_warnings"):
        for warning_msg in result["metadata_warnings"]:
            st.warning(warning_msg)

    st.markdown("---")

    # ── Tabs ───────────────────────────────────────────────────────────────────
    st.markdown('<div class="section-heading"><span class="section-index">01</span><span>Analysis results</span></div>', unsafe_allow_html=True)
    tabs = st.tabs([
        "Statistics",
        "Similarity",
        "Mutations",
        "Translation",
        "AI Interpretation",
        "Raw Sequence",
    ])

    # ── Tab 1: Statistics ──────────────────────────────────────────────────────
    with tabs[0]:
        if sequence_type == "protein":
            st.markdown("#### Amino Acid Composition")

            col1, col2 = st.columns([1, 1])
            with col1:
                st.plotly_chart(viz.plot_amino_acid_bar(dist), width='stretch')
            with col2:
                st.markdown("#### Protein Statistics")
                st.markdown(f"**Sequence Length:** {stats['length']} aa")
                st.markdown(f"**Unique residues:** {stats.get('unique_residues', 'N/A')}")
                st.markdown(f"**Residue diversity:** {len([v for v in dist['counts'].values() if v > 0])} / {len(dist['counts'])}")
                st.markdown(f"**Most abundant residue:** {max(dist['counts'], key=dist['counts'].get)}")
                st.markdown(f"**GRAVY:** {protein_stats.get('gravy', 'N/A')}")
                st.markdown(f"**Instability index:** {protein_stats.get('instability_index', 'N/A')}")
                st.markdown(f"**Aliphatic index:** {protein_stats.get('aliphatic_index', 'N/A')}")

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
                st.plotly_chart(viz.plot_nucleotide_pie(dist), width='stretch')
            with col2:
                st.plotly_chart(viz.plot_nucleotide_bar(dist), width='stretch')
            with col3:
                if gc_reference.get("available"):
                    st.plotly_chart(
                        viz.plot_gc_gauge(
                            stats["gc_content"],
                            reference_low=max(0.0, gc_reference["value"] - get_gc_content_stats_cached(organism).get("stdev_gc", 0.0)),
                            reference_high=min(100.0, gc_reference["value"] + get_gc_content_stats_cached(organism).get("stdev_gc", 0.0)),
                        ),
                        width='stretch',
                    )
                    st.caption(f"Based on {gc_reference['n']:,} sequences of {organism} in database")
                else:
                    st.plotly_chart(viz.plot_gc_gauge(stats["gc_content"]), width='stretch')
                    st.info(gc_reference["fallback_reason"])

            st.plotly_chart(
                viz.plot_gc_sliding_window(sequence, window=window_size),
                width='stretch',
            )
            skew_profile = bio.gc_skew_profile(sequence, window=window_size)
            if skew_profile:
                st.plotly_chart(viz.plot_gc_skew_profile(skew_profile), width='stretch')

            methylation = result.get("methylation_context") or {}
            st.markdown("#### Cytosine Methylation Context (plant-specific)")
            if methylation:
                methyl_cols = st.columns(3)
                for column, label, key in zip(methyl_cols, ("CG", "CHG", "CHH"), ("cg", "chg", "chh")):
                    column.metric(label, f"{methylation[key]['pct']:.2f}%", f"{methylation[key]['count']} cytosines")
                st.caption("Plant methylation contexts: CG, CHG (H = A/T/C), and CHH. Percentages use classifiable cytosines.")

            quality = result.get("quality_report", {})
            if quality.get("applicable", True):
                quality_status = "Yes" if quality.get("valid") else "No"
                reason = f" ({quality.get('reason')})" if quality.get("reason") else ""
                st.markdown(
                    f"**Passes collection quality filter:** {quality_status}{reason} "
                    f"({quality.get('n_pct', 0):.2f}% N, threshold {quality.get('threshold_pct', 5):.2f}%)"
                )

            if low_complexity.get("regions"):
                st.warning(
                    f"{low_complexity['coverage_pct']:.1f}% of this sequence is repetitive/low complexity; "
                    "similarity matches involving these regions may not reflect true homology."
                )

            if result.get("codon_usage"):
                st.markdown("#### Codon usage")
                query_usage = result["codon_usage"]
                query_total = sum(query_usage.values()) or 1
                divergent = []
                for codon, count in query_usage.items():
                    query_pct = count / query_total
                    row = {"Codon": codon, "Sequence %": round(query_pct * 100, 2)}
                    if codon_reference.get("available"):
                        species_pct = organism_codon_usage.get(codon, 0.0)
                        row.update({"Species %": round(species_pct * 100, 2), "Delta %": round((query_pct - species_pct) * 100, 2)})
                    divergent.append(row)
                if codon_reference.get("available"):
                    divergent.sort(key=lambda row: abs(row["Delta %"]), reverse=True)
                st.dataframe(pd.DataFrame(divergent[:10]), hide_index=True, width="stretch")
                if not codon_reference.get("available"):
                    st.info(codon_reference["fallback_reason"] + " Species comparison is omitted until the minimum sample size is reached.")
                if len(sequence) % 3:
                    st.caption("The sequence was truncated to complete codons; the trailing bases were excluded.")

                cai = bio.codon_adaptation_index(sequence, organism_codon_usage) if codon_reference.get("available") else None
                if cai is not None:
                    st.metric("Codon Adaptation Index (CAI)", f"{cai:.3f}")
                    st.caption("Approximation based on the organism-wide codon distribution, not a highly expressed-gene reference set.")

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
                if length_reference.get("available"):
                    mean_length = float(length_reference["value"])
                    delta = stats["length"] - mean_length
                    st.markdown(f"**Length vs species:** {delta:+.0f} bp versus the species mean ({mean_length:.0f} bp, n={length_reference['n']}).")
                else:
                    st.info(length_reference["fallback_reason"])
            with stat_col2:
                st.markdown(f"""
| Property | Value |
|---|---|
| Is coding length (×3) | `{'Yes' if stats['is_coding_length'] else 'No'}` |
| Contains ATG in any frame | `{'Yes' if stats['has_start_codon'] else 'No'}` |
| Contains stop codon in any frame | `{'Yes' if stats['has_stop_codon'] else 'No'}` |
| Complete ORF found (start→stop, same frame) | `{'Yes' if stats.get('has_complete_orf') else 'No'}` |
| A count | `{dist['counts']['A']}` |
| T count | `{dist['counts']['T']}` |
| G count | `{dist['counts']['G']}` |
| C count | `{dist['counts']['C']}` |
                """)
                # "Contains ATG/stop in any frame" above are independent
                # existence checks across all 6 reading frames — they don't
                # imply a start and stop belong to the same ORF. Only
                # "Complete ORF found" (from bioinformatics.find_orfs, which
                # actually pairs a start with its in-frame stop) supports a
                # "this sequence contains a real gene" claim; see
                # sequence_statistics()'s has_complete_orf docstring.
                if stats['has_start_codon'] and stats['has_stop_codon'] and not stats.get('has_complete_orf'):
                    st.caption(
                        "⚠️ An ATG and a stop codon both exist somewhere in the sequence, but not as "
                        "a matching start→stop pair in the same frame — see the ORFs tab/section for what was actually found."
                    )

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

            restriction_sites = result.get("restriction_sites", [])
            st.markdown("#### Restriction sites")
            if restriction_sites:
                st.dataframe(pd.DataFrame(restriction_sites), hide_index=True, width="stretch")
            else:
                st.caption("No sites for the common enzymes in the panel were detected.")

            primer_hints = result.get("primer_hints")
            if primer_hints:
                st.markdown("#### Primer design hints")
                primer_cols = st.columns(2)
                primer_cols[0].metric("5' primer Tm", f"{primer_hints['forward_tm']:.1f} °C", "GC clamp: yes" if primer_hints["forward_gc_clamp"] else "GC clamp: no")
                primer_cols[1].metric("3' primer Tm", f"{primer_hints['reverse_tm']:.1f} °C", "GC clamp: yes" if primer_hints["reverse_gc_clamp"] else "GC clamp: no")
                st.caption(f"Wallace estimate over 20 bp candidates. Forward: `{primer_hints['forward_sequence']}`; reverse: `{primer_hints['reverse_sequence']}`")

    # ── Tab 2: Similarity ──────────────────────────────────────────────────────
    with tabs[1]:
        st.markdown("#### Database Similarity Search")

        similarity_source = result.get("similarity_search_source", "local_database")
        similarity_candidate_count = result.get("similarity_candidate_count")
        similarity_candidate_pool_count = result.get("similarity_candidate_pool_count")
        similarity_elapsed_seconds = result.get("similarity_elapsed_seconds")
        similarity_prefiltered_count = result.get("similarity_prefiltered_count", 0)
        similarity_search_mode = result.get("similarity_search_mode", "Balanced")
        info_lines = [f"**Search mode:** `{similarity_search_mode}`"]
        if similarity_source:
            info_lines.append(f"**Source:** `{similarity_source}`")
        if similarity_candidate_count is not None:
            info_lines.append(f"**Candidates evaluated:** `{similarity_candidate_count}`")
        if similarity_candidate_pool_count is not None:
            info_lines.append(f"**Candidate pool:** `{similarity_candidate_pool_count}`")
        if similarity_elapsed_seconds is not None:
            info_lines.append(f"**Similarity workflow time:** `{similarity_elapsed_seconds:.3f} s`")
        if similarity_prefiltered_count:
            info_lines.append(f"**Skipped by prefilter:** `{similarity_prefiltered_count}`")
        if info_lines:
            st.markdown(" — ".join(info_lines))

        if low_complexity.get("regions"):
            st.warning(
                f"{low_complexity['coverage_pct']:.1f}% of this sequence is repetitive/low complexity; "
                "similarity matches involving these regions may not reflect true homology."
            )

        skipped_reason = result.get("similarity_skipped_reason")
        if skipped_reason == "sequence_too_long":
            st.warning(
                "Similarity search was not run because the sequence exceeds the alignment length threshold. "
                "This is not the same as no matches being found."
            )
        elif skipped_reason == "alignment_cost_too_high":
            st.warning(
                "Similarity search was not run because the estimated cost across all candidates is too high. "
                "This is not the same as no matches being found."
            )
        elif not similarity_results:
            st.warning("Similarity search completed but found no matches.")
        else:
            st.plotly_chart(
                viz.plot_similarity_scores(similarity_results),
                width='stretch',
            )

            if best_match:
                best_class = sim.classify_similarity(best_match["similarity_score"])
                st.markdown(
                    f"**Result confidence:** {best_class['emoji']} "
                    f"{best_class['label']} — {best_class['interpretation']}"
                )

            # Enhanced similarity analysis: top 3 comparison & confidence overview
            if len(similarity_results) >= 2:
                st.markdown("---")
                st.markdown("##### Top Matches Summary")
                top3_table = viz.build_top3_comparison_table(similarity_results, len(result.get("sequence", "")))
                if top3_table.get("rows"):
                    st.markdown("| Rank | Gene | Similarity | Trait | Organism | Coverage | Gaps |")
                    st.markdown("|------|------|-----------|-------|----------|----------|------|")
                    for row in top3_table["rows"]:
                        st.markdown(
                            f"| {row['rank']} | {row['gene']} | {row['similarity']} | "
                            f"{row['trait']} | {row['organism']} | {row['coverage']} | {row['gaps']} |"
                        )

            for i, match in enumerate(similarity_results):
                classification = sim.classify_similarity(match["similarity_score"])
                # Clean up gene name display: remove leading underscore-tag tokens
                # (e.g. _arr_, _arrow_) that can appear as artifact prefixes
                raw_name = match.get("gene_name", "") or ""
                clean_name = re.sub(r"^_[a-z]{2,20}[_-]", "", raw_name, flags=re.IGNORECASE)
                # Fallback to original if cleaning produced empty string
                display_name = clean_name if clean_name else raw_name
                with st.expander(
                    f"{classification['emoji']}  {display_name} — {match['similarity_score']:.1f}% similarity"
                ):
                    c1, c2 = st.columns(2)
                    with c1:
                        st.markdown(f"**Gene:** {display_name}")
                        st.markdown(f"**Trait:** {match['trait']}")
                        st.markdown(f"**Organism:** {match['organism']}")
                        st.markdown(f"**Accession:** {match['accession']}")
                    with c2:
                        st.markdown(f"**Similarity (aligned identity):** {match['similarity_score']:.1f}%")
                        st.markdown(f"**Alignment:** {match.get('alignment_method', 'global')}")
                        if match.get("alignment", {}).get("algorithm"):
                            st.markdown(f"**Algorithm:** {match['alignment']['algorithm']}")
                        st.markdown(f"**Level:** {classification['label']}")
                        st.markdown(f"**Interpretation:** {classification['interpretation']}")
                        st.markdown(f"**Description:** {match['description']}")

                    if match.get("alignment"):
                        st.markdown("**Alignment Map:**")
                        st.plotly_chart(
                            viz.plot_alignment(match["alignment"]["alignment_map"]),
                            width='stretch',
                            key=f"alignment_{i}",
                        )

                        # Enhanced visualizations (1-5)
                        query_len = len(result.get("sequence", ""))
                        metrics = viz.build_similarity_metrics_table(match, query_len)

                        # Alignment coverage heatmap (3)
                        if match.get("alignment", {}).get("seq1_aligned"):
                            st.plotly_chart(
                                viz.plot_alignment_coverage_heatmap(match, query_len),
                                width='stretch',
                                key=f"coverage_{i}",
                            )

                        # Confidence gauge (5)
                        if metrics:
                            col_conf, col_metrics = st.columns([1, 2])
                            with col_conf:
                                st.plotly_chart(
                                    viz.plot_confidence_gauge(metrics),
                                    width='stretch',
                                    key=f"confidence_{i}",
                                    use_container_width=True,
                                )
                            with col_metrics:
                                st.markdown("**Alignment Metrics (2):**")
                                if metrics:
                                    st.markdown(
                                        f"- **Aligned Length:** {metrics.get('alignment_length', 'N/A')} bp\n"
                                        f"- **Matches:** {metrics.get('matches', 'N/A')}\n"
                                        f"- **Mismatches:** {metrics.get('mismatches', 'N/A')}\n"
                                        f"- **Gaps:** {metrics.get('total_gaps', 'N/A')} ({metrics.get('gap_percent', 0):.1f}%)\n"
                                        f"- **Coverage:** {metrics.get('coverage_percent', 0):.1f}%\n"
                                        f"- **Identity:** {metrics.get('identity_percent', 0):.1f}%"
                                    )

                        # Gene context card (1)
                        context = viz.build_match_context_card(match)
                        st.markdown("**Gene Context (1):**")
                        st.markdown(
                            f"- **Description:** {context.get('description', 'No description')}\n"
                            f"- **Accession:** {context.get('accession', 'N/A')}\n"
                            f"- **Source:** {context.get('source', 'Unknown')}"
                        )

    # ── Tab 3: Mutations ───────────────────────────────────────────────────────
    with tabs[2]:
        st.markdown("#### Mutation Analysis")

        if not mutation_report:
            st.info("No mutation report — run analysis with a database match first.")
        else:
            mc1, mc2, mc3, mc4, mc5 = st.columns(5)
            mc1.metric("Substitutions", mutation_report["total_mutations"])
            mc2.metric("Indels", mutation_report.get("total_indels", 0))
            mc3.metric(
                "Identity (aligned bases only)",
                f"{mutation_report.get('non_gap_identity_percent', mutation_report['identity_percent'])}%",
                help="Matches ÷ compared positions only (gaps excluded) — matches the "
                     "'Compared positions' count shown below.",
            )
            mc4.metric(
                "Identity (full alignment)",
                f"{mutation_report['identity_percent']}%",
                help="Matches ÷ full alignment length, gaps included in the denominator "
                     "(BLAST-style) — will read lower than the aligned-bases-only identity "
                     "whenever there are indels, even with zero substitutions.",
            )
            mc5.metric("Compared positions (no gaps)", f"{mutation_report['compared_length']}")

            st.plotly_chart(
                viz.plot_mutation_map(mutation_report, mutation_report["compared_length"]),
                width='stretch',
            )

            if mutation_report.get("alignment"):
                aln_data = mutation_report["alignment"]
                match_line = "".join(
                    "|" if a == b and a != "-" else (" " if a == "-" or b == "-" else "X")
                    for a, b in zip(aln_data["query_aligned"], aln_data["reference_aligned"])
                )
                st.markdown(f"**Alignment method:** {aln_data.get('algorithm', 'Needleman-Wunsch')}")
                st.plotly_chart(
                    viz.plot_alignment({
                        "query": aln_data["query_aligned"],
                        "reference": aln_data["reference_aligned"],
                        "match_line": match_line,
                    }),
                    width='stretch',
                )

            mutations = mutation_report.get("mutations", [])
            indels = mutation_report.get("indels", [])
            if mutations:
                st.markdown("#### Substitutions")
                st.markdown(
                    "| Ref pos | Query pos | Reference | Query | Type |\n"
                    "|---|---|---|---|---|\n" +
                    "\n".join(
                        f"| {m['position_reference']} | {m['position_query']} | `{m['reference']}` | `{m['query']}` | {m['type'].capitalize()} |"
                        for m in mutations[:50]
                    )
                )
                if len(mutations) > 50:
                    st.info(f"Showing first 50 of {len(mutations)} substitutions.")
            if indels:
                st.markdown("#### Indels")
                st.markdown(
                    "| Ref pos | Query pos | Reference | Query | Type |\n"
                    "|---|---|---|---|---|\n" +
                    "\n".join(
                        f"| {m['position_reference']} | {m['position_query']} | `{m['reference']}` | `{m['query']}` | {m['type'].capitalize()} |"
                        for m in indels[:50]
                    )
                )
            if not mutations and not indels:
                st.success("No differences after global alignment — sequences are identical.")

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

            st.markdown("#### All Reading Frames (+ and - strands)")
            all_frames = bio.translate_all_frames(sequence, include_reverse=True)
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
            "Plant Gene Analyzer — Analysis Report",
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

else:
    # ── Welcome screen ──────────────────────────────────────────────────────────
    st.markdown(
        """
        <div class="welcome-panel">
            <p class="welcome-icon">🧬</p>
            <h3 class="welcome-title">Ready to Analyze</h3>
            <p class="welcome-text">
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
            | Database Similarity | Global alignment (Needleman-Wunsch) vs. reference genes |
            | Alignment | Needleman-Wunsch / Smith-Waterman / star MSA |
            | Mutation Detection | Substitutions and indels after global alignment |
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

    st.info(
        "🧪 Looking for alignments, distance matrices, phylogeny, or standalone "
        "protein analysis? See **Outils indépendants** in the sidebar navigation."
    )

