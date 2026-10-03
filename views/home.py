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
import alignment_engine as aln
from organism_reference import get_organism_reference
import similarityengine as sim
import visualization as viz
import export_utils as export_util
import re
import sequence_loader as loader
import config
import pipeline
import trait_research as tr
from i18n import translate, language_selector, current_lang, translate_input_type

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
    except Exception:
        logger.exception("Error loading CSS")
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
    except Exception:
        logger.exception("Error loading background video")


# analyze_sequence_record() and get_alignment_map() now live in pipeline.py
# (see import below) — this keeps the analysis orchestration testable and
# reusable independently of the Streamlit UI.


def _render_variant_table(items: list[dict], limit: int = 50) -> None:
    """Styled replacement for the old raw-markdown pipe table used for
    substitutions and indels. Same st.dataframe + column_config approach
    as the Similarity top-3 table, kept as a shared helper so both tables
    can't drift out of sync in appearance if one is tweaked later.
    """
    def translated_variant_term(value: object) -> str:
        raw = str(value or "unknown")
        normalized = re.sub(r"[^a-z0-9]+", "_", raw.lower()).strip("_")
        return translate(f"results.variant_{normalized}", default=raw.replace("_", " ").capitalize())

    columns = {
        "Ref pos": translate("ui.reference_position"),
        "Query pos": translate("ui.query_position"),
        "Reference": translate("results.reference_label"),
        "Query": translate("results.query_label"),
        "Type": translate("ui.type"),
        "Consequence": translate("ui.consequence"),
        "Codon": translate("results.codon"),
        "Amino acid": translate("results.amino_acid"),
        "Impact": translate("results.impact"),
        "BLOSUM62": "BLOSUM62",
    }
    rows = [
        {
            columns["Ref pos"]: m.get("position_reference", m.get("start_position_reference", "")),
            columns["Query pos"]: m.get("position_query", m.get("start_position_query", "")),
            columns["Reference"]: m.get("reference", m.get("bases", "-") if m.get("type") == "deletion" else "-"),
            columns["Query"]: m.get("query", m.get("bases", "-") if m.get("type") == "insertion" else "-"),
            columns["Type"]: translated_variant_term(m.get("type", "variant")),
            columns["Consequence"]: translated_variant_term(m.get("consequence", "Not classified")),
            columns["Codon"]: (
                f"{m['ref_codon']} → {m['query_codon']}"
                if m.get("ref_codon") and m.get("query_codon")
                else "—"
            ),
            columns["Amino acid"]: (
                f"{m['ref_amino_acid']} → {m['query_amino_acid']}"
                if m.get("ref_amino_acid") and m.get("query_amino_acid")
                else "—"
            ),
            columns["Impact"]: translated_variant_term(m.get("impact_class", m.get("consequence", "Not classified"))),
            columns["BLOSUM62"]: m.get("blosum62_score"),
        }
        for m in items[:limit]
    ]
    st.dataframe(
        pd.DataFrame(rows),
        hide_index=True,
        width='stretch',
        column_config={
            columns["Ref pos"]: st.column_config.NumberColumn(columns["Ref pos"], width="small"),
            columns["Query pos"]: st.column_config.NumberColumn(columns["Query pos"], width="small"),
            columns["Reference"]: st.column_config.TextColumn(columns["Reference"], width="small"),
            columns["Query"]: st.column_config.TextColumn(columns["Query"], width="small"),
            columns["Type"]: st.column_config.TextColumn(columns["Type"], width="medium"),
            columns["Consequence"]: st.column_config.TextColumn(columns["Consequence"], width="medium"),
            columns["Codon"]: st.column_config.TextColumn(columns["Codon"], width="medium"),
            columns["Amino acid"]: st.column_config.TextColumn(columns["Amino acid"], width="medium"),
            columns["Impact"]: st.column_config.TextColumn(columns["Impact"], width="medium"),
            columns["BLOSUM62"]: st.column_config.NumberColumn(columns["BLOSUM62"], width="small"),
        },
    )


def _is_admin_view() -> bool:
    """Gate for developer-only diagnostics (DB metadata count, gene search,
    metadata preview) in the sidebar -- these are debugging aids, not
    something a regular user analyzing a sequence needs to see, and having
    them always-on clutters the sidebar for everyone.

    Checked via a `?admin_key=...` URL query parameter matched against a
    secret configured in `.streamlit/secrets.toml` (key: ADMIN_KEY) or the
    ADMIN_KEY environment variable (fallback for hosts without Streamlit
    secrets support). Deliberately not a real login system -- just enough
    to keep this panel out of a regular user's way once this app is shared
    with anyone besides its developer.

    IMPORTANT: if no ADMIN_KEY is configured at all, this returns True
    (panel visible) to match today's always-on local-dev behavior and
    avoid silently hiding it before anyone's set a key. Set an ADMIN_KEY
    secret BEFORE deploying this app anywhere other people can reach it --
    otherwise this panel stays visible to every visitor.
    """
    configured_key = None
    try:
        configured_key = st.secrets.get("ADMIN_KEY")
    except Exception:
        pass  # no secrets.toml configured at all -- fall through to env var
    if not configured_key:
        configured_key = os.environ.get("ADMIN_KEY")
    if not configured_key:
        return True
    return st.query_params.get("admin_key", "") == configured_key


@st.cache_data(show_spinner=False)
def _cached_analyze(
    record_json: str,
    input_type: str,
    reading_frame: int,
    top_n_matches: int,
    similarity_deep_search: bool,
    reference_sequence: str,
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
        reference_sequence=reference_sequence or None,
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
    st.markdown(f"## 🧬 {translate('ui.app_title')}")
    st.markdown("---")

    with st.expander(translate('ui.about'), expanded=False):
        st.markdown(
            translate('ui.build_about') + "\n" + "\n".join(
                f"- {item}" for item in translate('ui.about_items', default=[])
            )
        )
    st.markdown("---")

    st.markdown(f"### {translate('ui.settings')}")
    st.markdown("#### 🧬 Recherche de similarité")
    top_n_matches = st.slider(
        translate('ui.top_matches'),
        min_value=1,
        max_value=config.MAX_TOP_N_MATCHES,
        value=config.DEFAULT_TOP_N_MATCHES,
        help=translate('ui.top_matches_help', default="Number of best-matching genes to display."),
    )
    similarity_deep_search = st.checkbox(
        translate('ui.deep_search'), value=False,
        help=translate('ui.deep_search_help', default="Disable the alignment length prefilter and evaluate more candidates. This is slower, but increases sensitivity for short or divergent queries."),
    )

    st.markdown("#### 📊 Statistiques")
    window_size = st.slider(
        translate('ui.window_size'),
        min_value=config.MIN_WINDOW_SIZE,
        max_value=config.MAX_WINDOW_SIZE,
        value=config.DEFAULT_WINDOW_SIZE,
        step=5,
        help=translate('ui.window_size_help', default="Window size (bp) for the GC content profile chart."),
    )
    skew_window_size = st.slider(
        translate(
            'ui.skew_window_size',
            default="Sliding window (GC / AT skew)",
        ),
        min_value=config.MIN_SKEW_WINDOW_SIZE,
        max_value=config.MAX_SKEW_WINDOW_SIZE,
        value=config.DEFAULT_SKEW_WINDOW_SIZE,
        step=config.SKEW_WINDOW_STEP,
        help=translate(
            'ui.skew_window_size_help',
            default="Larger windows smooth local noise and make broad compositional asymmetry easier to read.",
        ),
    )

    st.markdown("#### 🔤 Traduction")
    reading_frame = st.selectbox(
        translate('ui.reading_frame'),
        options=config.READING_FRAMES,
        format_func=lambda frame: f"{frame:+d}",
    )

    st.markdown("#### 🔧 Entrée")
    input_type_options = config.SUPPORTED_INPUT_TYPES
    sequence_input_type = st.selectbox(
        translate('ui.input_type'),
        options=input_type_options,
        format_func=lambda value: translate_input_type(value, lang=current_lang()),
        help=translate('ui.input_type_help', default="Choose the sequence type or let the app detect it automatically."),
    )
    st.markdown("---")

    admin_view = _is_admin_view() and st.checkbox(
        translate('ui.show_database_diagnostics', default="Show database diagnostics"),
        value=False,
        help=translate(
            'ui.show_database_diagnostics_help',
            default="Open database counts, search, and metadata preview."
        )
    )
    if admin_view:
        st.markdown(f"### {translate('ui.database')}")

    db = None
    metadata = None
    metadata_available = False

    if get_gene_count is not None and search_gene_metadata is not None:
        try:
            # A cheap COUNT(*) rather than materializing all ~56k rows just
            # to call len() on them. Cached for 5 minutes so a widget
            # interaction elsewhere on the page doesn't re-issue it.
            # Always computed regardless of admin_view -- metadata_available
            # gates real app behavior elsewhere (e.g. line ~629), it's only
            # the *display* of this diagnostic info that's admin-gated.
            total_genes = get_gene_count_cached()
            metadata_available = total_genes > 0

            if admin_view:
                st.success(f"✅ {total_genes} {translate('results.gene_records_available')}")
                st.markdown(translate('ui.metadata_load_help'))

                gene_search = st.text_input(
                    translate('ui.search_gene'),
                    value="",
                    help=translate('ui.metadata_filter_help'),
                )
                if gene_search:
                    query = gene_search.strip()
                    # Server-side ILIKE search (see postgres_utils.search_gene_metadata)
                    # -- only the ~20 rows actually shown ever leave Postgres,
                    # instead of pulling all ~56k rows into Python on every
                    # keystroke and filtering them in a list comprehension.
                    match_count = count_gene_metadata_matches_cached(query)
                    filtered = search_gene_metadata_cached(query, limit=20)
                    st.write(translate('ui.metadata_count_summary', count=len(filtered), match_count=match_count))
                else:
                    filtered = search_gene_metadata_cached("", limit=20)
                    st.info(translate('ui.showing_sample'))

                with st.expander(translate('ui.preview_metadata')):
                    for gene in filtered:
                        symbol = gene.get("symbol", "Unknown")
                        gene_id = gene.get("gene_id", "n/a")
                        trait = ", ".join(gene.get("traits", [])[:3]) or "No trait specified"
                        description = gene.get("description", "No description")
                        st.markdown(f"- **{symbol}** (`{gene_id}`) — {trait} — {description}")

                st.info(translate('ui.full_db_load'))

        except Exception:
            logger.exception("Lightweight gene metadata load failed")
            if admin_view:
                st.warning(translate('ui.metadata_load_error'))
            db = load_gene_database_cached(str(config.DATABASE_PATH))
    else:
        db = load_gene_database_cached(str(config.DATABASE_PATH))

    # Below: real failure states stay visible to everyone (a regular user
    # deserves to know analysis may not work), but the happy-path "X genes
    # loaded" confirmation is admin-only diagnostic noise for a working app.
    if db is not None:
        if not db:
            st.error(f"❌ {translate('ui.no_genes_available')}")
        elif isinstance(db, dict) and db and admin_view:
            fallback_note = " (local fallback; PostgreSQL unavailable)" if not metadata_available else ""
            st.success(f"✅ {len(db)} genes loaded{fallback_note}")
    elif metadata_available:
        if admin_view:
            st.info(translate('ui.metadata_load_help'))
    else:
        st.error(f"❌ {translate('ui.no_genes_available')}")


# ─── Top language control ─────────────────────────────────────────────────────
top_spacer, top_language = st.columns([5, 1])
with top_language:
    language_selector(key="top_language_selector")


# ─── Main header ───────────────────────────────────────────────────────────────
st.markdown(
    f"""
    <div class="hero-panel">
        <h1>🧬 {translate('ui.app_title')}</h1>
        <p class="hero-subtitle">{translate('ui.hero_subtitle')}</p>
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("---")


# ─── Input section ─────────────────────────────────────────────────────────────
col_input, col_demo = st.columns([2, 1])

with col_input:
    st.markdown(f"### {translate('ui.sequence_input')}")

    if sequence_input_type != "DNA":
        st.markdown(
            translate('ui.input_accepts', type=translate_input_type(sequence_input_type, lang=current_lang()))
        )

    uploaded_file = st.file_uploader(
        translate('ui.upload_file'),
        type=["fasta", "fa", "txt"],
        help=translate('ui.file_support_help'),
    )

    raw_sequence = st.text_area(
        translate('ui.paste_sequence'),
        height=140,
        placeholder=translate('ui.sequence_placeholder'),
    )
    reference_sequence = st.text_area(
        translate("ui.reference_sequence"),
        height=90,
        placeholder=translate("ui.reference_sequence_placeholder"),
        help=translate("ui.reference_sequence_help"),
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
            st.info(f"{translate('ui.file_loaded')} {uploaded_file.name}")
    elif raw_sequence:
        records = loader.parse_fasta(raw_sequence)
    
    if records:
        if len(records) > 1:
            record_options = [
                f"{idx + 1}. {r.get('metadata', {}).get('name', r['header'])}"
                for idx, r in enumerate(records)
            ]
            selected_index = st.selectbox(
                translate('ui.select_sequence'),
                options=list(range(len(records))),
                format_func=lambda i: record_options[i],
            )
            analyze_all = st.checkbox(
                translate('ui.analyze_all'),
                value=False,
                help="If checked, all parsed FASTA records will be analyzed in batch.",
            )
            if analyze_all:
                st.success(translate('ui.batch_summary_text', count=len(records)))
            else:
                raw_sequence = records[selected_index]["sequence"]
                st.info(f"{translate('ui.selected_sequence')} {record_options[selected_index]}")
        elif len(records) == 1:
            raw_sequence = records[0]["sequence"]
    
with col_demo:
    st.markdown(f"### {translate('ui.quick_demo')}")
    selected_demo = st.selectbox(
        translate('ui.choose_demo'),
        options=list(DEMO_SEQUENCES.keys()),
        label_visibility="collapsed",
    )
    if selected_demo != "Select a demo…":
        demo = DEMO_SEQUENCES[selected_demo]
        st.markdown(f"*{demo['desc']}*")
        if st.button(translate('ui.load_demo_sequence')):
            raw_sequence = demo["seq"]
            st.session_state["loaded_demo"] = demo["seq"]

    if "loaded_demo" in st.session_state and not raw_sequence:
        raw_sequence = st.session_state["loaded_demo"]

if raw_sequence and sequence_input_type != "Protein":
    preview_dna = bio.clean_sequence(raw_sequence, sequence_type="dna")
    if preview_dna:
        with st.expander(translate('ui.reading_frame_preview')):
            st.caption(translate('ui.reading_frame_caption'))
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

analyze_btn = st.button(f"🔬 {translate('ui.analyze_button')}", type="primary")

st.markdown("---")

st.info(translate('ui.independent_tools_notice'))


# ─── Analysis pipeline ─────────────────────────────────────────────────────────
if analyze_btn or (raw_sequence and "last_result" in st.session_state):

    if analyze_btn and not raw_sequence and not records:
        st.warning(f"⚠️ {translate('errors.no_sequence')}")
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
                            reference_sequence,
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
            logger.exception("Unexpected error during analysis")
            
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
    variant_report = result.get("variant_report") or {}
    from aiinterpreter import AIInterpreter

    interpretation = AIInterpreter(
        stats,
        similarity_results,
        mutation_report,
        lang=current_lang(),
    ).full_report()
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
        st.success(translate("results.batch_analysis_complete", count=len(last_results)))
        st.markdown(
            f"**{translate('ui.batch_summary')}:** "
            + translate("results.batch_summary", count=len(last_results), gc=average_gc)
        )
        if len(last_results) > 1:
            summary_rows = []
            for idx, item in enumerate(last_results, start=1):
                similarity_value = "—"
                if item["best_match"]:
                    similarity_score = item["best_match"]["similarity_score"]
                    similarity_value = f"{similarity_score:.1f}"
                summary_rows.append({
                    translate("ui.sequence"): item["header"],
                    translate("ui.input_type"): translate(f"ai.sequence_type_{item['sequence_type']}"),
                    translate("ui.length"): item["stats"]["length"],
                    translate("results.best_match_label"): item["best_match"]["gene_name"] if item["best_match"] else "—",
                    translate("results.similarity_percent"): similarity_value,
                })
            st.table(summary_rows)
            st.markdown(f"#### {translate('results.batch_details')}")
            for idx, item in enumerate(last_results, start=1):
                item_unit = "aa" if item["sequence_type"] == "protein" else "bp"
                item_type = translate(f"ai.sequence_type_{item['sequence_type']}")
                with st.expander(f"{idx}. {item['header']} — {item_type} ({item['stats']['length']} {item_unit})"):
                    st.markdown(f"- **{translate('results.best_match_label')}:** {item['best_match']['gene_name'] if item['best_match'] else '—'}")
                    st.markdown(f"- **{translate('ui.similarity')}:** {item['best_match']['similarity_score']:.1f}%" if item['best_match'] else f"- **{translate('ui.similarity')}:** —")
                    if item['sequence_type'] == 'dna':
                        st.markdown(f"- **{translate('results.orfs_found')}:** {len(item['orfs'])}")
                        if item['orfs']:
                            st.markdown(f"- **{translate('results.longest_orf')}:** {item['orfs'][0]['length']} bp, {translate('ui.reading_frame')} {item['orfs'][0]['frame']}")
                    else:
                        st.markdown(f"- **{translate('results.protein_weight')}:** {item['stats'].get('molecular_weight', 'N/A')} Da")
                        st.markdown(f"- **{translate('results.estimated_pi')}:** {item['stats'].get('isoelectric_point', 'N/A')}")
                        st.markdown(f"- **{translate('results.hydrophobicity_label')}:** {item['stats'].get('hydrophobicity', 'N/A')}")
    else:
        length_unit = "aa" if sequence_type == "protein" else "bp"
        st.success(translate("results.analysis_complete", length=f"{stats['length']:,}", unit=length_unit))
    
    # ── Export Options ─────────────────────────────────────────────────────────
    st.markdown("---")
    st.markdown(f"### 📥 {translate('ui.export_results')}")
    export_col1, export_col2, export_col3, export_col4, export_col5, export_col6, export_col7 = st.columns(7)
    
    with export_col1:
        if st.button(f"📄 {translate('ui.download_json')}"):
            try:
                json_path = export_util.export_results_json(result)
                with open(json_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        f"📥 {translate('results.json_report')}",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.json",
                        mime="application/json",
                    )
                logger.info(f"JSON export created: {json_path}")
                st.success(translate("results.export_success", format="JSON"))
            except Exception as e:
                logger.error(f"JSON export failed: {e}")
                st.error(translate("results.export_failed", error=e))
    
    with export_col2:
        if st.button(f"📊 {translate('ui.download_csv')}"):
            try:
                csv_path = export_util.export_results_csv(result)
                with open(csv_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        f"📥 {translate('results.csv_report')}",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.csv",
                        mime="text/csv",
                    )
                logger.info(f"CSV export created: {csv_path}")
                st.success(translate("results.export_success", format="CSV"))
            except Exception as e:
                logger.error(f"CSV export failed: {e}")
                st.error(translate("results.export_failed", error=e))
    
    with export_col3:
        if st.button(f"🌐 {translate('ui.download_html')}"):
            try:
                html_path = export_util.export_results_html(result)
                with open(html_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        f"📥 {translate('results.html_report')}",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.html",
                        mime="text/html",
                    )
                logger.info(f"HTML export created: {html_path}")
                st.success(translate("results.export_success", format="HTML"))
            except Exception as e:
                logger.error(f"HTML export failed: {e}")
                st.error(translate("results.export_failed", error=e))
    with export_col4:
        if st.button(f"📑 {translate('ui.download_xlsx')}"):
            try:
                xlsx_path = export_util.export_results_xlsx(result)
                with open(xlsx_path, "rb") as f:
                    st.download_button(
                        f"📥 {translate('results.xlsx_report')}",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    )
                logger.info(f"XLSX export created: {xlsx_path}")
                st.success(translate("results.export_success", format="XLSX"))
            except Exception as e:
                logger.error(f"XLSX export failed: {e}")
                st.error(translate("results.export_failed", error=e))

    with export_col5:
        if st.button(f"🧬 {translate('ui.download_fasta')}"):
            try:
                fasta_path = export_util.export_results_fasta(result)
                with open(fasta_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        f"📥 {translate('results.fasta_sequence')}",
                        f.read(),
                        file_name=f"analysis_{stats['length']}{'aa' if sequence_type == 'protein' else 'bp'}.fasta",
                        mime="text/plain",
                    )
                logger.info(f"FASTA export created: {fasta_path}")
                st.success(translate("results.export_success", format="FASTA"))
            except Exception as e:
                logger.error(f"FASTA export failed: {e}")
                st.error(translate("results.export_failed", error=e))

    with export_col6:
        if sequence_type == "protein":
            st.caption(translate("results.gff3_dna_only"))
        elif st.button(f"🧭 {translate('ui.download_gff3')}"):
            try:
                gff3_path = export_util.export_results_gff3(result)
                with open(gff3_path, "r", encoding="utf-8") as f:
                    st.download_button(
                        f"📥 {translate('results.gff3_annotations')}",
                        f.read(),
                        file_name=f"analysis_{stats['length']}bp.gff3",
                        mime="text/plain",
                    )
                logger.info(f"GFF3 export created: {gff3_path}")
                st.success(translate("results.export_success", format="GFF3"))
            except Exception as e:
                logger.error(f"GFF3 export failed: {e}")
                st.error(translate("results.export_failed", error=e))

    with export_col7:
        methods_paragraph = build_methods_paragraph(
            result, {"gc": gc_reference, "codon": codon_reference, "length": length_reference}
        )
        st.download_button(
            translate('ui.copy_methods_paragraph'),
            methods_paragraph,
            file_name="methods_paragraph.txt",
            mime="text/plain",
            help="Copy the generated methods sentence from the downloaded text.",
        )

    st.markdown("---")

    # ── KPI Metrics ────────────────────────────────────────────────────────────
    st.markdown(f"### {translate('ui.sequence_overview')}")
    if sequence_type == "protein":
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric(translate("results.length_aa"), f"{stats['length']:,}")
        m2.metric(translate("results.unique_residues"), f"{stats.get('unique_residues', 'N/A')}")
        m3.metric(translate("results.most_abundant"), f"{max(dist['counts'], key=dist['counts'].get)}")
        m4.metric(
            translate("results.best_match_label"),
            best_match["gene_name"] if best_match else "—",
            f"{best_match['similarity_score']:.1f}%" if best_match else None,
        )
        m5.metric(
            translate("results.mutations_label"),
            mutation_report["total_mutations"] if mutation_report else "—",
        )
    else:
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric(translate("results.length_bp"), f"{stats['length']:,}")
        m2.metric(translate("results.gc_content_label"), f"{stats['gc_content']}%")
        m3.metric(translate("results.at_content_label"), f"{stats['at_content']}%")
        m4.metric(
            translate("results.best_match_label"),
            best_match["gene_name"] if best_match else "—",
            f"{best_match['similarity_score']:.1f}%" if best_match else None,
        )
        m5.metric(
            translate("results.mutations_label"),
            mutation_report["total_mutations"] if mutation_report else "—",
        )

    if result.get("header_metadata"):
        header_meta = result["header_metadata"]
        header_notes = []
        if header_meta.get("gc"):
            header_notes.append(f"{translate('results.header_gc')}: {header_meta['gc']}")
        if header_meta.get("trait"):
            header_notes.append(f"{translate('results.header_trait')}: {header_meta['trait']}")
        if header_notes:
            st.info(f"**{translate('results.header_annotations')}:** " + ", ".join(header_notes))

    if result.get("metadata_warnings"):
        for warning_msg in result["metadata_warnings"]:
            st.warning(warning_msg)

    st.markdown("---")

    # ── Tabs ───────────────────────────────────────────────────────────────────
    st.markdown(f'<div class="section-heading"><span class="section-index">01</span><span>{translate("ui.analysis_results")}</span></div>', unsafe_allow_html=True)
    tabs = st.tabs([
        translate("ui.statistics"),
        translate("ui.similarity"),
        translate("ui.mutations_tab"),
        translate("ui.translation_tab"),
        translate("ui.ai_interpretation"),
        translate("ui.raw_sequence"),
    ])

    # ── Tab 1: Statistics ──────────────────────────────────────────────────────
    with tabs[0]:
        if sequence_type == "protein":
            st.markdown(f"#### {translate('results.amino_acid_composition')}")

            col1, col2 = st.columns([1, 1])
            with col1:
                st.plotly_chart(viz.plot_amino_acid_bar(dist), width='stretch')
            with col2:
                st.markdown(f"#### {translate('results.protein_statistics')}")
                st.markdown(f"**{translate('results.sequence_length')}:** {stats['length']} aa")
                st.markdown(f"**{translate('results.unique_residues')}:** {stats.get('unique_residues', 'N/A')}")
                st.markdown(f"**{translate('results.residue_diversity')}:** {len([v for v in dist['counts'].values() if v > 0])} / {len(dist['counts'])}")
                st.markdown(f"**{translate('results.most_abundant_residue')}:** {max(dist['counts'], key=dist['counts'].get)}")
                st.markdown(f"**GRAVY:** {protein_stats.get('gravy', 'N/A')}")
                st.markdown(f"**Guruprasad instability index:** {protein_stats.get('instability_index', 'N/A')}")
                st.caption(translate("ui.instability_index_help"))
                st.markdown(f"**Aliphatic index:** {protein_stats.get('aliphatic_index', 'N/A')}")

            if motifs:
                st.markdown(f"#### {translate('results.protein_motifs_found')}")
                for motif in motifs:
                    st.markdown(
                        f"- **{motif['name']}** (`{motif['motif']}`) — "
                        f"Position {motif['start']}–{motif['end']}  "
                        f"Match: `{motif['match']}`"
                    )
            else:
                st.info(translate("results.no_known_protein_motifs"))

        else:
            st.markdown(f"#### {translate('results.nucleotide_composition_gc')}")

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
                    st.caption(translate("results.gc_reference_basis", count=f"{gc_reference['n']:,}", organism=organism))
                else:
                    st.plotly_chart(viz.plot_gc_gauge(stats["gc_content"]), width='stretch')
                    st.info(gc_reference["fallback_reason"])

            st.plotly_chart(
                viz.plot_gc_sliding_window(sequence, window=window_size),
                width='stretch',
            )
            effective_skew_window = min(skew_window_size, len(sequence))
            skew_profile = bio.gc_skew_profile(sequence, window=effective_skew_window)
            if skew_profile:
                st.plotly_chart(
                    viz.plot_gc_skew_profile(skew_profile, window=effective_skew_window),
                    width='stretch',
                )

            methylation = result.get("methylation_context") or {}
            st.markdown(f"#### {translate('results.methylation_context')}")
            if methylation:
                methyl_cols = st.columns(3)
                for column, label, key in zip(methyl_cols, ("CG", "CHG", "CHH"), ("cg", "chg", "chh")):
                    column.metric(label, f"{methylation[key]['pct']:.2f}%", f"{methylation[key]['count']} cytosines")
                st.caption(translate("results.methylation_note"))

            quality = result.get("quality_report", {})
            if quality.get("applicable", True):
                quality_status = translate("results.yes") if quality.get("valid") else translate("results.no")
                reason = f" ({quality.get('reason')})" if quality.get("reason") else ""
                st.markdown(
                    f"**{translate('results.quality_filter')}:** {quality_status}{reason} "
                    f"({quality.get('n_pct', 0):.2f}% N, threshold {quality.get('threshold_pct', 5):.2f}%)"
                )

            if low_complexity.get("regions"):
                st.warning(
                    f"{low_complexity['coverage_pct']:.1f}% of this sequence is repetitive/low complexity; "
                    "similarity matches involving these regions may not reflect true homology."
                )

            if result.get("codon_usage"):
                st.markdown(f"#### {translate('results.codon_usage')}")
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
                    st.info(codon_reference["fallback_reason"] + " " + translate("results.species_comparison_unavailable"))
                if len(sequence) % 3:
                    st.caption(translate("results.trailing_bases_excluded"))

                cai = bio.codon_adaptation_index(sequence, organism_codon_usage) if codon_reference.get("available") else None
                if cai is not None:
                    st.metric(translate("results.cai_label"), f"{cai:.3f}")
                    st.caption(translate("results.cai_help"))

            st.markdown(f"#### {translate('results.detailed_statistics')}")
            stat_col1, stat_col2 = st.columns(2)
            with stat_col1:
                st.markdown(f"""
| {translate('results.property')} | {translate('results.value')} |
|---|---|
| {translate('results.sequence_length')} | `{stats['length']} bp` |
| {translate('results.gc_content_label')} | `{stats['gc_content']}%` |
| {translate('results.at_content_label')} | `{stats['at_content']}%` |
| {translate('results.gc_at_ratio')} | `{stats.get('gc_ratio', 'N/A')}` |
                """)
                if length_reference.get("available"):
                    mean_length = float(length_reference["value"])
                    delta = stats["length"] - mean_length
                    st.markdown(translate("results.length_vs_species", delta=f"{delta:+.0f}", mean=f"{mean_length:.0f}", count=length_reference["n"]))
                else:
                    st.info(length_reference["fallback_reason"])
            with stat_col2:
                st.markdown(f"""
| {translate('results.property')} | {translate('results.value')} |
|---|---|
| {translate('results.coding_length_multiple')} | `{translate('results.yes') if stats['is_coding_length'] else translate('results.no')}` |
| {translate('results.contains_atg')} | `{translate('results.yes') if stats['has_start_codon'] else translate('results.no')}` |
| {translate('results.contains_stop')} | `{translate('results.yes') if stats['has_stop_codon'] else translate('results.no')}` |
| {translate('results.complete_orf_found')} | `{translate('results.yes') if stats.get('has_complete_orf') else translate('results.no')}` |
| {translate('results.count_a')} | `{dist['counts']['A']}` |
| {translate('results.count_t')} | `{dist['counts']['T']}` |
| {translate('results.count_g')} | `{dist['counts']['G']}` |
| {translate('results.count_c')} | `{dist['counts']['C']}` |
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
                        translate("results.orf_complete_caution")
                    )

            if motifs:
                st.markdown(f"#### {translate('results.regulatory_motifs_found')}")
                for motif in motifs:
                    st.markdown(
                        f"- **{motif['name']}** (`{motif['motif']}`) — "
                        f"Position {motif['start']}–{motif['end']}  "
                        f"Match: `{motif['match']}`"
                    )
            else:
                st.info(translate("results.no_regulatory_motifs"))

            restriction_sites = result.get("restriction_sites", [])
            st.markdown(f"#### {translate('results.restriction_sites')}")
            if restriction_sites:
                st.dataframe(pd.DataFrame(restriction_sites), hide_index=True, width="stretch")
            else:
                st.caption(translate("results.no_restriction_sites"))

            primer_hints = result.get("primer_hints")
            if primer_hints:
                st.markdown(f"#### {translate('results.primer_design_hints')}")
                primer_cols = st.columns(2)
                primer_cols[0].metric(translate("results.forward_primer_tm"), f"{primer_hints['forward_tm']:.1f} °C", translate("results.gc_clamp_yes") if primer_hints["forward_gc_clamp"] else translate("results.gc_clamp_no"))
                primer_cols[1].metric(translate("results.reverse_primer_tm"), f"{primer_hints['reverse_tm']:.1f} °C", translate("results.gc_clamp_yes") if primer_hints["reverse_gc_clamp"] else translate("results.gc_clamp_no"))
                st.caption(translate("results.primer_wallace_estimate", forward=primer_hints["forward_sequence"], reverse=primer_hints["reverse_sequence"]))

    # ── Tab 2: Similarity ──────────────────────────────────────────────────────
    with tabs[1]:
        st.markdown(f"#### {translate('results.similarity_search_title')}")
        st.caption(
            translate("results.similarity_explanation")
        )

        skipped_reason = result.get("similarity_skipped_reason")
        similarity_source = result.get("similarity_search_source", "local_database")
        similarity_candidate_count = result.get("similarity_candidate_count")
        similarity_candidate_pool_count = result.get("similarity_candidate_pool_count")
        similarity_candidate_pool_requested = result.get("similarity_candidate_pool_requested")
        similarity_elapsed_seconds = result.get("similarity_elapsed_seconds")
        similarity_prefiltered_count = result.get("similarity_prefiltered_count", 0)
        similarity_search_mode = result.get("similarity_search_mode", "Balanced")
        info_lines = [f"**{translate('results.search_mode')}:** `{similarity_search_mode}`"]
        if similarity_source:
            info_lines.append(f"**{translate('results.source')}:** `{similarity_source}`")
        if similarity_candidate_count is not None and not skipped_reason:
            info_lines.append(f"**{translate('results.candidates_evaluated')}:** `{similarity_candidate_count}`")
        if similarity_candidate_pool_count is not None:
            if similarity_candidate_pool_requested is not None and similarity_candidate_pool_requested > similarity_candidate_pool_count:
                info_lines.append(translate("results.candidate_pool_reduced", count=similarity_candidate_pool_count, requested=similarity_candidate_pool_requested))
            else:
                info_lines.append(translate("results.candidate_pool", count=similarity_candidate_pool_count))
        if similarity_elapsed_seconds is not None:
            elapsed_label = translate("results.candidate_selection_time") if skipped_reason else translate("results.similarity_workflow_time")
            info_lines.append(f"**{elapsed_label}:** `{similarity_elapsed_seconds:.3f} s`")
        if similarity_prefiltered_count:
            info_lines.append(f"**{translate('results.skipped_prefilter')}:** `{similarity_prefiltered_count}`")
        if info_lines:
            st.markdown(" — ".join(info_lines))

        if low_complexity.get("regions"):
            st.warning(
                translate("results.low_complexity_warning", coverage=f"{low_complexity['coverage_pct']:.1f}")
            )

        if skipped_reason == "sequence_too_long":
            st.warning(translate("results.sequence_too_long_warning", length=f"{len(sequence):,}", limit=f"{config.MAX_ALIGNMENT_SEQUENCE_LENGTH:,}"))
        elif skipped_reason == "alignment_cost_too_high":
            st.warning(translate("results.alignment_cost_warning"))
        elif not similarity_results:
            st.warning(translate("results.no_similarity_matches"))
        else:
            st.plotly_chart(
                viz.plot_similarity_scores(similarity_results),
                width='stretch',
            )

            if best_match:
                best_class = sim.classify_similarity(best_match["similarity_score"])
                confidence_label = translate(f"results.{best_class['level']}_similarity")
                confidence_text = translate(f"results.{best_class['level']}_similarity_interpretation")
                st.markdown(
                    f"**{translate('results.result_confidence')}:** {best_class['emoji']} "
                    f"{confidence_label} — {confidence_text}"
                )

            # Enhanced similarity analysis: top 3 comparison & confidence overview
            if len(similarity_results) >= 2:
                st.markdown("---")
                st.markdown(f"##### {translate('ui.top_matches_summary')}")
                top3_table = viz.build_top3_comparison_table(similarity_results, len(result.get("sequence", "")))
                if top3_table.get("rows"):
                    def _pct_to_float(value) -> float:
                        # Rows from build_top3_comparison_table are pre-formatted
                        # display strings (e.g. "99.9%") — strip everything but
                        # the numeric part so st.dataframe can sort/bar them
                        # instead of treating them as opaque text.
                        try:
                            return float(str(value).replace("%", "").strip())
                        except (TypeError, ValueError):
                            return 0.0

                    table_rows = []
                    for row, match in zip(top3_table["rows"], similarity_results):
                        classification = sim.classify_similarity(match["similarity_score"])
                        table_rows.append({
                            translate("results.rank"): row["rank"],
                            translate("results.confidence"): classification["emoji"],
                            translate("glossary.gene").capitalize(): row["gene"],
                            translate("results.similarity_percent"): _pct_to_float(row["similarity"]),
                            translate("glossary.trait").capitalize(): row["trait"],
                            translate("glossary.organism").capitalize(): row["organism"],
                            translate("results.coverage"): _pct_to_float(row["coverage"]),
                            translate("results.gaps"): _pct_to_float(row["gaps"]),
                        })

                    df_top3 = pd.DataFrame(table_rows)
                    st.dataframe(
                        df_top3,
                        hide_index=True,
                        width='stretch',
                        column_config={
                            translate("results.rank"): st.column_config.NumberColumn(translate("results.rank"), width="small"),
                            translate("results.confidence"): st.column_config.TextColumn("", width="small"),
                            translate("glossary.gene").capitalize(): st.column_config.TextColumn(translate("glossary.gene").capitalize(), width="medium"),
                            translate("results.similarity_percent"): st.column_config.ProgressColumn(
                                translate("results.similarity_percent"),
                                help=translate("results.global_similarity_help"),
                                format="%.1f%%",
                                min_value=0,
                                max_value=100,
                            ),
                            translate("glossary.trait").capitalize(): st.column_config.TextColumn(translate("glossary.trait").capitalize(), width="large"),
                            translate("glossary.organism").capitalize(): st.column_config.TextColumn(translate("glossary.organism").capitalize(), width="medium"),
                            translate("results.coverage"): st.column_config.ProgressColumn(
                                translate("results.global_coverage"),
                                format="%.1f%%",
                                min_value=0,
                                max_value=100,
                            ),
                            translate("results.gaps"): st.column_config.NumberColumn(translate("results.gaps"), format="%.1f%%"),
                        },
                    )

                    # Conservation view across the top matches: an independent
                    # multiple alignment (star_alignment, reference = query),
                    # NOT the pairwise alignment gaps already computed for
                    # ranking -- those differ per-candidate and can't be
                    # stacked into one consistent column-by-column view.
                    # Reuses star_alignment() (alignment_engine.py) and
                    # plot_msa_table() (visualization.py), both already
                    # built for the separate Phylogeny section -- nothing
                    # new to maintain, just wired into Similarity too.
                    conservation_candidates = similarity_results[:5]
                    msa_sequences, msa_labels = [], []
                    query_raw = result.get("sequence", "")
                    if query_raw:
                        msa_sequences.append(query_raw)
                        msa_labels.append("Query")
                    for match in conservation_candidates:
                        ref_aligned = (
                            match.get("alignment", {})
                            .get("alignment_map", {})
                            .get("reference", "")
                        )
                        ref_raw = ref_aligned.replace("-", "")
                        if ref_raw:
                            raw_name = match.get("gene_name", "") or ""
                            clean_name = re.sub(r"^_[a-z]{2,20}[_-]", "", raw_name, flags=re.IGNORECASE)
                            msa_sequences.append(ref_raw)
                            msa_labels.append(clean_name if clean_name else raw_name)

                    if len(msa_sequences) >= 2:
                        st.markdown("---")
                        st.markdown(translate("results.conservation_top_matches", count=len(msa_sequences) - 1))
                        st.caption(translate("results.conservation_alignment_explanation"))
                        try:
                            msa_result = aln.star_alignment(msa_sequences, seq_type=sequence_type)
                        except Exception:
                            logger.exception("Conservation MSA failed for top matches")
                            msa_result = {"error": "MSA failed"}

                        if msa_result.get("aligned_sequences"):
                            aligned = msa_result["aligned_sequences"]
                            full_width = len(aligned[0]) if aligned else 0
                            # Keep the embedded view readable. A wider window is
                            # available through Plotly fullscreen, but 40
                            # columns keeps bases and position labels legible
                            # without requiring that extra interaction.
                            max_cols = 40
                            window = [seq[:max_cols] for seq in aligned]
                            if full_width > max_cols:
                                st.info(translate("results.showing_aligned_columns", shown=max_cols, total=full_width))
                            st.plotly_chart(
                                viz.plot_msa_table(window, labels=msa_labels),
                                width='stretch',
                                key="conservation_msa",
                            )
                            st.caption(translate("results.conservation_score", score=f"{msa_result.get('conservation_score', 0):.1f}"))

            for i, match in enumerate(similarity_results):
                classification = sim.classify_similarity(match["similarity_score"])
                # Clean up gene name display: remove leading underscore-tag tokens
                # (e.g. _arr_, _arrow_) that can appear as artifact prefixes
                raw_name = match.get("gene_name", "") or ""
                clean_name = re.sub(r"^_[a-z]{2,20}[_-]", "", raw_name, flags=re.IGNORECASE)
                # Fallback to original if cleaning produced empty string
                display_name = clean_name if clean_name else raw_name
                localized_level = translate(f"results.{classification['level']}_similarity")
                localized_interpretation = translate(f"results.{classification['level']}_similarity_interpretation")
                with st.expander(
                    f"{classification['emoji']}  {display_name} — {match['similarity_score']:.1f}% {translate('ui.similarity').lower()}"
                ):
                    c1, c2 = st.columns(2)
                    with c1:
                        st.markdown(f"**{translate('glossary.gene').capitalize()}:** {display_name}")
                        st.markdown(f"**{translate('glossary.trait').capitalize()}:** {match['trait']}")
                        st.markdown(f"**{translate('glossary.organism').capitalize()}:** {match['organism']}")
                        st.markdown(f"**{translate('results.accession')}:** {match['accession']}")
                    with c2:
                        st.metric(
                            translate("results.global_similarity"),
                            f"{match['similarity_score']:.2f}%",
                            help=translate("results.global_similarity_help"),
                        )
                        local_coverage = match.get("local_coverage_percent")
                        if local_coverage is not None:
                            st.metric(
                                translate("results.local_coverage"),
                                f"{local_coverage:.2f}%",
                                help=translate("results.local_coverage_help"),
                            )
                            local_identity = match.get("local_identity")
                            if local_identity is not None:
                                st.caption(translate("results.local_identity", identity=f"{local_identity:.2f}"))
                            if local_coverage < 90.0:
                                st.warning(translate("results.partial_match_warning"))
                        st.markdown(f"**{translate('results.alignment_method')}:** {match.get('alignment_method', 'global')}")
                        if match.get("alignment", {}).get("algorithm"):
                            st.markdown(f"**{translate('results.algorithm')}:** {match['alignment']['algorithm']}")
                        st.markdown(f"**{translate('results.level')}:** {localized_level}")
                        st.markdown(f"**{translate('results.interpretation')}:** {localized_interpretation}")
                        st.markdown(f"**{translate('results.description')}:** {match['description']}")

                    if match.get("alignment"):
                        alignment_map = match["alignment"]["alignment_map"]
                        st.markdown(f"**{translate('ui.alignment_map')}:**")
                        st.caption(translate("results.showing_alignment_window"))
                        st.plotly_chart(
                            viz.plot_alignment(alignment_map),
                            width='stretch',
                            key=f"alignment_{i}",
                        )
                        st.plotly_chart(
                            viz.plot_alignment_overview(alignment_map),
                            width='stretch',
                            key=f"alignment_overview_{i}",
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
                                st.markdown(f"**{translate('ui.alignment_metrics')}:**")
                                if metrics:
                                    st.markdown(
                                        f"- **{translate('results.aligned_length')}:** {metrics.get('alignment_length', 'N/A')} bp\n"
                                        f"- **{translate('results.matches')}:** {metrics.get('matches', 'N/A')}\n"
                                        f"- **{translate('results.mismatches')}:** {metrics.get('mismatches', 'N/A')}\n"
                                        f"- **{translate('results.gaps')}:** {metrics.get('total_gaps', 'N/A')} ({metrics.get('gap_percent', 0):.1f}%)\n"
                                        f"- **{translate('results.coverage')}:** {metrics.get('coverage_percent', 0):.1f}%\n"
                                        f"- **{translate('results.identity')}:** {metrics.get('identity_percent', 0):.1f}%"
                                    )

                        # Gene context card (1)
                        context = viz.build_match_context_card(match)
                        st.markdown(f"**{translate('ui.gene_context')}:**")
                        st.markdown(
                            f"- **{translate('results.description')}:** {context.get('description', 'No description')}\n"
                            f"- **{translate('results.accession')}:** {context.get('accession', 'N/A')}\n"
                            f"- **{translate('results.source')}:** {context.get('source', 'Unknown')}"
                        )

    # ── Tab 3: Mutations ───────────────────────────────────────────────────────
    with tabs[2]:
        st.markdown(f"#### {translate('ui.mutation_analysis')}")

        if not mutation_report:
            warnings = result.get("metadata_warnings", [])
            reference_warning = next(
                (warning for warning in warnings if "No close reference found" in warning),
                None,
            )
            if reference_warning:
                st.info(translate("results.no_close_reference"))
            else:
                st.info(translate("results.no_mutation_report"))
        else:
            raw_substitutions = mutation_report.get("mutations", [])
            classified_substitutions = variant_report.get("substitutions") or raw_substitutions
            raw_by_position = {
                (item.get("position_reference"), item.get("position_query")): item
                for item in raw_substitutions
            }
            substitutions = []
            for item in classified_substitutions:
                enriched = dict(item)
                raw_item = raw_by_position.get(
                    (item.get("position_reference"), item.get("position_query")),
                    {},
                )
                enriched.setdefault("type", raw_item.get("type", "variant"))
                substitutions.append(enriched)
            indel_blocks = variant_report.get("indel_blocks") or mutation_report.get("indels", [])
            transitions = sum(1 for item in substitutions if item.get("type") == "transition")
            transversions = sum(1 for item in substitutions if item.get("type") == "transversion")
            mutation_rate = mutation_report.get("mutation_rate_percent", 0)

            st.info(
                translate(
                    "ui.mutation_summary",
                    substitutions=len(substitutions),
                    indels=len(indel_blocks),
                    rate=mutation_rate,
                )
            )
            if result.get("mutation_reference_source") == "explicit_reference":
                st.caption(translate("ui.mutation_reference_explicit"))

            with st.expander(translate("ui.understand_results")):
                st.markdown(
                    f"- **{translate('ui.transitions')}** : {translate('ui.transition_definition')}\n"
                    f"- **{translate('ui.transversions')}** : {translate('ui.transversion_definition')}\n"
                    f"- **{translate('ui.indels')}** : {translate('ui.indel_definition')}\n"
                    f"- **Identity** : {translate('ui.identity_definition')}"
                )

            summary_cols = st.columns(5)
            summary_cols[0].metric(translate('ui.substitutions'), len(substitutions))
            summary_cols[1].metric(translate('ui.indels'), len(indel_blocks))
            summary_cols[2].metric(translate("ui.transitions"), transitions)
            summary_cols[3].metric(translate("ui.transversions"), transversions)
            summary_cols[4].metric(translate("ui.substitution_rate"), f"{mutation_rate:.2f}%")

            frameshift_count = sum(1 for item in indel_blocks if item.get("frameshift"))
            consequences = {}
            for item in substitutions:
                consequence = item.get("consequence", "unknown")
                consequences[consequence] = consequences.get(consequence, 0) + 1
            if frameshift_count:
                st.warning(translate("results.frameshift_warning", count=frameshift_count))
            elif consequences:
                important = sum(
                    count for name, count in consequences.items()
                    if name in {"missense", "nonsense", "readthrough", "radical", "downstream_of_frameshift"}
                )
                if important:
                    st.warning(translate("results.functional_mutation_warning", count=important))
                elif all(name == "silent" for name in consequences):
                    st.success(translate("results.silent_mutations_notice"))
                else:
                    st.info(translate("results.uncertain_mutation_impact"))

            important_only = st.checkbox(
                translate("ui.important_mutations_only"),
                value=False,
                help=translate("ui.important_mutations_help"),
            )
            important_consequences = {"missense", "nonsense", "readthrough", "radical", "downstream_of_frameshift"}
            displayed_substitutions = [
                item for item in substitutions
                if not important_only or item.get("consequence") in important_consequences
            ]
            displayed_indels = [
                {**item, "consequence": "frameshift" if item.get("frameshift") else "in frame"}
                for item in indel_blocks
                if not important_only or item.get("frameshift")
            ]

            identity_cols = st.columns(3)
            identity_cols[0].metric(
                translate("ui.identity_aligned"),
                f"{mutation_report.get('non_gap_identity_percent', mutation_report['identity_percent'])}%",
                help=translate("ui.identity_aligned_help"),
            )
            identity_cols[1].metric(
                translate("ui.identity_full"),
                f"{mutation_report['identity_percent']}%",
                help=translate("ui.identity_full_help"),
            )
            identity_cols[2].metric(translate('ui.compared_positions'), f"{mutation_report['compared_length']}")

            export_cols = st.columns(2)
            with export_cols[0]:
                mutation_csv_path = export_util.export_mutations_csv(result)
                with open(mutation_csv_path, "r", encoding="utf-8") as handle:
                    st.download_button(
                        translate("ui.download_mutations_csv"),
                        handle.read(),
                        file_name="mutations.csv",
                        mime="text/csv",
                    )
            with export_cols[1]:
                mutation_vcf_path = export_util.export_mutations_vcf(result)
                with open(mutation_vcf_path, "r", encoding="utf-8") as handle:
                    st.download_button(
                        translate("ui.download_substitutions_vcf"),
                        handle.read(),
                        file_name="mutations.vcf",
                        mime="text/plain",
                        help=translate("ui.vcf_help"),
                    )

            st.plotly_chart(
                viz.plot_mutation_map(
                    {**mutation_report, "indel_blocks": indel_blocks},
                    max(mutation_report["query_length"], mutation_report["reference_length"]),
                    labels={
                        "title": translate("ui.mutation_map_title"),
                        "transition": f"{translate('ui.transitions')} ({translate('ui.yellow')})",
                        "transversion": f"{translate('ui.transversions')} ({translate('ui.red')})",
                        "indel": f"{translate('ui.indels')} ({translate('ui.gray')})",
                        "legend_title": translate("ui.substitution_type"),
                        "xaxis": f"{translate('ui.position')} (bp)",
                        "position": translate("ui.position"),
                        "reference_position": translate("ui.reference_position"),
                        "query_position": translate("ui.query_position"),
                        "change": translate("ui.change"),
                        "type": translate("ui.type"),
                        "event": translate("ui.event"),
                        "bases": translate("ui.bases"),
                        "consequence": translate("ui.consequence"),
                        "frameshift": translate("ui.frameshift"),
                        "in_frame": translate("ui.in_frame"),
                    },
                ),
                width='stretch',
            )

            window_size = st.slider(
                translate("ui.variant_window"),
                min_value=50,
                max_value=2000,
                value=500,
                step=50,
                help=translate("ui.variant_window_help"),
            )
            sequence_length = max(mutation_report["query_length"], mutation_report["reference_length"])
            frequency_rows = []
            for start in range(1, sequence_length + 1, window_size):
                end = min(sequence_length, start + window_size - 1)
                substitution_count = sum(
                    start <= (item.get("position_query") or item.get("position_reference", 0)) <= end
                    for item in substitutions
                )
                indel_count = sum(
                    start <= (item.get("start_position_query") or item.get("start_position_reference", 0)) <= end
                    for item in indel_blocks
                )
                window_consequences = sorted({
                    str(item.get("consequence", "unknown")).replace("_", " ")
                    for item in substitutions
                    if start <= (item.get("position_query") or item.get("position_reference", 0)) <= end
                })
                frequency_rows.append({
                    "Region (bp)": f"{start}-{end}",
                    "Substitutions": substitution_count,
                    "Indels": indel_count,
                    "Variants": substitution_count + indel_count,
                    "Consequences": ", ".join(window_consequences) or "—",
                })
            if frequency_rows:
                with st.expander(translate("ui.variant_frequency")):
                    st.dataframe(pd.DataFrame(frequency_rows), hide_index=True, width="stretch")

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

            if displayed_substitutions:
                st.markdown(f"#### {translate('ui.substitutions')}")
                _render_variant_table(displayed_substitutions)
                if len(displayed_substitutions) > 50:
                    st.info(translate("results.showing_substitutions", count=len(displayed_substitutions)))
            if displayed_indels:
                st.markdown(f"#### {translate('ui.indels')}")
                _render_variant_table(displayed_indels)
                if len(displayed_indels) > 50:
                    st.info(translate("results.showing_indels", count=len(displayed_indels)))
            if not displayed_substitutions and not displayed_indels:
                st.success(translate("results.no_differences_after_alignment"))

    # ── Tab 4: Translation ─────────────────────────────────────────────────────
    with tabs[3]:
        if sequence_type == "protein":
            st.markdown(f"#### {translate('results.sequence_type_protein')}")
            st.info(
                translate("results.sequence_type_protein_details")
            )
            st.markdown(f"#### {translate('results.protein_properties')}")
            st.markdown(f"**{translate('results.sequence_length')}:** {stats['length']} aa")
            st.markdown(f"**{translate('results.unique_residues')}:** {stats.get('unique_residues', 'N/A')}")
            st.markdown(f"**{translate('results.most_abundant_residue')}:** {max(dist['counts'], key=dist['counts'].get)}")
        else:
            st.markdown(f"#### {translate('results.selected_frame', frame=f'{reading_frame:+d}')}")
            tl = translation
            codon_count = len(tl.get("codons", []))
            status_label = translate("results.complete_orf_status") if tl["status"] == "complete" else translate("results.open_translation_status")
            m1, m2, m3, m4 = st.columns(4)
            m1.metric(translate("results.protein_length_short"), f"{tl['length']} aa")
            m2.metric(translate("results.complete_codons"), codon_count)
            m3.metric(translate("results.translation_status"), translate("results.complete") if tl["status"] == "complete" else translate("results.open"))
            m4.metric(translate("results.remaining_bases"), tl.get("remainder_nucleotides", 0))
            st.caption(status_label)
            ambiguous_count = sum(1 for base in sequence.upper() if base not in {"A", "T", "G", "C"})
            if ambiguous_count:
                st.warning(
                    translate("results.ambiguous_nucleotides", count=ambiguous_count)
                )

            if tl["protein"]:
                st.markdown(f"**{translate('results.protein_sequence')}**")
                st.code(tl.get("protein_with_stop", tl["protein"]), language=None)
            else:
                st.warning(translate("results.no_translation"))

            codon_rows = bio.translation_codon_rows(sequence, frame=reading_frame)
            if codon_rows:
                st.markdown(f"**{translate('results.codon_map')}**")
                st.dataframe(
                    pd.DataFrame(codon_rows).rename(columns={
                        "codon_index": translate("results.codon_number"), "start": translate("results.start_nt"), "end": translate("results.end_nt"),
                        "codon": translate("results.codon"), "amino_acid": translate("results.amino_acid"), "is_stop": translate("results.stop"),
                    }),
                    hide_index=True,
                    width="stretch",
                )
                st.download_button(
                    translate("results.download_translated_fasta"),
                    f">translated_frame_{reading_frame:+d}\n{tl.get('protein_with_stop', tl['protein'])}\n",
                    file_name=f"translated_frame_{reading_frame:+d}.fasta",
                    mime="text/plain",
                    key="translation_fasta",
                )

            st.markdown(f"#### {translate('results.six_frame_comparison')}")
            all_frames = bio.translate_all_frames(sequence, include_reverse=True)
            recommended_frame = max(
                all_frames.items(),
                key=lambda item: (item[1]["status"] == "complete", item[1]["length"]),
            )[0]
            st.info(
                translate("results.recommended_frame", frame=recommended_frame)
            )
            frame_rows = []
            for frame_name, frame_result in all_frames.items():
                frame_rows.append({
                    translate("results.frame"): frame_name,
                    translate("results.strand"): translate(f"results.strand_{frame_result.get('strand', 'forward')}").title(),
                    translate("results.protein_aa"): frame_result["length"],
                    translate("results.stop"): translate("results.yes") if frame_result["status"] == "complete" else translate("results.no"),
                    translate("results.complete_codons"): len(frame_result.get("codons", [])),
                    translate("results.remaining_bases"): frame_result.get("remainder_nucleotides", 0),
                })
            st.dataframe(pd.DataFrame(frame_rows), hide_index=True, width="stretch")
            for frame_name, frame_result in all_frames.items():
                expanded = frame_name == f"Frame {reading_frame:+d}"
                with st.expander(f"{frame_name} — {frame_result['length']} aa", expanded=expanded):
                    st.code(frame_result.get("protein_with_stop", frame_result["protein"]) or "(empty)", language=None)
                    st.caption(translate("results.complete_orf_status") if frame_result["status"] == "complete" else translate("results.open_translation_status"))

            st.markdown(f"#### {translate('results.predicted_orfs')}")
            orf_rows = [
                {
                    translate("results.strand_frame"): orf["frame"],
                    translate("results.start_nt"): orf["start"],
                    translate("results.end_nt"): orf["end"],
                    translate("results.length_nt"): orf["length"],
                    translate("results.protein_aa"): len(str(orf.get("protein", "")).replace("...[truncated]", "")),
                    translate("results.complete"): translate("results.yes") if orf["complete"] else translate("results.no"),
                }
                for orf in result.get("orfs", [])[:100]
            ]
            if orf_rows:
                st.dataframe(pd.DataFrame(orf_rows), hide_index=True, width="stretch")
                st.markdown(f"**{translate('results.orf_map')}**")
                sequence_length = max(len(sequence), 1)
                for index, orf in enumerate(result.get("orfs", [])[:20], start=1):
                    left = max(0.0, (int(orf["start"]) - 1) / sequence_length * 100)
                    width = max(1.0, int(orf["length"]) / sequence_length * 100)
                    strand_color = "#7A8B5C" if str(orf["frame"]).startswith("+") else "#B8873B"
                    st.markdown(
                        f"<div style='margin:0.25rem 0; color:#EDEAE0; font-size:0.85rem;'>"
                        f"<span style='display:inline-block;width:7rem;'>{orf['frame']} · {orf['start']}-{orf['end']}</span>"
                        f"<span style='display:inline-block;position:relative;width:calc(100% - 7rem);height:1.2rem;background:rgba(255,255,255,.08);'>"
                        f"<span style='position:absolute;left:{left:.2f}%;width:{width:.2f}%;min-width:8px;height:100%;background:{strand_color};' title='ORF {index}'></span>"
                        f"</span></div>",
                        unsafe_allow_html=True,
                    )
                gff_rows = [
                        f"{result.get('header', 'sequence')}\tPlantGeneAnalyzer\tORF\t{orf['start']}\t{orf['end']}\t.\t{ '+' if str(orf['frame']).startswith('+') else '-' }\t.\tID=orf_{idx + 1};frame={orf['frame']}"
                    for idx, orf in enumerate(result.get("orfs", [])[:100])
                ]
                st.download_button(
                    translate("results.download_orfs_gff3"),
                    "##gff-version 3\n" + "\n".join(gff_rows) + "\n",
                    file_name="predicted_orfs.gff3",
                    mime="text/plain",
                    key="translation_gff3",
                )
            else:
                st.info(translate("results.no_orfs"))

            st.markdown(f"#### {translate('results.complementary_sequences')}")
            comp_col1, comp_col2 = st.columns(2)
            with comp_col1:
                st.markdown(f"**{translate('results.complement_orientation')}:**")
                st.code(bio.complement(sequence[:80]) + ("…" if len(sequence) > 80 else ""), language=None)
            with comp_col2:
                st.markdown(f"**{translate('results.reverse_complement')}:**")
                st.code(bio.reverse_complement(sequence[:80]) + ("…" if len(sequence) > 80 else ""), language=None)
            st.caption(translate("results.translation_disclaimer"))

    # ── Tab 5: AI Interpretation ───────────────────────────────────────────────
    with tabs[4]:
        st.markdown(f"#### {translate('ui.ai_biological_interpretation')}")

        interp = interpretation

        # Overall summary
        st.info(f"**{translate('ai.summary_label')}:** {interp['overall_summary']}")

        # Confidence badge
        conf = interp["confidence_level"]
        conf_colors = {"High": "🟢", "Medium": "🟡", "Low": "🔴"}
        confidence_key = str(conf["level"]).lower()
        localized_confidence = translate(f"ai.confidence_level_{confidence_key}")
        st.markdown(
            f"**{translate('ai.confidence_label')}:** {conf_colors.get(conf['level'], '⚪')} {localized_confidence} — {conf['note']}"
        )

        st.markdown("---")

        # Two-column layout
        left, right = st.columns(2)

        with left:
            st.markdown(f"##### {translate('ai.sequence_profile_title')}")
            profile = interp["sequence_profile"]
            for note in profile["notes"]:
                st.markdown(f"- {note}")
            st.markdown(f"*{translate('ai.coding_potential_label')}: **{profile['coding_potential'].upper()}***")

            st.markdown(f"##### {translate('ai.gc_analysis_title')}")
            gc_interp = interp["gc_interpretation"]
            for line in gc_interp["interpretation"]:
                st.markdown(f"- {line}")
            st.markdown(f"*{gc_interp['stress_implication']}*")

            st.markdown(f"##### {translate('ai.functional_prediction_title')}")
            func = interp["functional_prediction"]
            for p in func["predictions"]:
                st.markdown(f"- {p}")

        with right:
            st.markdown(f"##### {translate('ai.similarity_title')}")
            sim_interp = interp["similarity_interpretation"]
            for line in sim_interp.get("interpretation", ["—"]):
                st.markdown(f"- {line}")

            st.markdown(f"##### {translate('ai.mutation_title')}")
            mut_interp = interp["mutation_interpretation"]
            for line in mut_interp.get("interpretation", ["—"]):
                st.markdown(f"- {line}")

            st.markdown(f"##### {translate('ai.stress_title')}")
            stress = interp["stress_resistance"]
            detected = stress.get("detected_resistance", {})
            if detected:
                for stress_type, detail in detected.items():
                    localized_stress = translate(f"ai.trait_{stress_type}", default=stress_type)
                    st.markdown(f"- **{localized_stress.upper()}:** {detail}")
            else:
                st.markdown(f"- {translate('ai.no_specific_stress')}")

        st.markdown("---")
        st.markdown(f"#### {translate('ai.recommendations_title')}")

        recs = interp["agricultural_recommendations"]
        priority_colors = {
            translate("ai.priority_high"): "🔴",
            translate("ai.priority_medium"): "🟡",
            translate("ai.priority_low"): "🔵",
            "HIGH": "🔴", "MEDIUM": "🟡", "LOW": "🔵",
        }
        for rec in recs:
            priority_icon = priority_colors.get(rec["priority"], "⚪")
            with st.expander(f"{priority_icon} [{rec['priority']}] {rec['category']}"):
                st.markdown(rec["recommendation"])

    # ── Tab 6: Raw Sequence ────────────────────────────────────────────────────
    with tabs[5]:
        st.markdown(f"#### {translate('ui.cleaned_sequence')}")
        if sequence_type == "protein":
            st.markdown(
                f"**{translate('results.raw_sequence_length', length=len(sequence), unit='aa')}**  |  "
                f"**{translate('results.valid_amino_acids')}**"
            )
        else:
            st.markdown(
                f"**{translate('results.raw_sequence_length', length=len(sequence), unit='bp')}**  |  "
                f"**{translate('results.gc_short')}:** {stats['gc_content']}%  |  "
                f"**{translate('results.valid_nucleotides')}**"
            )
        st.code(sequence, language=None)

        st.markdown(f"#### {translate('ui.download_report')}")
        if sequence_type == "protein":
            fasta_content = f">Query_sequence | length={len(sequence)}aa\n{sequence}\n"
        else:
            fasta_content = f">Query_sequence | length={len(sequence)}bp | GC={stats['gc_content']}%\n{sequence}\n"
        st.download_button(
            label=translate("results.download_as_fasta"),
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
        if st.button(translate('ui.run_annotation'), key="run_annotation"):
            from core_engines.annotation_engine import annotate_sequence
            try:
                anns = annotate_sequence(sequence, db=db)
                st.success(translate('ui.annotation_complete'))
                st.json(anns)
            except Exception as e:
                logger.error(f"Annotation failed: {e}")
                st.error(f"Annotation failed: {e}")

else:
    # ── Welcome screen ──────────────────────────────────────────────────────────
    st.markdown(
        f"""
        <div class="welcome-panel">
            <p class="welcome-icon">🧬</p>
            <h3 class="welcome-title">{translate('ui.welcome_title')}</h3>
            <p class="welcome-text">{translate('ui.welcome_message')}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.expander(translate('ui.what_does_app_analyze')):
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

    with st.expander(translate('ui.supported_input_formats')):
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
        "protein analysis? See **Independent Tools** in the sidebar navigation."
    )