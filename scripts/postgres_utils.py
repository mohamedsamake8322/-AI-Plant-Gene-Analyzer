#!/usr/bin/env python3
"""
PostgreSQL helper utilities for the plant gene project.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
from urllib.parse import quote_plus

from dotenv import load_dotenv
import psycopg
from psycopg import sql

from quality_rules import MAX_N_RATIO, MIN_SEQUENCE_LENGTH, validate_sequence_quality

try:
    from psycopg_pool import ConnectionPool
except ImportError:  # pragma: no cover
    ConnectionPool = None

ROOT = Path(__file__).resolve().parents[1]
load_dotenv(ROOT / ".env")

# k-mer size for the Postgres-side inverted index (see gene_kmers table /
# populate_kmer_index / find_candidate_genes_by_kmer below). 12 is a
# reasonable default for DNA: 4^12 ≈ 16.7M possible k-mers, long enough
# that unrelated sequences rarely share many by chance, short enough that
# a real homolog still shares plenty. Protein sequences use the same k
# over a 20-letter alphabet, which is far more specific per position.
KMER_K = 12
_DNA_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_PROTEIN_CODE = {c: i for i, c in enumerate("ACDEFGHIKLMNPQRSTVWY")}


def _kmer_hashes(sequence: str, k: int = KMER_K, seq_type: str = "dna") -> set[int]:
    """Deterministic integer encoding of every k-mer in `sequence`.

    Each k-mer is encoded as a base-B number (B=4 for DNA's ACGT, B=20 for
    the standard amino acids), so the same k-mer always hashes to the same
    integer across processes and machines. This matters because the index
    is built once (populate_kmer_index) and looked up later from a
    different process (the Streamlit app, possibly a different worker) —
    Python's built-in hash() is randomized per-process (PYTHONHASHSEED)
    and would silently break every lookup if used here instead.

    Windows containing a character outside the core alphabet (N, ambiguity
    codes, non-standard residues) are skipped rather than mapped to a
    fallback slot, so those low-information k-mers don't create spurious
    matches between otherwise unrelated sequences.
    """
    seq = sequence.upper()
    code = _PROTEIN_CODE if seq_type == "protein" else _DNA_CODE
    base = 20 if seq_type == "protein" else 4
    n = len(seq)
    if n < k:
        return set()

    hashes: set[int] = set()
    for i in range(n - k + 1):
        value = 0
        valid = True
        for ch in seq[i:i + k]:
            c = code.get(ch)
            if c is None:
                valid = False
                break
            value = value * base + c
        if valid:
            hashes.add(value)
    return hashes


def _load_streamlit_secret(name: str) -> str | None:
    try:
        import streamlit as st

        value = st.secrets.get(name)
        return str(value) if value is not None else None
    except Exception:
        return None


def _env_or_secret(name: str) -> str | None:
    return os.getenv(name) or _load_streamlit_secret(name)


DATABASE_URL = _env_or_secret("DATABASE_URL")
DB_HOST = _env_or_secret("DB_HOST")
DB_PORT = _env_or_secret("DB_PORT")
DB_NAME = _env_or_secret("DB_NAME")
DB_USER = _env_or_secret("DB_USER")
DB_PASSWORD = _env_or_secret("DB_PASSWORD")


def _resolve_database_url() -> str:
    database_url = _env_or_secret("DATABASE_URL")
    if database_url:
        return database_url

    db_host = _env_or_secret("DB_HOST")
    db_port = _env_or_secret("DB_PORT")
    db_name = _env_or_secret("DB_NAME")
    db_user = _env_or_secret("DB_USER")
    db_password = _env_or_secret("DB_PASSWORD")
    if db_host and db_name and db_user and db_password:
        return (
            f"postgresql://{quote_plus(db_user)}:{quote_plus(db_password)}@{db_host}"
            f":{db_port or 5432}/{quote_plus(db_name)}"
        )

    raise RuntimeError(
        "DATABASE_URL not set and DB_HOST/DB_NAME/DB_USER/DB_PASSWORD are not all configured"
    )


_pool: ConnectionPool | None = None

# Without these, a network that silently black-holes packets (drops them
# without sending RST/FIN -- common behavior of restrictive institutional
# firewalls/DPI) leaves psycopg waiting forever for a reply that will
# never come: no error, no timeout, just a permanent hang. connect_timeout
# bounds the initial handshake; the keepalive settings make the OS notice
# a dead connection (no ACK after repeated probes) and surface an error
# within roughly keepalives_idle + keepalives_interval * keepalives_count
# seconds instead of hanging indefinitely mid-query.
_CONNECT_KWARGS = {
    "connect_timeout": 10,
    "keepalives": 1,
    "keepalives_idle": 20,
    "keepalives_interval": 10,
    "keepalives_count": 3,
}


class _DirectConnectionPool:
    """Minimal fallback pool wrapper used when psycopg_pool is unavailable."""

    def __init__(self, dsn: str) -> None:
        self._dsn = dsn

    def connection(self) -> psycopg.Connection:
        return psycopg.connect(self._dsn, autocommit=True, **_CONNECT_KWARGS)


def _reset_connection(conn: psycopg.Connection) -> None:
    """Called by the pool whenever a connection is returned, so state left
    behind by one caller can't leak into the next caller that borrows this
    same physical connection. In particular: load_gene_database_from_postgres
    and populate_kmer_index temporarily set autocommit=False (required for
    named/server-side cursors) -- without this reset, that setting would
    otherwise silently persist onto the next unrelated call that borrows the
    same connection from the pool.
    """
    if not conn.closed:
        conn.rollback()
        conn.autocommit = True


def _get_pool() -> ConnectionPool:
    global _pool
    if _pool is None:
        if ConnectionPool is None:
            _pool = _DirectConnectionPool(_resolve_database_url())
        else:
            # min_size/max_size are deliberately modest: Neon/Supabase already
            # pool at the proxy level (e.g. pgbouncer on port 6543), so this
            # app-level pool only needs to be big enough to cover concurrent
            # requests within one Streamlit process, not to be the primary
            # pooling layer. Override via PG_POOL_MAX_SIZE if a deployment
            # genuinely needs more (e.g. many concurrent Streamlit sessions
            # sharing one process).
            max_size = int(os.getenv("PG_POOL_MAX_SIZE", "5"))
            _pool = ConnectionPool(
                _resolve_database_url(),
                min_size=1,
                max_size=max_size,
                kwargs={"autocommit": True, **_CONNECT_KWARGS},
                reset=_reset_connection,
                open=True,
            )
    return _pool


def get_connection():
    """Borrow a pooled connection. Use exactly as before:

        with get_connection() as conn:
            ...

    Returns the connection to the pool on exit instead of closing the
    underlying TCP connection, so repeated calls (common in a Streamlit app,
    which reruns its script on every widget interaction) reuse an existing
    connection instead of paying a fresh connect+auth handshake each time.
    """
    return _get_pool().connection()


def close_pool() -> None:
    """Close the connection pool. Not required for a long-running Streamlit
    process (the pool lives for the process's lifetime), but useful for
    scripts/tests that want a clean shutdown."""
    global _pool
    if _pool is not None:
        _pool.close()
        _pool = None


def create_tables() -> None:
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS genes (
                    id SERIAL PRIMARY KEY,
                    gene_id TEXT,
                    symbol TEXT,
                    organism TEXT,
                    sequence TEXT,
                    sequence_type TEXT,
                    description TEXT,
                    source TEXT,
                    source_url TEXT,
                    external_links JSONB,
                    expression_profiles JSONB,
                    pathways JSONB,
                    publications JSONB,
                    annotations JSONB,
                    traits JSONB,
                    length INTEGER,
                    date_added TIMESTAMPTZ,
                    sequence_hash TEXT,
                    -- "sequence_backed" (has/had a real NCBI/UniProt/KEGG/
                    -- etc. sequence) vs "annotation_only" (a real source
                    -- record exists -- e.g. KEGG pathway, UniProt GO terms,
                    -- PlantTFDB TF family -- but no sequence was ever
                    -- attached to it) vs "plaza_only" (PLAZA enrichment
                    -- data with no sequence attached -- see
                    -- collect_all_sources.py PLAZA block). Similarity/
                    -- BLAST-type queries should filter
                    -- WHERE origin = 'sequence_backed'; relations/orthologs
                    -- queries can use any of the three.
                    origin TEXT DEFAULT 'sequence_backed',
                    -- Orthologs and PLAZA family IDs live here, in their
                    -- own queryable column (previously would only have
                    -- lived inside the now-removed `record` JSONB blob).
                    relations JSONB DEFAULT '{}'::jsonb
                );
                """
            )
            # Both ADD COLUMN IF NOT EXISTS, so re-running create_tables()
            # against an existing (pre-refactor) table backfills these new
            # columns instead of requiring a fresh table.
            cur.execute("ALTER TABLE genes ADD COLUMN IF NOT EXISTS origin TEXT DEFAULT 'sequence_backed';")
            cur.execute("ALTER TABLE genes ADD COLUMN IF NOT EXISTS relations JSONB DEFAULT '{}'::jsonb;")
            cur.execute("CREATE INDEX IF NOT EXISTS genes_origin_idx ON genes (origin);")
            # NULL != NULL in Postgres, so a plain UNIQUE(gene_id) does not
            # prevent duplicates for records identified only by `symbol`
            # (common for GEO/Expression Atlas entries). A unique index on
            # COALESCE(gene_id, symbol) closes that gap and gives us a
            # single, reliable ON CONFLICT target for upserts.
            cur.execute(
                """
                CREATE UNIQUE INDEX IF NOT EXISTS genes_dedup_key_idx
                ON genes (COALESCE(gene_id, symbol));
                """
            )
            # Two different gene_id/symbol values (e.g. the same gene fetched
            # from NCBI under an accession and from UniProt under a different
            # one) can carry the EXACT SAME sequence. The index above cannot
            # catch that, because it dedups on identifier, not on biological
            # content. sequence_hash + find_gene_key_by_sequence_hash() below
            # close that gap: before insert, callers check whether this exact
            # sequence already exists under another key, and if so, redirect
            # the insert to that key so it merges via the upsert above
            # instead of creating a redundant row. Not UNIQUE on purpose —
            # enforcement happens app-side (see dedupe_by_sequence) so a rare
            # race between parallel --workers doesn't hard-fail an insert;
            # a non-unique index here is just for fast lookup.
            cur.execute("CREATE INDEX IF NOT EXISTS genes_sequence_hash_idx ON genes (sequence_hash);")
            cur.execute("CREATE INDEX IF NOT EXISTS genes_organism_idx ON genes (organism);")
            cur.execute("CREATE INDEX IF NOT EXISTS genes_sequence_type_idx ON genes (sequence_type);")
            cur.execute("CREATE INDEX IF NOT EXISTS genes_source_idx ON genes (source);")
            # Supports the length-ratio prefilter in similarityengine.compare_with_database:
            # that filter is currently applied in Python after loading metadata, but this
            # index lets a future range-query version (WHERE length BETWEEN ...) push the
            # same filter down to Postgres instead of scanning all 55k+ rows in the app.
            cur.execute("CREATE INDEX IF NOT EXISTS genes_length_idx ON genes (length);")

            # Trigram indexes so search_gene_metadata's ILIKE '%...%' queries
            # (symbol/gene_id/description search from the sidebar) use an
            # index scan instead of a sequential scan over the whole table —
            # matters once the table is in the tens/hundreds of thousands of
            # rows. Wrapped in try/except: CREATE EXTENSION requires a
            # privilege some managed Postgres roles don't have, and search
            # still works (just slower) without it.
            try:
                cur.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm;")
                cur.execute("CREATE INDEX IF NOT EXISTS genes_symbol_trgm_idx ON genes USING GIN (symbol gin_trgm_ops);")
                cur.execute("CREATE INDEX IF NOT EXISTS genes_gene_id_trgm_idx ON genes USING GIN (gene_id gin_trgm_ops);")
                cur.execute("CREATE INDEX IF NOT EXISTS genes_description_trgm_idx ON genes USING GIN (description gin_trgm_ops);")
                cur.execute("CREATE INDEX IF NOT EXISTS genes_sequence_trgm_idx ON genes USING GIN (sequence gin_trgm_ops);")
            except Exception as e:
                logging.getLogger("postgres_utils").warning(f"pg_trgm indexes not created (missing privilege?): {e}")

            # Inverted k-mer index: (kmer -> gene_key), populated separately
            # by populate_kmer_index() since hashing every sequence is a
            # one-time-ish batch job, not something to redo on every table
            # creation. This is the table that lets similarity search find
            # candidates with a single indexed query instead of ever
            # comparing a submitted sequence against all rows in `genes`.
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS gene_kmers (
                    kmer BIGINT NOT NULL,
                    gene_key TEXT NOT NULL
                );
                """
            )
            cur.execute("CREATE INDEX IF NOT EXISTS gene_kmers_kmer_idx ON gene_kmers (kmer);")
            cur.execute(
                "CREATE UNIQUE INDEX IF NOT EXISTS gene_kmers_unique_idx ON gene_kmers (kmer, gene_key);"
            )
            cur.execute("ALTER TABLE genes ADD COLUMN IF NOT EXISTS kmer_indexed BOOLEAN DEFAULT FALSE;")


_UPSERT_SQL = sql.SQL(
    """
    INSERT INTO genes (
        gene_id, symbol, organism, sequence, sequence_type,
        description, source, source_url, external_links,
        expression_profiles, pathways, publications,
        annotations, traits, length, date_added, sequence_hash,
        origin, relations
    ) VALUES (
        %(gene_id)s, %(symbol)s, %(organism)s, %(sequence)s, %(sequence_type)s,
        %(description)s, %(source)s, %(source_url)s, %(external_links)s,
        %(expression_profiles)s, %(pathways)s, %(publications)s,
        %(annotations)s, %(traits)s, %(length)s, %(date_added)s, %(sequence_hash)s,
        %(origin)s, %(relations)s
    )
    ON CONFLICT (COALESCE(gene_id, symbol))
    DO UPDATE SET
        gene_id = EXCLUDED.gene_id,
        symbol = EXCLUDED.symbol,
        organism = EXCLUDED.organism,
        -- Don't let a source with no sequence (e.g. PlantTFDB's download
        -- fallback, which only has TF metadata) blank out a real sequence
        -- already collected from another source (e.g. NCBI) for the same
        -- gene_id. Only overwrite when the incoming value is non-empty.
        sequence = COALESCE(NULLIF(EXCLUDED.sequence, ''), genes.sequence),
        sequence_type = COALESCE(NULLIF(EXCLUDED.sequence_type, ''), genes.sequence_type),
        description = COALESCE(NULLIF(EXCLUDED.description, ''), genes.description),
        -- Merge instead of overwrite: if this gene_id appeared more than
        -- once in the raw JSON (once per source, before full merging),
        -- each upsert should ADD to the known source list, not replace it
        -- with whatever this particular incoming record happened to carry.
        source = (
            SELECT string_agg(DISTINCT s, ',' ORDER BY s)
            FROM unnest(
                string_to_array(COALESCE(genes.source, ''), ',')
                || string_to_array(COALESCE(EXCLUDED.source, ''), ',')
            ) AS s
            WHERE s != ''
        ),
        source_url = COALESCE(NULLIF(EXCLUDED.source_url, ''), genes.source_url),
        -- JSONB objects: shallow-merge instead of replace, so new keys from
        -- a later source add to what's already known rather than erasing it.
        external_links = genes.external_links || EXCLUDED.external_links,
        annotations = genes.annotations || EXCLUDED.annotations,
        relations = genes.relations || EXCLUDED.relations,
        -- JSONB arrays: keep the existing array if the incoming one is empty.
        expression_profiles = CASE WHEN EXCLUDED.expression_profiles = '[]'::jsonb
            THEN genes.expression_profiles ELSE EXCLUDED.expression_profiles END,
        pathways = CASE WHEN EXCLUDED.pathways = '[]'::jsonb
            THEN genes.pathways ELSE EXCLUDED.pathways END,
        publications = CASE WHEN EXCLUDED.publications = '[]'::jsonb
            THEN genes.publications ELSE EXCLUDED.publications END,
        -- Merge instead of replace: same reasoning as source above -- a
        -- later duplicate upsert for the same gene_id shouldn't blow away
        -- traits collected from an earlier one just because its own
        -- traits list happens to be shorter.
        traits = (
            SELECT COALESCE(jsonb_agg(DISTINCT elem), '[]'::jsonb)
            FROM jsonb_array_elements(genes.traits || EXCLUDED.traits) AS elem
        ),
        -- Don't let a 0/NULL length (no sequence fetched) overwrite a real
        -- known length from a previous insert.
        length = COALESCE(NULLIF(EXCLUDED.length, 0), genes.length),
        date_added = EXCLUDED.date_added,
        sequence_hash = COALESCE(NULLIF(EXCLUDED.sequence_hash, ''), genes.sequence_hash),
        -- Never downgrade a gene's origin when a later pass touches the
        -- same key (shouldn't normally happen given how plaza_only keys
        -- are prefixed, but defensive). Priority order, best to worst:
        -- sequence_backed (real sequence) > annotation_only (real source
        -- record, no sequence) > plaza_only (pure PLAZA ortholog stub).
        origin = CASE
            WHEN genes.origin = 'sequence_backed' OR EXCLUDED.origin = 'sequence_backed'
                THEN 'sequence_backed'
            WHEN genes.origin = 'annotation_only' OR EXCLUDED.origin = 'annotation_only'
                THEN 'annotation_only'
            ELSE EXCLUDED.origin
        END;
    """
)


def sequence_hash(sequence: str | None) -> str | None:
    """Deterministic content hash of a sequence, used to detect the same
    biological sequence arriving under two different gene_id/symbol values
    (e.g. the same gene from NCBI vs UniProt). Normalizes case/whitespace
    first so trivial formatting differences don't produce different hashes.
    Uses plain md5 (via Postgres's built-in md5(), no extension required —
    same reasoning as the pg_trgm privilege fallback above) since this is
    for deduplication, not security.
    """
    if not sequence:
        return None
    normalized = sequence.upper().replace(" ", "").replace("\n", "").strip()
    if not normalized:
        return None
    return hashlib.md5(normalized.encode("utf-8")).hexdigest()


def find_gene_key_by_sequence_hash(seq_hash: str, conn: psycopg.Connection | None = None) -> str | None:
    """Return the canonical gene_key (COALESCE(gene_id, symbol)) of an
    existing row whose sequence hashes to `seq_hash`, or None if no such
    row exists yet. Used by dedupe_by_sequence() before insert.
    """
    if not seq_hash:
        return None

    def _query(c: psycopg.Connection) -> str | None:
        with c.cursor() as cur:
            cur.execute(
                "SELECT COALESCE(gene_id, symbol) FROM genes WHERE sequence_hash = %s LIMIT 1;",
                (seq_hash,),
            )
            row = cur.fetchone()
            return row[0] if row else None

    if conn is not None:
        return _query(conn)
    with get_connection() as conn:
        return _query(conn)


def dedupe_by_sequence(record: dict, conn: psycopg.Connection | None = None) -> dict:
    """Prepare a record for insert, redirecting it to an existing gene's key
    if its exact sequence is already in the database under a different
    gene_id/symbol.

    Call this BEFORE insert_gene_record(). It mutates and returns `record`:
      - always sets record["sequence_hash"]
      - if the sequence already exists under a different key, overwrites
        record["gene_id"] with that existing key (so the upsert's
        ON CONFLICT COALESCE(gene_id, symbol) target merges into the
        existing row instead of creating a duplicate one), and preserves
        the original incoming id as an alternate accession under
        external_links so no traceability is lost.

    This deliberately reuses the existing merge-on-conflict upsert rather
    than adding a second, separate merge path — the sequence-hash match
    just changes *which row* the same upsert logic targets.
    """
    seq, _seq_type = extract_primary_sequence(record)
    seq_hash = sequence_hash(seq)
    record["sequence_hash"] = seq_hash
    if not seq_hash:
        return record

    incoming_key = record.get("gene_id") or record.get("symbol")
    existing_key = find_gene_key_by_sequence_hash(seq_hash, conn=conn)

    if existing_key and existing_key != incoming_key:
        links = dict(record.get("external_links") or {})
        if incoming_key:
            # record here is the raw pre-restructure dict, which carries
            # "sources_summary" (a list, e.g. ["ncbi", "uniprot"]) rather
            # than the flat "source" string that only exists after
            # restructure_to_schema() runs. Mirror that same fallback
            # logic here so this alt_id key reflects the real source(s)
            # instead of always falling through to "unknown_source".
            sources_summary = record.get("sources_summary")
            if sources_summary:
                src = ",".join(sources_summary)
            else:
                src = record.get("source") or "unknown_source"
            links[f"alt_id_{src}"] = incoming_key
        record["external_links"] = links
        record["gene_id"] = existing_key
        # Leave symbol as-is; COALESCE(gene_id, symbol) resolves to gene_id
        # regardless, and keeping symbol lets it still be searched on.

    return record


def extract_primary_sequence(record: dict) -> tuple[str | None, str | None]:
    """
    Returns (sequence_string, sequence_type) for whichever ONE sequence
    this gene's row should carry for search/similarity purposes.

    WHY THIS EXISTS: the pipeline's schema changed (see collect_all_sources.py
    restructure_to_schema()) from one flat sequence string + sequence_type
    per record to a nested {"dna": ..., "rna": ..., "protein": ...} dict,
    since one gene can now legitimately have all three at once. But this
    Postgres table still stores ONE sequence per row (same reasoning as
    before: pg_trgm/kmer similarity search operates on a single sequence
    column) -- so a priority order is needed to pick which one. DNA is
    preferred first since that's what the app's primary search/analysis
    flow (Statistics, Similarity, BLAST) is built around; protein next
    (still useful for similarity search); RNA last since NCBI mRNA/cDNA
    entries in this pipeline are stored with T not U anyway (see
    _VALID_CHARS above) and mostly overlap with DNA content.

    Nothing is lost by only keeping one here: the full nested dict is
    still stored as-is inside the `record` JSONB column (see
    _record_to_params), so any caller that needs all three can still get
    them from there.

    Backward compatible: if `record["sequence"]` is already a plain string
    (old flat format, in case an older JSON file is loaded), it's used
    directly instead of assuming the new dict shape.
    """
    seq_field = record.get("sequence")

    if isinstance(seq_field, str):
        # Old flat format.
        return (seq_field or None), record.get("sequence_type")

    if isinstance(seq_field, dict):
        for seq_type in ("dna", "protein", "rna"):
            value = seq_field.get(seq_type)
            if value:
                return value, seq_type

    return None, None


def _record_to_params(record: dict) -> dict:
    if not record.get("gene_id") and not record.get("symbol"):
        raise ValueError("Record must contain gene_id or symbol")

    sequence, sequence_type = extract_primary_sequence(record)

    # annotation (singular) is the new nested schema's key; annotations
    # (plural) was the old flat one. Support both so an older JSON file
    # doesn't silently lose its GO terms / KEGG data.
    annotation = record.get("annotation") or record.get("annotations") or {}
    kegg_pathways = annotation.get("kegg_pathways") or record.get("pathways") or []
    literature = record.get("literature") or {}
    publications = literature.get("publications") or record.get("publications") or []
    relations = record.get("relations") or {}

    # sources_summary (new, a list) replaces the old single "source" string
    # column's meaning -- join it so the existing TEXT column still gets a
    # readable value ("ncbi,uniprot,plaza") instead of needing a schema
    # change just for this. Falls back to the old flat "source" field.
    sources_summary = record.get("sources_summary")
    if sources_summary:
        source = ",".join(sources_summary)
    else:
        source = record.get("source")

    return {
        "gene_id": record.get("gene_id"),
        "symbol": record.get("symbol"),
        "organism": record.get("organism"),
        "sequence": sequence,
        "sequence_type": sequence_type,
        # New schema has no top-level "description" -- common_name is the
        # closest equivalent (populated from PLAZA's gene_description.csv,
        # see collect_all_sources.py PLAZA block).
        "description": record.get("description") or record.get("common_name"),
        "source": source,
        "source_url": record.get("source_url"),
        "external_links": json.dumps(record.get("external_links", {})),
        # NOTE: the new nested schema (restructure_to_schema) does not
        # currently emit an "expression" field at all -- GEO/Expression
        # Atlas data collected upstream is silently dropped during
        # restructuring. Not fixed here; flagged as a separate follow-up.
        "expression_profiles": json.dumps(record.get("expression_profiles", [])),
        "pathways": json.dumps(kegg_pathways),
        "publications": json.dumps(publications),
        "annotations": json.dumps(annotation),
        "traits": json.dumps(record.get("traits", [])),
        # Orthologs / gene-family IDs from PLAZA -- new column, see
        # create_tables(). Without this, all the PLAZA orthology work
        # would be stored in `record` JSONB only and not easily queryable.
        "relations": json.dumps(relations),
        # Defensive fallback only -- restructure_to_schema() now always
        # sets origin explicitly (sequence_backed / annotation_only /
        # plaza_only), so this default should rarely if ever trigger. It
        # deliberately does NOT default to "sequence_backed": doing so was
        # exactly the bug that silently mislabeled UniProt/KEGG/PlantTFDB
        # records with no real sequence (see collect_all_sources.py fix,
        # 2026-08-22). "annotation_only" is the safe, non-optimistic default
        # when origin is unexpectedly absent.
        "origin": record.get("origin", "annotation_only"),
        "length": record.get("length") or (len(sequence) if sequence else None),
        "date_added": record.get("date_added"),
        "sequence_hash": record.get("sequence_hash") or sequence_hash(sequence),
    }


def is_valid_sequence(
    sequence: str | None, sequence_type: str | None = "dna",
    min_length: int = MIN_SEQUENCE_LENGTH, max_n_ratio: float = MAX_N_RATIO,
) -> tuple[bool, str]:
    """Reject records before they ever reach the insert path: too short to
    be useful for downstream fine-tuning, too many ambiguous bases (N), or
    containing characters that don't belong to the declared sequence type.
    Returns (is_valid, reason) — reason is empty when valid, so callers can
    log/count why a record was skipped.
    """
    return validate_sequence_quality(sequence, sequence_type, min_length, max_n_ratio)


def backfill_sequence_hashes() -> int:
    """One-time migration for genes already in the table before this column
    existed (your current ~56k rows). Computes the hash server-side with a
    single UPDATE (Postgres's built-in md5(), matching sequence_hash()
    above) instead of looping row-by-row from Python. Safe to re-run — only
    touches rows where sequence_hash IS NULL.
    """
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("ALTER TABLE genes ADD COLUMN IF NOT EXISTS sequence_hash TEXT;")
            cur.execute(
                """
                UPDATE genes
                SET sequence_hash = md5(UPPER(REPLACE(REPLACE(sequence, ' ', ''), E'\\n', '')))
                WHERE sequence_hash IS NULL AND sequence IS NOT NULL AND sequence <> '';
                """
            )
            n = cur.rowcount
            cur.execute("CREATE INDEX IF NOT EXISTS genes_sequence_hash_idx ON genes (sequence_hash);")
    return n


def insert_gene_record(record: dict, conn: psycopg.Connection | None = None) -> None:
    """Insert or update a single gene record.

    Pass an existing `conn` (e.g. from `get_connection()`) when inserting
    many records in a loop, to avoid opening a new connection per row.
    """
    params = _record_to_params(record)
    if conn is not None:
        with conn.cursor() as cur:
            cur.execute(_UPSERT_SQL, params)
        return
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(_UPSERT_SQL, params)


def insert_gene_records(records: list[dict], batch_size: int = 500) -> int:
    """Insert or update many gene records over a single reused connection.

    Uses `executemany` in batches instead of opening a fresh connection per
    record, which matters especially against a pooled endpoint like
    Supabase's pgbouncer (port 6543).

    Returns the number of records successfully upserted; records missing
    both `gene_id` and `symbol` are skipped with a warning.
    """
    valid_params = []
    for record in records:
        try:
            valid_params.append(_record_to_params(record))
        except ValueError:
            logging.getLogger("postgres_utils").warning(
                "Skipping record without gene_id or symbol: %r", record.get("description")
            )

    inserted = 0
    with get_connection() as conn:
        with conn.cursor() as cur:
            for start in range(0, len(valid_params), batch_size):
                batch = valid_params[start : start + batch_size]
                cur.executemany(_UPSERT_SQL, batch)
                inserted += len(batch)
    return inserted


def load_gene_database_from_postgres(batch_size: int = 1000) -> dict:
    """Load all gene records from Postgres.

    Paginates with keyset pagination (WHERE id > last_id ORDER BY id LIMIT
    batch_size) instead of a single fetchall() or a long-lived named
    server-side cursor. Two reasons:

    1. Memory: with tens of thousands of rows, an unbounded fetchall() can
       exhaust available memory (this bit us with "out of memory for query
       result" once the dataset passed ~49k rows -- at the time, each row
       also duplicated its data into a since-removed `record` JSONB blob,
       which made it worse).
    2. Network resilience: a named server-side cursor is one long-lived
       connection held open for the entire fetch. On a restrictive network
       (institutional firewall/DPI that silently drops packets rather than
       sending RST), that single long-lived stream can stall forever with
       no error. Short, independent per-batch queries bound the damage: a
       stalled batch hits statement_timeout and raises clearly, instead of
       hanging the whole load with no feedback.
    """
    records: dict[str, dict] = {}
    last_id = 0
    batch_num = 0
    with get_connection() as conn:
        with conn.cursor() as cur:
            # Per-connection safety net: if a single batch query stalls
            # server-side for any reason, fail after 60s with a clear
            # error instead of hanging silently.
            cur.execute("SET statement_timeout = '60s';")
        try:
            while True:
                batch_num += 1
                with conn.cursor() as cur:
                    cur.execute(
                        """
                        SELECT id, gene_id, symbol, organism, sequence, sequence_type,
                               description, source, source_url, external_links,
                               expression_profiles, pathways, publications,
                               annotations, traits, length, date_added,
                               origin, relations
                        FROM genes
                        WHERE id > %s
                        ORDER BY id
                        LIMIT %s;
                        """,
                        (last_id, batch_size),
                    )
                    rows = cur.fetchall()

                if not rows:
                    break

                for row in rows:
                    (
                        row_id,
                        gene_id,
                        symbol,
                        organism,
                        sequence,
                        sequence_type,
                        description,
                        source,
                        source_url,
                        external_links,
                        expression_profiles,
                        pathways,
                        publications,
                        annotations,
                        traits,
                        length,
                        date_added,
                        origin,
                        relations,
                    ) = row

                    if isinstance(external_links, str):
                        external_links = json.loads(external_links)
                    if isinstance(expression_profiles, str):
                        expression_profiles = json.loads(expression_profiles)
                    if isinstance(pathways, str):
                        pathways = json.loads(pathways)
                    if isinstance(publications, str):
                        publications = json.loads(publications)
                    if isinstance(annotations, str):
                        annotations = json.loads(annotations)
                    if isinstance(traits, str):
                        traits = json.loads(traits)
                    if isinstance(relations, str):
                        relations = json.loads(relations)

                    key = gene_id or symbol
                    last_id = row_id
                    if not key:
                        continue

                    records[key] = {
                        "gene_id": gene_id,
                        "symbol": symbol,
                        "organism": organism,
                        "sequence": sequence,
                        "sequence_type": sequence_type,
                        "description": description,
                        "source": source,
                        "source_url": source_url,
                        "external_links": external_links or {},
                        "expression_profiles": expression_profiles or [],
                        "pathways": pathways or [],
                        "publications": publications or [],
                        "annotations": annotations or {},
                        "traits": traits or [],
                        "length": length,
                        "date_added": date_added,
                        "origin": origin,
                        "relations": relations or {},
                    }

                print(f"  ... lot {batch_num} charge ({len(records)} enregistrements au total)")

                if len(rows) < batch_size:
                    break  # dernier lot, pas la peine de refaire un aller-retour
        finally:
            with conn.cursor() as cur:
                cur.execute("RESET statement_timeout;")
            conn.commit()  # read-only, just closes the transaction cleanly


    return records


def load_gene_sequences_by_keys(keys: list[str]) -> dict:
    """Fetch full sequence data for a specific, already-known-small set of
    gene_id/symbol keys.

    This is the counterpart to load_gene_database_metadata_from_postgres():
    the UI loads lightweight metadata for all ~56k genes up front, the
    similarity search prefilters candidates against that metadata (length
    ratio, then k-mer/composition), and only the surviving handful of keys
    (typically tens to low hundreds, not tens of thousands) are looked up
    here with a single indexed query. This avoids ever pulling every full
    sequence over the network/into memory just to score one query sequence.

    A single round trip via `= ANY(%s)` against the same
    `COALESCE(gene_id, symbol)` expression used for upserts (so it hits the
    genes_dedup_key_idx unique index) is used rather than one query per key.
    """
    if not keys:
        return {}

    records: dict[str, dict] = {}
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT gene_id, symbol, organism, sequence, sequence_type,
                       description, source, source_url, traits, length
                FROM genes
                WHERE COALESCE(gene_id, symbol) = ANY(%s);
                """,
                (list(keys),),
            )
            for row in cur.fetchall():
                (
                    gene_id, symbol, organism, sequence, sequence_type,
                    description, source, source_url, traits, length,
                ) = row
                if isinstance(traits, str):
                    traits = json.loads(traits)
                key = gene_id or symbol
                if not key:
                    continue
                records[key] = {
                    "gene_id": gene_id,
                    "symbol": symbol,
                    "organism": organism,
                    "sequence": sequence,
                    "sequence_type": sequence_type,
                    "description": description,
                    "source": source,
                    "source_url": source_url,
                    "traits": traits or [],
                    "length": length,
                }
    return records


def get_gene_count() -> int:
    """Cheap COUNT(*) for the sidebar header — never materializes rows."""
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM genes;")
            return cur.fetchone()[0]


def get_gc_content_stats_for_organism(organism: str) -> dict[str, float | int]:
    """Return mean and standard deviation of GC% for sequence-backed genes."""
    if not organism:
        return {"mean_gc": 0.0, "stdev_gc": 0.0, "n_sequences": 0}
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT AVG(gc_value), STDDEV_SAMP(gc_value), COUNT(*)
                FROM (
                    SELECT 100.0 * (
                        LENGTH(sequence) - LENGTH(REPLACE(sequence, 'G', ''))
                        + LENGTH(sequence) - LENGTH(REPLACE(sequence, 'C', ''))
                    ) / NULLIF(LENGTH(sequence), 0) AS gc_value
                    FROM genes
                    WHERE organism = %s AND sequence IS NOT NULL AND sequence <> ''
                ) AS values;
                """,
                (organism,),
            )
            mean_gc, stdev_gc, n_sequences = cur.fetchone()
    return {
        "mean_gc": round(float(mean_gc or 0.0), 2),
        "stdev_gc": round(float(stdev_gc or 0.0), 2),
        "n_sequences": int(n_sequences or 0),
    }


def get_codon_usage_for_organism(organism: str) -> dict[str, float]:
    """Return per-codon frequencies aggregated server-side for an organism."""
    if not organism:
        return {}
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                WITH codons AS (
                    SELECT SUBSTRING(sequence FROM pos FOR 3) AS codon
                    FROM genes, generate_series(1, LENGTH(sequence) - 2, 3) AS pos
                    WHERE organism = %s AND sequence IS NOT NULL AND sequence <> ''
                ), counts AS (
                    SELECT codon, COUNT(*)::float AS count
                    FROM codons
                    WHERE codon ~ '^[ATGC]{3}$'
                    GROUP BY codon
                )
                SELECT codon, count / NULLIF(SUM(count) OVER (), 0)
                FROM counts;
                """,
                (organism,),
            )
            return {codon: float(frequency) for codon, frequency in cur.fetchall()}


def search_gene_metadata(query: str | None, limit: int = 20, offset: int = 0) -> list[dict]:
    """Server-side metadata search: only `limit` rows are ever pulled into
    the app, regardless of how large the genes table is. Runs ILIKE across
    the indexed text columns (see the pg_trgm GIN indexes in create_tables)
    instead of loading every row into Python and filtering there.
    """
    q = (query or "").strip()
    where_sql = ""
    params: dict = {"limit": limit, "offset": offset}
    if q:
        where_sql = (
            "WHERE gene_id ILIKE %(pattern)s OR symbol ILIKE %(pattern)s "
            "OR description ILIKE %(pattern)s OR source ILIKE %(pattern)s "
            "OR traits::text ILIKE %(pattern)s"
        )
        params["pattern"] = f"%{q}%"

    records: list[dict] = []
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                f"""
                SELECT gene_id, symbol, organism, description, source, traits, length
                FROM genes
                {where_sql}
                ORDER BY id
                LIMIT %(limit)s OFFSET %(offset)s;
                """,
                params,
            )
            for row in cur.fetchall():
                gene_id, symbol, organism, description, source, traits, length = row
                if isinstance(traits, str):
                    traits = json.loads(traits)
                records.append({
                    "gene_id": gene_id,
                    "symbol": symbol,
                    "organism": organism,
                    "description": description,
                    "source": source,
                    "traits": traits or [],
                    "length": length,
                })
    return records


def count_gene_metadata_matches(query: str | None) -> int:
    """Matching row count for the same filter search_gene_metadata applies —
    lets the UI show "showing 20 of N" without pulling all N rows."""
    q = (query or "").strip()
    if not q:
        return get_gene_count()
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT COUNT(*) FROM genes
                WHERE gene_id ILIKE %(pattern)s OR symbol ILIKE %(pattern)s
                   OR description ILIKE %(pattern)s OR source ILIKE %(pattern)s
                   OR traits::text ILIKE %(pattern)s;
                """,
                {"pattern": f"%{q}%"},
            )
            return cur.fetchone()[0]


def populate_kmer_index(
    k: int = KMER_K,
    batch_size: int = 2_000,
    flush_every: int = 200_000,
    rebuild: bool = False,
) -> int:
    """Populate a compact, storage-safe candidate index for the Neon-backed app.

    The previous implementation wrote one row per k-mer into `gene_kmers` for
    every gene, which can exceed Neon's storage quota for a 55k+ gene catalog.
    To keep the app functional within those limits, this function now performs
    a lightweight marking pass: it only updates the `kmer_indexed` flag on the
    genes table and leaves the actual candidate selection to the trigram-based
    lookup implemented in find_candidate_genes_by_kmer().

    Returns the number of genes marked as indexed in this run.
    """
    with get_connection() as conn:
        with conn.cursor() as cur:
            if rebuild:
                cur.execute("UPDATE genes SET kmer_indexed = FALSE;")
                try:
                    cur.execute("TRUNCATE gene_kmers;")
                except Exception:
                    pass

    indexed_genes = 0
    pending_keys: list[str] = []
    progress_every = max(1_000, min(10_000, batch_size * 5))

    def _flush(write_cur) -> None:
        if pending_keys:
            write_cur.execute(
                "UPDATE genes SET kmer_indexed = TRUE WHERE COALESCE(gene_id, symbol) = ANY(%s);",
                (pending_keys,),
            )
        pending_keys.clear()

    with get_connection() as conn:
        conn.autocommit = False
        try:
            with conn.cursor(name="kmer_populate_cursor") as read_cur, conn.cursor() as write_cur:
                read_cur.itersize = batch_size
                read_cur.execute(
                    """
                    SELECT COALESCE(gene_id, symbol) AS gene_key
                    FROM genes
                    WHERE (kmer_indexed IS NOT TRUE) AND sequence IS NOT NULL AND sequence <> '';
                    """
                )
                for (gene_key,) in read_cur:
                    pending_keys.append(gene_key)
                    indexed_genes += 1
                    if indexed_genes % progress_every == 0:
                        print(f"Marked {indexed_genes} gene(s) as indexed (compact trigram prefilter active)...")
                    if len(pending_keys) >= flush_every:
                        _flush(write_cur)
                        conn.commit()
                _flush(write_cur)
                conn.commit()
        except Exception:
            conn.rollback()
            raise

    print("Compact trigram-based candidate prefilter is active; no full k-mer row table was populated.")
    return indexed_genes


def find_candidate_genes_by_kmer(
    query: str | set[int] | list[int],
    limit: int = 200,
    length_ratio: float = 3.0,
) -> list[tuple[str, float]]:
    """Return up to `limit` candidate genes for a query sequence.

    The preferred path uses a compact trigram similarity search on the
    `genes.sequence` column (backed by a GIN trigram index) rather than the
    storage-heavy `gene_kmers` table. This keeps candidate selection working
    on Neon without hitting the storage quota while still avoiding a full
    O(n*m) scan over every gene sequence.

    length_ratio: same length-ratio prefilter concept as
    similarityengine.compare_with_database's max_length_ratio -- candidates
    outside [query_len/length_ratio, query_len*length_ratio] are excluded
    server-side before scoring. Kept as a plain parameter (default 3.0)
    rather than importing config.py directly, so this module stays runnable
    standalone (`python scripts/postgres_utils.py ...`) without depending on
    the project root being on sys.path. Callers that already have config in
    scope should pass config.LENGTH_RATIO_PREFILTER explicitly to keep this
    in sync with similarityengine's own prefilter.

    A legacy list of k-mer hashes is still accepted for backward
    compatibility; in that case it falls back to the old `gene_kmers` table
    if present, otherwise returns an empty list.
    """
    if isinstance(query, str):
        sequence = query.upper().replace(" ", "")
        if not sequence:
            return []
        query_len = len(sequence)
        min_len = max(1, int(query_len / length_ratio))
        max_len = int(query_len * length_ratio)
        with get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    SELECT COALESCE(gene_id, symbol) AS gene_key,
                           similarity(sequence, %(query)s) AS score
                    FROM genes
                    WHERE sequence IS NOT NULL
                      AND sequence <> ''
                      AND sequence %% %(query)s
                      AND length BETWEEN %(min_len)s AND %(max_len)s
                    ORDER BY score DESC
                    LIMIT %(limit)s;
                    """,
                    {"query": sequence, "min_len": min_len, "max_len": max_len, "limit": limit},
                )
                return [(row[0], float(row[1])) for row in cur.fetchall()]

    kmer_list = list(query)
    if not kmer_list:
        return []
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM gene_kmers;")
            if cur.fetchone()[0] == 0:
                return []
            cur.execute(
                """
                SELECT gene_key, COUNT(*) AS shared
                FROM gene_kmers
                WHERE kmer = ANY(%s)
                GROUP BY gene_key
                ORDER BY shared DESC
                LIMIT %s;
                """,
                (kmer_list, limit),
            )
            return [(row[0], float(row[1])) for row in cur.fetchall()]


def find_candidate_genes_by_kmer_exhaustive(
    query: str,
    limit: int = 500,
) -> list[tuple[str, float]]:
    """
    **Étage 1 (Deep Search):** Exhaustive trigram scan WITHOUT length prefilter.
    
    Returns up to `limit` candidates ranked by trigram similarity, considering
    ALL genes in the database regardless of length. This is the intensive stage
    that justifies the "Deep" label — it scans the full 56k gene table (but still
    uses the fast pg_trgm GIN index, not a brute-force O(n*m) alignment).
    
    Used by find_similar_genes_deep() in similarityengine.py as the first stage.
    The second stage (Needleman-Wunsch alignment) is applied only to the top
    candidates returned here.
    
    Args:
        query: DNA/RNA/protein sequence (string)
        limit: max candidates to return (default 500)
    
    Returns:
        List of (gene_key, trigram_similarity_score) sorted descending by score
    """
    sequence = query.upper().replace(" ", "")
    if not sequence:
        return []
    
    with get_connection() as conn:
        with conn.cursor() as cur:
            # No length filter — query against all 56k genes
            cur.execute(
                """
                SELECT COALESCE(gene_id, symbol) AS gene_key,
                       similarity(sequence, %(query)s) AS score
                FROM genes
                WHERE sequence IS NOT NULL
                  AND sequence <> ''
                  AND sequence %% %(query)s
                ORDER BY score DESC
                LIMIT %(limit)s;
                """,
                {"query": sequence, "limit": limit},
            )
            return [(row[0], float(row[1])) for row in cur.fetchall()]


def find_gene_keys_by_length_range(
    min_length: int, max_length: int, sequence_type: str | None = None, limit: int = 500
) -> list[str]:
    """Server-side length-range candidate lookup — fallback path for when
    gene_kmers hasn't been populated yet (see populate_kmer_index). Uses
    the genes_length_idx index so this is a fast range scan on Postgres's
    side instead of requiring the app to hold a full metadata dict in
    memory just to filter by length in Python.
    """
    where = "WHERE length BETWEEN %(min_len)s AND %(max_len)s"
    params: dict = {"min_len": min_length, "max_len": max_length, "limit": limit}
    if sequence_type:
        where += " AND sequence_type = %(seq_type)s"
        params["seq_type"] = sequence_type

    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                f"SELECT COALESCE(gene_id, symbol) FROM genes {where} LIMIT %(limit)s;",
                params,
            )
            return [row[0] for row in cur.fetchall() if row[0]]


def load_gene_database_metadata_from_postgres(batch_size: int = 1000) -> dict:
    """Load only metadata needed for UI search without reading full sequences."""
    records: dict[str, dict] = {}
    with get_connection() as conn:
        conn.autocommit = False
        try:
            with conn.cursor(name="gene_metadata_cursor") as cur:
                cur.itersize = batch_size
                cur.execute(
                    """
                    SELECT gene_id, symbol, organism, description, source, traits, length
                    FROM genes;
                    """
                )
                for row in cur:
                    gene_id, symbol, organism, description, source, traits, length = row
                    if isinstance(traits, str):
                        traits = json.loads(traits)
                    key = gene_id or symbol
                    if not key:
                        continue
                    records[key] = {
                        "gene_id": gene_id,
                        "symbol": symbol,
                        "organism": organism,
                        "description": description,
                        "source": source,
                        "traits": traits or [],
                        "length": length,
                    }
        finally:
            conn.commit()
    return records

def load_json_records(path: Path) -> list[dict]:
    import json

    raw = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(raw, dict) and "genes" in raw:
        return raw["genes"]
    if isinstance(raw, list):
        return raw
    if isinstance(raw, dict) and raw.get("gene_id"):
        return [raw]
    raise ValueError("Unsupported JSON format for gene records")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Postgres maintenance for the plant gene project. "
        "Run 'create-tables' once after provisioning a new database, and "
        "'populate-kmer-index' once after the initial data collection (and "
        "again after ingesting new batches of genes)."
    )
    parser.add_argument(
        "command",
        choices=["create-tables", "populate-kmer-index", "backfill-sequence-hashes"],
    )
    parser.add_argument("--k", type=int, default=KMER_K, help=f"k-mer length (default {KMER_K})")
    parser.add_argument(
        "--batch-size", type=int, default=2_000,
        help="Rows to read per cursor fetch before flushing (default 2000)",
    )
    parser.add_argument(
        "--flush-every", type=int, default=200_000,
        help="How many k-mer rows to accumulate before flushing to Postgres (default 200000)",
    )
    parser.add_argument(
        "--rebuild", action="store_true",
        help="Re-index every gene from scratch (needed after changing --k)",
    )
    args = parser.parse_args()

    if args.command == "create-tables":
        create_tables()
        print("Tables and indexes ensured.")
    elif args.command == "populate-kmer-index":
        n = populate_kmer_index(
            k=args.k,
            batch_size=args.batch_size,
            flush_every=args.flush_every,
            rebuild=args.rebuild,
        )
        print(f"Indexed {n} gene(s) with k={args.k}.")
    elif args.command == "backfill-sequence-hashes":
        n = backfill_sequence_hashes()
        print(f"Backfilled sequence_hash for {n} existing gene(s).")