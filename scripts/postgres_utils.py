#!/usr/bin/env python3
"""
PostgreSQL helper utilities for the plant gene project.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from urllib.parse import quote_plus

from dotenv import load_dotenv
import psycopg
from psycopg import sql

try:
    # Optional: pooled connections avoid paying a fresh TCP+TLS+auth
    # handshake (and, on Neon, a possible compute cold-start) on every
    # single query. Falls back to plain psycopg.connect() below if the
    # package isn't installed -- nothing else changes.
    from psycopg_pool import ConnectionPool
except ImportError:  # pragma: no cover
    ConnectionPool = None

ROOT = Path(__file__).resolve().parents[1]
load_dotenv(ROOT / ".env")


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


def get_connection() -> psycopg.Connection:
    return psycopg.connect(_resolve_database_url(), autocommit=True)


_POOL: "ConnectionPool | None" = None


def get_pool() -> "ConnectionPool":
    """Return a lazily-created, process-wide connection pool.

    Use this (via `pooled_connection()` below) for short, frequent reads
    such as the sidebar metadata search, where opening a brand-new
    connection per Streamlit rerun/keystroke is the dominant cost against
    a remote/serverless Postgres (e.g. Neon). Prefer Neon's pooled
    connection string (port 6543 / pgbouncer) as DATABASE_URL when using
    this, since it stacks with server-side pooling.
    """
    global _POOL
    if ConnectionPool is None:
        raise RuntimeError(
            "psycopg_pool is not installed. Run `pip install psycopg[pool]` "
            "to enable pooled connections, or use get_connection() instead."
        )
    if _POOL is None:
        _POOL = ConnectionPool(
            _resolve_database_url(),
            min_size=1,
            max_size=5,
            kwargs={"autocommit": True},
        )
    return _POOL


class _PooledOrPlainConnection:
    """Context manager returning a pooled connection when available,
    otherwise a regular one-off connection. Keeps call sites simple and
    keeps this module usable even without psycopg_pool installed."""

    def __enter__(self) -> psycopg.Connection:
        if ConnectionPool is not None:
            self._ctx = get_pool().connection()
            return self._ctx.__enter__()
        self._ctx = get_connection()
        return self._ctx.__enter__() if hasattr(self._ctx, "__enter__") else self._ctx

    def __exit__(self, *exc):
        return self._ctx.__exit__(*exc)


def pooled_connection() -> _PooledOrPlainConnection:
    return _PooledOrPlainConnection()


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
                    kmer_signature JSONB,
                    expression_profiles JSONB,
                    pathways JSONB,
                    publications JSONB,
                    annotations JSONB,
                    traits JSONB,
                    length INTEGER,
                    date_added TIMESTAMPTZ,
                    record JSONB
                );
                """
            )
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
            cur.execute("CREATE INDEX IF NOT EXISTS genes_organism_idx ON genes (organism);")
            cur.execute("CREATE INDEX IF NOT EXISTS genes_sequence_type_idx ON genes (sequence_type);")
            cur.execute("CREATE INDEX IF NOT EXISTS genes_source_idx ON genes (source);")
            # Index on kmer_signature (JSONB) can speed up candidate selection
            # using containment or expression indices in future queries.
            cur.execute("CREATE INDEX IF NOT EXISTS genes_kmer_signature_idx ON genes USING gin (kmer_signature);")

            # Trigram indexes let the sidebar's "search gene ID, symbol, or
            # trait" box run as a fast server-side ILIKE query instead of
            # pulling all ~56k metadata rows into Python and filtering with
            # a list comprehension on every keystroke/rerun (see
            # search_gene_metadata below).
            cur.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm;")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS genes_gene_id_trgm_idx "
                "ON genes USING gin (gene_id gin_trgm_ops);"
            )
            cur.execute(
                "CREATE INDEX IF NOT EXISTS genes_symbol_trgm_idx "
                "ON genes USING gin (symbol gin_trgm_ops);"
            )
            cur.execute(
                "CREATE INDEX IF NOT EXISTS genes_description_trgm_idx "
                "ON genes USING gin (description gin_trgm_ops);"
            )
            # GIN index on the traits JSONB array for containment/search.
            cur.execute("CREATE INDEX IF NOT EXISTS genes_traits_gin_idx ON genes USING gin (traits);")


_UPSERT_SQL = sql.SQL(
    """
    INSERT INTO genes (
        gene_id, symbol, organism, sequence, sequence_type,
        description, source, source_url, external_links, kmer_signature,
        expression_profiles, pathways, publications,
        annotations, traits, length, date_added, record
    ) VALUES (
        %(gene_id)s, %(symbol)s, %(organism)s, %(sequence)s, %(sequence_type)s,
        %(description)s, %(source)s, %(source_url)s, %(external_links)s, %(kmer_signature)s,
        %(expression_profiles)s, %(pathways)s, %(publications)s,
        %(annotations)s, %(traits)s, %(length)s, %(date_added)s, %(record)s
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
        source = EXCLUDED.source,
        source_url = COALESCE(NULLIF(EXCLUDED.source_url, ''), genes.source_url),
        -- JSONB objects: shallow-merge instead of replace, so new keys from
        -- a later source add to what's already known rather than erasing it.
        external_links = genes.external_links || EXCLUDED.external_links,
        annotations = genes.annotations || EXCLUDED.annotations,
        -- Keep existing kmer_signature when incoming one is null/empty.
        kmer_signature = COALESCE(EXCLUDED.kmer_signature, genes.kmer_signature),
        -- JSONB arrays: keep the existing array if the incoming one is empty.
        expression_profiles = CASE WHEN EXCLUDED.expression_profiles = '[]'::jsonb
            THEN genes.expression_profiles ELSE EXCLUDED.expression_profiles END,
        pathways = CASE WHEN EXCLUDED.pathways = '[]'::jsonb
            THEN genes.pathways ELSE EXCLUDED.pathways END,
        publications = CASE WHEN EXCLUDED.publications = '[]'::jsonb
            THEN genes.publications ELSE EXCLUDED.publications END,
        traits = CASE WHEN EXCLUDED.traits = '[]'::jsonb
            THEN genes.traits ELSE EXCLUDED.traits END,
        -- Don't let a 0/NULL length (no sequence fetched) overwrite a real
        -- known length from a previous insert.
        length = COALESCE(NULLIF(EXCLUDED.length, 0), genes.length),
        date_added = EXCLUDED.date_added,
        record = EXCLUDED.record;
    """
)


def _record_to_params(record: dict) -> dict:
    if not record.get("gene_id") and not record.get("symbol"):
        raise ValueError("Record must contain gene_id or symbol")
    return {
        "gene_id": record.get("gene_id"),
        "symbol": record.get("symbol"),
        "organism": record.get("organism"),
        "sequence": record.get("sequence"),
        "sequence_type": record.get("sequence_type"),
        "description": record.get("description"),
        "source": record.get("source"),
        "source_url": record.get("source_url"),
        "external_links": json.dumps(record.get("external_links", {})),
        "expression_profiles": json.dumps(record.get("expression_profiles", [])),
        "pathways": json.dumps(record.get("pathways", [])),
        "publications": json.dumps(record.get("publications", [])),
        "annotations": json.dumps(record.get("annotations", {})),
        "traits": json.dumps(record.get("traits", [])),
        "kmer_signature": json.dumps(record.get("kmer_signature", [])),
        "length": record.get("length") or (len(record.get("sequence", "")) if record.get("sequence") else None),
        "date_added": record.get("date_added"),
        "record": json.dumps(record),
    }


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
    """Load all gene records from Postgres, streaming in batches instead of
    a single fetchall(). With tens of thousands of rows -- each duplicating
    its data once in individual columns and again in the `record` JSONB
    blob -- an unbounded fetchall() can exhaust available memory (this bit
    us with "out of memory for query result" once the dataset passed ~49k
    rows). A named (server-side) cursor keeps only `batch_size` rows in
    memory at a time.
    """
    records: dict[str, dict] = {}
    with get_connection() as conn:
        conn.autocommit = False
        try:
            with conn.cursor(name="gene_export_cursor") as cur:
                cur.itersize = batch_size
                cur.execute(
                    """
                    SELECT gene_id, symbol, organism, sequence, sequence_type,
                           description, source, source_url, external_links,
                           expression_profiles, pathways, publications,
                           annotations, traits, length, date_added, record
                    FROM genes;
                    """
                )
                for row in cur:
                    (
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
                        record,
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
                    if isinstance(record, str):
                        record = json.loads(record)

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
                        "external_links": external_links or {},
                        "expression_profiles": expression_profiles or [],
                        "pathways": pathways or [],
                        "publications": publications or [],
                        "annotations": annotations or {},
                        "traits": traits or [],
                        "length": length,
                        "date_added": date_added,
                        "record": record,
                    }
        finally:
            conn.commit()  # read-only, just closes the transaction cleanly

    return records


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

def get_gene_count() -> int:
    """Cheap COUNT(*) for sidebar display -- doesn't touch row data at all."""
    with pooled_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM genes;")
            (count,) = cur.fetchone()
            return int(count)


def _metadata_row_to_dict(row: tuple) -> dict:
    gene_id, symbol, organism, description, source, traits, length = row
    if isinstance(traits, str):
        traits = json.loads(traits)
    return {
        "gene_id": gene_id,
        "symbol": symbol,
        "organism": organism,
        "description": description,
        "source": source,
        "traits": traits or [],
        "length": length,
    }


def search_gene_metadata(query: str | None = None, limit: int = 20, offset: int = 0) -> list[dict]:
    """Search/paginate gene metadata server-side.

    Replaces the pattern of loading all ~56k metadata rows into Python and
    filtering them with a list comprehension on every sidebar keystroke.
    Only `limit` rows ever cross the wire. `query` is matched against
    gene_id, symbol, description, source, and traits using the trigram/GIN
    indexes created in create_tables(); pass `query=None` (or "") for a
    plain paginated browse (e.g. the "sample of 20" preview).
    """
    base_select = """
        SELECT gene_id, symbol, organism, description, source, traits, length
        FROM genes
    """
    with pooled_connection() as conn:
        with conn.cursor() as cur:
            if query:
                like = f"%{query}%"
                cur.execute(
                    base_select
                    + """
                    WHERE gene_id ILIKE %(like)s
                       OR symbol ILIKE %(like)s
                       OR description ILIKE %(like)s
                       OR source ILIKE %(like)s
                       OR traits::text ILIKE %(like)s
                    ORDER BY gene_id NULLS LAST
                    LIMIT %(limit)s OFFSET %(offset)s;
                    """,
                    {"like": like, "limit": limit, "offset": offset},
                )
            else:
                cur.execute(
                    base_select + " ORDER BY gene_id NULLS LAST LIMIT %(limit)s OFFSET %(offset)s;",
                    {"limit": limit, "offset": offset},
                )
            return [_metadata_row_to_dict(row) for row in cur.fetchall()]


def count_gene_metadata_matches(query: str) -> int:
    """Count matches for `query` without fetching the matching rows themselves
    (used to show "Showing N matching records" without loading N of them)."""
    like = f"%{query}%"
    with pooled_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT COUNT(*) FROM genes
                WHERE gene_id ILIKE %(like)s
                   OR symbol ILIKE %(like)s
                   OR description ILIKE %(like)s
                   OR source ILIKE %(like)s
                   OR traits::text ILIKE %(like)s;
                """,
                {"like": like},
            )
            (count,) = cur.fetchone()
            return int(count)


def find_candidate_gene_ids_by_kmers(kmers: list[str], limit: int = 500) -> list[str]:
    """Return a list of gene_id (or symbol) candidates ordered by k-mer overlap.

    This executes a server-side aggregation that counts matching k-mers in
    `kmer_signature` JSONB arrays and returns the top `limit` gene identifiers.
    """
    if not kmers:
        return []
    with pooled_connection() as conn:
        with conn.cursor() as cur:
            # Unnest the JSONB kmer_signature array and count overlaps with the
            # provided kmers list. Use = ANY(%s) to pass the Python list safely.
            cur.execute(
                """
                SELECT COALESCE(gene_id, symbol) AS gid, COUNT(*) AS common
                FROM genes, jsonb_array_elements_text(kmer_signature) AS k(kmer)
                WHERE k.kmer = ANY(%s)
                GROUP BY gid
                ORDER BY common DESC
                LIMIT %s;
                """,
                (kmers, limit),
            )
            rows = cur.fetchall()
            return [row[0] for row in rows]


def find_candidate_gene_ids_by_kmers_weighted(kmers: list[str], limit: int = 500, min_jaccard: float = 0.01) -> list[str]:
    """Return gene ids ordered by Jaccard similarity between provided `kmers` and stored `kmer_signature`.

    Jaccard = |intersection| / (|query_kmers| + jsonb_array_length(kmer_signature) - |intersection|)
    The query computes the intersection count and filters by minimum Jaccard.
    """
    if not kmers:
        return []
    with pooled_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT gid FROM (
                  SELECT COALESCE(gene_id, symbol) AS gid,
                         COUNT(*) AS common,
                         jsonb_array_length(kmer_signature) AS ref_count
                  FROM genes, jsonb_array_elements_text(kmer_signature) AS k(kmer)
                  WHERE k.kmer = ANY(%s)
                  GROUP BY gid, ref_count
                ) t
                WHERE (common::float / (%s + ref_count - common)) >= %s
                ORDER BY (common::float / (%s + ref_count - common)) DESC
                LIMIT %s;
                """,
                (kmers, len(kmers), min_jaccard, len(kmers), limit),
            )
            rows = cur.fetchall()
            return [row[0] for row in rows]


def load_gene_records_by_ids(ids: list[str]) -> dict:
    """Load full gene records for the provided list of gene_id/symbol identifiers.

    Returns a dict keyed by gene_id/symbol like `load_gene_database_from_postgres`.
    """
    if not ids:
        return {}
    placeholders = ','.join(['%s'] * len(ids))
    query = (
        f"SELECT gene_id, symbol, organism, sequence, sequence_type, description, source, source_url, "
        f"external_links, expression_profiles, pathways, publications, annotations, traits, length, date_added, record "
        f"FROM genes WHERE COALESCE(gene_id, symbol) IN ({placeholders});"
    )
    with pooled_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(query, ids)
            records: dict[str, dict] = {}
            for row in cur.fetchall():
                (
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
                    record,
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
                if isinstance(record, str):
                    record = json.loads(record)

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
                    "external_links": external_links or {},
                    "expression_profiles": expression_profiles or [],
                    "pathways": pathways or [],
                    "publications": publications or [],
                    "annotations": annotations or {},
                    "traits": traits or [],
                    "length": length,
                    "date_added": date_added,
                    "record": record,
                }
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