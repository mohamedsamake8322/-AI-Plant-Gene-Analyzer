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


_UPSERT_SQL = sql.SQL(
    """
    INSERT INTO genes (
        gene_id, symbol, organism, sequence, sequence_type,
        description, source, source_url, external_links,
        expression_profiles, pathways, publications,
        annotations, traits, length, date_added, record
    ) VALUES (
        %(gene_id)s, %(symbol)s, %(organism)s, %(sequence)s, %(sequence_type)s,
        %(description)s, %(source)s, %(source_url)s, %(external_links)s,
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
        # Named (server-side) cursors need a real transaction -- autocommit=True
        # (the default from get_connection()) would close it prematurely.
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
