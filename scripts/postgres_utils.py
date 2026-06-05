#!/usr/bin/env python3
"""
PostgreSQL helper utilities for the plant gene project.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from urllib.parse import quote_plus

from dotenv import load_dotenv
import psycopg
from psycopg import sql

ROOT = Path(__file__).resolve().parents[1]
load_dotenv(ROOT / ".env")

DATABASE_URL = os.getenv("DATABASE_URL")
DB_HOST = os.getenv("DB_HOST")
DB_PORT = os.getenv("DB_PORT")
DB_NAME = os.getenv("DB_NAME")
DB_USER = os.getenv("DB_USER")
DB_PASSWORD = os.getenv("DB_PASSWORD")

if DB_HOST and DB_NAME and DB_USER and DB_PASSWORD:
    DATABASE_URL = (
        f"postgresql://{quote_plus(DB_USER)}:{quote_plus(DB_PASSWORD)}@{DB_HOST}"
        f":{DB_PORT or 5432}/{quote_plus(DB_NAME)}"
    )
elif not DATABASE_URL:
    raise RuntimeError(
        "DATABASE_URL not set and DB_HOST/DB_NAME/DB_USER/DB_PASSWORD are not all configured"
    )


def get_connection() -> psycopg.Connection:
    return psycopg.connect(DATABASE_URL, autocommit=True)


def create_tables() -> None:
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS genes (
                    id SERIAL PRIMARY KEY,
                    gene_id TEXT UNIQUE,
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


def insert_gene_record(record: dict) -> None:
    if not record.get("gene_id") and not record.get("symbol"):
        raise ValueError("Record must contain gene_id or symbol")

    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL(
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
                    ON CONFLICT (gene_id)
                    DO UPDATE SET
                        symbol = EXCLUDED.symbol,
                        organism = EXCLUDED.organism,
                        sequence = EXCLUDED.sequence,
                        sequence_type = EXCLUDED.sequence_type,
                        description = EXCLUDED.description,
                        source = EXCLUDED.source,
                        source_url = EXCLUDED.source_url,
                        external_links = EXCLUDED.external_links,
                        expression_profiles = EXCLUDED.expression_profiles,
                        pathways = EXCLUDED.pathways,
                        publications = EXCLUDED.publications,
                        annotations = EXCLUDED.annotations,
                        traits = EXCLUDED.traits,
                        length = EXCLUDED.length,
                        date_added = EXCLUDED.date_added,
                        record = EXCLUDED.record;
                    """
                ),
                {
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
                },
            )


def load_gene_database_from_postgres() -> dict:
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT gene_id, symbol, organism, sequence, sequence_type,
                       description, source, source_url, external_links,
                       expression_profiles, pathways, publications,
                       annotations, traits, length, date_added, record
                FROM genes;
                """
            )
            rows = cur.fetchall()

    records: dict[str, dict] = {}
    for row in rows:
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
