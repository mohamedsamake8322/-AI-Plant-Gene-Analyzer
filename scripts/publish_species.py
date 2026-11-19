#!/usr/bin/env python3
"""S7: publish an assembly only after the complete S6 report is green."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validate_bioinfo_pipeline import connect, validate


def _progress(message: str) -> None:
    print(f"[S7] {message}", file=sys.stderr, flush=True)


def publish(species: str) -> dict:
    _progress(f"Validating {species}...")
    report = validate(species)
    _progress("Validation complete.")
    item = report["species"].get(species)
    if not report["ready_for_publication"] or not item or not item["ready"]:
        raise RuntimeError(
            "S7 stopped: S6 validation failed. Inspect the JSON report before retrying."
        )

    _progress("Updating registry and refreshing v_gene_text; this may take a while...")
    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            UPDATE species_registry
            SET enabled = TRUE, status = 'published', completed_phase = 7,
                updated_at = now()
            WHERE species = %s
            """,
            (species,),
        )
        if cursor.rowcount != 1:
            raise RuntimeError(f"No unique registry row found for {species}.")
        cursor.execute("SET LOCAL work_mem = '64MB'")
        cursor.execute("SET LOCAL max_parallel_workers_per_gather = 0")
        cursor.execute("SET LOCAL enable_sort = off")
        cursor.execute("REFRESH MATERIALIZED VIEW v_gene_text")
        _progress("v_gene_text refresh complete; writing publication audit...")
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS species_publish_audit (
                id BIGSERIAL PRIMARY KEY,
                species TEXT NOT NULL,
                published_at TIMESTAMPTZ NOT NULL DEFAULT now(),
                validation_report JSONB NOT NULL
            )
            """
        )
        cursor.execute(
            """
            INSERT INTO species_publish_audit (species, validation_report)
            VALUES (%s, %s::jsonb)
            """,
            (species, json.dumps(report, ensure_ascii=False)),
        )

    _progress("Database transaction committed; invalidating ESM-2 embedding cache...")
    from core_engines.esm2_similarity import reload_embedding_index
    reload_embedding_index()
    _progress("Publication complete.")
    return {
        "species": species,
        "status": "published",
        "published_at": datetime.now(timezone.utc).isoformat(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="S7: publish one validated species")
    parser.add_argument("--species", required=True)
    args = parser.parse_args()
    print(json.dumps(publish(args.species), indent=2))


if __name__ == "__main__":
    main()
