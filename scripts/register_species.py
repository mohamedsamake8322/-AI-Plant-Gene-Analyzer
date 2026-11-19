#!/usr/bin/env python3
"""S0: register assemblies from the manifest in the bioinfo_pipeline DB.

The manifest is the only species-specific input. Registration is idempotent:
re-running it updates the same registry row and never inserts duplicates.
Publication remains a separate S7 operation after QC has passed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core_engines.bioinfo_db import connect as connect_bioinfo_pipeline

VALID_ACCESSION = re.compile(r"^GC[AF]_\d{9}\.\d+$")


def _connection():
    return connect_bioinfo_pipeline()
    load_dotenv(ROOT / ".env")
    url = __import__("os").environ.get("BIOINFO_PIPELINE_URL")
    if url:
        return psycopg2.connect(url)
    import os

    return psycopg2.connect(
        host=os.getenv("BIOINFO_PIPELINE_HOST", "localhost"),
        port=int(os.getenv("BIOINFO_PIPELINE_PORT", "5432")),
        user=os.getenv("BIOINFO_PIPELINE_USER", "postgres"),
        password=os.getenv("BIOINFO_PIPELINE_PASSWORD"),
        dbname=os.getenv("BIOINFO_PIPELINE_DBNAME", "bioinfo_pipeline"),
    )


def _manifest_hash(entry: dict) -> str:
    payload = json.dumps(entry, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _assemblies(species_key: str, entry: dict) -> list[tuple[str, str, str, str]]:
    assemblies = entry.get("assemblies")
    if not assemblies:
        assemblies = [{
            "name": entry.get("reference_assembly"),
            "accession": entry.get("refseq_accession"),
        }]
    rows = []
    for assembly in assemblies:
        name = str(assembly.get("name") or entry.get("reference_assembly") or "").strip()
        configured_accession = (
            assembly.get("accession")
            or assembly.get("assembly_accession")
            or (entry.get("refseq_accession") if len(assemblies) == 1 else None)
        )
        accession = str(
            configured_accession
            or ""
        ).strip()
        if not name or not VALID_ACCESSION.match(accession):
            raise ValueError(
                f"{species_key}: each assembly needs a valid GCA/GCF accession and name"
            )
        default_key = (
            species_key
            if len(assemblies) == 1
            else f"{species_key}_{name.lower().replace('-', '_').replace('.', '_')}"
        )
        key = str(assembly.get("species") or default_key).strip()
        assemblage = str(assembly.get("assemblage") or key).strip()
        rows.append((key, name, accession, assemblage))
    return rows


def register(manifest_path: Path) -> int:
    document = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    species_entries = document.get("species") or {}
    if not isinstance(species_entries, dict) or not species_entries:
        raise ValueError("Manifest must contain a non-empty species mapping.")

    rows = []
    for species_key, entry in species_entries.items():
        if not isinstance(entry, dict):
            raise ValueError(f"{species_key}: manifest entry must be a mapping")
        for registry_key, assembly_name, accession, assemblage in _assemblies(species_key, entry):
            assembly_entry = next(
                (
                    candidate
                    for candidate in entry.get("assemblies", [])
                    if candidate.get("name") == assembly_name
                ),
                {},
            )
            source_entry = {**entry, **assembly_entry}
            download_manifest = (
                ROOT / "pipeline_data" / "ncbi" / registry_key / "download_manifest.json"
            )
            downloaded = {}
            if download_manifest.exists():
                downloaded = json.loads(download_manifest.read_text(encoding="utf-8"))
            rows.append((
                registry_key,
                f"{entry.get('common_name', species_key)} ({assembly_name})",
                assemblage,
                accession,
                _manifest_hash({
                    "species": species_key,
                    "registry_key": registry_key,
                    "entry": entry,
                    "assembly": assembly_name,
                    "assemblage": assemblage,
                    "accession": accession,
                }),
                json.dumps({
                    "roles": entry.get("roles", []),
                    "modules": entry.get("modules", {}),
                    "source_species": species_key,
                }),
                json.dumps({
                    "source": source_entry.get("source", "NCBI RefSeq"),
                    "download_date": downloaded.get(
                        "download_date", source_entry.get("download_date")
                    ),
                    "files": entry.get(
                        "source_files",
                        downloaded.get("files", source_entry.get("files", {})),
                    ),
                    "checksums": source_entry.get("checksums", {}),
                }),
                json.dumps(entry.get("tool_versions", {})),
                entry.get("source_protein_count"),
            ))

    with _connection() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS species_registry (
                species TEXT PRIMARY KEY,
                display_name TEXT NOT NULL,
                assemblage TEXT NOT NULL,
                assembly_accession TEXT NOT NULL UNIQUE,
                enabled BOOLEAN NOT NULL DEFAULT FALSE,
                status TEXT NOT NULL DEFAULT 'registered',
                completed_phase SMALLINT NOT NULL DEFAULT 0,
                manifest_hash TEXT,
                modules JSONB NOT NULL DEFAULT '{}'::jsonb,
                source_metadata JSONB NOT NULL DEFAULT '{}'::jsonb,
                tool_versions JSONB NOT NULL DEFAULT '{}'::jsonb,
                source_protein_count BIGINT,
                updated_at TIMESTAMPTZ NOT NULL DEFAULT now()
            )
            """
        )
        cursor.execute(
            """
            ALTER TABLE species_registry
                ADD COLUMN IF NOT EXISTS source_metadata JSONB NOT NULL DEFAULT '{}'::jsonb,
                ADD COLUMN IF NOT EXISTS tool_versions JSONB NOT NULL DEFAULT '{}'::jsonb,
                ADD COLUMN IF NOT EXISTS source_protein_count BIGINT
            """
        )
        for row in rows:
            cursor.execute(
                """
                UPDATE species_registry
                SET species = %s
                WHERE assembly_accession = %s AND species <> %s
                """,
                (row[0], row[3], row[0]),
            )
        cursor.executemany(
            """
            INSERT INTO species_registry
                (species, display_name, assemblage, assembly_accession,
                enabled, status, completed_phase, manifest_hash, modules,
                source_metadata, tool_versions, source_protein_count, updated_at)
            VALUES (%s, %s, %s, %s, FALSE, 'registered', 0, %s, %s::jsonb,
                   %s::jsonb, %s::jsonb, %s, now())
            ON CONFLICT (species) DO UPDATE SET
                display_name = EXCLUDED.display_name,
                assemblage = EXCLUDED.assemblage,
                assembly_accession = EXCLUDED.assembly_accession,
                manifest_hash = EXCLUDED.manifest_hash,
                modules = EXCLUDED.modules,
                source_metadata = EXCLUDED.source_metadata,
                tool_versions = EXCLUDED.tool_versions,
                source_protein_count = EXCLUDED.source_protein_count,
                updated_at = now()
            """,
            rows,
        )
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="S0: register species assemblies from species.yaml")
    parser.add_argument("--manifest", default=str(ROOT / "species.yaml"))
    args = parser.parse_args()
    count = register(Path(args.manifest))
    print(f"Registered {count} assembly record(s); none were published.")


if __name__ == "__main__":
    main()
