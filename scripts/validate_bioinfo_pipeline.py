#!/usr/bin/env python3
"""S6: reproducible quantitative QC for bioinfo_pipeline.

The validator deliberately reports missing evidence separately from zero
coverage. This prevents a table containing only successful annotation hits
from being mistaken for a completed scan.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import psycopg2
from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core_engines.bioinfo_db import connect as connect_bioinfo_pipeline
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
EMBED_DIM = 1280


def connect():
    return connect_bioinfo_pipeline()
    load_dotenv(ROOT / ".env")
    if os.getenv("BIOINFO_PIPELINE_URL"):
        return psycopg2.connect(os.environ["BIOINFO_PIPELINE_URL"])
    return psycopg2.connect(
        host=os.getenv("BIOINFO_PIPELINE_HOST", "localhost"),
        port=int(os.getenv("BIOINFO_PIPELINE_PORT", "5432")),
        user=os.getenv("BIOINFO_PIPELINE_USER", "postgres"),
        password=os.getenv("BIOINFO_PIPELINE_PASSWORD"),
        dbname=os.getenv("BIOINFO_PIPELINE_DBNAME", "bioinfo_pipeline"),
    )


def _table_exists(cursor, table: str) -> bool:
    cursor.execute("SELECT to_regclass(%s) IS NOT NULL", (table,))
    return bool(cursor.fetchone()[0])


def _columns(cursor, table: str) -> set[str]:
    cursor.execute(
        """
        SELECT column_name FROM information_schema.columns
        WHERE table_name = %s AND table_schema = ANY(current_schemas(false))
        """,
        (table,),
    )
    return {row[0] for row in cursor.fetchall()}


def _decode_vector(value) -> np.ndarray:
    if isinstance(value, (bytes, bytearray, memoryview)):
        return np.frombuffer(value, dtype=np.float32)
    return np.asarray(value, dtype=np.float32)


def _fasta_stats(path: Path) -> dict:
    count = 0
    nonempty = 0
    invalid = 0
    lengths: list[int] = []
    current: list[str] = []

    def finish() -> None:
        nonlocal nonempty, invalid
        sequence = "".join(current).replace(" ", "").replace("\r", "").upper()
        if sequence:
            nonempty += 1
            lengths.append(len(sequence))
            invalid += len(set(sequence) - VALID_AA)

    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    finish()
                count += 1
                current = []
            else:
                current.append(line.strip())
        if current:
            finish()
    return {
        "records": count,
        "nonempty": nonempty,
        "invalid_residues": invalid,
        "max_length": max(lengths, default=0),
        "over_1022": sum(length > 1022 for length in lengths),
    }


def _declared_fasta_paths(source_metadata: object) -> list[Path]:
    if not isinstance(source_metadata, dict):
        return []
    files = source_metadata.get("files", {})
    values = files.values() if isinstance(files, dict) else files if isinstance(files, list) else []
    paths = []
    for value in values:
        if isinstance(value, dict):
            value = value.get("path")
        if isinstance(value, str) and value.lower().endswith((".faa", ".fasta", ".fa", ".faa.gz")):
            path = Path(value)
            if path.exists():
                paths.append(path)
    return paths


def _declared_file_record(source_metadata: object, path: Path) -> dict:
    files = source_metadata.get("files", {}) if isinstance(source_metadata, dict) else {}
    if not isinstance(files, dict):
        return {}
    for record in files.values():
        if isinstance(record, dict) and Path(record.get("path", "")).name == path.name:
            return record
    return {}


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _percent(n: int, total: int) -> float:
    return round(100 * n / total, 2) if total else 0.0


def validate(species: str | None = None) -> dict:
    report = {"species": {}, "global": {}, "ready_for_publication": True}
    with connect() as connection, connection.cursor() as cursor:
        if not _table_exists(cursor, "species_registry"):
            raise RuntimeError("species_registry is missing; run S0 first.")
        registry_columns = _columns(cursor, "species_registry")
        selected = """
            SELECT species, assemblage, assembly_accession, status,
                   modules, source_metadata, source_protein_count
            FROM species_registry
            WHERE (%s IS NULL OR species = %s)
            ORDER BY species
        """
        cursor.execute(selected, (species, species))
        registries = cursor.fetchall()
        if not registries:
            raise RuntimeError("No matching species registry row.")

        required_tables = (
            "dim_protein", "dim_sequence", "esm2_embeddings",
            "fact_functional_annotation",
        )
        for table in required_tables:
            report["global"][f"{table}_exists"] = _table_exists(cursor, table)
        report["global"]["annotation_scan_evidence_exists"] = any(
            _table_exists(cursor, name)
            for name in ("annotation_scan_evidence", "annotation_runs", "pipeline_runs")
        )
        report["ready_for_publication"] = all(
            report["global"][f"{table}_exists"] for table in required_tables
        )

        for (
            species_key, assemblage, assembly_accession, status,
            modules, source_metadata, source_protein_count,
        ) in registries:
            modules = modules or {}
            item = {"status": status, "checks": {}, "ready": True}
            checks = item["checks"]
            protein_columns = _columns(cursor, "dim_protein")
            seq_columns = _columns(cursor, "dim_sequence")
            embedding_columns = _columns(cursor, "esm2_embeddings")

            cursor.execute(
                "SELECT COUNT(*) FROM dim_protein WHERE assemblage = %s",
                (assemblage,),
            )
            protein_count = cursor.fetchone()[0]
            checks["dim_protein_count"] = protein_count
            gene_id_columns = [
                column for column in ("gene_id", "plaza_gene_id")
                if column in protein_columns
            ]
            if gene_id_columns:
                cursor.execute(
                    f"""SELECT COUNT(*) FROM dim_protein
                       WHERE assemblage = %s AND COALESCE(
                       {", ".join(gene_id_columns)}) IS NOT NULL""",
                    (assemblage,),
                )
                checks["dim_protein_gene_id_pct"] = _percent(cursor.fetchone()[0], protein_count)
            if {"chr", "start_pos", "stop_pos"} <= protein_columns:
                cursor.execute(
                    """SELECT COUNT(*) FROM dim_protein
                       WHERE assemblage = %s AND chr IS NOT NULL
                       AND start_pos IS NOT NULL AND stop_pos IS NOT NULL""",
                    (assemblage,),
                )
                checks["dim_protein_coordinate_pct"] = _percent(cursor.fetchone()[0], protein_count)

            fasta_stats = {}
            for path in _declared_fasta_paths(source_metadata):
                stats = _fasta_stats(path)
                expected_checksums = (
                    source_metadata.get("checksums", {})
                    if isinstance(source_metadata, dict)
                    else {}
                )
                expected = expected_checksums.get(str(path), expected_checksums.get(path.name))
                expected = expected or _declared_file_record(source_metadata, path).get("md5")
                if expected:
                    stats["md5"] = _md5(path)
                    stats["md5_matches_manifest"] = stats["md5"].lower() == str(expected).lower()
                fasta_stats[str(path)] = stats
            checks["source_fasta"] = fasta_stats
            source_count = source_protein_count or sum(
                stats["records"] for stats in fasta_stats.values()
            )
            checks["source_protein_count"] = source_count or None
            checks["source_count_matches_dim_protein"] = (
                source_count == protein_count if source_count else False
            )
            checks["source_checksum_failures"] = sum(
                1 for stats in fasta_stats.values()
                if stats.get("md5_matches_manifest") is False
            )

            sequence_count = 0
            sequence_md5_by_accession = {}
            over_1022 = 0
            checks["sequence_nonempty_pct"] = 0.0
            checks["sequence_alphabet_invalid_residues"] = 0
            checks["sequence_length_mismatch_count"] = 0
            checks["cds_matches_protein_status_counts"] = {}
            if report["global"]["dim_sequence_exists"]:
                sequence_column = next(
                    (name for name in ("protein_sequence", "sequence") if name in seq_columns),
                    None,
                )
                if sequence_column:
                    cursor.execute(
                        f"""SELECT COUNT(*) FROM dim_sequence
                            WHERE assembly_accession = %s
                            AND {sequence_column} IS NOT NULL
                            AND {sequence_column} <> ''""",
                        (assembly_accession,),
                    )
                    sequence_count = cursor.fetchone()[0]
                    checks["sequence_nonempty_pct"] = _percent(sequence_count, protein_count)
                    checks["sequence_alphabet_invalid_residues"] = 0
                    checks["sequence_length_mismatch_count"] = 0
                    sequence_select = ["protein_accession", sequence_column]
                    if "sequence_md5" in seq_columns:
                        sequence_select.append("sequence_md5")
                    if "cds_matches_protein" in seq_columns:
                        sequence_select.append("cds_matches_protein")
                    md5_index = (
                        sequence_select.index("sequence_md5")
                        if "sequence_md5" in sequence_select else None
                    )
                    status_index = (
                        sequence_select.index("cds_matches_protein")
                        if "cds_matches_protein" in sequence_select else None
                    )
                    cursor.execute(
                        f"""SELECT {", ".join(sequence_select)} FROM dim_sequence
                            WHERE assembly_accession = %s""",
                        (assembly_accession,),
                    )
                    for row in cursor.fetchall():
                        accession, sequence = row[:2]
                        value = str(sequence or "").upper()
                        checks["sequence_alphabet_invalid_residues"] += len(
                            set(value) - VALID_AA
                        )
                        over_1022 += len(value) > 1022
                        sequence_md5_by_accession[accession] = (
                            row[md5_index] if md5_index is not None
                            and row[md5_index]
                            else hashlib.md5(value.encode("ascii", "ignore")).hexdigest()
                        )
                        if status_index is not None and row[status_index]:
                            status = row[status_index]
                            checks["cds_matches_protein_status_counts"][status] = (
                                checks["cds_matches_protein_status_counts"].get(status, 0)
                                + 1
                            )
                else:
                    checks["sequence_schema_error"] = "No protein_sequence/sequence column."
            checks["dim_sequence_protein_count"] = sequence_count

            embedding_count = 0
            embedding_bad_dimension = 0
            embedding_nan = 0
            embedding_zero = 0
            md5_mismatch = 0
            md5_missing = 0
            embedding_accessions = set()
            embedding_md5_status_counts = {}
            if report["global"]["esm2_embeddings_exists"]:
                cursor.execute(
                    """SELECT COUNT(*) FROM esm2_embeddings
                       WHERE assembly_accession = %s""",
                    (assembly_accession,),
                )
                embedding_count = cursor.fetchone()[0]
                select_columns = ["protein_accession", "embedding"]
                if "sequence_md5" in embedding_columns:
                    select_columns.append("sequence_md5")
                if "md5_verification_status" in embedding_columns:
                    select_columns.append("md5_verification_status")
                with connection.cursor(name="embedding_qc_cursor") as embedding_cursor:
                    embedding_cursor.itersize = 2048
                    embedding_cursor.execute(
                        f"""SELECT {", ".join(select_columns)} FROM esm2_embeddings
                            WHERE assembly_accession = %s""",
                        (assembly_accession,),
                    )
                    for row in embedding_cursor:
                        embedding_accessions.add(row[0])
                        vector = _decode_vector(row[1])
                        embedding_bad_dimension += vector.shape != (EMBED_DIM,)
                        embedding_nan += int(np.isnan(vector).any())
                        embedding_zero += int(
                            vector.shape == (EMBED_DIM,) and np.linalg.norm(vector) == 0
                        )
                        md5_index = (
                            select_columns.index("sequence_md5")
                            if "sequence_md5" in select_columns
                            else None
                        )
                        if md5_index is not None and row[md5_index]:
                            md5_mismatch += int(
                                sequence_md5_by_accession.get(row[0]) != row[md5_index]
                            )
                        else:
                            md5_missing += 1
                        status_index = (
                            select_columns.index("md5_verification_status")
                            if "md5_verification_status" in select_columns
                            else None
                        )
                        if status_index is not None:
                            status = row[status_index] or "unset"
                            embedding_md5_status_counts[status] = (
                                embedding_md5_status_counts.get(status, 0) + 1
                            )
            cursor.execute(
                """SELECT COUNT(*) FROM dim_protein
                   WHERE assemblage = %s AND protein_accession = ANY(%s)""",
                (assemblage, list(embedding_accessions)),
            )
            embedding_in_dim = cursor.fetchone()[0] if embedding_accessions else 0
            checks.update({
                "embedding_count": embedding_count,
                "embedding_protein_coverage_pct": _percent(embedding_in_dim, protein_count),
                "embedding_orphan_count": embedding_count - embedding_in_dim,
                "embedding_bad_dimension": embedding_bad_dimension,
                "embedding_nan_vectors": embedding_nan,
                "embedding_zero_vectors": embedding_zero,
                "embedding_sequence_md5_mismatches": md5_mismatch,
                "embedding_sequence_md5_missing": md5_missing,
                "embedding_sequence_md5_check": (
                    "sequence_md5" in embedding_columns
                    and "sequence_md5" in seq_columns
                    and md5_missing == 0
                    and md5_mismatch == 0
                ),
                "embedding_md5_verification_status_counts": embedding_md5_status_counts,
                "proteins_over_1022aa": over_1022,
            })

            source_names = {
                "interpro": ("interpro",),
                "go": ("go", "go_plaza", "go_swissprot"),
                "swissprot": ("swissprot", "go_swissprot"),
            }
            annotation_sources = {}
            if _table_exists(cursor, "annotation_pfam_domain"):
                cursor.execute(
                    """SELECT COUNT(DISTINCT protein_accession)
                       FROM annotation_pfam_domain
                       WHERE assembly_accession = %s""",
                    (assembly_accession,),
                )
                annotation_sources["pfam"] = _percent(
                    cursor.fetchone()[0], protein_count
                )
            else:
                cursor.execute(
                    """SELECT COUNT(DISTINCT protein_accession)
                       FROM fact_functional_annotation
                       WHERE assemblage = %s AND source = 'pfam'""",
                    (assemblage,),
                )
                annotation_sources["pfam"] = _percent(
                    cursor.fetchone()[0], protein_count
                )
            for source_name, aliases in source_names.items():
                cursor.execute(
                    """SELECT COUNT(DISTINCT protein_accession)
                       FROM fact_functional_annotation
                       WHERE assemblage = %s AND source = ANY(%s)""",
                    (assemblage, list(aliases)),
                )
                annotation_sources[source_name] = _percent(cursor.fetchone()[0], protein_count)
            checks["annotation_coverage_pct_by_source"] = annotation_sources
            checks["pfam_below_sanity_threshold"] = (
                "pfam" in annotation_sources and annotation_sources["pfam"] < 60
            )
            checks["pfam_above_sanity_threshold"] = (
                "pfam" in annotation_sources and annotation_sources["pfam"] > 90
            )
            if _table_exists(cursor, "annotation_runs"):
                cursor.execute(
                    """
                    SELECT tool_name
                    FROM annotation_runs
                    WHERE assembly_accession=%s
                      AND run_completed IS TRUE
                      AND (
                          (tool_name='pfam_hmmscan' AND database_release='38.2')
                          OR
                          (tool_name='swissprot_diamond' AND database_release='2026_03')
                      )
                    """,
                    (assembly_accession,),
                )
                evidence_tools = {row[0] for row in cursor.fetchall()}
                required_tools = {"pfam_hmmscan", "swissprot_diamond"}
                checks["annotation_scan_evidence_tools"] = sorted(evidence_tools)
                checks["annotation_scan_evidence_missing_tools"] = sorted(
                    required_tools - evidence_tools
                )
                checks["annotation_scan_evidence"] = not (
                    required_tools - evidence_tools
                )
            elif _table_exists(cursor, "annotation_scan_evidence"):
                cursor.execute(
                    """
                    SELECT tool_name
                    FROM annotation_scan_evidence
                    WHERE assembly_accession=%s
                      AND evidence_status='completed'
                    """,
                    (assembly_accession,),
                )
                evidence_tools = {row[0] for row in cursor.fetchall()}
                required_tools = {"pfam_hmmscan", "swissprot_diamond"}
                checks["annotation_scan_evidence_tools"] = sorted(evidence_tools)
                checks["annotation_scan_evidence_missing_tools"] = sorted(
                    required_tools - evidence_tools
                )
                checks["annotation_scan_evidence"] = not (
                    required_tools - evidence_tools
                )
            else:
                checks["annotation_scan_evidence_tools"] = []
                checks["annotation_scan_evidence_missing_tools"] = [
                    "pfam_hmmscan",
                    "swissprot_diamond",
                ]
                checks["annotation_scan_evidence"] = False

            orphan_checks = {}
            cursor.execute(
                """SELECT COUNT(*) FROM esm2_embeddings e
                   LEFT JOIN dim_protein p ON p.assemblage = %s
                   AND p.protein_accession = e.protein_accession
                   WHERE e.assembly_accession = %s AND p.protein_accession IS NULL""",
                (assemblage, assembly_accession),
            )
            orphan_checks["embeddings"] = cursor.fetchone()[0]
            if _table_exists(cursor, "fact_functional_annotation"):
                cursor.execute(
                    """SELECT COUNT(*) FROM fact_functional_annotation a
                       LEFT JOIN dim_protein p ON p.assemblage = a.assemblage
                       AND p.protein_accession = a.protein_accession
                       WHERE a.assemblage = %s AND p.protein_accession IS NULL""",
                    (assemblage,),
                )
                orphan_checks["annotations"] = cursor.fetchone()[0]
            checks["orphan_rows"] = orphan_checks
            checks["plaza_coordinate_concordance_pct"] = None
            checks["plaza_coordinate_discordant_links"] = None
            checks["plaza_coordinate_concordance_counts"] = None
            if (
                _table_exists(cursor, "dim_gene_protein")
                and "coordinate_concordance" in _columns(cursor, "dim_gene_protein")
            ):
                cursor.execute(
                    """
                    SELECT
                        COUNT(*) FILTER (WHERE gene_source='plaza'
                                         AND coordinate_concordance='concordant'),
                        COUNT(*) FILTER (WHERE gene_source='plaza'),
                        COUNT(*) FILTER (WHERE gene_source='plaza'
                                         AND coordinate_concordance='discordant'),
                        COUNT(*) FILTER (WHERE gene_source='plaza'
                                         AND coordinate_concordance='non_verifiable')
                    FROM dim_gene_protein
                    WHERE assembly_name=%s
                    """,
                    (assemblage,),
                )
                concordant, total_links, discordant, non_verifiable = cursor.fetchone()
                checks["plaza_coordinate_concordance_pct"] = _percent(
                    concordant, total_links
                )
                checks["plaza_coordinate_discordant_links"] = discordant
                checks["plaza_coordinate_concordance_counts"] = {
                    "concordant": concordant,
                    "discordant": discordant,
                    "non_verifiable": non_verifiable,
                    "total": total_links,
                }
                checks["plaza_coordinate_concordance_alert"] = (
                    total_links > 0 and checks["plaza_coordinate_concordance_pct"] < 95
                )

            module_failures = []
            if modules.get("gene_id") and checks.get("dim_protein_gene_id_pct", 0) < 100:
                module_failures.append("gene_id")
            if modules.get("coordinates") and checks.get("dim_protein_coordinate_pct", 0) < 100:
                module_failures.append("coordinates")
            if modules.get("orthology") and not _table_exists(cursor, "fact_orthology"):
                module_failures.append("orthology")
            checks["declared_module_failures"] = module_failures
            failures = []
            if not checks["source_count_matches_dim_protein"]:
                failures.append("source count does not match dim_protein")
            if checks["source_checksum_failures"]:
                failures.append("source checksum mismatch")
            if checks["sequence_nonempty_pct"] != 100:
                failures.append("sequences are not 100% non-empty")
            if checks["embedding_protein_coverage_pct"] != 100:
                failures.append("embedding coverage is below 100%")
            if checks["embedding_orphan_count"]:
                failures.append("embedding orphans exist")
            if not checks["embedding_sequence_md5_check"]:
                failures.append("embedding sequence MD5 verification is incomplete")
            if checks["pfam_below_sanity_threshold"]:
                failures.append("Pfam coverage is below 60%")
            if not checks["annotation_scan_evidence"]:
                failures.append("annotation scan evidence is missing")
            failures.extend(module_failures)
            checks["failed_criteria"] = failures
            item["ready"] = (
                protein_count > 0
                and checks["source_count_matches_dim_protein"]
                and checks["source_checksum_failures"] == 0
                and checks["embedding_sequence_md5_check"]
                and checks["sequence_nonempty_pct"] == 100
                and checks["embedding_protein_coverage_pct"] == 100
                and embedding_bad_dimension == embedding_nan == embedding_zero == 0
                and checks["annotation_scan_evidence"]
                and not module_failures
                and all(value == 0 for value in orphan_checks.values())
            )
            report["species"][species_key] = item
            report["ready_for_publication"] &= item["ready"]
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="S6: validate registered bioinfo_pipeline assemblies")
    parser.add_argument("--species")
    parser.add_argument("--out")
    args = parser.parse_args()
    report = validate(args.species)
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    print(rendered)
    if args.out:
        Path(args.out).write_text(rendered + "\n", encoding="utf-8")
    if not report["ready_for_publication"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
