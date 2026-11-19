#!/usr/bin/env python3
"""Validate and load frozen Pfam/DIAMOND evidence into bioinfo_pipeline.

The command is a dry run unless ``--apply`` is supplied.  It only reads
``annotation/<species>/pfam_hmmscan/38.2`` and
``annotation/<species>/swissprot_diamond/2026_03``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

try:
    from psycopg2.extras import execute_values
except ImportError:  # pragma: no cover - reported when the command is run
    execute_values = None

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SPECIES = {
    "oryza_sativa": ("GCF_001433935.1", 42580),
    "zea_mays": ("GCF_902167145.1", 57345),
    "nicotiana_tabacum_k326": ("GCF_000715075.1", 91379),
    "nicotiana_tabacum_tn90": ("GCF_000715135.1", 84350),
}
PFAM_RELEASE = "38.2"
DIAMOND_RELEASE = "2026_03"
DIAMOND_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart",
    "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp",
    "stitle",
]
PFAM_COLUMNS = [
    "protein_accession", "pfam_accession", "pfam_name", "i_evalue",
    "dom_score", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from",
    "env_to", "description",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_fingerprint(sequences: Iterator[str]) -> str:
    md5s = sorted(hashlib.md5(sequence.encode()).hexdigest() for sequence in sequences)
    return hashlib.sha256("\n".join(md5s).encode()).hexdigest()


def require_evidence(path: Path, tool: str, release: str) -> dict[str, Any]:
    if not path.exists():
        raise RuntimeError(f"Evidence absente: {path}")
    evidence = json.loads(path.read_text(encoding="utf-8"))
    if evidence.get("tool_name") != tool:
        raise RuntimeError(f"{path}: outil inattendu")
    if evidence.get("database_release") != release:
        raise RuntimeError(f"{path}: release inattendue")
    if evidence.get("run_completed") is not True:
        raise RuntimeError(f"{path}: run_completed n'est pas true")
    return evidence


def iter_pfam(path: Path) -> Iterator[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != PFAM_COLUMNS:
            raise RuntimeError(f"{path}: en-tête Pfam inattendu")
        for row in reader:
            yield row


def swissprot_accession(sseqid: str) -> str:
    parts = sseqid.split("|")
    if len(parts) >= 3 and parts[0] == "sp":
        return parts[1]
    return parts[0].split()[0]


def parse_float(value: str, field: str, path: Path) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise RuntimeError(f"{path}: {field} non numérique: {value!r}") from exc


def iter_diamond(path: Path, thresholds: dict[str, Any]) -> Iterator[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        for line_number, fields in enumerate(reader, 1):
            if len(fields) != len(DIAMOND_COLUMNS):
                raise RuntimeError(
                    f"{path}:{line_number}: {len(fields)} colonnes, "
                    f"{len(DIAMOND_COLUMNS)} attendues"
                )
            values = dict(zip(DIAMOND_COLUMNS, fields))
            numeric = {
                key: parse_float(values[key], key, path)
                for key in ("pident", "evalue", "bitscore", "qcovhsp", "scovhsp")
            }
            grade = (
                numeric["evalue"] <= float(thresholds["evalue"])
                and numeric["pident"] >= float(thresholds["min_identity"])
                and numeric["qcovhsp"] >= float(thresholds["min_query_cov"])
                and numeric["scovhsp"] >= float(thresholds["min_subject_cov"])
            )
            yield {**values, **numeric, "swissprot_accession": swissprot_accession(values["sseqid"]),
                   "is_transfer_grade": grade}


def validate_annotation_files(data_root: Path, species: str) -> dict[str, Any]:
    assembly, expected = SPECIES[species]
    base = data_root / "annotation" / species
    pfam_dir = base / "pfam_hmmscan" / PFAM_RELEASE
    diamond_dir = base / "swissprot_diamond" / DIAMOND_RELEASE
    pfam_ev = require_evidence(pfam_dir / "evidence.json", "pfam_hmmscan", PFAM_RELEASE)
    diamond_ev = require_evidence(
        diamond_dir / "evidence.json", "swissprot_diamond", DIAMOND_RELEASE
    )
    if pfam_ev["assembly_accession"] != assembly or diamond_ev["assembly_accession"] != assembly:
        raise RuntimeError(f"{species}: assemblage incohérent dans evidence.json")
    if pfam_ev["sequence_set_fingerprint"] != diamond_ev["sequence_set_fingerprint"]:
        raise RuntimeError(f"{species}: fingerprints Pfam/DIAMOND différents")
    pfam_path = pfam_dir / "domains.tsv"
    diamond_path = diamond_dir / "hits.tsv"
    for path, evidence in ((pfam_path, pfam_ev), (diamond_path, diamond_ev)):
        if not path.exists() or sha256_file(path) != evidence["output_tsv_sha256"]:
            raise RuntimeError(f"{species}: empreinte invalide pour {path}")

    pfam_rows = 0
    pfam_proteins: set[str] = set()
    for row in iter_pfam(pfam_path):
        pfam_rows += 1
        pfam_proteins.add(row["protein_accession"])
    diamond_rows = 0
    diamond_proteins: set[str] = set()
    transfer_proteins: set[str] = set()
    thresholds = diamond_ev.get("go_transfer_thresholds")
    if not isinstance(thresholds, dict):
        raise RuntimeError(f"{species}: go_transfer_thresholds absent")
    for row in iter_diamond(diamond_path, thresholds):
        diamond_rows += 1
        diamond_proteins.add(row["qseqid"])
        if row["is_transfer_grade"]:
            transfer_proteins.add(row["qseqid"])

    expected_counts = {
        "pfam_rows": pfam_ev["domain_rows"],
        "pfam_proteins": pfam_ev["proteins_with_domain"],
        "diamond_rows": diamond_ev["hit_rows"],
        "diamond_proteins": diamond_ev["proteins_with_hit"],
        "transfer_proteins": diamond_ev["proteins_with_transfer_grade_hit"],
        "transfer_rows": diamond_ev["transfer_grade_hit_rows"],
    }
    observed = {
        "pfam_rows": pfam_rows,
        "pfam_proteins": len(pfam_proteins),
        "diamond_rows": diamond_rows,
        "diamond_proteins": len(diamond_proteins),
        "transfer_proteins": len(transfer_proteins),
        "transfer_rows": sum(
            1 for row in iter_diamond(diamond_path, thresholds) if row["is_transfer_grade"]
        ),
    }
    if expected_counts != observed:
        raise RuntimeError(f"{species}: comptages != evidence: {observed} != {expected_counts}")
    return {
        "species": species,
        "assembly": assembly,
        "expected_proteins": expected,
        "pfam_evidence": pfam_ev,
        "diamond_evidence": diamond_ev,
        "pfam_path": pfam_path,
        "diamond_path": diamond_path,
        "counts": observed,
    }


def database_fingerprint(cursor: Any, assembly: str, expected: int) -> str:
    cursor.execute(
        "SELECT protein_sequence FROM dim_sequence "
        "WHERE assembly_accession=%s AND protein_sequence IS NOT NULL",
        (assembly,),
    )
    sequences = [row[0].upper() for row in cursor.fetchall()]
    if len(sequences) != expected:
        raise RuntimeError(f"{assembly}: {len(sequences)} séquences en base, attendu {expected}")
    return sequence_fingerprint(iter(sequences))


def validate_database_accessions(
    cursor: Any, assembly: str, item: dict[str, Any]
) -> None:
    cursor.execute(
        "SELECT protein_accession FROM dim_protein WHERE assemblage=%s",
        (item["species"],),
    )
    known = {row[0] for row in cursor.fetchall()}
    for path, rows in (
        (item["pfam_path"], iter_pfam(item["pfam_path"])),
        (
            item["diamond_path"],
            iter_diamond(
                item["diamond_path"],
                item["diamond_evidence"]["go_transfer_thresholds"],
            ),
        ),
    ):
        missing: set[str] = set()
        for row in rows:
            accession = (
                row["protein_accession"]
                if path == item["pfam_path"]
                else row["qseqid"]
            )
            if accession not in known:
                missing.add(accession)
                if len(missing) == 10:
                    break
        if missing:
            raise RuntimeError(
                f"{assembly}: accessions absentes de dim_protein dans {path}: "
                f"{sorted(missing)}"
            )


def run_id(evidence: dict[str, Any]) -> str:
    value = "|".join(
        str(evidence.get(key, ""))
        for key in ("assembly_accession", "tool_name", "database_release",
                    "params_hash", "output_tsv_sha256")
    )
    return hashlib.sha256(value.encode()).hexdigest()


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS annotation_runs (
    run_id TEXT PRIMARY KEY, assembly_accession TEXT NOT NULL,
    species_key TEXT NOT NULL, tool_name TEXT NOT NULL, tool_version TEXT NOT NULL,
    database_release TEXT NOT NULL, params_hash TEXT NOT NULL,
    sequence_set_fingerprint TEXT NOT NULL, input_sha256 TEXT,
    output_sha256 TEXT NOT NULL, started_at TIMESTAMPTZ NOT NULL,
    finished_at TIMESTAMPTZ NOT NULL, run_environment TEXT,
    scanned_proteins INTEGER NOT NULL, run_completed BOOLEAN NOT NULL,
    UNIQUE (assembly_accession, tool_name, database_release)
);
CREATE TABLE IF NOT EXISTS annotation_pfam_domain (
    run_id TEXT NOT NULL REFERENCES annotation_runs(run_id),
    assembly_accession TEXT NOT NULL, protein_accession TEXT NOT NULL,
    pfam_accession TEXT NOT NULL, pfam_accession_base TEXT NOT NULL,
    pfam_version INTEGER, pfam_name TEXT NOT NULL, i_evalue DOUBLE PRECISION,
    dom_score DOUBLE PRECISION, hmm_from INTEGER, hmm_to INTEGER,
    ali_from INTEGER, ali_to INTEGER, env_from INTEGER, env_to INTEGER,
    description TEXT,
    PRIMARY KEY (assembly_accession, protein_accession, pfam_accession, ali_from, ali_to)
);
CREATE TABLE IF NOT EXISTS annotation_diamond_hit (
    run_id TEXT NOT NULL REFERENCES annotation_runs(run_id),
    assembly_accession TEXT NOT NULL, qseqid TEXT NOT NULL, sseqid TEXT NOT NULL,
    swissprot_accession TEXT NOT NULL, pident DOUBLE PRECISION, length INTEGER,
    qlen INTEGER, slen INTEGER, qstart INTEGER, qend INTEGER, sstart INTEGER,
    send INTEGER, evalue DOUBLE PRECISION, bitscore DOUBLE PRECISION,
    qcovhsp DOUBLE PRECISION, scovhsp DOUBLE PRECISION, stitle TEXT NOT NULL,
    is_transfer_grade BOOLEAN NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_annotation_pfam_assembly_protein
    ON annotation_pfam_domain(assembly_accession, protein_accession);
CREATE INDEX IF NOT EXISTS idx_annotation_diamond_assembly_protein
    ON annotation_diamond_hit(assembly_accession, qseqid);
"""


def ensure_tables(cursor: Any) -> None:
    cursor.execute(SCHEMA_SQL)


def load_species(connection: Any, item: dict[str, Any], apply: bool, replace: bool) -> None:
    assembly = item["assembly"]
    tools = [
        ("pfam_hmmscan", PFAM_RELEASE, item["pfam_evidence"], item["pfam_path"]),
        ("swissprot_diamond", DIAMOND_RELEASE, item["diamond_evidence"], item["diamond_path"]),
    ]
    try:
        with connection.cursor() as cursor:
            ensure_tables(cursor)
            for tool, release, evidence, output_path in tools:
                rid = run_id(evidence)
                cursor.execute(
                    "SELECT output_sha256 FROM annotation_runs "
                    "WHERE assembly_accession=%s AND tool_name=%s AND database_release=%s",
                    (assembly, tool, release),
                )
                existing = cursor.fetchone()
                if existing and existing[0] != evidence["output_tsv_sha256"] and not replace:
                    raise RuntimeError(
                        f"{assembly}/{tool}/{release}: run différent déjà chargé; "
                        "utiliser --replace explicitement"
                    )
                if existing and existing[0] == evidence["output_tsv_sha256"] and not replace:
                    continue
                if existing:
                    cursor.execute(
                        "DELETE FROM annotation_pfam_domain WHERE assembly_accession=%s AND run_id IN "
                        "(SELECT run_id FROM annotation_runs WHERE tool_name=%s AND database_release=%s)",
                        (assembly, tool, release),
                    )
                    cursor.execute(
                        "DELETE FROM annotation_diamond_hit WHERE assembly_accession=%s AND run_id IN "
                        "(SELECT run_id FROM annotation_runs WHERE tool_name=%s AND database_release=%s)",
                        (assembly, tool, release),
                    )
                    cursor.execute(
                        "DELETE FROM annotation_runs WHERE assembly_accession=%s AND tool_name=%s AND database_release=%s",
                        (assembly, tool, release),
                    )
                now = datetime.now(timezone.utc)
                cursor.execute(
                    """INSERT INTO annotation_runs
                    (run_id, assembly_accession, species_key, tool_name, tool_version,
                     database_release, params_hash, sequence_set_fingerprint, input_sha256,
                     output_sha256, started_at, finished_at, run_environment,
                     scanned_proteins, run_completed)
                    VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)""",
                    (rid, assembly, item["species"], tool, evidence["tool_version"],
                     release, evidence["params_hash"], evidence["sequence_set_fingerprint"],
                     evidence.get("query_file_sha256"), evidence["output_tsv_sha256"],
                     evidence["started_at"], evidence["finished_at"],
                     evidence.get("run_environment"), evidence["scanned_proteins"], True),
                )
                if tool == "pfam_hmmscan":
                    rows = []
                    for row in iter_pfam(output_path):
                        accession, _, version = row["pfam_accession"].partition(".")
                        rows.append((rid, assembly, row["protein_accession"],
                                     row["pfam_accession"], accession,
                                     int(version) if version.isdigit() else None,
                                     row["pfam_name"], float(row["i_evalue"]),
                                     float(row["dom_score"]), int(row["hmm_from"]),
                                     int(row["hmm_to"]), int(row["ali_from"]),
                                     int(row["ali_to"]), int(row["env_from"]),
                                     int(row["env_to"]), row["description"]))
                    execute_values(cursor, "INSERT INTO annotation_pfam_domain VALUES %s", rows, page_size=5000)
                else:
                    rows = []
                    for row in iter_diamond(output_path, evidence["go_transfer_thresholds"]):
                        rows.append((rid, assembly, row["qseqid"], row["sseqid"],
                                     row["swissprot_accession"], row["pident"],
                                     int(row["length"]), int(row["qlen"]), int(row["slen"]),
                                     int(row["qstart"]), int(row["qend"]), int(row["sstart"]),
                                     int(row["send"]), row["evalue"], row["bitscore"],
                                     row["qcovhsp"], row["scovhsp"], row["stitle"],
                                     row["is_transfer_grade"]))
                    execute_values(cursor, "INSERT INTO annotation_diamond_hit VALUES %s", rows, page_size=5000)
        connection.commit()
    except Exception:
        connection.rollback()
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=ROOT / "pipeline_data")
    parser.add_argument("--only", nargs="*", choices=sorted(SPECIES))
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--report", type=Path, default=ROOT / "results" / "annotation_load_report.json")
    args = parser.parse_args(argv)
    if args.replace and not args.apply:
        parser.error("--replace requires --apply")
    ready = ROOT / "check_annotation_ready.py"
    check = subprocess.run(
        [sys.executable, str(ready), "--data-root", str(args.data_root)],
        check=False,
    )
    if check.returncode:
        return check.returncode
    targets = args.only or list(SPECIES)
    items = [validate_annotation_files(args.data_root, species) for species in targets]
    from core_engines.bioinfo_db import connect
    report: dict[str, Any] = {"dry_run": not args.apply, "species": {}}
    with connect() as connection:
        with connection.cursor() as cursor:
            for item in items:
                actual = database_fingerprint(cursor, item["assembly"], item["expected_proteins"])
                expected = item["pfam_evidence"]["sequence_set_fingerprint"]
                if actual != expected:
                    raise RuntimeError(f"{item['species']}: fingerprint dim_sequence différente")
                validate_database_accessions(cursor, item["assembly"], item)
                report["species"][item["species"]] = {
                    "assembly": item["assembly"], "sequence_set_fingerprint": actual,
                    "counts": item["counts"], "status": "validated",
                }
        if args.apply:
            for item in items:
                load_species(connection, item, True, args.replace)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
