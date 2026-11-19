#!/usr/bin/env python3
"""Rebuild sequence/protein dimensions from NCBI FASTA and GFF3.

Default mode is a dry run. ``--apply`` writes only after the preview has been
inspected; ``--prune`` archives rows absent from the source before removing
them from dim_protein.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import warnings
import re
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import unquote

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from psycopg2.extras import execute_values

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect


def opener(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("rt", encoding="utf-8", errors="replace")


def header_fields(header: str) -> dict[str, str]:
    fields = {}
    for key, value in re.findall(r"\[([^=\]]+)=([^\]]*)\]", header):
        fields[key] = value
    return fields


def read_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    current = None
    with opener(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                header = line[1:]
                accession = header.split(None, 1)[0].split("|")[-1]
                fields = header_fields(header)
                protein_id = fields.get("protein_id") or accession
                current = protein_id
                records[current] = {
                    "accession": protein_id,
                    "description": header,
                    "sequence": "",
                    "gene_id": fields.get("GeneID") or fields.get("gene"),
                    "gene_symbol": fields.get("gene") or fields.get("locus_tag"),
                }
            elif current:
                records[current]["sequence"] += line.strip().upper()
    return records


def read_cds(path: Path) -> dict[str, str]:
    return {
        key: record["sequence"]
        for key, record in read_fasta(path).items()
    }


def read_assembly_report(path: Path) -> dict[str, str]:
    mapping = {}
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 7:
                mapping[fields[6]] = fields[2]
                mapping[fields[4]] = fields[2]
    return mapping


def read_gff(
    path: Path, chromosome_names: dict[str, str] | None = None
) -> dict[str, dict[str, object]]:
    features: dict[str, dict[str, object]] = {}
    result: dict[str, dict[str, object]] = {}
    with opener(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"CDS", "gene", "mRNA"}:
                continue
            attrs = {
                key: unquote(value)
                for key, value in re.findall(
                    r"(?:^|;)([^=;]+)=([^;]*)", fields[8]
                )
            }
            if attrs.get("ID"):
                features[attrs["ID"]] = {
                    "parent": attrs.get("Parent"),
                    "attrs": attrs,
                }
            protein = attrs.get("protein_id") or attrs.get("Name")
            if not protein:
                continue
            if fields[2] != "CDS":
                continue
            merged = dict(attrs)
            parent_id = attrs.get("Parent", "").split(",", 1)[0]
            visited = set()
            while parent_id and parent_id not in visited:
                visited.add(parent_id)
                parent = features.get(parent_id, {})
                parent_attrs = parent.get("attrs", {})
                for key, value in parent_attrs.items():
                    merged.setdefault(key, value)
                parent_id = str(parent.get("parent") or "").split(",", 1)[0]
            gene_id = None
            for reference in merged.get("Dbxref", "").split(","):
                if reference.startswith("GeneID:"):
                    gene_id = reference.removeprefix("GeneID:")
                    break
            row = result.setdefault(protein, {
                "chr": (chromosome_names or {}).get(fields[0], fields[0]),
                "start_pos": int(fields[3]),
                "stop_pos": int(fields[4]), "strand": fields[6],
                "gene_id": gene_id,
                "gene_symbol": merged.get("gene") or merged.get("gene_symbol"),
                "locus_tag": merged.get("locus_tag"),
                "product": merged.get("product"),
                "transl_table": merged.get("transl_table"),
                "gff_notes": [],
                "gff_exception": [],
                "gff_inference": [],
                "partial": bool(merged.get("partial")),
                "transl_except": bool(merged.get("transl_except")),
                "phase": fields[7],
            })
            row["start_pos"] = min(row["start_pos"], int(fields[3]))
            row["stop_pos"] = max(row["stop_pos"], int(fields[4]))
            row["gene_id"] = row["gene_id"] or gene_id
            row["gene_symbol"] = row["gene_symbol"] or merged.get("gene")
            row["locus_tag"] = row["locus_tag"] or merged.get("locus_tag")
            row["transl_table"] = row["transl_table"] or merged.get("transl_table")
            row["partial"] = row["partial"] or bool(merged.get("partial"))
            row["transl_except"] = row["transl_except"] or bool(merged.get("transl_except"))
            row["phase"] = row["phase"] or fields[7]
            for key, target in (
                ("Note", "gff_notes"),
                ("exception", "gff_exception"),
                ("inference", "gff_inference"),
            ):
                value = merged.get(key)
                if value:
                    row[target].extend(
                        item.strip() for item in value.split(",") if item.strip()
                    )
    return result


def sequence_md5(sequence: str) -> str:
    return hashlib.md5(sequence.encode("ascii", "ignore")).hexdigest()


def _translation_alignment(expected: str, observed: str) -> dict[str, int | float]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(expected, observed)[0]
    coordinates = alignment.coordinates
    matches = mismatches = gaps = aligned_length = 0
    for index in range(coordinates.shape[1] - 1):
        start_a, end_a = coordinates[0, index:index + 2]
        start_b, end_b = coordinates[1, index:index + 2]
        length_a = int(end_a - start_a)
        length_b = int(end_b - start_b)
        if length_a and length_b:
            compared = min(length_a, length_b)
            left = expected[int(start_a):int(start_a) + compared]
            right = observed[int(start_b):int(start_b) + compared]
            matches += sum(a == b for a, b in zip(left, right))
            mismatches += compared - sum(a == b for a, b in zip(left, right))
            aligned_length += compared
            if length_a != length_b:
                gaps += abs(length_a - length_b)
                aligned_length += abs(length_a - length_b)
        else:
            gaps += max(length_a, length_b)
            aligned_length += max(length_a, length_b)
    return {
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "alignment_length": aligned_length,
        "global_identity": matches / max(aligned_length, 1),
        "min_length_identity": matches / max(min(len(expected), len(observed)), 1),
    }


def translation_check(proteins: dict, cds: dict, coordinates: dict) -> dict:
    compared = 0
    exceptions = []
    status_counts = Counter()
    legacy_counts = Counter()
    status_by_accession = {}
    for accession, record in proteins.items():
        sequence = cds.get(accession)
        if not sequence:
            continue
        compared += 1
        gff_row = coordinates.get(accession, {})
        table = int(gff_row.get("transl_table") or 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            translated = str(Seq(sequence).translate(table=table, to_stop=False))
        expected = record["sequence"].rstrip("*")
        observed = translated.rstrip("*")
        matches = sum(a == b for a, b in zip(expected, observed))
        legacy_identity = matches / max(len(expected), len(observed), 1)
        alignment = _translation_alignment(expected, observed)
        if expected == observed:
            status = "exact"
        elif (
            alignment["global_identity"] >= 0.99
            and abs(len(expected) - len(observed)) <= 3
        ):
            status = "near"
        else:
            status = "divergent"
        status_counts[status] += 1
        status_by_accession[accession] = status
        if legacy_identity < 0.99 or len(expected) != len(observed):
            legacy_counts["below_99_or_length"] += 1
            exceptions.append({
                "protein_accession": accession,
                "legacy_identity": round(legacy_identity, 6),
                **{key: round(value, 6) if isinstance(value, float) else value
                   for key, value in alignment.items()},
                "cds_matches_protein": status,
                "protein_length": len(expected),
                "translated_length": len(observed),
                "cds_length": len(sequence),
                "cds_length_mod_3": len(sequence) % 3,
                "transl_table": table,
                "gff_partial": bool(gff_row.get("partial")),
                "gff_transl_except": bool(gff_row.get("transl_except")),
                "gff_exception": list(dict.fromkeys(gff_row.get("gff_exception", []))),
                "gff_inference": list(dict.fromkeys(gff_row.get("gff_inference", []))),
                "gff_explained": bool(
                    gff_row.get("partial")
                    or gff_row.get("transl_except")
                    or gff_row.get("gff_exception")
                ),
                "gff_modified_relative_to_genomic": any(
                    "modified relative to" in text.lower()
                    and "genomic" in text.lower()
                    for text in (
                        gff_row.get("gff_notes", [])
                        + gff_row.get("gff_exception", [])
                        + gff_row.get("gff_inference", [])
                    )
                ),
            })
    return {
        "compared": compared,
        "status_counts": dict(status_counts),
        "status_by_accession": status_by_accession,
        "legacy_exception_count": len(exceptions),
        "legacy_counts": dict(legacy_counts),
        "exception_prefix_counts": dict(
            Counter(item["protein_accession"].split("_", 1)[0] for item in exceptions)
        ),
        "exception_global_identity_buckets": dict(
            Counter(
                "<50%" if item["global_identity"] < 0.5
                else "50-90%" if item["global_identity"] < 0.9
                else "90-99%" if item["global_identity"] < 0.99
                else ">=99%"
                for item in exceptions
            )
        ),
        "exception_legacy_identity_buckets": dict(
            Counter(
                "<50%" if item["legacy_identity"] < 0.5
                else "50-90%" if item["legacy_identity"] < 0.9
                else "90-99%" if item["legacy_identity"] < 0.99
                else ">=99%"
                for item in exceptions
            )
        ),
        "above_99_pct": compared - len(exceptions),
        "modified_relative_to_genomic": sum(
            item["gff_modified_relative_to_genomic"] for item in exceptions
        ),
        "gff_explained_count": sum(item["gff_explained"] for item in exceptions),
        "exceptions": exceptions,
    }


def canonical_chromosome(value: object) -> str:
    match = re.search(r"(?:chr|chromosome)?0*(\d+)$", str(value or ""), re.IGNORECASE)
    return str(int(match.group(1))) if match else str(value or "").lower()


def coordinate_status(plaza: tuple, gff_row: dict | None) -> str:
    if not plaza or not gff_row:
        return "non_verifiable"
    chromosome, start, stop, strand = plaza
    if start is None or stop is None or strand is None:
        return "non_verifiable"
    if canonical_chromosome(chromosome) != canonical_chromosome(gff_row["chr"]):
        return "discordant"
    if strand != gff_row["strand"]:
        return "discordant"
    overlap = max(
        0,
        min(stop, gff_row["stop_pos"]) - max(start, gff_row["start_pos"]) + 1,
    )
    cds_length = gff_row["stop_pos"] - gff_row["start_pos"] + 1
    return "concordant" if overlap / max(cds_length, 1) >= 0.5 else "discordant"


def _prune_absent_proteins(cursor, assemblage: str, accessions: set[str]) -> None:
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS dim_protein_archive AS
        SELECT p.*, now() AS archived_at, NULL::text AS successor_accession,
               NULL::text AS archive_reason
        FROM dim_protein p WHERE false
    """)
    cursor.execute("""
        ALTER TABLE dim_protein_archive
            ADD COLUMN IF NOT EXISTS successor_accession TEXT,
            ADD COLUMN IF NOT EXISTS archive_reason TEXT
    """)
    cursor.execute("""
        CREATE TEMP TABLE rebuild_source_accessions (
            protein_accession TEXT PRIMARY KEY
        ) ON COMMIT DROP
    """)
    execute_values(
        cursor,
        "INSERT INTO rebuild_source_accessions VALUES %s",
        [(accession,) for accession in accessions],
    )
    cursor.execute("""
        INSERT INTO dim_protein_archive
        SELECT p.*, now(), NULL::text, 'absent_from_source'
        FROM dim_protein p
        WHERE p.assemblage=%s AND NOT EXISTS (
          SELECT 1 FROM rebuild_source_accessions s
          WHERE s.protein_accession=p.protein_accession)
          AND NOT EXISTS (
              SELECT 1 FROM dim_protein_archive a
              WHERE a.assemblage=p.assemblage
                AND a.protein_accession=p.protein_accession)
    """, (assemblage,))
    cursor.execute("""
        CREATE TEMP TABLE rebuild_old_accessions (
            protein_accession TEXT PRIMARY KEY
        ) ON COMMIT DROP
    """)
    cursor.execute("""
        INSERT INTO rebuild_old_accessions
        SELECT p.protein_accession
        FROM dim_protein p
        WHERE p.assemblage=%s
          AND NOT EXISTS (
              SELECT 1 FROM rebuild_source_accessions s
              WHERE s.protein_accession=p.protein_accession)
        UNION
        SELECT p.protein_accession
        FROM dim_protein_archive p
        WHERE p.assemblage=%s
          AND NOT EXISTS (
              SELECT 1 FROM rebuild_source_accessions s
              WHERE s.protein_accession=p.protein_accession)
    """, (assemblage, assemblage))
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS fact_functional_annotation_archive
        AS SELECT a.*, now() AS archived_at, NULL::text AS archive_reason
        FROM fact_functional_annotation a WHERE false
    """)
    cursor.execute("""
        INSERT INTO fact_functional_annotation_archive
        SELECT a.*, now(), 'protein_absent_from_source'
        FROM fact_functional_annotation a
        JOIN rebuild_old_accessions old
          ON old.protein_accession=a.protein_accession
        WHERE a.assemblage=%s
          AND NOT EXISTS (
              SELECT 1 FROM fact_functional_annotation_archive old
              WHERE old.assemblage=a.assemblage
                AND old.protein_accession=a.protein_accession
                AND old.id=a.id)
    """, (assemblage,))
    cursor.execute("""
        DELETE FROM fact_functional_annotation a
        USING rebuild_old_accessions old
        WHERE a.assemblage=%s
          AND a.protein_accession=old.protein_accession
    """, (assemblage,))
    cursor.execute("""
        DELETE FROM dim_protein p
        WHERE p.assemblage=%s AND NOT EXISTS (
          SELECT 1 FROM rebuild_source_accessions s
          WHERE s.protein_accession=p.protein_accession)
    """, (assemblage,))


def rebuild(assemblage: str, assembly_accession: str, protein_fasta: Path,
            cds_fasta: Path, gff: Path, assembly_report: Path | None = None,
            apply: bool = False, prune: bool = False,
            confirm_production: bool = False) -> dict:
    proteins = read_fasta(protein_fasta)
    cds = read_cds(cds_fasta)
    chromosome_names = read_assembly_report(assembly_report) if assembly_report else {}
    coordinates = read_gff(gff, chromosome_names)
    translation = translation_check(proteins, cds, coordinates)
    accessions = set(proteins)
    with connect() as connection, connection.cursor() as cursor:
        if apply and connection.info.dbname == "bioinfo_pipeline" and not confirm_production:
            raise RuntimeError(
                "Refus d'écrire dans bioinfo_pipeline sans "
                "--confirm-production."
            )
        cursor.execute("SELECT protein_accession FROM dim_protein WHERE assemblage=%s", (assemblage,))
        existing = {row[0] for row in cursor.fetchall()}
        cursor.execute(
            """
            SELECT protein_accession, gene_id, gene_symbol, description,
                   chr, start_pos, stop_pos, strand
            FROM dim_protein WHERE assemblage=%s
            """,
            (assemblage,),
        )
        existing_protein_rows = {
            row[0]: row[1:] for row in cursor.fetchall()
        }
        cursor.execute(
            """
            SELECT protein_accession, term_desc
            FROM fact_functional_annotation
            WHERE assemblage=%s AND source='description'
            """,
            (assemblage,),
        )
        existing_descriptions = {}
        for accession, description in cursor.fetchall():
            existing_descriptions.setdefault(accession, set()).add(description or "")
        description_counts = {"new": 0, "same": 0, "different": 0}
        for accession, record in proteins.items():
            values = existing_descriptions.get(accession)
            if not values:
                description_counts["new"] += 1
            elif record["description"] in values:
                description_counts["same"] += 1
            else:
                description_counts["different"] += 1
        desired_protein_rows = {
            key: (
                coordinates.get(key, {}).get("gene_id") or value.get("gene_id"),
                coordinates.get(key, {}).get("gene_symbol") or value.get("gene_symbol"),
                value["description"],
                coordinates.get(key, {}).get("chr"),
                coordinates.get(key, {}).get("start_pos"),
                coordinates.get(key, {}).get("stop_pos"),
                coordinates.get(key, {}).get("strand"),
            )
            for key, value in proteins.items()
        }
        changed_proteins = sum(
            accession in existing_protein_rows
            and existing_protein_rows[accession] != desired_protein_rows[accession]
            for accession in accessions
        )
        cursor.execute(
            """
            SELECT gene_id, protein_accession, chr, start_pos, stop_pos, strand
            FROM dim_gene_protein
            WHERE assembly_name=%s AND gene_source='plaza'
            """,
            (assemblage,),
        )
        concordance = {"concordant": 0, "discordant": 0, "non_verifiable": 0}
        distance_histogram = {}
        plaza_links = cursor.fetchall()
        for gene_id, accession, chromosome, start, stop, strand in plaza_links:
            status = coordinate_status(
                (chromosome, start, stop, strand), coordinates.get(accession)
            )
            concordance[status] += 1
            if status == "discordant" and accession in coordinates:
                gff_row = coordinates[accession]
                distance = max(
                    0,
                    start - gff_row["stop_pos"],
                    gff_row["start_pos"] - stop,
                )
                bucket = "<10" if distance < 10 else "10-99" if distance < 100 else "100-999" if distance < 1000 else ">=1000"
                distance_histogram[bucket] = distance_histogram.get(bucket, 0) + 1
        cursor.execute(
            "SELECT to_regclass('public.dim_sequence') IS NOT NULL"
        )
        sequence_exists = cursor.fetchone()[0]
        if sequence_exists:
            cursor.execute(
                "SELECT protein_accession, sequence_md5 FROM dim_sequence WHERE assembly_accession=%s",
                (assembly_accession,),
            )
            existing_sequences = {row[0]: row[1] for row in cursor.fetchall()}
        else:
            existing_sequences = {}
        preview = {
            "source_proteins": len(accessions),
            "existing_dim_protein": len(existing),
            "new_proteins": len(accessions - existing),
            "updated_proteins": len(accessions & existing),
            "changed_proteins": changed_proteins,
            "absent_from_source": len(existing - accessions),
            "cds_matched": len(accessions & set(cds)),
            "coordinates_matched": len(accessions & set(coordinates)),
            "geneid_links": sum(bool(row.get("gene_id")) for row in coordinates.values()),
            "description_comparison": description_counts,
            "plaza_coordinate_concordance": concordance,
            "plaza_coordinate_distance_histogram": distance_histogram,
            "translation": translation,
            "sequences_new": len(accessions - set(existing_sequences)),
            "sequences_updated": len(accessions & set(existing_sequences)),
            "sequences_changed": sum(
                accession in existing_sequences
                and existing_sequences[accession] != sequence_md5(record["sequence"])
                for accession, record in proteins.items()
            ),
            "dry_run": not apply,
        }
        if not apply:
            return preview
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS dim_sequence (
              assembly_accession TEXT NOT NULL, protein_accession TEXT NOT NULL,
              protein_sequence TEXT, cds_sequence TEXT, sequence_md5 TEXT,
              cds_matches_protein TEXT,
              PRIMARY KEY (assembly_accession, protein_accession))
        """)
        cursor.execute(
            "ALTER TABLE dim_sequence ADD COLUMN IF NOT EXISTS "
            "cds_matches_protein TEXT"
        )
        for column, definition in (
            ("gene_id", "TEXT"), ("gene_symbol", "TEXT"), ("locus_tag", "TEXT"),
            ("description", "TEXT")
        ):
            cursor.execute(f"ALTER TABLE dim_protein ADD COLUMN IF NOT EXISTS {column} {definition}")
        for column in ("gene_symbol", "locus_tag"):
            cursor.execute(
                f"ALTER TABLE dim_gene_protein ADD COLUMN IF NOT EXISTS {column} TEXT"
            )
        cursor.execute(
            "ALTER TABLE dim_gene_protein ADD COLUMN IF NOT EXISTS "
            "coordinate_concordance TEXT"
        )
        execute_values(cursor, """
            INSERT INTO dim_protein
              (protein_accession, assemblage, gene_id, gene_symbol, description,
               chr, start_pos, stop_pos, strand)
            VALUES %s
            ON CONFLICT (assemblage, protein_accession) DO UPDATE SET
              gene_id=COALESCE(dim_protein.gene_id, EXCLUDED.gene_id),
              gene_symbol=COALESCE(dim_protein.gene_symbol, EXCLUDED.gene_symbol),
              description=COALESCE(dim_protein.description, EXCLUDED.description),
              chr=COALESCE(EXCLUDED.chr, dim_protein.chr),
              start_pos=COALESCE(EXCLUDED.start_pos, dim_protein.start_pos),
              stop_pos=COALESCE(EXCLUDED.stop_pos, dim_protein.stop_pos),
              strand=COALESCE(EXCLUDED.strand, dim_protein.strand)
        """, [
            (key, assemblage,
             coordinates.get(key, {}).get("gene_id") or proteins[key].get("gene_id"),
             coordinates.get(key, {}).get("gene_symbol") or proteins[key].get("gene_symbol"),
             proteins[key]["description"],
             coordinates.get(key, {}).get("chr"),
             coordinates.get(key, {}).get("start_pos"),
             coordinates.get(key, {}).get("stop_pos"),
             coordinates.get(key, {}).get("strand"))
            for key in proteins
        ], page_size=500)
        execute_values(cursor, """
            INSERT INTO dim_gene_protein
                (assembly_name, gene_source, gene_id, protein_accession,
                 chr, start_pos, stop_pos, strand, gene_symbol, locus_tag)
            VALUES %s
            ON CONFLICT DO NOTHING
        """, [
            (assemblage, "ncbi_geneid", coordinates[key]["gene_id"], key,
             coordinates[key].get("chr"), coordinates[key].get("start_pos"),
             coordinates[key].get("stop_pos"), coordinates[key].get("strand"),
             coordinates[key].get("gene_symbol"), coordinates[key].get("locus_tag"))
            for key in proteins
            if coordinates.get(key, {}).get("gene_id")
        ], page_size=500)
        execute_values(cursor, """
            INSERT INTO dim_sequence
              (assembly_accession, protein_accession, protein_sequence, cds_sequence,
               sequence_md5, cds_matches_protein)
            VALUES %s
            ON CONFLICT (assembly_accession, protein_accession) DO UPDATE SET
              protein_sequence=EXCLUDED.protein_sequence, cds_sequence=EXCLUDED.cds_sequence,
              sequence_md5=EXCLUDED.sequence_md5,
              cds_matches_protein=EXCLUDED.cds_matches_protein
        """, [
            (assembly_accession, key, value["sequence"], cds.get(key),
             sequence_md5(value["sequence"]),
             translation["status_by_accession"].get(key))
            for key, value in proteins.items()
        ], page_size=500)
        cursor.execute("""
            CREATE TEMP TABLE rebuild_products (
                assemblage TEXT,
                protein_accession TEXT,
                term_desc TEXT
            ) ON COMMIT DROP
        """)
        execute_values(
            cursor,
            "INSERT INTO rebuild_products VALUES %s",
            [
                (
                    assemblage,
                    key,
                    coordinates.get(key, {}).get("product")
                    or value["description"],
                )
                for key, value in proteins.items()
            ],
            page_size=500,
        )
        cursor.execute("""
            INSERT INTO fact_functional_annotation
                (assemblage, protein_accession, source, term_id, term_desc, evidence)
            SELECT d.assemblage, d.protein_accession, 'refseq_product',
                   d.protein_accession, d.term_desc, 'NCBI FASTA'
            FROM rebuild_products d
            WHERE NOT EXISTS (
                SELECT 1 FROM fact_functional_annotation a
                WHERE a.assemblage=d.assemblage
                  AND a.protein_accession=d.protein_accession
                  AND a.source='refseq_product'
                  AND a.term_desc=d.term_desc
            )
        """)
        execute_values(
            cursor,
            """
            UPDATE dim_gene_protein AS link
            SET coordinate_concordance = data.status
            FROM (VALUES %s) AS data(
                assembly_name, gene_source, gene_id, protein_accession, status
            )
            WHERE link.assembly_name=data.assembly_name
              AND link.gene_source=data.gene_source
              AND link.gene_id=data.gene_id
              AND link.protein_accession=data.protein_accession
            """,
            [
                (
                    assemblage,
                    "plaza",
                    gene_id,
                    accession,
                    coordinate_status(
                        (chromosome, start, stop, strand),
                        coordinates.get(accession),
                    ),
                )
                for gene_id, accession, chromosome, start, stop, strand in plaza_links
            ],
            page_size=500,
        )
        if prune:
            _prune_absent_proteins(cursor, assemblage, accessions)
        preview["dry_run"] = False
        return preview


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblage", required=True)
    parser.add_argument("--assembly-accession", required=True)
    parser.add_argument("--protein-fasta", type=Path, required=True)
    parser.add_argument("--cds-fasta", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--prune", action="store_true")
    parser.add_argument("--confirm-production", action="store_true")
    args = parser.parse_args()
    if args.prune and not args.apply:
        parser.error("--prune requires --apply")
    import json
    print(json.dumps(rebuild(
        args.assemblage, args.assembly_accession, args.protein_fasta,
        args.cds_fasta, args.gff, args.assembly_report,
        args.apply, args.prune, args.confirm_production,
    ), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
