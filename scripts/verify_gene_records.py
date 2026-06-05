#!/usr/bin/env python3
"""
List and verify gene records imported into the project (Postgres or cleaned JSON).

Usage examples:
  python scripts/verify_gene_records.py --db
  python scripts/verify_gene_records.py --json-file data/clean/plant_data_clean.json --csv-out report.csv --split-out-dir reports
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any


def try_import_postgres_utils() -> Any:
    try:
        import postgres_utils as pu

        return pu
    except Exception:
        try:
            from scripts import postgres_utils as pu

            return pu
        except Exception:
            return None


def load_from_json(path: Path) -> list[dict]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(raw, dict) and "genes" in raw:
        return raw["genes"]
    if isinstance(raw, list):
        return raw
    if isinstance(raw, dict) and raw.get("gene_id"):
        return [raw]
    if isinstance(raw, dict):
        for k in ("genes", "ensembl", "ncbi", "geo", "expression_atlas"):
            if k in raw and isinstance(raw[k], list):
                return raw[k]
    return []


def resolve_json_path(path: Path) -> Path:
    if path.exists():
        return path
    alt = Path("data") / "clean" / path.name
    if alt.exists():
        return alt
    raise FileNotFoundError(path)


def classify_record(rec: dict) -> str:
    seq = rec.get("sequence") or rec.get("seq") or ""
    seq = str(seq).strip()
    seq_type = str(rec.get("sequence_type") or rec.get("type") or rec.get("bio_type") or "").lower()
    if seq:
        if seq_type in ("protein", "aa", "amino_acid", "peptide"):
            return "protein"
        if seq_type in ("rna", "mrna", "trna", "rrna", "transcript"):
            return "rna"
        if seq_type in ("dna", "genomic", "cdna", "cds", "nucleotide"):
            return "dna"
        if any(char in seq.upper() for char in "BFJLOUZX"):
            return "protein"
        if all(char in "ACGTUN" for char in seq.upper()):
            return "dna"
        if all(char in "ACGTURYSWKMBDHVN" for char in seq.upper()):
            return "rna"
        return "unknown"
    source = str(rec.get("source") or "").lower()
    if "geo" in source or rec.get("accession", "").startswith("GSE"):
        return "metadata_geo"
    if "atlas" in source or rec.get("experiment_id"):
        return "metadata_atlas"
    return "metadata"


def normalize_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=False)


def create_report_row(rec: dict, preview_length: int, show_seq: bool) -> dict:
    gene_id = normalize_value(rec.get("gene_id") or rec.get("accession") or rec.get("experiment_id") or "")
    symbol = normalize_value(rec.get("symbol") or rec.get("gene_query") or "")
    organism = normalize_value(rec.get("organism") or rec.get("species") or "")
    seq = normalize_value(rec.get("sequence") or rec.get("seq") or "")
    seq_type = normalize_value(rec.get("sequence_type") or rec.get("type") or rec.get("bio_type") or "")
    record_type = classify_record(rec)
    length = rec.get("length") or (len(seq) if seq else "")
    source = normalize_value(rec.get("source") or "")
    date_added = normalize_value(rec.get("date_added") or "")
    url = normalize_value(rec.get("url") or rec.get("api_url") or "")
    title = normalize_value(rec.get("title") or rec.get("record", {}).get("description") or "")
    preview = seq[:preview_length] + ("..." if len(seq) > preview_length else "") if seq else ""
    row = {
        "gene_id": gene_id,
        "symbol": symbol,
        "organism": organism,
        "type": record_type,
        "sequence_type": seq_type,
        "length": length,
        "source": source,
        "date_added": date_added,
        "title": title,
        "url": url,
        "sequence_preview": preview,
    }
    if show_seq:
        row["full_sequence"] = seq
    return row


def write_split_json(rows: list[dict], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    groups: dict[str, list[dict]] = {}
    for row in rows:
        groups.setdefault(row["type"], []).append(row)
    for rec_type, items in groups.items():
        path = out_dir / f"{rec_type}.json"
        path.write_text(json.dumps(items, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Wrote {len(items)} records to {path}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Verify gene records from Postgres or cleaned JSON")
    p.add_argument("--db", action="store_true", help="Load records from PostgreSQL using scripts/postgres_utils.py")
    p.add_argument("--json-file", help="Path to cleaned JSON file (default: data/clean/plant_data_clean.json)")
    p.add_argument("--limit", type=int, default=0, help="Limit number of records displayed (0 = all)")
    p.add_argument("--csv-out", help="Write CSV report to this file")
    p.add_argument("--json-out", help="Write JSON report to this file")
    p.add_argument("--split-out-dir", help="Write separate JSON files by record type")
    p.add_argument("--show-seq", action="store_true", help="Show sequence (preview length controlled by --preview) if present")
    p.add_argument("--preview", type=int, default=80, help="Sequence preview length")
    args = p.parse_args(argv)

    records: list[dict] = []

    if args.db:
        pu = try_import_postgres_utils()
        if not pu:
            print("Cannot import postgres_utils. Ensure you run from project root and have dependencies installed.")
            sys.exit(2)
        try:
            dbmap = pu.load_gene_database_from_postgres()
            records = list(dbmap.values())
        except Exception as e:
            print("Failed to load from Postgres:", e)
            sys.exit(3)
    else:
        path = Path(args.json_file or Path("data") / "clean" / "plant_data_clean.json")
        try:
            path = resolve_json_path(path)
        except FileNotFoundError:
            print("JSON file not found:", path)
            sys.exit(2)
        records = load_from_json(path)

    total = len(records)
    print(f"Loaded {total} record(s)")
    if args.limit and args.limit > 0:
        records = records[: args.limit]

    rows = [create_report_row(rec, args.preview, args.show_seq) for rec in records]
    type_counts: dict[str, int] = {}
    for row in rows:
        type_counts[row["type"]] = type_counts.get(row["type"], 0) + 1

    print("\nType counts:")
    for rec_type, count in sorted(type_counts.items()):
        print(f"  {rec_type}: {count}")

    print("\nSummary:")
    print(f"{'gene_id':<20}{'symbol':<15}{'organism':<30}{'type':<15}{'len':>6}{'source':>15}")
    print("-" * 110)
    for r in rows:
        print(f"{r['gene_id'][:20]:<20}{r['symbol'][:15]:<15}{r['organism'][:30]:<30}{r['type'][:15]:<15}{str(r['length']):>6}{str(r['source'])[:15]:>15}")
        if args.show_seq and r.get("full_sequence"):
            print(f"  SEQ: {r['full_sequence']}\n")
        elif r['sequence_preview']:
            print(f"  PREVIEW: {r['sequence_preview']}\n")

    if args.csv_out:
        out_path = Path(args.csv_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["gene_id"])
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        print(f"Wrote CSV report to {out_path}")

    if args.json_out:
        out_path = Path(args.json_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Wrote JSON report to {out_path}")

    if args.split_out_dir:
        write_split_json(rows, Path(args.split_out_dir))


if __name__ == "__main__":
    main()
