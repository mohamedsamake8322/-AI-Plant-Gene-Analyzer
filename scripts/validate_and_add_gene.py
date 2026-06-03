#!/usr/bin/env python3
"""
Validate and add gene records to genes_database.json using genes_schema.json.
Usage:
  python scripts/validate_and_add_gene.py validate --db genes_database.json
  python scripts/validate_and_add_gene.py add --db genes_database.json --file new_gene.json
"""

import sys
import json
from pathlib import Path
import datetime

ROOT = Path(__file__).resolve().parents[1]
DB_PATH_DEFAULT = ROOT / "genes_database.json"
SCHEMA_PATH = ROOT / "genes_schema.json"


def clean_sequence(s: str) -> str:
    return "".join(s.split()).upper()


def load_db(path: Path) -> dict:
    if not path.exists():
        return {"metadata": {"version": "1.0", "last_update": datetime.datetime.utcnow().isoformat()+"Z", "sources": []}, "genes": []}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_db(path: Path, data: dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# Try to use jsonschema if available, otherwise fallback to simple checks
try:
    import jsonschema
    _HAS_JSONSCHEMA = True
except Exception:
    _HAS_JSONSCHEMA = False


def validate_with_schema(record: dict) -> tuple[bool, str]:
    if _HAS_JSONSCHEMA:
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        try:
            jsonschema.validate(instance={"metadata": {"version":"1.0","last_update":datetime.datetime.utcnow().isoformat()+"Z","sources":[]}, "genes": [record]}, schema=schema)
            return True, ""
        except jsonschema.ValidationError as e:
            return False, str(e)
    # fallback simple validation
    return fallback_validate(record)


def fallback_validate(record: dict) -> tuple[bool, str]:
    if not isinstance(record, dict):
        return False, "Record must be a JSON object"
    for k in ("gene_id", "symbol", "sequence", "sequence_type"):
        if k not in record or not record[k]:
            return False, f"Missing required field: {k}"
    if record["sequence_type"] not in ("dna", "protein"):
        return False, "sequence_type must be 'dna' or 'protein'"
    seq = clean_sequence(record["sequence"])
    if not seq:
        return False, "Empty sequence"
    # basic alphabet check
    if record["sequence_type"] == "dna":
        for c in seq:
            if c not in set("ACGTN"):
                return False, f"Invalid character in DNA sequence: {c}"
    return True, ""


def add_record_to_db(record: dict, db_path: Path) -> tuple[bool, str]:
    db = load_db(db_path)
    rec = dict(record)  # copy
    rec["sequence"] = clean_sequence(rec["sequence"])
    rec["length"] = len(rec["sequence"])
    rec.setdefault("date_added", datetime.datetime.utcnow().isoformat() + "Z")
    # ensure unique gene_id or symbol
    for existing in db.get("genes", []):
        if existing.get("gene_id") == rec.get("gene_id") or existing.get("symbol") == rec.get("symbol"):
            return False, f"Record with same gene_id or symbol already exists: {rec.get('gene_id') or rec.get('symbol')}"
    valid, msg = validate_with_schema(rec)
    if not valid:
        return False, msg
    db.setdefault("genes", []).append(rec)
    db["metadata"] = db.get("metadata", {})
    db["metadata"]["last_update"] = datetime.datetime.utcnow().isoformat() + "Z"
    save_db(db_path, db)
    return True, "Added"


def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")

    p_val = sub.add_parser("validate")
    p_val.add_argument("--db", default=str(DB_PATH_DEFAULT))

    p_add = sub.add_parser("add")
    p_add.add_argument("--db", default=str(DB_PATH_DEFAULT))
    p_add.add_argument("--file", required=True, help="JSON file with a single gene record")

    args = p.parse_args(argv)
    if args.cmd == "validate":
        db = load_db(Path(args.db))
        print("DB has", len(db.get("genes", [])), "genes")
        # validate all entries
        bad = 0
        for i, rec in enumerate(db.get("genes", []), 1):
            ok, msg = validate_with_schema(rec)
            if not ok:
                print(f"[{i}] INVALID: {msg}")
                bad += 1
        if bad:
            print(f"Validation failed: {bad} bad records")
            sys.exit(2)
        print("All records valid")
        return

    if args.cmd == "add":
        path = Path(args.file)
        if not path.exists():
            print("Input file not found", args.file)
            return
        rec = json.loads(path.read_text(encoding="utf-8"))
        ok, msg = add_record_to_db(rec, Path(args.db))
        if not ok:
            print("Failed:", msg)
            sys.exit(3)
        print("Success:", msg)
        return

    p.print_help()


if __name__ == "__main__":
    main(sys.argv[1:])
