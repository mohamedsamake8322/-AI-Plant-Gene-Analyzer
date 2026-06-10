#!/usr/bin/env python3
"""
Collect multiple sequence types (DNA, RNA, protein) for a plant species.
Merges all results into a single clean JSON and imports to PostgreSQL.

Usage:
  python scripts/collect_multi_type.py \
    --plant "Oryza sativa" \
    --out-clean data/clean/rice_all_types.json \
    --retmax 300
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_pipeline as pipeline_module
import clean_data as clean_data


def collect_and_clean_type(
    plant_term: str,
    seq_type: str,
    retmax: int,
    temp_raw: Path,
    temp_clean: Path,
) -> list[dict]:
    """Collect one sequence type and return cleaned records."""
    print(f"\n[→] Collecting {seq_type.upper()} sequences for '{plant_term}'...")

    ncbi_db_map = {
        "dna": "nucleotide",
        "protein": "protein",
    }
    ncbi_db = ncbi_db_map.get(seq_type, "nucleotide")

    pipeline_argv = [
        "--ncbi-term", plant_term,
        "--ncbi-db", ncbi_db,
        "--retmax", str(retmax),
        "--out-raw", str(temp_raw),
        "--out-clean", str(temp_clean),
        "--skip-load",
    ]

    pipeline_module.main(pipeline_argv)

    if temp_clean.exists():
        raw = json.loads(temp_clean.read_text(encoding="utf-8"))
        records = raw.get("genes", [])
        print(f"  ✓ Collected {len(records)} {seq_type.upper()} records")
        return records
    return []


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Collect multiple sequence types (DNA, protein) for a plant species."
    )
    parser.add_argument("--plant", required=True, help="Plant species term (e.g., 'Oryza sativa')")
    parser.add_argument("--retmax", type=int, default=300, help="Max records per type (default: 300)")
    parser.add_argument(
        "--types",
        default="dna,protein",
        help="Comma-separated sequence types to collect (default: dna,protein)",
    )
    parser.add_argument("--out-clean", default=None, help="Output clean JSON file")
    parser.add_argument("--out-raw", default=None, help="Output raw JSON file")
    parser.add_argument("--create-tables", action="store_true", help="Create PostgreSQL tables")
    parser.add_argument("--skip-load", action="store_true", help="Skip PostgreSQL import")
    args = parser.parse_args(argv)

    plant_clean_name = args.plant.lower().replace(" ", "_")
    out_clean = Path(args.out_clean) if args.out_clean else (
        ROOT / "data" / "clean" / f"{plant_clean_name}_all_types.json"
    )
    out_raw = Path(args.out_raw) if args.out_raw else (
        ROOT / f"raw_collect_{plant_clean_name}.json"
    )

    out_clean.parent.mkdir(parents=True, exist_ok=True)
    out_raw.parent.mkdir(parents=True, exist_ok=True)

    seq_types = [t.strip().lower() for t in args.types.split(",")]
    print(f"🌱 Collecting {', '.join(seq_types).upper()} for '{args.plant}'...\n")

    all_records: dict[str, dict] = {}
    temp_dir = ROOT / ".temp_collect"
    temp_dir.mkdir(exist_ok=True)

    for seq_type in seq_types:
        temp_raw = temp_dir / f"raw_{seq_type}.json"
        temp_clean = temp_dir / f"clean_{seq_type}.json"

        records = collect_and_clean_type(
            args.plant, seq_type, args.retmax, temp_raw, temp_clean
        )

        for record in records:
            gene_id = record.get("gene_id") or record.get("symbol")
            if gene_id:
                if gene_id not in all_records:
                    all_records[gene_id] = record
                else:
                    # Merge metadata if gene already exists (e.g., from another type)
                    existing = all_records[gene_id]
                    if not existing.get("sequence") and record.get("sequence"):
                        all_records[gene_id]["sequence"] = record["sequence"]
                    if not existing.get("sequence_type") and record.get("sequence_type"):
                        all_records[gene_id]["sequence_type"] = record["sequence_type"]

    print(f"\n✓ Merged {len(all_records)} unique genes\n")

    # Write combined clean JSON
    out_data = {
        "metadata": {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "plant": args.plant,
            "sequence_types": seq_types,
            "count": len(all_records),
        },
        "genes": list(all_records.values()),
    }

    out_clean.write_text(json.dumps(out_data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"✓ Wrote {len(all_records)} records to {out_clean}\n")

    # Optional: Import to PostgreSQL
    if not args.skip_load:
        import load_to_postgres as load_to_postgres

        load_argv = ["--json-file", str(out_clean)]
        if args.create_tables:
            load_argv.insert(0, "--create-tables")

        print("[3/3] Loading to PostgreSQL...\n")
        load_to_postgres.main(load_argv)

    # Cleanup temp directory
    import shutil
    if temp_dir.exists():
        shutil.rmtree(temp_dir)

    print("\n✓ Pipeline complete!")


if __name__ == "__main__":
    main()
