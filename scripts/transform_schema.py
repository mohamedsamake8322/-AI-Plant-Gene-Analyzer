#!/usr/bin/env python3
"""
transform_schema.py
-------------------
Transform legacy JSON schema to professional schema.
Updates all existing JSON files and integrates with app.

Usage:
  # Transform a single file
  python scripts/transform_schema.py --json data/clean/rice_genes.json

  # Transform all files in a directory
  python scripts/transform_schema.py --dir data/clean

  # Transform and save to new location
  python scripts/transform_schema.py --json rice_genes.json --out rice_genes_pro.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import professional_schema as pschema


def transform_json_file(
    json_path: Path,
    out_path: Optional[Path] = None,
    backup: bool = True,
) -> Path:
    """Transform a single JSON file to professional schema."""
    
    if not json_path.exists():
        raise FileNotFoundError(f"JSON file not found: {json_path}")
    
    print(f"📖 Reading: {json_path}")
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    # Extract genes from various formats
    genes = []
    if isinstance(data, dict):
        if "genes" in data:
            genes = data.get("genes", [])
        elif all(isinstance(v, dict) for v in data.values()):
            genes = list(data.values())
        else:
            genes = [data]
    elif isinstance(data, list):
        genes = data
    
    print(f"🔄 Transforming {len(genes)} records...")
    
    # Transform to professional schema
    pro_records = []
    organisms = set()
    seq_types = set()
    
    for legacy_gene in genes:
        if not isinstance(legacy_gene, dict):
            continue
        
        try:
            pro_record = pschema.transform_legacy_gene(legacy_gene)
            pro_records.append(pro_record)
            
            if pro_record.core.organism:
                organisms.add(pro_record.core.organism.scientific_name)
            for seq in pro_record.sequences:
                seq_types.add(seq.type)
                
        except Exception as e:
            print(f"  ⚠️  Failed to transform {legacy_gene.get('gene_id', '?')}: {e}")
    
    # Create output structure
    schema_meta = pschema.create_schema_metadata(
        total_records=len(pro_records),
        unique_organisms=len(organisms),
        sequence_types=sorted(list(seq_types)),
        validated_count=len(pro_records),
    )
    
    output_data = {
        **schema_meta,
        "genes": [r.to_dict() for r in pro_records],
    }
    
    # Save to file
    if out_path is None:
        out_path = json_path.with_stem(json_path.stem + "_pro")
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Backup original
    if backup and json_path != out_path and json_path.exists():
        backup_path = json_path.with_suffix(".bak")
        if not backup_path.exists():
            import shutil
            shutil.copy(json_path, backup_path)
            print(f"💾 Backed up to: {backup_path}")
    
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output_data, f, ensure_ascii=False, indent=2)
    
    print(f"✅ Saved {len(pro_records)} professional records to: {out_path}")
    print(f"   - Organisms: {len(organisms)}")
    print(f"   - Sequence types: {', '.join(sorted(seq_types))}")
    
    return out_path


def transform_directory(
    dir_path: Path,
    pattern: str = "*.json",
    out_suffix: str = "_pro",
) -> list[Path]:
    """Transform all JSON files in a directory."""
    
    json_files = list(dir_path.glob(pattern))
    if not json_files:
        print(f"⚠️  No JSON files found in: {dir_path}")
        return []
    
    print(f"📂 Found {len(json_files)} JSON files in {dir_path}\n")
    
    transformed = []
    for json_file in json_files:
        out_file = json_file.with_stem(json_file.stem + out_suffix)
        try:
            result = transform_json_file(json_file, out_file, backup=True)
            transformed.append(result)
            print()
        except Exception as e:
            print(f"❌ Error transforming {json_file}: {e}\n")
    
    return transformed


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Transform legacy JSON schema to professional schema",
    )
    parser.add_argument(
        "--json",
        type=Path,
        help="Path to JSON file to transform",
    )
    parser.add_argument(
        "--dir",
        type=Path,
        help="Path to directory containing JSON files",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Output file path (if --json specified)",
    )
    parser.add_argument(
        "--pattern",
        default="*.json",
        help="Glob pattern for files in directory (default: *.json)",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Don't create backup files",
    )
    parser.add_argument(
        "--suffix",
        default="_pro",
        help="Suffix for output files (default: _pro)",
    )
    
    args = parser.parse_args(argv)
    
    if args.json:
        transform_json_file(args.json, args.out, backup=not args.no_backup)
    elif args.dir:
        transform_directory(args.dir, args.pattern, args.suffix)
    else:
        print("❌ Specify either --json or --dir")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
