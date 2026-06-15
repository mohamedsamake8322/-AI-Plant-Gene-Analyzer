#!/usr/bin/env python3
"""
Rebuild master_plant_db.json from per-species *_all_sources.json files.

Usage:
  python scripts/rebuild_master.py
  python scripts/rebuild_master.py --out data/clean/master_plant_db.json

This script deduplicates by gene_id + organism and writes a nicely
formatted JSON master file.
"""
from pathlib import Path
import json
import argparse


def main():
    parser = argparse.ArgumentParser(description="Rebuild master_plant_db.json from species files")
    parser.add_argument("--species-dir", default="data/clean/species", help="Directory with per-species JSON files")
    parser.add_argument("--out", default="data/clean/master_plant_db.json", help="Output master JSON file")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    sp_dir = root / args.species_dir
    out_path = root / args.out

    if not sp_dir.exists():
        raise SystemExit(f"Species directory not found: {sp_dir}")

    master_genes = []
    seen = set()

    files = sorted(sp_dir.glob("*_all_sources.json"))
    for p in files:
        try:
            if p.stat().st_size == 0:
                print(f"skip empty: {p.name}")
                continue
            data = json.loads(p.read_text(encoding="utf-8"))
            genes = data.get("genes", [])
            for g in genes:
                gid = g.get("gene_id") or g.get("symbol") or ""
                org = g.get("organism", "")
                key = f"{gid}::{org}"
                if not gid:
                    # keep records with no gene_id but avoid duplicates by full record
                    key = json.dumps(g, sort_keys=True)
                if key in seen:
                    continue
                seen.add(key)
                master_genes.append(g)
        except Exception as e:
            print(f"error reading {p}: {e}")

    master = {
        "metadata": {
            "generated_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
            "description": "Master plant genomics database — rebuilt from species files",
            "total_genes": len(master_genes),
            "total_species": len({g.get("organism", "Unknown") for g in master_genes}),
        },
        "genes": master_genes,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(master, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"WROTE {out_path} ({len(master_genes)} genes)")


if __name__ == "__main__":
    main()
