#!/usr/bin/env python3
"""
Stream rebuild of master_plant_db.json to avoid loading all records in memory.

Writes a temporary JSONL file of gene records, collects simple stats,
then writes the final `master_plant_db.json` by streaming the JSONL lines
into a JSON array. Removes the temp file afterwards.

Usage: python scripts/rebuild_master_stream.py
"""
from pathlib import Path
import json
import argparse


def stream_rebuild(species_dir: Path, out_path: Path, tmp_path: Path):
    total = 0
    organisms = set()
    seq_type_counts = {}
    source_counts = {}

    tmp_path.parent.mkdir(parents=True, exist_ok=True)
    with tmp_path.open("w", encoding="utf-8") as tmpf:
        for p in sorted(species_dir.glob("*_all_sources.json")):
            try:
                if p.stat().st_size == 0:
                    print(f"skip empty: {p.name}")
                    continue
                data = json.loads(p.read_text(encoding="utf-8"))
                for g in data.get("genes", []):
                    # write single line JSON record
                    tmpf.write(json.dumps(g, ensure_ascii=False))
                    tmpf.write("\n")
                    total += 1
                    org = g.get("organism", "Unknown")
                    organisms.add(org)
                    st = g.get("sequence_type", "unknown")
                    seq_type_counts[st] = seq_type_counts.get(st, 0) + 1
                    src = g.get("source", "unknown")
                    source_counts[src] = source_counts.get(src, 0) + 1
            except Exception as e:
                print(f"error reading {p}: {e}")

    # build metadata
    metadata = {
        "generated_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        "description": "Master plant genomics database — rebuilt from species files (stream)",
        "total_genes": total,
        "total_species": len(organisms),
        "organism_counts": {o: None for o in organisms},
        "sequence_type_counts": seq_type_counts,
        "source_counts": source_counts,
    }

    # write final JSON by streaming lines from tmp
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as outf:
        outf.write(json.dumps({"metadata": metadata}, ensure_ascii=False, indent=2)[:-2])
        outf.write(",\n  " + '"genes": [\n')
        first = True
        with tmp_path.open("r", encoding="utf-8") as tmpf:
            for line in tmpf:
                line = line.rstrip("\n")
                if not line:
                    continue
                if first:
                    outf.write("    ")
                    outf.write(line)
                    first = False
                else:
                    outf.write(",\n    ")
                    outf.write(line)
        outf.write("\n  ]\n}\n")

    # remove tmp
    try:
        tmp_path.unlink()
    except Exception:
        pass


def main():
    parser = argparse.ArgumentParser(description="Stream rebuild master JSON from species files")
    parser.add_argument("--species-dir", default="data/clean/species")
    parser.add_argument("--out", default="data/clean/master_plant_db.json")
    parser.add_argument("--tmp", default="data/clean/.master_tmp.jsonl")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    sp_dir = root / args.species_dir
    out_path = root / args.out
    tmp_path = root / args.tmp

    if not sp_dir.exists():
        raise SystemExit(f"Species dir not found: {sp_dir}")

    print("stream rebuilding from", sp_dir)
    stream_rebuild(sp_dir, out_path, tmp_path)
    print("wrote", out_path)


if __name__ == "__main__":
    main()
