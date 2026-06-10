#!/usr/bin/env python3
"""
Normalize and clean collected plant data into a canonical gene list.

Usage:
  python scripts/clean_data.py --in raw_collect.json --out data/clean/plant_data_clean.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data" / "clean" / "plant_data_clean.json"
DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)


def normalize_record(item: dict) -> dict | None:
    if not isinstance(item, dict):
        return None

    gene_id = item.get("gene_id") or item.get("accession") or item.get("symbol")
    symbol = item.get("symbol") or item.get("gene_id") or item.get("accession")
    seq = item.get("sequence") or item.get("seq") or ""
    seq = str(seq).replace("\n", "").strip()
    seq_type = item.get("sequence_type") or item.get("type") or None

    if not seq:
        return None

    if not seq_type:
        seq_type = "dna" if all(c in "ACGTN" for c in seq.upper()) else "protein"
    if seq_type not in ("dna", "protein"):
        seq_type = "dna"

    external_links = dict(item.get("external_links") or {})
    if item.get("url"):
        external_links.setdefault("source", item["url"])
    if item.get("source_url"):
        external_links.setdefault("source_api", item["source_url"])
    if item.get("accession"):
        external_links.setdefault("accession", item["accession"])

    rec = {
        "gene_id": gene_id,
        "symbol": symbol,
        "organism": item.get("organism") or item.get("species") or "Unknown",
        "sequence": seq.upper(),
        "sequence_type": seq_type,
        "description": item.get("description") or item.get("title") or "",
        "traits": item.get("traits") or item.get("trait") or [],
        "external_links": external_links,
        "expression_profiles": item.get("expression_profiles") or [],
        "pathways": item.get("pathways") or [],
        "publications": item.get("publications") or [],
        "annotations": item.get("annotations") or {},
        "length": len(seq),
        "date_added": item.get("date_added") or datetime.utcnow().isoformat() + "Z",
        "source": item.get("source") or "collected",
    }
    return rec


def clean(input_path: Path, output_path: Path) -> int:
    raw = json.loads(input_path.read_text(encoding="utf-8"))
    candidates: list[dict] = []
    if isinstance(raw, dict):
        for k in ("genes", "geo", "ensembl", "expression_atlas", "ncbi", "records"):
            if k in raw and isinstance(raw[k], list):
                candidates.extend(raw[k])
        if not candidates:
            if all(isinstance(v, dict) for v in raw.values()):
                candidates.extend(raw.values())
            elif raw.get("gene_id"):
                candidates.append(raw)
    elif isinstance(raw, list):
        candidates.extend(raw)

    cleaned: dict[str, dict] = {}
    skipped = 0
    for item in candidates:
        rec = normalize_record(item)
        if not rec:
            skipped += 1
            continue
        key = rec.get("gene_id") or rec.get("symbol")
        if not key:
            skipped += 1
            continue
        existing = cleaned.get(key)
        if existing:
            if (not existing.get("sequence")) and rec.get("sequence"):
                cleaned[key] = rec
        else:
            cleaned[key] = rec

    # Output as a JSON object with metadata and genes list
    out = {
        "metadata": {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "source_file": str(input_path),
            "count": len(cleaned),
        },
        "genes": list(cleaned.values()),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    if skipped:
        print(f"Skipped {skipped} non-sequence metadata record(s)")
    return len(cleaned)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Clean and normalize collected plant data")
    p.add_argument("--in", dest="input", required=True, help="Input JSON file from collectors")
    p.add_argument("--out", dest="out", default=str(DEFAULT_OUT), help="Output cleaned JSON file")
    args = p.parse_args(argv)

    input_path = Path(args.input)
    out_path = Path(args.out)
    if not input_path.exists():
        print("Input file not found:", input_path)
        raise SystemExit(2)

    n = clean(input_path, out_path)
    print(f"Cleaned {n} gene record(s) -> {out_path}")


if __name__ == "__main__":
    main()
