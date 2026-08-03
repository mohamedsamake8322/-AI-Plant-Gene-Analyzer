#!/usr/bin/env python3
"""
Fetch plant gene and protein data from the Ensembl REST API (no API key required).
Docs: https://rest.ensembl.org
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"
ENSEMBL_BASE = "https://rest.ensembl.org"
DEFAULT_SPECIES = "arabidopsis_thaliana"


def rest_get(path: str, accept: str = "application/json") -> str | dict:
    url = f"{ENSEMBL_BASE}{path}"
    req = urllib.request.Request(url, headers={"Content-Type": accept, "Accept": accept})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            body = resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        detail = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Ensembl HTTP {e.code} for {path}: {detail[:300]}") from e
    if accept == "application/json":
        return json.loads(body)
    return body


def lookup_symbol(species: str, symbol: str) -> dict:
    path = f"/lookup/symbol/{urllib.parse.quote(species)}/{urllib.parse.quote(symbol)}?expand=1"
    return rest_get(f"{path}&content-type=application/json")


def lookup_id(feature_id: str) -> dict:
    path = f"/lookup/id/{urllib.parse.quote(feature_id)}?expand=1"
    return rest_get(f"{path}&content-type=application/json")


def fetch_sequence(feature_id: str, seq_type: str = "cdna") -> str:
    path = f"/sequence/id/{urllib.parse.quote(feature_id)}?type={seq_type}"
    seq = rest_get(f"{path}&content-type=text/plain", accept="text/plain")
    return str(seq).upper().replace("\n", "")


def pick_transcript_id(info: dict) -> str | None:
    transcripts = info.get("Transcript") or []
    if not transcripts:
        return None
    canonical = next((t for t in transcripts if t.get("is_canonical")), transcripts[0])
    return canonical.get("id")


def make_gene_record(info: dict, sequence: str, species: str, sequence_type: str) -> dict:
    symbol = info.get("display_name") or info.get("name") or info.get("id")
    gene_id = info.get("id") or symbol
    desc = info.get("description") or info.get("biotype") or ""
    return {
        "gene_id": gene_id,
        "symbol": symbol,
        "organism": species.replace("_", " ").title(),
        "traits": [],
        "sequence": sequence,
        "sequence_type": sequence_type,
        "description": desc,
        "external_links": {
            "ensembl": f"https://plants.ensembl.org/{species}/Gene/Summary?g={gene_id}",
        },
        "expression_profiles": [],
        "pathways": [],
        "publications": [],
        "source": "Ensembl",
        "source_url": f"https://rest.ensembl.org/lookup/id/{gene_id}",
        "annotations": {
            "biotype": info.get("biotype"),
            "assembly": info.get("assembly_name"),
            "start": info.get("start"),
            "end": info.get("end"),
            "strand": info.get("strand"),
        },
    }


def fetch_gene(species: str, symbol: str | None, feature_id: str | None, seq_type: str) -> dict:
    if feature_id:
        info = lookup_id(feature_id)
    elif symbol:
        info = lookup_symbol(species, symbol)
    else:
        raise ValueError("Provide --symbol or --id")

    transcript_id = pick_transcript_id(info) or info.get("id")
    if not transcript_id:
        raise RuntimeError("No transcript or gene ID returned by Ensembl")

    sequence = fetch_sequence(transcript_id, seq_type=seq_type)
    if not sequence:
        raise RuntimeError(f"Empty sequence for {transcript_id}")
    sequence_type = "protein" if seq_type == "protein" else "dna"
    return make_gene_record(info, sequence, species, sequence_type)


def main(argv: list[str]) -> None:
    p = argparse.ArgumentParser(description="Ensembl Plants gene collector")
    p.add_argument("--symbol", "-s", help="Gene symbol, e.g. DREB1A")
    p.add_argument("--id", help="Ensembl stable ID (e.g. AT1G01010 or ENSG...)")
    p.add_argument("--species", default=DEFAULT_SPECIES, help=f"Ensembl species (default: {DEFAULT_SPECIES})")
    p.add_argument(
        "--seq-type",
        default="cdna",
        choices=["cdna", "cds", "genomic", "protein"],
        help="Sequence type (cdna, cds, genomic, or protein)",
    )
    p.add_argument("--add", action="store_true", help="Insert into genes_database.json")
    p.add_argument("--out", help="Write gene JSON to file")
    p.add_argument("--dbpath", default=str(DEFAULT_DB))
    args = p.parse_args(argv)

    if not args.symbol and not args.id:
        p.error("Provide --symbol or --id")

    try:
        record = fetch_gene(args.species, args.symbol, args.id, args.seq_type)
    except (RuntimeError, ValueError) as e:
        print(f"Failed: {e}")
        sys.exit(1)

    time.sleep(0.1)

    if args.out:
        Path(args.out).write_text(json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Wrote record to {args.out}")

    if args.add:
        sys.path.insert(0, str(ROOT))
        import scripts.validate_and_add_gene as validator

        ok, msg = validator.add_record_to_db(record, Path(args.dbpath))
        print(f"{record['gene_id']}: {msg}")
        if not ok:
            sys.exit(1)
        return

    print(json.dumps(record, indent=2, ensure_ascii=False))
    print(f"\nFetched {record['symbol']} ({len(record['sequence'])} bp). Use --add to insert into DB.")


if __name__ == "__main__":
    main(sys.argv[1:])
