#!/usr/bin/env python3
"""Annotate protein sequences in species JSON files with Swiss-Prot DIAMOND."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

if __package__:
    from .protein_json_annotations import annotate_species_file, find_species_files, output_path_for
else:
    from protein_json_annotations import annotate_species_file, find_species_files, output_path_for

REFERENCE_DIR = ROOT / "data" / "reference"
DEFAULT_DB = REFERENCE_DIR / "swissprot_diamond.dmnd"


def _tool_command(name: str, args: list[str]) -> list[str]:
    if shutil.which(name):
        return [name, *args]
    if os.name == "nt" and shutil.which("wsl.exe"):
        return ["wsl.exe", "--exec", name, *(_wsl_path(arg) for arg in args)]
    raise RuntimeError(f"{name} is unavailable; install it or run this script from WSL.")


def _wsl_path(value: str) -> str:
    path = Path(value)
    if not path.is_absolute() or len(value) < 3 or value[1:3] != ":\\":
        return value
    return f"/mnt/{value[0].lower()}/{value[3:].replace(chr(92), '/') }"


def parse_hits(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = {}
    if not path.exists():
        return hits
    with path.open(encoding="utf-8") as output:
        for line in output:
            fields = line.rstrip("\n").split("\t", 10)
            if len(fields) != 11:
                continue
            target_id = fields[1]
            target_parts = target_id.split("|")
            accession = target_parts[1] if len(target_parts) > 1 else target_id
            hits.setdefault(fields[0], []).append({
                "accession": accession,
                "target_id": target_id,
                "identity": float(fields[2]),
                "alignment_length": int(fields[3]),
                "query_start": int(fields[4]),
                "query_end": int(fields[5]),
                "target_start": int(fields[6]),
                "target_end": int(fields[7]),
                "evalue": float(fields[8]),
                "bit_score": float(fields[9]),
                "description": fields[10],
            })
    return hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Species JSON file or directory containing *_all_sources.json files.")
    parser.add_argument("--output-dir", type=Path,
                        help="Write annotated copies here; default updates each input atomically.")
    parser.add_argument("--reference-db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--limit", type=int,
                        help="Process at most N unique protein sequences per species (partial run).")
    parser.add_argument("--max-targets", type=int, default=5)
    parser.add_argument("--evalue", type=float, default=1e-5)
    args = parser.parse_args()

    if args.threads < 1 or args.batch_size < 1 or args.max_targets < 1:
        parser.error("threads, batch-size, and max-targets must be positive")
    if args.limit is not None and args.limit < 1:
        parser.error("limit must be positive")
    database = args.reference_db.resolve()
    if not database.is_file():
        raise FileNotFoundError(f"DIAMOND database not found: {database}")

    def run_batch(query: Path, output: Path, threads: int) -> None:
        output.unlink(missing_ok=True)
        command = [
            "blastp", "--query", str(query), "--db", str(database),
            "--out", str(output), "--outfmt", "6",
            "qseqid", "sseqid", "pident", "length", "qstart", "qend",
            "sstart", "send", "evalue", "bitscore", "stitle",
            "--max-target-seqs", str(args.max_targets), "--evalue", str(args.evalue),
            "--threads", str(threads), "--sensitive",
        ]
        subprocess.run(_tool_command("diamond", command), check=True)

    files = find_species_files(args.input)
    if not files:
        parser.error(f"No *_all_sources.json files found in {args.input}")
    for species_file in files:
        destination = output_path_for(species_file, args.output_dir)
        annotate_species_file(
            input_path=species_file,
            output_path=destination,
            tool_name="diamond",
            annotation_key="diamond_hits",
            reference_files=[database],
            options={"evalue": args.evalue, "max_targets": args.max_targets},
            batch_size=args.batch_size,
            threads=args.threads,
            limit=args.limit,
            run_batch=run_batch,
            parse_output=parse_hits,
        )


if __name__ == "__main__":
    main()