#!/usr/bin/env python3
"""Annotate protein sequences in species JSON files with Pfam-A domains."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

if __package__:
    from .protein_json_annotations import annotate_species_file, find_species_files, output_path_for
else:
    from protein_json_annotations import annotate_species_file, find_species_files, output_path_for

REFERENCE_DIR = ROOT / "data" / "reference"
DEFAULT_HMM = REFERENCE_DIR / "Pfam-A.hmm"
INDEX_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")


def _tool_command(name: str, args: list[str]) -> list[str]:
    if shutil.which(name):
        return [name, *args]
    if os.name == "nt" and shutil.which("wsl.exe"):
        return ["wsl.exe", "--exec", name, *(_wsl_path(arg) for arg in args)]
    raise RuntimeError(f"{name} is unavailable; install it or run this script from WSL.")


_WSL_UNC_RE = re.compile(r"^\\\\wsl(?:\$|\.localhost)\\[^\\]+\\(.*)$")


def _wsl_path(value: str) -> str:
    # UNC path into WSL's own native filesystem, e.g.
    # \\wsl$\Ubuntu\home\user\reference\Pfam-A.hmm -> /home/user/reference/Pfam-A.hmm
    # Native paths here avoid the slow Windows<->WSL DrvFs bridge entirely.
    unc_match = _WSL_UNC_RE.match(value)
    if unc_match:
        return "/" + unc_match.group(1).replace("\\", "/")
    if len(value) < 3 or value[1:3] != ":\\":
        return value
    return f"/mnt/{value[0].lower()}/{value[3:].replace(chr(92), '/') }"


def ensure_hmm_database(database: Path) -> list[Path]:
    database = database.resolve()
    archive = database.with_suffix(database.suffix + ".gz")
    if not database.is_file():
        if not archive.is_file():
            raise FileNotFoundError(f"Pfam HMM file not found: {database} or {archive}")
        temporary = database.with_name(f".{database.name}.tmp")
        try:
            with gzip.open(archive, "rb") as source, temporary.open("wb") as target:
                shutil.copyfileobj(source, target, length=1024 * 1024)
            temporary.replace(database)
        except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(
                f"Cannot decompress {archive}: the gzip archive is invalid or incomplete. "
                "Re-download it and verify with `gzip -t` before retrying."
            ) from exc

    index_files = hmm_index_files(database)
    if not all(path.is_file() and path.stat().st_mtime >= database.stat().st_mtime for path in index_files):
        subprocess.run(_tool_command("hmmpress", ["-f", str(database)]), check=True)
    return [database, *index_files]


def hmm_index_files(database: Path) -> list[Path]:
    return [Path(f"{database}{suffix}") for suffix in INDEX_SUFFIXES]


def parse_domains(path: Path) -> dict[str, list[dict]]:
    domains: dict[str, list[dict]] = {}
    if not path.exists():
        return domains
    with path.open(encoding="utf-8") as output:
        for line in output:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                continue
            domains.setdefault(fields[3], []).append({
                "accession": fields[1] if fields[1] != "-" else fields[0],
                "name": fields[0],
                "independent_evalue": float(fields[12]),
                "score": float(fields[13]),
                "hmm_start": int(fields[15]),
                "hmm_end": int(fields[16]),
                "start": int(fields[17]),
                "end": int(fields[18]),
                "description": fields[22] if len(fields) > 22 else "",
            })
    return domains


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Species JSON file or directory containing *_all_sources.json files.")
    parser.add_argument("--output-dir", type=Path,
                        help="Write annotated copies here; default updates each input atomically.")
    parser.add_argument("--reference-hmm", type=Path, default=DEFAULT_HMM)
    parser.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--limit", type=int,
                        help="Process at most N unique protein sequences per species (partial run).")
    args = parser.parse_args()

    if args.threads < 1 or args.batch_size < 1:
        parser.error("threads and batch-size must be positive")
    if args.limit is not None and args.limit < 1:
        parser.error("limit must be positive")
    reference_files = ensure_hmm_database(args.reference_hmm)

    def run_batch(query: Path, output: Path, threads: int) -> None:
        output.unlink(missing_ok=True)
        command = [
            "--cut_ga", "--noali", "--cpu", str(threads),
            "--domtblout", str(output), str(reference_files[0]), str(query),
        ]
        subprocess.run(_tool_command("hmmscan", command), check=True)

    files = find_species_files(args.input)
    if not files:
        parser.error(f"No *_all_sources.json files found in {args.input}")
    for species_file in files:
        destination = output_path_for(species_file, args.output_dir)
        annotate_species_file(
            input_path=species_file,
            output_path=destination,
            tool_name="pfam",
            annotation_key="pfam_domains",
            reference_files=reference_files,
            options={"cutoff": "gathering"},
            batch_size=args.batch_size,
            threads=args.threads,
            limit=args.limit,
            run_batch=run_batch,
            parse_output=parse_domains,
        )


if __name__ == "__main__":
    main()