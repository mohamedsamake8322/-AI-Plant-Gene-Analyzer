#!/usr/bin/env python3
"""Run `collect_all_plants.ps1` via PowerShell and capture a timestamped log.

Usage (from repository root):
  python collect/run_collect_with_log.py -Workers 4 -Retmax 300

This creates `logs/collect_run_YYYYMMDD_HHMMSS.log` and returns the PS exit code.
"""
from __future__ import annotations

import argparse
import subprocess
from datetime import datetime
from pathlib import Path
import sys


def build_ps_args(args: argparse.Namespace) -> list[str]:
    ps_args: list[str] = []
    ps_args += ["-Workers", str(args.Workers)]
    ps_args += ["-Retmax", str(args.Retmax)]
    if args.Sources:
        ps_args += ["-Sources", args.Sources]
    if args.Category:
        ps_args += ["-Category", args.Category]
    if args.LoadDB:
        ps_args.append("-LoadDB")
    if args.CreateTables:
        ps_args.append("-CreateTables")
    if args.Force:
        ps_args.append("-Force")
    if args.NoMerge:
        ps_args.append("-NoMerge")
    if args.DryRun:
        ps_args.append("-DryRun")
    return ps_args


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Run collect_all_plants.ps1 and capture a log file")
    p.add_argument("-Workers", type=int, default=4)
    p.add_argument("-Retmax", type=int, default=300)
    p.add_argument("-Sources", type=str, default="ncbi,uniprot,kegg,planttfdb,pubmed")
    p.add_argument("-Category", type=str, default="")
    p.add_argument("-LoadDB", action="store_true")
    p.add_argument("-CreateTables", action="store_true")
    p.add_argument("-Force", action="store_true")
    p.add_argument("-NoMerge", action="store_true")
    p.add_argument("-DryRun", action="store_true")
    p.add_argument("-Sequential", action="store_true", help="Force sequential run (Workers=1)")

    args = p.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[1]
    script_dir = Path(__file__).resolve().parent
    target = script_dir / "collect_all_plants.ps1"
    if not target.exists():
        print(f"Error: target script not found: {target}")
        return 2

    workers = 1 if args.Sequential else args.Workers
    args.Workers = workers

    log_dir = repo_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"collect_run_{ts}.log"

    # Build PowerShell command
    ps_invocation = [
        "powershell",
        "-ExecutionPolicy",
        "Bypass",
        "-NoProfile",
        "-File",
        str(target),
    ]
    ps_invocation += build_ps_args(args)

    print(f"Running: {' '.join(ps_invocation)}")
    print(f"Log file: {log_file}")
    if args.DryRun:
        print("Dry run, not executing")
        return 0

    # Open log file and run
    with open(log_file, "wb") as lf:
        try:
            proc = subprocess.run(ps_invocation, stdout=lf, stderr=subprocess.STDOUT)
            rc = proc.returncode
            if rc == 0:
                print(f"Collection finished successfully. See log: {log_file}")
            else:
                print(f"Collection exited with code {rc}. See log: {log_file}")
            return rc
        except FileNotFoundError:
            print("Error: 'powershell' executable not found in PATH.")
            return 3
        except Exception as e:
            print(f"Fatal error while running collection: {e}")
            return 4


if __name__ == "__main__":
    rc = main()
    sys.exit(rc)
