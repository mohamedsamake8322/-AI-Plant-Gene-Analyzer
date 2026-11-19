#!/usr/bin/env python
"""Diagnostic LECTURE SEULE : DIAMOND et Pfam (HMMER) sont-ils branches ?

Verifie :
  1. les executables (diamond, hmmscan, hmmpress) : PATH, dossiers du projet, WSL ;
  2. les bases (Swiss-Prot .dmnd, Pfam-A.hmm + index .h3*) ;
  3. les preuves de scan enregistrees (versions, commandes, bases) ;
  4. le code Python (hors .venv) qui appelle ces outils.

Usage (depuis C:\\Downloads\\IA, venv active) :
    python scripts\\diagnose_annotation_tools.py
    python scripts\\diagnose_annotation_tools.py --root C:\\Downloads\\IA --out results\\diag_outils.json

Ne modifie rien : ni base, ni fichiers (sauf le rapport JSON demande).
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

TOOLS = {
    "diamond": ["diamond", "diamond.exe"],
    "hmmscan": ["hmmscan", "hmmscan.exe"],
    "hmmpress": ["hmmpress", "hmmpress.exe"],
}
VERSION_ARGS = {"diamond": ["version"], "hmmscan": ["-h"], "hmmpress": ["-h"]}
SKIP_DIRS = {".venv", "venv", "__pycache__", ".git", "node_modules", "backups"}
CODE_PATTERN = re.compile(r"diamond|hmmscan|hmmsearch|hmmpress|pfam_scan", re.I)
EVIDENCE_KEYS = re.compile(
    r"version|database|db_|command|tool|sha256|release|e_?value|evalue|date|"
    r"species|assembly|sensitiv|threshold",
    re.I,
)


def run(cmd, timeout=20):
    try:
        out = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, shell=False
        )
        text = (out.stdout or "") + (out.stderr or "")
        return out.returncode, text.strip()
    except FileNotFoundError:
        return None, "introuvable"
    except subprocess.TimeoutExpired:
        return None, "timeout"
    except Exception as exc:  # noqa: BLE001
        return None, f"erreur: {exc}"


def walk(root: Path):
    for p in root.rglob("*"):
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        yield p


def check_executables(root: Path):
    report = {}
    local_dirs = [d for d in [root / "tools", root / "bin", root / "pipeline_data" / "tools"] if d.exists()]
    for name, candidates in TOOLS.items():
        entry = {"path_env": None, "project_dirs": [], "version": None, "wsl": None}
        for c in candidates:
            found = shutil.which(c)
            if found:
                entry["path_env"] = found
                break
        for d in local_dirs:
            for c in candidates:
                for hit in d.rglob(c):
                    entry["project_dirs"].append(str(hit))
        exe = entry["path_env"] or (entry["project_dirs"][0] if entry["project_dirs"] else None)
        if exe:
            rc, txt = run([exe] + VERSION_ARGS[name])
            first = next((ln for ln in txt.splitlines() if ln.strip()), "")
            entry["version"] = first[:200]
        report[name] = entry

    # HMMER n'existe pas en natif sur Windows : verifier WSL.
    if shutil.which("wsl"):
        for name in ("hmmscan", "diamond"):
            rc, txt = run(["wsl", "which", name], timeout=40)
            report[name]["wsl"] = txt if rc == 0 else "absent dans WSL"
    return report


def check_databases(root: Path):
    patterns = {
        "swissprot_dmnd": ["*.dmnd"],
        "swissprot_fasta": ["*swiss*prot*.fa*", "uniprot_sprot*"],
        "pfam_hmm": ["Pfam-A.hmm", "Pfam-A*.hmm*"],
        "pfam_index": ["*.h3m", "*.h3i", "*.h3f", "*.h3p"],
    }
    found = {k: [] for k in patterns}
    for p in walk(root):
        if not p.is_file():
            continue
        for key, pats in patterns.items():
            if any(p.match(pat) for pat in pats):
                st = p.stat()
                found[key].append(
                    {
                        "path": str(p),
                        "size_mb": round(st.st_size / 1e6, 1),
                        "modified": datetime.fromtimestamp(st.st_mtime).isoformat(timespec="minutes"),
                    }
                )
    return found


def check_evidence(root: Path):
    out = []
    for p in walk(root):
        if p.suffix.lower() != ".json" or not p.is_file():
            continue
        if p.stat().st_size > 5_000_000:
            continue
        try:
            raw = p.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if not re.search(r"hmmscan|diamond", raw, re.I):
            continue
        try:
            data = json.loads(raw)
        except json.JSONDecodeError:
            continue
        picked = {}

        def collect(obj, prefix=""):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    key = f"{prefix}{k}"
                    if isinstance(v, (dict, list)):
                        collect(v, key + ".")
                    elif EVIDENCE_KEYS.search(str(k)) and len(picked) < 25:
                        picked[key] = str(v)[:160]
            elif isinstance(obj, list):
                for i, v in enumerate(obj[:3]):
                    collect(v, f"{prefix}{i}.")

        collect(data)
        out.append({"file": str(p), "extract": picked})
        if len(out) >= 12:
            break
    return out


def check_code(root: Path):
    hits = {}
    for p in walk(root):
        if p.suffix.lower() not in {".py", ".ipynb"} or not p.is_file():
            continue
        if p.stat().st_size > 2_000_000:
            continue
        try:
            lines = p.read_text(encoding="utf-8", errors="ignore").splitlines()
        except OSError:
            continue
        matches = [
            (i + 1, ln.strip()[:140])
            for i, ln in enumerate(lines)
            if CODE_PATTERN.search(ln)
        ]
        if matches:
            uses_subprocess = any("subprocess" in ln for ln in lines)
            hits[str(p.relative_to(root))] = {
                "n_matches": len(matches),
                "uses_subprocess": uses_subprocess,
                "examples": matches[:4],
            }
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    root = Path(args.root).resolve()

    report = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "root": str(root),
        "python": sys.version.split()[0],
        "executables": check_executables(root),
        "databases": check_databases(root),
        "scan_evidence": check_evidence(root / "pipeline_data") if (root / "pipeline_data").exists() else [],
        "code_references": check_code(root),
    }

    ex = report["executables"]
    summary = {
        "diamond_available": bool(ex["diamond"]["path_env"] or ex["diamond"]["project_dirs"] or (ex["diamond"]["wsl"] not in (None, "absent dans WSL"))),
        "hmmscan_available": bool(ex["hmmscan"]["path_env"] or ex["hmmscan"]["project_dirs"] or (ex["hmmscan"]["wsl"] not in (None, "absent dans WSL"))),
        "swissprot_dmnd_found": bool(report["databases"]["swissprot_dmnd"]),
        "pfam_hmm_found": bool(report["databases"]["pfam_hmm"]),
        "pfam_index_found": bool(report["databases"]["pfam_index"]),
        "files_calling_tools": sorted(report["code_references"]),
    }
    report["summary"] = summary

    text = json.dumps(report, indent=2, ensure_ascii=False)
    print(text)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
        print(f"\nRapport ecrit : {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
