#!/usr/bin/env python3
"""
Validate that core bioinformatics engines produce alignment-based results.
Run: python scripts/validate_bioinformatics.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import alignment_engine as aln
import bioinformatics as bio
import distance_engine as dist
import phylogeny_engine as phylo
import similarityengine as sim
import numpy as np


def main() -> int:
    report: dict = {"checks": [], "passed": 0, "failed": 0}

    def check(name: str, condition: bool, detail: str = "") -> None:
        status = "PASS" if condition else "FAIL"
        report["checks"].append({"name": name, "status": status, "detail": detail})
        if condition:
            report["passed"] += 1
        else:
            report["failed"] += 1
        print(f"[{status}] {name}" + (f" — {detail}" if detail else ""))

    seq = "ATGCGATCGATCG"
    aln_res = aln.needleman_wunsch(seq, seq)
    check("Needleman-Wunsch self-identity", aln_res["identity_percent"] == 100.0)

    mut = bio.detect_mutations("ATGCGATCA", "ATGCGATCG")
    check("Mutations include alignment block", "alignment" in mut)
    check("Mutations report substitutions", mut["total_mutations"] >= 1)

    msa = aln.star_alignment(["ATGCGATCG", "ATGCGATCA", "ATGCGATCC"])
    check("Star MSA row count", msa["num_sequences"] == 3)

    dm = dist.distance_matrix(
        [
            {"name": "A", "sequence": "ATGCGATCG"},
            {"name": "B", "sequence": "ATGCGATCA"},
        ],
        method="kimura",
    )
    check("Distance matrix uses alignment", bool(dm.get("aligned_sequences")))

    db = {
        "G1": {
            "sequence": seq,
            "trait": "t",
            "description": "d",
            "organism": "o",
            "accession": "a",
        }
    }
    match = sim.compare_with_database("ATGCGATCA", db, top_n=1)[0]
    check("Similarity uses NW algorithm", "Needleman" in str(match["alignment"].get("algorithm", "")))

    tree = phylo.upgma(np.array([[0, 0.1], [0.1, 0]]), ["A", "B"])
    check("UPGMA Newick", tree.get("newick", "").endswith(";"))

    out = ROOT / "logs" / "bioinformatics_validation.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nReport saved to {out}")
    print(f"Passed: {report['passed']} | Failed: {report['failed']}")
    return 0 if report["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
