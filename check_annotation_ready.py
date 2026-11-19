#!/usr/bin/env python3
"""Contrôle en lecture seule des preuves d'annotation (DIAMOND + Pfam) des 4 espèces, avant tout chargement en base.
Usage (WSL ou Windows) :  python3 check_annotation_ready.py --data-root /mnt/c/Downloads/IA/pipeline_data
Écrit seulement annotation/_annotation_readiness.json."""
import argparse, gzip, hashlib, json, sys
from pathlib import Path

SPECIES = [("oryza_sativa", "GCF_001433935.1", 42580), ("zea_mays", "GCF_902167145.1", 57345),
           ("nicotiana_tabacum_k326", "GCF_000715075.1", 91379), ("nicotiana_tabacum_tn90", "GCF_000715135.1", 84350)]

def sha256_file(p):
    d = hashlib.sha256()
    with open(p, "rb") as h:
        for c in iter(lambda: h.read(1 << 20), b""):
            d.update(c)
    return d.hexdigest()

def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    rec, acc, parts = {}, None, []
    with op(path, "rt") as h:
        for line in h:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    rec[acc] = "".join(parts).upper()
                acc, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if acc is not None:
        rec[acc] = "".join(parts).upper()
    return rec

def fingerprint(prot):
    return hashlib.sha256("\n".join(sorted(hashlib.md5(s.encode()).hexdigest() for s in prot.values())).encode()).hexdigest()

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data-root", required=True, type=Path)
    root = ap.parse_args().data_root
    annot = root / "annotation"
    report, problems = {}, []
    pf_gz = root / "references" / "pfam" / "38.2" / "Pfam-A.hmm.gz"
    pf_sha = sha256_file(pf_gz) if pf_gz.exists() else None
    print(f"Pfam 38.2 de référence : {'présent' if pf_sha else 'ABSENT (contrôle de l empreinte sautée)'}")
    for key, asm, expected in SPECIES:
        r = report[key] = {}
        found = sorted(root.rglob(f"*{asm}*protein.faa.gz")) or sorted(root.rglob(f"*{asm}*protein.faa"))
        if not found:
            problems.append(f"{key}: FASTA introuvable"); r["fasta"] = None; fp = None
        else:
            prot = read_fasta(found[0]); fp = fingerprint(prot)
            r["fasta_proteins"] = len(prot)
            if len(prot) != expected:
                problems.append(f"{key}: {len(prot)} protéines dans le FASTA, attendu {expected}")
        # --- Pfam ---
        pdir = annot / key / "pfam_hmmscan" / "38.2"
        pe = pdir / "evidence.json"
        if not pe.exists():
            problems.append(f"{key}: Pfam evidence.json absent"); r["pfam"] = None
        else:
            ev = json.loads(pe.read_text())
            tsv = pdir / "domains.tsv"
            chunks = list((pdir / "chunks").glob("chunk_*.domtbl"))
            checks = {"run_completed": ev.get("run_completed") is True,
                      "scanned_equals_expected": ev.get("scanned_proteins") == expected,
                      "fingerprint_matches_fasta": fp is not None and ev.get("sequence_set_fingerprint") == fp,
                      "release_38_2": ev.get("database_release") == "38.2",
                      "pfam_sha_matches_reference": pf_sha is None or ev.get("database_hmm_gz_sha256") == pf_sha,
                      "tsv_present_and_sha_ok": tsv.exists() and ev.get("output_tsv_sha256") == sha256_file(tsv),
                      "chunks_match_evidence": len(chunks) == ev.get("chunks")}
            r["pfam"] = {"checks": checks, "scanned": ev.get("scanned_proteins"), "with_domain": ev.get("proteins_with_domain"),
                         "without_domain": ev.get("proteins_without_domain"), "domain_rows": ev.get("domain_rows"),
                         "hmmer": ev.get("tool_version"), "environment": ev.get("run_environment", "colab")}
            problems += [f"{key}: Pfam {k} NON" for k, v in checks.items() if not v]
        # --- DIAMOND ---
        ddirs = sorted((annot / key / "swissprot_diamond").glob("*/evidence.json"))
        if not ddirs:
            problems.append(f"{key}: DIAMOND evidence.json absent"); r["diamond"] = None
        else:
            de = json.loads(ddirs[-1].read_text())
            checks = {"run_completed": de.get("run_completed") is True,
                      "scanned_equals_expected": de.get("scanned_proteins") == expected,
                      "fingerprint_matches_fasta": fp is not None and de.get("sequence_set_fingerprint") == fp}
            r["diamond"] = {"checks": checks, "release": ddirs[-1].parent.name, "tool_version": de.get("tool_version"),
                            "scanned": de.get("scanned_proteins"), "with_hit": de.get("proteins_with_hit"),
                            "with_transfer_grade_hit": de.get("proteins_with_transfer_grade_hit")}
            if len(ddirs) > 1:
                problems.append(f"{key}: plusieurs versions DIAMOND ({[d.parent.name for d in ddirs]}) : en garder une seule")
            problems += [f"{key}: DIAMOND {k} NON" for k, v in checks.items() if not v]
    print()
    for key, r in report.items():
        print(f"=== {key} (FASTA {r.get('fasta_proteins')}) ===")
        for tool in ("pfam", "diamond"):
            t = r.get(tool)
            if not t:
                print(f"  {tool}: ABSENT"); continue
            ok = all(t["checks"].values())
            print(f"  {tool}: {'OK' if ok else 'PROBLEME'} | " + ", ".join(f"{k}={v}" for k, v in t.items() if k != "checks"))
            for k, v in t["checks"].items():
                if not v:
                    print(f"      NON : {k}")
    (annot / "_annotation_readiness.json").write_text(json.dumps({"report": report, "problems": problems}, indent=2), encoding="utf-8")
    print("\nVERDICT :", "PRÊT pour le chargement en base" if not problems else f"{len(problems)} problème(s) :")
    for p in problems:
        print("  -", p)
    sys.exit(1 if problems else 0)

main()
