#!/usr/bin/env python
"""Test d'equivalence : scans DIAMOND + Pfam EN DIRECT (WSL) contre les resultats GELES de K326.

Prend N proteines de K326 (graine fixe), les scanne dans WSL avec les memes outils et
parametres que le calcul par lots, puis compare aux fichiers geles :
    pipeline_data\\annotation\\nicotiana_tabacum_k326\\swissprot_diamond\\2026_03\\hits.tsv
    pipeline_data\\annotation\\nicotiana_tabacum_k326\\pfam_hmmscan\\38.2\\domains.tsv

Lecture seule cote projet (rien n'est ecrit en base ni dans pipeline_data). Les fichiers de travail
sont dans /tmp/live_eq de WSL ; seul le rapport JSON est ecrit (--out).

Usage (depuis C:\\Downloads\\IA, venv active) :
    python scripts\\test_live_vs_frozen.py --out results\\test_live_vs_frozen_k326.json
    python scripts\\test_live_vs_frozen.py --n 20 --seed 17 --skip-pfam      # DIAMOND seul
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
import subprocess
import sys
import time
from pathlib import Path

SPECIES = "nicotiana_tabacum_k326"
FAA = Path("pipeline_data/ncbi/nicotiana_tabacum_k326/GCF_000715075.1_ASM71507v2_protein.faa.gz")
ANN = Path("pipeline_data/annotation") / SPECIES
HITS = ANN / "swissprot_diamond" / "2026_03" / "hits.tsv"
DIAMOND_EVID = ANN / "swissprot_diamond" / "2026_03" / "evidence.json"
DOMAINS = ANN / "pfam_hmmscan" / "38.2" / "domains.tsv"

OUTFMT = ("qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore "
          "qcovhsp scovhsp stitle")
HIT_COLS = OUTFMT.split()
WORK = "/tmp/live_eq"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def read_fasta_gz(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, parts = None, []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name:
        seqs[name] = "".join(parts)
    return seqs


def load_frozen_hits(path: Path, ids: set[str]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {i: [] for i in ids}
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t", len(HIT_COLS) - 1)
            if p[0] in out:
                out[p[0]].append(dict(zip(HIT_COLS, p)))
    return out


def load_frozen_domains(path: Path, ids: set[str]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {i: [] for i in ids}
    with open(path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["protein_accession"] in out:
                out[row["protein_accession"]].append(row)
    return out


def frozen_protein_sets(hits_path: Path, dom_path: Path) -> tuple[set[str], set[str]]:
    with_hit, with_dom = set(), set()
    with open(hits_path, encoding="utf-8") as f:
        for line in f:
            with_hit.add(line.split("\t", 1)[0])
    with open(dom_path, encoding="utf-8") as f:
        next(f)
        for line in f:
            with_dom.add(line.split("\t", 1)[0])
    return with_hit, with_dom


def wsl(script: str, stdin: bytes | None = None, timeout: int = 3600) -> tuple[int, str, str]:
    r = subprocess.run(["wsl.exe", "--exec", "bash", "-c", script], input=stdin,
                       capture_output=True, timeout=timeout)
    return r.returncode, r.stdout.decode("utf-8", "replace"), r.stderr.decode("utf-8", "replace")


def close(a: float, b: float, rel=1e-3, abs_=1e-300) -> bool:
    return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))


def compare_diamond(frozen: list[dict], live: list[dict]) -> dict:
    def key(r):
        return (r["sseqid"], int(r["qstart"]), int(r["sstart"]))

    fz, lv = {key(r): r for r in frozen}, {key(r): r for r in live}
    common = set(fz) & set(lv)
    only_f, only_l = set(fz) - set(lv), set(lv) - set(fz)
    field_diffs = []
    for k in common:
        a, b = fz[k], lv[k]
        for c in ("length", "qend", "send"):
            if int(a[c]) != int(b[c]):
                field_diffs.append((k[0], c, a[c], b[c]))
        if abs(float(a["pident"]) - float(b["pident"])) > 0.05:
            field_diffs.append((k[0], "pident", a["pident"], b["pident"]))
        if abs(float(a["bitscore"]) - float(b["bitscore"])) > 0.5:
            field_diffs.append((k[0], "bitscore", a["bitscore"], b["bitscore"]))
        if not close(float(a["evalue"]), float(b["evalue"])):
            field_diffs.append((k[0], "evalue", a["evalue"], b["evalue"]))
    # Difference acceptable : egalite de score a la limite max-target-seqs (5e rang)
    tie_only = False
    if (only_f or only_l) and fz and lv:
        scores = [float(fz[k]["bitscore"]) for k in only_f] + [float(lv[k]["bitscore"]) for k in only_l]
        min_common = min((float(fz[k]["bitscore"]) for k in common), default=None)
        tie_only = (min_common is not None and len(fz) >= 5 and len(lv) >= 5
                    and all(abs(s - min_common) <= 0.5 for s in scores))
    return {
        "n_frozen": len(fz), "n_live": len(lv), "n_common": len(common),
        "only_frozen": sorted(k[0] for k in only_f), "only_live": sorted(k[0] for k in only_l),
        "field_diffs": field_diffs[:10],
        "top_hit_same": bool(frozen and live and frozen[0]["sseqid"] == live[0]["sseqid"]) or (not frozen and not live),
        "exact": not (only_f or only_l or field_diffs),
        "tie_at_limit_only": bool(tie_only and not field_diffs),
    }


def parse_domtbl(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split(None, 22)
        if len(p) < 21:
            continue
        rows.append({
            "protein_accession": p[3], "pfam_accession": p[1], "pfam_name": p[0],
            "i_evalue": p[12], "dom_score": p[13], "hmm_from": p[15], "hmm_to": p[16],
            "ali_from": p[17], "ali_to": p[18], "env_from": p[19], "env_to": p[20],
        })
    return rows


def compare_pfam(frozen: list[dict], live: list[dict]) -> dict:
    def key(r):
        return (r["pfam_accession"], int(r["ali_from"]), int(r["ali_to"]),
                int(r["env_from"]), int(r["env_to"]))

    fz, lv = {key(r): r for r in frozen}, {key(r): r for r in live}
    common = set(fz) & set(lv)
    diffs = []
    for k in common:
        a, b = fz[k], lv[k]
        if abs(float(a["dom_score"]) - float(b["dom_score"])) > 0.15:
            diffs.append((k[0], "dom_score", a["dom_score"], b["dom_score"]))
        if not close(float(a["i_evalue"]), float(b["i_evalue"]), rel=0.02):
            diffs.append((k[0], "i_evalue", a["i_evalue"], b["i_evalue"]))
    return {
        "n_frozen": len(fz), "n_live": len(lv), "n_common": len(common),
        "only_frozen": sorted(str(k) for k in set(fz) - set(lv))[:5],
        "only_live": sorted(str(k) for k in set(lv) - set(fz))[:5],
        "field_diffs": diffs[:10],
        "exact": not (set(fz) ^ set(lv)) and not diffs,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--skip-pfam", action="store_true")
    ap.add_argument("--skip-diamond", action="store_true")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    root = Path(args.root).resolve()
    t0 = time.time()
    checks: list[tuple[str, str, str]] = []

    def rec(name, status, detail=""):
        checks.append((name, status, detail))
        print(f"[{status:4}] {name}" + (f" : {detail}" if detail else ""))

    # 1. Pre-requis : meme fichier de requetes que le calcul gele
    faa = root / FAA
    ev = json.loads((root / DIAMOND_EVID).read_text(encoding="utf-8"))
    h = sha256(faa)
    rec("FASTA K326 = fichier du calcul gele", "PASS" if h == ev["query_file_sha256"] else "FAIL", h[:16])
    if h != ev["query_file_sha256"]:
        sys.exit(1)
    rc, out, err = wsl("command -v diamond hmmscan; diamond version; "
                       "ls -l $HOME/reference/swissprot_diamond.dmnd $HOME/reference/Pfam-A.hmm.h3m 2>&1")
    print(out.strip())
    ok_tools = rc == 0 and "diamond version 2.2.8" in out
    rec("WSL : diamond 2.2.8 et bases presentes", "PASS" if ok_tools else "FAIL", err.strip()[:120])
    if not ok_tools:
        sys.exit(1)

    # 2. Echantillon stratifie
    seqs = read_fasta_gz(faa)
    with_hit, with_dom = frozen_protein_sets(root / HITS, root / DOMAINS)
    allp = sorted(seqs)
    both = sorted(with_hit & with_dom)
    only_dom = sorted(with_dom - with_hit)
    only_hit = sorted(with_hit - with_dom)
    neither = sorted(set(allp) - with_hit - with_dom)
    rng = random.Random(args.seed)
    k = max(1, args.n // 5)
    pick = []
    for pool, m in ((both, args.n - 3 * k), (only_dom, k), (only_hit, k), (neither, k)):
        pick += rng.sample(pool, min(m, len(pool)))
    for forced in ("NP_001311551.2",):
        if forced in seqs and forced not in pick:
            pick[0] = forced
    pick = list(dict.fromkeys(pick))
    longs = sorted([p for p in allp if len(seqs[p]) > 1022])
    if longs:
        pick.append(rng.choice(longs))  # une proteine > 1022 aa
    pick = list(dict.fromkeys(pick))
    fasta = "".join(f">{p}\n{seqs[p]}\n" for p in pick)
    rec("Echantillon", "INFO", f"{len(pick)} proteines (avec hit+domaine, domaine seul, hit seul, aucun, 1 longue)")

    fz_hits = load_frozen_hits(root / HITS, set(pick))
    fz_dom = load_frozen_domains(root / DOMAINS, set(pick))

    # 3. Envoi du FASTA dans WSL (par stdin : pas de /mnt/c)
    rc, out, err = wsl(f"rm -rf {WORK} && mkdir -p {WORK} && cat > {WORK}/q.faa && sync", stdin=fasta.encode())
    if rc != 0:
        rec("Envoi du FASTA dans WSL", "FAIL", err[:150])
        sys.exit(1)

    report = {"species": SPECIES, "n": len(pick), "seed": args.seed, "proteins": pick,
              "diamond": {}, "pfam": {}}

    # 4. DIAMOND en direct
    if not args.skip_diamond:
        td = time.time()
        cmd = (f"cd {WORK} && diamond blastp -d $HOME/reference/swissprot_diamond -q q.faa "
               f"--sensitive --evalue 1e-5 --max-target-seqs 5 --threads {args.threads} "
               f"--outfmt 6 {OUTFMT} -o hits_live.tsv 2> diamond.err; echo rc=$?; cat hits_live.tsv")
        rc, out, err = wsl(cmd)
        report["diamond_seconds"] = round(time.time() - td, 1)
        if not out.startswith("rc=0"):
            rec("DIAMOND en direct", "FAIL", (out[:80] + err[:150]))
        else:
            live: dict[str, list[dict]] = {p: [] for p in pick}
            for line in out.splitlines()[1:]:
                p = line.split("\t", len(HIT_COLS) - 1)
                if p[0] in live:
                    live[p[0]].append(dict(zip(HIT_COLS, p)))
            n_exact = n_tie = n_bad = 0
            for p in pick:
                c = compare_diamond(fz_hits[p], live[p])
                report["diamond"][p] = c
                if c["exact"]:
                    n_exact += 1
                elif c["tie_at_limit_only"]:
                    n_tie += 1
                else:
                    n_bad += 1
            rec("DIAMOND live = gele", "PASS" if n_bad == 0 else "FAIL",
                f"identiques {n_exact}, egalite au 5e rang {n_tie}, differents {n_bad} "
                f"({report['diamond_seconds']} s)")
            for p in pick:
                c = report["diamond"][p]
                if not c["exact"] and not c["tie_at_limit_only"]:
                    print(f"       ! {p}: gele {c['n_frozen']} / live {c['n_live']} ; "
                          f"seulement gele {c['only_frozen'][:3]} ; seulement live {c['only_live'][:3]} ; "
                          f"champs {c['field_diffs'][:2]}")

    # 5. Pfam en direct
    if not args.skip_pfam:
        tp = time.time()
        cmd = (f"cd {WORK} && hmmscan --cut_ga --noali --cpu {args.threads} --domtblout dom_live.domtbl "
               f"$HOME/reference/Pfam-A.hmm q.faa > hmmscan.out 2> hmmscan.err; echo rc=$?; cat dom_live.domtbl")
        rc, out, err = wsl(cmd, timeout=7200)
        report["pfam_seconds"] = round(time.time() - tp, 1)
        if not out.startswith("rc=0"):
            rec("Pfam en direct", "FAIL", (out[:80] + err[:150]))
        else:
            rows = parse_domtbl(out.split("\n", 1)[1] if "\n" in out else "")
            live_dom: dict[str, list[dict]] = {p: [] for p in pick}
            for r in rows:
                if r["protein_accession"] in live_dom:
                    live_dom[r["protein_accession"]].append(r)
            n_exact = n_bad = 0
            for p in pick:
                c = compare_pfam(fz_dom[p], live_dom[p])
                report["pfam"][p] = c
                n_exact += c["exact"]
                n_bad += not c["exact"]
            rec("Pfam live = gele", "PASS" if n_bad == 0 else "FAIL",
                f"identiques {n_exact}, differents {n_bad} ({report['pfam_seconds']} s)")
            for p in pick:
                c = report["pfam"][p]
                if not c["exact"]:
                    print(f"       ! {p}: gele {c['n_frozen']} / live {c['n_live']} ; "
                          f"seulement gele {c['only_frozen'][:2]} ; seulement live {c['only_live'][:2]} ; "
                          f"champs {c['field_diffs'][:2]}")

    wsl(f"rm -rf {WORK}")
    report["checks"] = checks
    report["verdict"] = "OK" if all(s != "FAIL" for _, s, _ in checks) else "ECHEC"
    report["duration_seconds"] = round(time.time() - t0, 1)
    print(f"\nVerdict : {report['verdict']}  ({report['duration_seconds']} s)")
    if args.out:
        o = Path(args.out)
        o.parent.mkdir(parents=True, exist_ok=True)
        o.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Rapport : {o}")


if __name__ == "__main__":
    main()
