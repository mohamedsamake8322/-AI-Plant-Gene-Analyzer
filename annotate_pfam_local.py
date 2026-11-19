#!/usr/bin/env python3
"""Pfam (hmmscan --cut_ga) par lots : version LOCALE (WSL/Linux) de la cellule 5 du notebook 01.

Même fichiers de sortie, même preuve de passage (evidence.json) : le chargeur lit indifféremment
l'un ou l'autre. Ne pas traiter la même espèce à la fois ici et sur Colab.

Exemple (dans WSL, Drive monté en /mnt/g) :
  python3 annotate_pfam_local.py --data-root /mnt/g/MyDrive/files/data --pfam-release 38.2 \
      --only zea_mays nicotiana_tabacum_k326 --threads 8
"""
import argparse, gzip, hashlib, json, os, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SPECIES = [
    {"key": "oryza_sativa",           "assembly": "GCF_001433935.1", "expected": 42580},
    {"key": "zea_mays",               "assembly": "GCF_902167145.1", "expected": 57345},
    {"key": "nicotiana_tabacum_k326", "assembly": "GCF_000715075.1", "expected": 91379},
    {"key": "nicotiana_tabacum_tn90", "assembly": "GCF_000715135.1", "expected": 84350},
]
PFAM_PARAMS = {"tool": "hmmscan", "cutoff": "--cut_ga", "alignments": False}   # NE PAS MODIFIER (identique au notebook)
DOMTBL_COLUMNS = ["target_name", "target_accession", "tlen", "query_name", "query_accession", "qlen",
                  "full_evalue", "full_score", "full_bias", "domain_n", "domain_of", "c_evalue",
                  "i_evalue", "dom_score", "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to",
                  "env_from", "env_to", "acc", "description"]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def read_fasta(path):
    opener = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    records, accession, parts = {}, None, []
    with opener as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession is not None:
                    records[accession] = "".join(parts).upper()
                accession = line[1:].split()[0]
                assert accession not in records, f"Accession dupliquée : {accession}"
                parts = []
            else:
                parts.append(line.replace(" ", ""))
    if accession is not None:
        records[accession] = "".join(parts).upper()
    return records


def sequence_set_fingerprint(proteins):
    md5s = sorted(hashlib.md5(s.encode()).hexdigest() for s in proteins.values())
    return sha256_text("\n".join(md5s))


def find_query_file(data_root, assembly):
    found = sorted(data_root.rglob(f"*{assembly}*protein.faa.gz")) or sorted(data_root.rglob(f"*{assembly}*protein.faa"))
    if len(found) > 1:
        raise RuntimeError(f"Plusieurs fichiers pour {assembly} : {found}")
    return found[0] if found else None


def parse_domtbl(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split(None, 22)
            if len(parts) == 22:
                parts.append("")
            rows.append(parts)
    return pd.DataFrame(rows, columns=DOMTBL_COLUMNS)


def domtbl_complete(path):
    path = Path(path)
    if not path.exists():
        return False
    last = ""
    with open(path) as handle:
        for line in handle:
            if line.strip():
                last = line.strip()
    return last == "# [ok]"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", required=True, type=Path, help="dossier 'files/data' du Drive")
    ap.add_argument("--pfam-release", required=True, help="ex. 38.2 : dossier references/pfam/<version> du Drive")
    ap.add_argument("--only", nargs="*", help="clés d'espèces à traiter (défaut : toutes)")
    ap.add_argument("--jobs", type=int, default=1, help="lots traités en parallèle (processus hmmscan simultanés)")
    ap.add_argument("--threads", type=int, default=None, help="threads par processus hmmscan (défaut : cœurs / jobs)")
    ap.add_argument("--chunk", type=int, default=2000)
    ap.add_argument("--benchmark-n", type=int, default=0, help="N protéines seulement : mesure du temps, jamais une preuve")
    ap.add_argument("--workdir", type=Path, default=Path.home() / "pfam_work")
    args = ap.parse_args()
    if args.threads is None:
        args.threads = max(1, (os.cpu_count() or 2) // max(1, args.jobs))
    assert args.threads >= 1 and args.chunk >= 1 and args.jobs >= 1

    if shutil.which("hmmscan") is None or shutil.which("hmmpress") is None:
        sys.exit("hmmscan/hmmpress introuvables : sudo apt-get install hmmer (dans WSL).")
    help_text = subprocess.run(["hmmscan", "-h"], capture_output=True, text=True).stdout
    found = re.search(r"HMMER\s+([^\n]+)", help_text)
    hmmer_version = found.group(1).strip() if found else "inconnue"

    pf_gz = args.data_root / "references" / "pfam" / args.pfam_release / "Pfam-A.hmm.gz"
    if not pf_gz.exists():
        sys.exit(f"Fichier Pfam introuvable : {pf_gz} (celui téléchargé par le notebook 01)")
    pfam_sha = sha256_file(pf_gz)
    work = args.workdir / args.pfam_release; work.mkdir(parents=True, exist_ok=True)
    hmm = work / "Pfam-A.hmm"
    marker = work / "Pfam-A.pressed.ok"     # écrit seulement après un hmmpress réussi
    if marker.exists() and marker.read_text().strip() != pfam_sha:
        marker.unlink()
    if not marker.exists():
        # état inconnu ou interrompu : on repart proprement (décompression puis indexation)
        for stale in work.glob("Pfam-A.hmm*"):
            stale.unlink()
        print("Préparation de Pfam (décompression + hmmpress, quelques minutes, une seule fois)...", flush=True)
        tmp_hmm = work / "Pfam-A.hmm.tmp"
        with gzip.open(pf_gz, "rb") as src, tmp_hmm.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 20)
        tmp_hmm.replace(hmm)
        with open(hmm, "rb") as handle:
            handle.seek(max(0, hmm.stat().st_size - 64))
            if not handle.read().rstrip().endswith(b"//"):
                sys.exit("Pfam-A.hmm tronqué (ne finit pas par '//') : relancer.")
        proc = subprocess.run(["hmmpress", "-f", str(hmm)], capture_output=True, text=True)
        if proc.returncode != 0:
            sys.exit("hmmpress a échoué : " + proc.stderr[-500:])
        suffixes = (".h3f", ".h3i", ".h3m", ".h3p")
        assert all(hmm.with_suffix(".hmm" + s_).exists() for s_ in suffixes), "Index hmmpress incomplet"
        marker.write_text(pfam_sha)
    models = sum(1 for line in open(hmm) if line.startswith("NAME "))
    print(f"HMMER {hmmer_version} | Pfam {args.pfam_release} : {models:,} modèles | {args.jobs} lot(s) en parallèle x {args.threads} thread(s)")

    phash = sha256_text(json.dumps(PFAM_PARAMS, sort_keys=True))
    targets = [s for s in SPECIES if not args.only or s["key"] in args.only]
    annot = args.data_root / "annotation"
    for sp in targets:
        key, assembly = sp["key"], sp["assembly"]
        path = find_query_file(args.data_root, assembly)
        if path is None:
            print(f"[ABSENT] {key} : aucun protein.faa(.gz) pour {assembly}"); continue
        proteins = read_fasta(path)
        assert all(s and "*" not in s for s in proteins.values()), f"{key} : séquence vide ou contenant '*'"
        if sp["expected"] is not None:
            assert len(proteins) == sp["expected"], f"{key} : {len(proteins)} protéines, attendu {sp['expected']}"
        fingerprint = sequence_set_fingerprint(proteins)
        ids = list(proteins)[:args.benchmark_n] if args.benchmark_n else list(proteins)

        folder = annot / key / "pfam_hmmscan" / (args.pfam_release + ("_BENCHMARK" if args.benchmark_n else ""))
        chunks_dir = folder / "chunks"; chunks_dir.mkdir(parents=True, exist_ok=True)
        evidence_path = folder / "evidence.json"
        if not args.benchmark_n and evidence_path.exists():
            ev = json.loads(evidence_path.read_text())
            if (ev.get("run_completed") and ev.get("sequence_set_fingerprint") == fingerprint
                    and ev.get("params_hash") == phash and ev.get("database_release") == args.pfam_release):
                print(f"[{key}] Pfam déjà fait pour {args.pfam_release} : ignoré"); continue

        # Les lots sont numérotés selon la taille de lot : toute reprise doit utiliser exactement les mêmes réglages.
        manifest = {"chunk": args.chunk, "n_ids": len(ids), "ids_sha256": sha256_text("\n".join(ids))}
        manifest_path = chunks_dir / "manifest.json"
        if args.benchmark_n:                      # benchmark : toujours repartir de zéro
            shutil.rmtree(chunks_dir); chunks_dir.mkdir(parents=True)
        if manifest_path.exists() and chunks_dir.exists():
            if json.loads(manifest_path.read_text()) != manifest:
                sys.exit(f"[{key}] Des lots existent déjà avec d'autres réglages (taille de lot ou protéines différentes) : "
                         f"relancer avec les mêmes options, ou supprimer {chunks_dir}")
        elif any(chunks_dir.glob("chunk_*.domtbl")):
            sys.exit(f"[{key}] Lots existants sans manifeste (origine inconnue) : supprimer {chunks_dir} pour repartir proprement")
        manifest_path.write_text(json.dumps(manifest))

        started = datetime.now(timezone.utc); t0 = time.time()
        n_chunks = (len(ids) + args.chunk - 1) // args.chunk
        pending = [i for i in range(n_chunks) if not domtbl_complete(chunks_dir / f"chunk_{i:05d}.domtbl")]

        def scan_chunk(index):
            out = chunks_dir / f"chunk_{index:05d}.domtbl"
            part = ids[index * args.chunk:(index + 1) * args.chunk]
            fasta = work / f"{key}_chunk_{index:05d}.faa"
            with fasta.open("w") as handle:
                for accession in part:
                    handle.write(f">{accession}\n{proteins[accession]}\n")
            tmp_out = work / f"{key}_chunk_{index:05d}.domtbl"
            proc = subprocess.run(["hmmscan", "--cpu", str(args.threads), PFAM_PARAMS["cutoff"], "--noali",
                                   "--domtblout", str(tmp_out), "-o", "/dev/null", str(hmm), str(fasta)],
                                  capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"hmmscan a échoué ({key}, lot {index}) : {proc.stderr[-500:]}")
            if not domtbl_complete(tmp_out):
                raise RuntimeError(f"Lot {index} incomplet (pas de '# [ok]')")
            shutil.copy(tmp_out, out)
            fasta.unlink(missing_ok=True); tmp_out.unlink(missing_ok=True)
            return index

        finished = n_chunks - len(pending)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(scan_chunk, i) for i in pending]
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    for f in futures:
                        f.cancel()
                    sys.exit(f"Échec : {exc}")
                finished += 1
                print(f"[{key}] lots terminés {finished}/{n_chunks} | {(time.time() - t0) / 60:.1f} min", flush=True)

        chunk_files = sorted(chunks_dir.glob("chunk_*.domtbl"))
        completed = len(chunk_files) == n_chunks and all(domtbl_complete(p) for p in chunk_files)
        domains = pd.concat([parse_domtbl(p) for p in chunk_files], ignore_index=True) if chunk_files else pd.DataFrame(columns=DOMTBL_COLUMNS)
        assert set(domains["query_name"]) <= set(ids), "Domaines pour des protéines absentes de la requête"
        domains_out = domains.rename(columns={"query_name": "protein_accession", "target_accession": "pfam_accession",
                                              "target_name": "pfam_name"})[
            ["protein_accession", "pfam_accession", "pfam_name", "i_evalue", "dom_score",
             "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "description"]]
        domains_out.to_csv(folder / "domains.tsv", sep="\t", index=False)
        scanned = len(ids) if completed else None
        with_domain = int(domains["query_name"].nunique())
        evidence = {
            "assembly_accession": assembly, "species_key": key,
            "tool_name": "pfam_hmmscan", "tool_version": hmmer_version,
            "database_name": f"pfam_{args.pfam_release}", "database_release": args.pfam_release,
            "database_version_text": f"Pfam {args.pfam_release} (fichier de référence partagé du Drive)",
            "database_hmm_gz_sha256": pfam_sha, "database_models": models,
            "parameters": PFAM_PARAMS, "params_hash": phash,
            "query_file": path.name, "query_file_sha256": sha256_file(path),
            "sequence_set_fingerprint": fingerprint, "input_proteins": len(proteins),
            "queried_proteins": len(ids), "chunks": n_chunks, "run_completed": completed,
            "scanned_proteins": scanned, "proteins_with_domain": with_domain,
            "proteins_without_domain": (scanned - with_domain) if scanned is not None else None,
            "domain_rows": int(len(domains)), "output_tsv": "domains.tsv",
            "output_tsv_sha256": sha256_file(folder / "domains.tsv"),
            "started_at": started.isoformat(), "finished_at": datetime.now(timezone.utc).isoformat(),
            "duration_seconds_this_session": round(time.time() - t0, 1), "cpu_count": os.cpu_count(),
            "benchmark_only": bool(args.benchmark_n), "run_environment": "local",
        }
        evidence_path.write_text(json.dumps(evidence, indent=2), encoding="utf-8")
        print(f"[{key}] preuve écrite : {evidence_path}")
        if args.benchmark_n:
            per = (time.time() - t0) / len(ids)
            print(f"[{key}] BENCHMARK : {per:.2f} s/protéine -> espèce complète ~{per * len(proteins) / 3600:.1f} h")


if __name__ == "__main__":
    main()
