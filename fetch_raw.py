"""fetch_raw.py : telechargement des fichiers bruts depuis sources.yaml.

Regles appliquees :
  - une source n'est telechargee que si license_checked: true et status != blocked
  - la verification SSL est toujours active (sauf `verify_ssl: false` explicite dans la source)
  - les fichiers sont ecrits dans data/raw/<espece>/<source>/<version>/ et jamais ecrases
  - chaque fichier est enregistre dans data/manifest.csv (url, taille, md5, sha256, date, release)
  - reprise automatique des telechargements interrompus (fichier .part)
  - controle du disque avant de commencer
  - version "mouvante" (rolling/release/courante/...) -> dossier rolling_AAAA-MM-JJ

Usage :
    python fetch_raw.py sources.yaml --source ncbi_refseq_nitab_tn90 --dry-run
    python fetch_raw.py sources.yaml --source ncbi_refseq_nitab_tn90
    python fetch_raw.py sources.yaml --source maizegdb_b73_v5 --include-optional
"""
import argparse
import csv
import hashlib
import re
import shutil
import ssl
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import yaml

UA = {"User-Agent": "Mozilla/5.0 (fetch_raw)"}
MANIFEST_FIELDS = ["source_id", "species", "version", "file", "url", "bytes", "md5", "sha256",
                   "downloaded_at_utc", "last_modified", "etag", "remote_release",
                   "checksum_verified"]
TIMEOUT = 30            # secondes sans aucune donnee -> abandon de la tentative
MIN_SPEED = 20_000      # octets/s minimum sur une fenetre de 30 s, sinon la tentative est abandonnee
ROLLING_PREFIXES = ("rolling", "release", "courante", "current", "continu", "date")


def slug(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")[:60] or "na"


def gate_reason(src):
    if not src.get("license_checked"):
        return "licence non verifiee (license_checked: false)"
    if src.get("status") == "blocked":
        return "source bloquee (status: blocked)"
    return None


def planned_files(src, include_optional):
    """None si la source n'a ni liste `files` ni `single_file` dans sources.yaml :
    la source est alors sautee proprement (un batch --source A --source B continue).

    `single_file: true` : `url` pointe directement sur UN fichier (ou une requete d'API) ;
    `file_name` est obligatoire si l'URL contient ?, & ou : (interdits dans un nom Windows),
    sinon le nom est deduit de la fin de l'URL."""
    if src.get("single_file"):
        name = src.get("file_name") or src["url"].rstrip("/").rsplit("/", 1)[-1]
        return [{"name": name, "size_mb": src.get("size_mb", "?")}]
    files = src.get("files")
    if not files:
        return None
    return [f for f in files if include_optional or not f.get("optional")]


def hash_file(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def download(url, dest, verify_ssl=True, retries=4, timeout=None, min_speed=None):
    """Telecharge url vers dest (reprise possible) avec progression et detection de blocage.

    Une tentative est abandonnee (puis reprise la ou elle s'est arretee) si aucune donnee
    n'arrive pendant `timeout` s, ou si le debit moyen sur 30 s tombe sous `min_speed` o/s."""
    timeout = TIMEOUT if timeout is None else timeout
    min_speed = MIN_SPEED if min_speed is None else min_speed
    part = dest.with_name(dest.name + ".part")
    ctx = None if verify_ssl else ssl._create_unverified_context()
    last_exc = None
    for attempt in range(retries):
        pos = part.stat().st_size if part.exists() else 0
        headers = dict(UA)
        if pos:
            headers["Range"] = f"bytes={pos}-"
            print(f"(reprise a {pos / 1e6:.1f} Mo) ", end="", flush=True)
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                resumed = bool(pos and r.status == 206)
                mode = "ab" if resumed else "wb"
                length = r.headers.get("Content-Length")
                total = (int(length) + (pos if resumed else 0)) if length else None
                release = " ".join(x for x in (r.headers.get("X-UniProt-Release", ""),
                                               r.headers.get("X-UniProt-Release-Date", "")) if x)
                meta = {"last_modified": r.headers.get("Last-Modified", ""),
                        "etag": r.headers.get("ETag", ""), "release": release}
                got = pos if resumed else 0
                t0 = last_print = win_t = time.time()
                win_bytes = 0
                with open(part, mode) as f:
                    while True:
                        chunk = r.read(256 * 1024)
                        if not chunk:
                            break
                        f.write(chunk)
                        got += len(chunk)
                        win_bytes += len(chunk)
                        now = time.time()
                        if now - last_print >= 5:
                            tot = f"/{total / 1e6:.1f}" if total else ""
                            print(f"[{got / 1e6:.1f}{tot} Mo, {win_bytes / max(now - win_t, 1e-6) / 1e3:.0f} Ko/s] ",
                                  end="", flush=True)
                            last_print = now
                        if now - win_t >= 30:
                            if win_bytes / (now - win_t) < min_speed:
                                raise TimeoutError(f"debit trop faible "
                                                   f"({win_bytes / (now - win_t) / 1e3:.1f} Ko/s)")
                            win_t, win_bytes = now, 0
            part.replace(dest)
            return meta
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and part.exists():      # deja complet
                part.replace(dest)
                return {"last_modified": "", "etag": "", "release": ""}
            last_exc = exc
            if exc.code in (401, 403, 404):
                break                                   # inutile de reessayer
        except Exception as exc:                        # reseau, timeout, blocage, SSL...
            last_exc = exc
            print(f"(tentative {attempt + 1}/{retries} abandonnee : {exc}) ", end="", flush=True)
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"echec {url} : {last_exc}")


def parse_md5_file(path):
    """Lit un fichier de type md5sum : '<md5>  ./nom' -> {nom: md5}."""
    sums = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        m = re.match(r"^([0-9a-fA-F]{32})\s+\*?(?:\./)?(.+)$", line.strip())
        if m:
            sums[Path(m.group(2)).name] = m.group(1).lower()
    return sums


def append_manifest(manifest, row):
    manifest.parent.mkdir(parents=True, exist_ok=True)
    new = not manifest.exists()
    if not new:
        with open(manifest, encoding="utf-8", newline="") as f:
            header = next(csv.reader(f), [])
        if header != MANIFEST_FIELDS:
            raise SystemExit(f"{manifest} a un ancien format de colonnes : "
                             "renomme-le (ex. manifest_old.csv) puis relance.")
    with open(manifest, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        if new:
            w.writeheader()
        w.writerow(row)


def run_source(src, args):
    sid = src["id"]
    species = src["species"]
    if isinstance(species, list):
        print(f"[{sid}] source multi-especes : non geree par ce script, a decouper par espece.")
        return False
    reason = gate_reason(src)
    if reason:
        print(f"\n[{sid}]  BLOQUE : {reason}")
        return False

    files = planned_files(src, args.include_optional)
    if files is not None and args.only:
        files = [f for f in files if any(k.lower() in f["name"].lower() for k in args.only)]
        if not files:
            print(f"\n[{sid}]  Aucun fichier ne correspond a --only {args.only}")
            return False
    if files is None:
        print(f"\n[{sid}]  IGNOREE : pas encore de liste `files` dans sources.yaml "
              "(cette source n'est pas prete pour ce script -- a completer plus tard).")
        return False
    def _mb(f):
        try:
            return float(f.get("size_mb", 0))
        except (TypeError, ValueError):
            return 0.0            # taille inconnue ("?") : ne pas planter

    total_mb = sum(_mb(f) for f in files)
    print(f"\n[{sid}] {len(files)} fichiers, ~{total_mb:.0f} Mo declares")
    for f in files:
        opt = " (optionnel)" if f.get("optional") else ""
        print(f"    - {f['name']}  ~{f.get('size_mb', '?')} Mo{opt}")
    if args.dry_run:
        print("  (dry-run : rien n'est telecharge)")
        return True

    version = str(src.get("version", "na"))
    if version.lower().startswith(ROLLING_PREFIXES):
        version = "rolling_" + datetime.now(timezone.utc).strftime("%Y-%m-%d")
    dest_dir = Path(args.raw_dir) / species / sid / slug(version)
    dest_dir.mkdir(parents=True, exist_ok=True)
    free_gb = shutil.disk_usage(dest_dir).free / 1e9
    if free_gb - total_mb * 1.2 / 1000 < args.reserve_gb:
        print(f"  ARRET : disque insuffisant ({free_gb:.1f} Go libres, reserve {args.reserve_gb} Go).")
        return False

    if src.get("single_file"):
        base = None
        single_url = src["url"]
    else:
        base = src["url"].rstrip("/") + "/"
        single_url = None
    verify = src.get("verify_ssl", True)
    if not verify:
        print("  ATTENTION : verification SSL desactivee pour cette source (verify_ssl: false).")

    sums = {}
    if src.get("checksum_file") and base:
        cs_path = dest_dir / src["checksum_file"]
        try:
            if not cs_path.exists():
                download(base + src["checksum_file"], cs_path, verify)
            sums = parse_md5_file(cs_path)
        except RuntimeError as exc:
            print(f"  ATTENTION : fichier de checksums indisponible ({exc}). "
                  "Les fichiers seront telecharges SANS verification md5.")

    ok = True
    for f in files:
        name = f["name"]
        dest = dest_dir / name
        if dest.exists():
            print(f"  = {name} deja present (non touche)")
            continue
        print(f"  > {name} ...", end=" ", flush=True)
        file_url = single_url if single_url else base + name
        try:
            meta = download(file_url, dest, verify)
        except RuntimeError as exc:
            print("ECHEC")
            print(f"    {exc}")
            ok = False
            continue
        md5, sha = hash_file(dest)
        verified = ""
        if name in sums:
            verified = "oui" if sums[name] == md5 else "NON"
            if verified == "NON":
                dest.rename(dest.with_name(name + ".BAD"))
                print("CHECKSUM INCORRECT (fichier renomme .BAD)")
                ok = False
                continue
        append_manifest(Path(args.manifest), {
            "source_id": sid, "species": species, "version": version,
            "file": str(dest), "url": file_url, "bytes": dest.stat().st_size,
            "md5": md5, "sha256": sha,
            "downloaded_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "last_modified": meta["last_modified"], "etag": meta["etag"],
            "remote_release": meta["release"],
            "checksum_verified": verified})
        print(f"ok ({dest.stat().st_size / 1e6:.1f} Mo{', checksum ok' if verified == 'oui' else ''})")
    return ok


def main():
    global TIMEOUT, MIN_SPEED
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("yaml_path")
    ap.add_argument("--source", action="append", required=True, help="id de source (repetable)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--include-optional", action="store_true")
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--manifest", default="data/manifest.csv")
    ap.add_argument("--reserve-gb", type=float, default=5.0, help="espace disque a garder libre")
    ap.add_argument("--only", action="append", default=[],
                    help="ne telecharger que les fichiers dont le nom contient ce texte (repetable), "
                         "ex. --only genomic.gff --only gene_ontology")
    ap.add_argument("--timeout", type=int, default=TIMEOUT, help="secondes sans donnee avant abandon")
    ap.add_argument("--min-speed-kbs", type=float, default=MIN_SPEED / 1000,
                    help="debit minimal (Ko/s) sur 30 s avant abandon de la tentative")
    args = ap.parse_args()
    TIMEOUT, MIN_SPEED = args.timeout, int(args.min_speed_kbs * 1000)

    cfg = yaml.safe_load(open(args.yaml_path, encoding="utf-8"))
    by_id = {s["id"]: s for s in cfg["sources"]}
    all_ok = True
    for sid in args.source:
        if sid not in by_id:
            print(f"Source inconnue : {sid}")
            all_ok = False
            continue
        all_ok &= run_source(by_id[sid], args)
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()