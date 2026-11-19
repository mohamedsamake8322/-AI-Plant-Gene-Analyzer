#!/usr/bin/env python
"""Prepare les bases DIAMOND (Swiss-Prot 2026_03) et Pfam 38.2 dans WSL pour le
scan EN DIRECT. Simulation par defaut ; --apply pour construire.

Ce que fait le script
  Simulation (defaut, lecture seule, aucune ecriture) :
    - verifie WSL, les versions de diamond / HMMER, l'empreinte du binaire diamond ;
    - verifie l'empreinte SHA-256 de Pfam-A.hmm.gz et du FASTA Swiss-Prot dans le
      tar.gz, contre les valeurs des preuves (evidence.json) ;
    - verifie l'espace libre sur C: (le disque virtuel de WSL y vit).
  --apply :
    - decompresse Pfam-A.hmm dans WSL (~/reference), hmmpress, verifie 30134 modeles ;
    - extrait uniprot_sprot.fasta.gz, reverifie son SHA-256, diamond makedb ;
    - verifie 575 748 sequences ;
    - ecrit un journal JSON dans results\\.

Les bases sont creees dans le systeme de fichiers NATIF de WSL (rapide), pas
dans /mnt/c. Rien n'est modifie cote Windows ni dans PostgreSQL.

Usage (depuis C:\\Downloads\\IA, venv active) :
    python scripts\\prepare_live_annotation_assets.py
    python scripts\\prepare_live_annotation_assets.py --apply
    python scripts\\prepare_live_annotation_assets.py --apply --stage-via-windows   (si /mnt/c echoue)
    python scripts\\prepare_live_annotation_assets.py --verify-only   (controle sans rien construire)
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]

# Valeurs issues des evidence.json (identiques pour les 4 especes).
EXPECTED = {
    "pfam_gz_sha256": "2d82087b6c5c60d762cc767f98e8260b273134c215ab7efccf7440614a4e5dab",
    "swissprot_fasta_gz_sha256": "a9c3496aa727e61d25d9e1291c9ba28d2c2ad19c0e26569dc083d53ec4169536",
    "diamond_binary_sha256": "d3d72578fa97e87cc0846a06d2f9d7dc481a2106265f9546bc91512972fcd1f1",
    "diamond_version": "2.2.8",
    "hmmer_version": "3.4",
    "pfam_models": 30134,
    "swissprot_entries": 575748,
}
# Tailles (octets) des fichiers Pfam de la construction validee : Pfam-A.hmm decompresse
# depuis le .gz dont l'empreinte est conforme, puis hmmpress = "Pressed and indexed 30134 HMMs".
PFAM_SIZES = {
    "Pfam-A.hmm": 2246909846,
    "Pfam-A.hmm.h3f": 531685301,
    "Pfam-A.hmm.h3i": 2983393,
    "Pfam-A.hmm.h3m": 928801003,
    "Pfam-A.hmm.h3p": 1094551483,
}
NEEDED_GB = 4.5          # estimation : Pfam decompresse + index + Swiss-Prot + base
MIN_MARGIN_GB = 3.0      # marge a garder libre sur C: apres construction

DEFAULT_PFAM_GZ = ROOT / "pipeline_data" / "references" / "pfam" / "38.2" / "Pfam-A.hmm.gz"
DEFAULT_TAR = (
    ROOT / "pipeline_data" / "references" / "swissprot" / "2026_03"
    / "uniprot_sprot-only2026_03.tar.gz"
)

checks: list[dict] = []


def record(name: str, status: str, detail: str = "", blocking: bool = False) -> None:
    checks.append({"check": name, "status": status, "detail": detail, "blocking": blocking})
    icon = {"PASS": "OK  ", "WARN": "WARN", "FAIL": "FAIL"}[status]
    print(f"[{icon}] {name}" + (f" : {detail}" if detail else ""))


def wsl_bash(script: str, timeout: int = 120, stream: bool = False):
    cmd = ["wsl.exe", "--exec", "bash", "-lc", script]
    if stream:
        return subprocess.run(cmd, timeout=timeout).returncode, ""
    try:
        out = subprocess.run(
            cmd, capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=timeout,
        )
        return out.returncode, ((out.stdout or "") + (out.stderr or "")).strip()
    except FileNotFoundError:
        return None, "wsl.exe introuvable"
    except subprocess.TimeoutExpired:
        return None, "timeout"


def to_wsl_path(path: Path) -> str:
    s = str(path.resolve())
    if len(s) > 2 and s[1] == ":":
        return "/mnt/" + s[0].lower() + "/" + s[3:].replace("\\", "/")
    return s


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def find_swissprot_member(tar_path: Path) -> tuple[str | None, str | None]:
    """Nom du membre uniprot_sprot.fasta.gz dans le tar.gz et son SHA-256."""
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar:
            if Path(member.name).name == "uniprot_sprot.fasta.gz" and member.isfile():
                h = hashlib.sha256()
                src = tar.extractfile(member)
                assert src is not None
                for chunk in iter(lambda: src.read(8 * 1024 * 1024), b""):
                    h.update(chunk)
                return member.name, h.hexdigest()
    return None, None


def preflight(args) -> dict:
    info: dict = {}

    # 1. Espace disque
    free_gb = shutil.disk_usage(ROOT).free / 1024**3
    info["free_gb_before"] = round(free_gb, 2)
    margin = free_gb - NEEDED_GB
    if margin < MIN_MARGIN_GB:
        record(
            "Espace libre sur le disque du projet", "FAIL",
            f"{free_gb:.1f} Go libres, besoin ~{NEEDED_GB} Go + marge {MIN_MARGIN_GB} Go. "
            "Liberez de l'espace (ex. deplacer les dumps vers la cle USB).",
            blocking=True,
        )
    else:
        record("Espace libre sur le disque du projet", "PASS", f"{free_gb:.1f} Go libres")

    # 2. WSL
    rc, out = wsl_bash("echo $HOME; uname -r")
    if rc != 0:
        record("WSL accessible", "FAIL", out, blocking=True)
        return info
    home = out.splitlines()[0].strip()
    info["wsl_home"] = home
    record("WSL accessible", "PASS", f"HOME={home}")

    # 3. Outils et versions
    rc, out = wsl_bash("diamond version 2>&1 | head -1")
    info["diamond_version_text"] = out
    if rc != 0 or "diamond" not in out.lower():
        record("diamond installe dans WSL", "FAIL", out, blocking=True)
    elif EXPECTED["diamond_version"] not in out:
        record("Version diamond", "WARN",
               f"{out} (preuves : {EXPECTED['diamond_version']}) - resultats potentiellement differents")
    else:
        record("Version diamond", "PASS", out)

    rc, out = wsl_bash("hmmscan -h 2>&1 | sed -n 2p")
    info["hmmer_version_text"] = out
    if rc != 0 or "HMMER" not in out:
        record("hmmscan installe dans WSL", "FAIL", out, blocking=True)
    elif f"HMMER {EXPECTED['hmmer_version']}" not in out:
        record("Version HMMER", "WARN",
               f"{out} (preuves : {EXPECTED['hmmer_version']}) - resultats potentiellement differents")
    else:
        record("Version HMMER", "PASS", out)

    rc, out = wsl_bash("command -v hmmpress && echo present")
    record("hmmpress installe dans WSL", "PASS" if rc == 0 else "FAIL",
           "" if rc == 0 else "installer le paquet hmmer", blocking=(rc != 0))

    rc, out = wsl_bash("sha256sum \"$(command -v diamond)\" | cut -d' ' -f1")
    info["diamond_binary_sha256"] = out
    if out == EXPECTED["diamond_binary_sha256"]:
        record("Binaire diamond = celui des preuves", "PASS")
    else:
        record("Binaire diamond = celui des preuves", "WARN",
               f"{out[:16]}... different de {EXPECTED['diamond_binary_sha256'][:16]}... "
               "(meme version possible, binaire different) - a documenter")

    # 4. Fichiers sources
    if not args.pfam_gz.is_file():
        record("Pfam-A.hmm.gz present", "FAIL", str(args.pfam_gz), blocking=True)
    else:
        print("      (calcul SHA-256 de Pfam-A.hmm.gz ...)")
        sha = sha256_file(args.pfam_gz)
        info["pfam_gz_sha256"] = sha
        ok = sha == EXPECTED["pfam_gz_sha256"]
        record("Pfam-A.hmm.gz = fichier des preuves", "PASS" if ok else "FAIL",
               "" if ok else f"{sha} attendu {EXPECTED['pfam_gz_sha256']}", blocking=not ok)

    if not args.tar.is_file():
        record("Archive Swiss-Prot presente", "FAIL", str(args.tar), blocking=True)
    else:
        print("      (lecture de l'archive Swiss-Prot, ~1-2 min ...)")
        try:
            member, sha = find_swissprot_member(args.tar)
        except (tarfile.TarError, OSError, EOFError) as exc:
            member, sha = None, None
            record("Archive Swiss-Prot lisible", "FAIL", str(exc), blocking=True)
        else:
            info["swissprot_member"] = member
            info["swissprot_fasta_gz_sha256"] = sha
            if member is None:
                record("uniprot_sprot.fasta.gz dans l'archive", "FAIL", "membre introuvable", blocking=True)
            else:
                ok = sha == EXPECTED["swissprot_fasta_gz_sha256"]
                record("FASTA Swiss-Prot = fichier des preuves", "PASS" if ok else "FAIL",
                       member if ok else f"{sha} attendu {EXPECTED['swissprot_fasta_gz_sha256']}",
                       blocking=not ok)

    # 5. Etat existant
    rc, out = wsl_bash(f"ls -l {shlex.quote(home)}/{shlex.quote(args.wsl_dir)} 2>&1 | head -20")
    info["existing_assets"] = out
    if rc == 0 and "swissprot_diamond.dmnd" in out:
        record("Bases deja presentes dans WSL", "WARN",
               "swissprot_diamond.dmnd existe deja ; --apply la reconstruira")
    return info


def stage_inputs(args, info: dict, dest: Path, skip_pfam: bool = False) -> None:
    """Cote Windows : ecrit Pfam-A.hmm (decompresse) et uniprot_sprot.fasta.gz
    directement dans le dossier natif de WSL (via \\\\wsl.localhost\\...), sans que
    WSL ait a lire /mnt/c (pont DrvFs parfois instable)."""
    dest.mkdir(parents=True, exist_ok=True)
    chunk = 8 * 1024 * 1024

    if not skip_pfam:
        print("      Pfam : decompression Windows -> WSL (1,5 Go, quelques minutes) ...")
        part = dest / "Pfam-A.hmm.part"
        with gzip.open(args.pfam_gz, "rb") as src, part.open("wb") as dst:
            shutil.copyfileobj(src, dst, chunk)
        part.replace(dest / "Pfam-A.hmm")

    print("      Swiss-Prot : extraction du FASTA Windows -> WSL ...")
    part = dest / "uniprot_sprot.fasta.gz.part"
    # Iteration + arret au premier membre trouve : ne lit PAS la fin de l'archive
    # (qui peut etre tronquee sans que le FASTA soit touche).
    found = False
    with tarfile.open(args.tar, "r:gz") as tar:
        for member in tar:
            if member.name == info["swissprot_member"] and member.isfile():
                src = tar.extractfile(member)
                assert src is not None
                with part.open("wb") as dst:
                    shutil.copyfileobj(src, dst, chunk)
                found = True
                break
    if not found:
        raise RuntimeError(f"membre introuvable dans l'archive : {info['swissprot_member']}")
    part.replace(dest / "uniprot_sprot.fasta.gz")


def verify_assets(ref: str, deep: bool = False) -> dict:
    """Controles separes, chacun avec son propre delai ; retourne le detail."""
    q = shlex.quote
    res: dict = {}
    rc, out = wsl_bash(f"stat -c %s {q(ref)}/swissprot_diamond.dmnd 2>&1", timeout=60)
    size = int(out.strip()) if rc == 0 and out.strip().isdigit() else 0
    res["dmnd_bytes"] = size

    rc, out = wsl_bash(
        f"diamond dbinfo -d {q(ref)}/swissprot_diamond.dmnd 2>&1 | grep -E 'Sequences|Letters|Error'",
        timeout=180,
    )
    res["dbinfo"] = out
    seq_ok = size > 0 and re.search(rf"Sequences\s+{EXPECTED['swissprot_entries']}\b", out) is not None
    record("Base DIAMOND : 575 748 sequences", "PASS" if seq_ok else "FAIL",
           f"{size/1e6:.1f} Mo ; {out}" if not seq_ok else f"{size/1e6:.1f} Mo", blocking=not seq_ok)

    if deep:
        rc, out = wsl_bash(f"grep -c '^NAME ' {q(ref)}/Pfam-A.hmm", timeout=900)
        res["pfam_models"] = out
        models_ok = out.strip() == str(EXPECTED["pfam_models"])
        record("Pfam : 30 134 modeles (comptage complet)", "PASS" if models_ok else "FAIL",
               "" if models_ok else f"obtenu : {out!r}", blocking=not models_ok)
    else:
        models_ok = True   # controle leger : voir les tailles ci-dessous

    names = " ".join(q(f"{ref}/{n}") for n in PFAM_SIZES)
    rc, out = wsl_bash(f"stat -c '%n %s' {names}", timeout=120)
    res["pfam_files"] = out
    sizes = {}
    for line in out.splitlines():
        parts = line.rsplit(" ", 1)
        if len(parts) == 2 and parts[1].isdigit():
            sizes[PurePosixPath(parts[0]).name] = int(parts[1])
    idx_ok = all(sizes.get(n) == v for n, v in PFAM_SIZES.items())
    record("Pfam : 5 fichiers aux tailles de la construction validee", "PASS" if idx_ok else "FAIL",
           "" if idx_ok else f"obtenu : {sizes or out!r}", blocking=not idx_ok)
    res["ok"] = bool(seq_ok and models_ok and idx_ok)
    return res


def apply_build(args, info: dict) -> dict:
    home = info["wsl_home"]
    ref = f"{home}/{args.wsl_dir}"
    tmp = f"{ref}/_tmp"
    q = shlex.quote
    head = f"""set -euo pipefail
R={q(ref)}; T={q(tmp)}
mkdir -p "$R" "$T"
rm -f "$R"/*.part
"""
    common_verify = f"""
sync
echo '== verification'
diamond dbinfo -d "$R/swissprot_diamond.dmnd" | head -20
ls -l "$R"
"""
    if args.stage_via_windows:
        rc, distro = wsl_bash("echo $WSL_DISTRO_NAME")
        if rc != 0 or not distro:
            record("Nom de la distribution WSL", "FAIL", distro, blocking=True)
            return {"build_returncode": rc, "ok": False}
        rc, out = wsl_bash(head + "echo dossiers prets")
        if rc != 0:
            record("Preparation du dossier WSL", "FAIL", out, blocking=True)
            return {"build_returncode": rc, "ok": False}
        unc = Path("\\\\wsl.localhost\\" + distro + ref.replace("/", "\\"))
        print(f"\n--- Staging Windows -> {unc} ---")
        stage_inputs(args, info, unc, skip_pfam=args.skip_pfam)
        pfam_block = "" if args.skip_pfam else f"""
echo '== Pfam : hmmpress'
hmmpress -f "$R/Pfam-A.hmm"
n=$(grep -c '^NAME ' "$R/Pfam-A.hmm")
echo "modeles Pfam : $n"
test "$n" -eq {EXPECTED['pfam_models']}
"""
        body = pfam_block + f"""
echo '== Swiss-Prot : verification SHA-256'
F="$R/uniprot_sprot.fasta.gz"
s=$(sha256sum "$F" | cut -d' ' -f1)
test "$s" = {EXPECTED['swissprot_fasta_gz_sha256']}
echo '== DIAMOND : makedb'
diamond makedb --in "$F" -d "$R/swissprot_diamond" --threads 2
sync
rm -f "$F"
rm -rf "$T"
"""
        script = head + body + common_verify
    else:
        pfam_mnt = "" if args.skip_pfam else f"""echo '== Pfam : decompression'
gunzip -c {q(to_wsl_path(args.pfam_gz))} > "$R/Pfam-A.hmm.part"
mv "$R/Pfam-A.hmm.part" "$R/Pfam-A.hmm"
echo '== Pfam : hmmpress'
hmmpress -f "$R/Pfam-A.hmm"
n=$(grep -c '^NAME ' "$R/Pfam-A.hmm")
echo "modeles Pfam : $n"
test "$n" -eq {EXPECTED['pfam_models']}
"""
        script = head + pfam_mnt + f"""echo '== Swiss-Prot : extraction'
tar -xzf {q(to_wsl_path(args.tar))} --occurrence=1 -C "$T" {q(info['swissprot_member'])}
F="$T/{info['swissprot_member']}"
s=$(sha256sum "$F" | cut -d' ' -f1)
test "$s" = {EXPECTED['swissprot_fasta_gz_sha256']}
echo '== DIAMOND : makedb'
diamond makedb --in "$F" -d "$R/swissprot_diamond" --threads 2
sync
rm -rf "$T"
""" + common_verify

    print("\n--- Construction dans WSL (plusieurs minutes) ---")
    rc, _ = wsl_bash(script, timeout=3 * 3600, stream=True)
    result = {"build_returncode": rc, "mode": "stage_via_windows" if args.stage_via_windows else "wsl_reads_mnt_c"}

    v = verify_assets(ref)
    result["verification"] = v
    result["ok"] = rc == 0 and v["ok"]
    free_after = shutil.disk_usage(ROOT).free / 1024**3
    result["free_gb_after"] = round(free_after, 2)
    print(f"Espace libre apres construction : {free_after:.1f} Go")
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="construire les bases dans WSL")
    ap.add_argument("--force", action="store_true", help="ignorer l'avertissement d'espace disque")
    ap.add_argument("--deep", action="store_true",
                    help="avec --verify-only : compte aussi les 30 134 modeles (lit 2,2 Go, lourd)")
    ap.add_argument("--verify-only", action="store_true",
                    help="ne rien construire : verifier les bases deja presentes dans WSL")
    ap.add_argument("--skip-pfam", action="store_true",
                    help="ne pas refaire Pfam (deja presse et verifie dans WSL)")
    ap.add_argument("--stage-via-windows", action="store_true",
                    help="Windows ecrit les fichiers dans WSL (evite la lecture de /mnt/c par WSL)")
    ap.add_argument("--pfam-gz", type=Path, default=DEFAULT_PFAM_GZ)
    ap.add_argument("--tar", type=Path, default=DEFAULT_TAR)
    ap.add_argument("--wsl-dir", default="reference", help="dossier sous $HOME dans WSL")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    if args.verify_only:
        rc, home = wsl_bash("echo $HOME")
        if rc != 0:
            print("WSL inaccessible :", home)
            return 1
        v = verify_assets(f"{home}/{args.wsl_dir}", deep=args.deep)
        print("\nVerification des bases :", "OK" if v["ok"] else "ECHEC")
        return 0 if v["ok"] else 2

    print(f"Mode : {'APPLY' if args.apply else 'SIMULATION (aucune ecriture)'}\n")
    info = preflight(args)

    blocking = [c for c in checks if c["status"] == "FAIL" and c["blocking"]]
    if args.force:
        blocking = [c for c in blocking if not c["check"].startswith("Espace libre")]
    journal = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "mode": "apply" if args.apply else "simulation",
        "expected": EXPECTED,
        "info": info,
        "checks": checks,
    }

    if blocking:
        print(f"\nBLOQUE : {len(blocking)} controle(s) en echec. Rien n'est construit.")
        journal["status"] = "blocked"
        code = 1
    elif not args.apply:
        warns = sum(1 for c in checks if c["status"] == "WARN")
        print(f"\nPreflight OK ({warns} avertissement(s)). Relancez avec --apply pour construire.")
        journal["status"] = "preflight_ok_not_applied"
        code = 0
    else:
        journal["build"] = apply_build(args, info)
        journal["status"] = "built" if journal["build"]["ok"] else "build_failed"
        code = 0 if journal["build"]["ok"] else 2
        print(f"\nStatut : {journal['status']}")

    out = args.out or ROOT / "results" / (
        f"live_assets_{journal['mode']}_{datetime.now():%Y%m%d_%H%M%S}.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(journal, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Journal : {out}")
    return code


if __name__ == "__main__":
    sys.exit(main())
