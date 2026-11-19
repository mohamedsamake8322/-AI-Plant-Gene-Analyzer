#!/usr/bin/env python
"""Branche le panneau d'annotation en direct dans views/home.py (avant le bloc ESM-2).

Usage (depuis C:\\Downloads\\IA) :
    python scripts\\patch_home_live_annotation.py            # simulation
    python scripts\\patch_home_live_annotation.py --apply    # ecrit views\\home.py (sauvegarde .bak)

Prerequis : core_engines\\live_annotation.py copie dans le projet.
Idempotent. Pour desactiver sans toucher au code : set LIVE_ANNOTATION=0 avant de lancer Streamlit.
"""
import argparse
import shutil
import sys
from pathlib import Path

MARK = "# LIVE_ANNOTATION_PANEL"
ANCHOR = ("        esm2_results = None\n"
          "        esm2_error = None\n"
          '        if sequence_type == "protein":\n')
BLOCK = (f"        {MARK} (DIAMOND + Pfam en direct, preuve separee du score de similarite)\n"
         "        try:\n"
         "            from core_engines.live_annotation import render_live_annotation_panel\n"
         "            render_live_annotation_panel(\n"
         "                st, translate, pd, sequence, sequence_type, current_analysis_signature, logger, bio\n"
         "            )\n"
         "        except Exception:\n"
         '            logger.exception("Live annotation panel failed")\n\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--file", default="views/home.py")
    ap.add_argument("--engine", default="core_engines/live_annotation.py")
    a = ap.parse_args()
    path = Path(a.file)
    if not Path(a.engine).exists():
        sys.exit(f"ERREUR : {a.engine} absent. Copiez live_annotation.py dans core_engines\\ d'abord.")
    raw = path.read_bytes().decode("utf-8")
    crlf = "\r\n" in raw
    text = raw.replace("\r\n", "\n")
    if MARK in text:
        print("Deja applique : rien a faire.")
        return
    n = text.count(ANCHOR)
    if n != 1:
        sys.exit(f"ERREUR : ancre trouvee {n} fois (attendu 1). Rien n'est modifie.")
    new = text.replace(ANCHOR, BLOCK + ANCHOR)
    for needed in ("current_analysis_signature", "logger", "pd.", "bio."):
        if needed not in text:
            sys.exit(f"ERREUR : '{needed}' absent de {path} ; le panneau ne peut pas etre branche tel quel.")
    if not a.apply:
        print(f"Simulation OK : 1 insertion de {BLOCK.count(chr(10))} lignes avant le bloc ESM-2 de {path}. "
              "Relancez avec --apply.")
        return
    shutil.copy2(path, str(path) + ".bak_live_annotation")
    out = new.replace("\n", "\r\n") if crlf else new
    path.write_bytes(out.encode("utf-8"))
    print(f"Ecrit : {path} (sauvegarde : {path}.bak_live_annotation)")


if __name__ == "__main__":
    main()
