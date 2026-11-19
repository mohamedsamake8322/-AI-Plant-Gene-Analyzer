#!/usr/bin/env python
"""Branche le panneau d'annotation en direct dans l'onglet Protein Analysis de views/independent_tools.py.

Le panneau est place HORS du bloc `if st.button("Analyze protein")` (sinon son propre bouton ferait
disparaitre les resultats) : il lit la zone de texte `protein_input` et garde son resultat en session.

Usage (depuis C:\\Downloads\\IA) :
    python scripts\\patch_independent_tools_live_annotation.py            # simulation
    python scripts\\patch_independent_tools_live_annotation.py --apply    # ecrit (sauvegarde .bak)
Prerequis : core_engines\\live_annotation.py present. Idempotent. LIVE_ANNOTATION=0 pour masquer.
"""
import argparse
import shutil
import sys
from pathlib import Path

MARK = "# LIVE_ANNOTATION_PANEL"
ANCHOR = ('                    key="independent_protein_export_csv",\n'
          '                )\n'
          '\n'
          'with tool_tabs[4]:\n')
BLOCK = (
    f"    {MARK} (DIAMOND + Pfam en direct ; hors du bloc 'Analyze protein')\n"
    "    try:\n"
    "        import logging as _logging_live\n"
    "        import pandas as _pd_live\n"
    "        from core_engines.live_annotation import render_live_annotation_panel\n"
    '        _live_header, _live_seq = bio.parse_fasta_input((protein_input or "").strip())\n'
    '        _live_clean = bio.clean_sequence(_live_seq, sequence_type="protein")\n'
    '        _live_ok, _live_msg = bio.validate_sequence(_live_clean, sequence_type="protein")\n'
    "        if _live_clean and _live_ok:\n"
    "            render_live_annotation_panel(\n"
    '                st, translate, _pd_live, _live_clean, "protein",\n'
    '                ("independent_tools", len(_live_clean), _live_clean[:40]),\n'
    "                _logging_live.getLogger(__name__), bio,\n"
    "            )\n"
    "    except Exception as _live_error:\n"
    '        st.warning(f"Live annotation unavailable: {_live_error}")\n'
    "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--file", default="views/independent_tools.py")
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
    new = text.replace(ANCHOR, ANCHOR.replace("\nwith tool_tabs[4]:\n", "\n") + BLOCK + "with tool_tabs[4]:\n")
    if not a.apply:
        print(f"Simulation OK : 1 insertion de {BLOCK.count(chr(10))} lignes avant 'with tool_tabs[4]' de {path}. "
              "Relancez avec --apply.")
        return
    shutil.copy2(path, str(path) + ".bak_live_annotation")
    path.write_bytes((new.replace("\n", "\r\n") if crlf else new).encode("utf-8"))
    print(f"Ecrit : {path} (sauvegarde : {path}.bak_live_annotation)")


if __name__ == "__main__":
    main()
