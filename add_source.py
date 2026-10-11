"""add_source.py : ajoute un bloc de sources (fichier .yaml) a la fin de sources.yaml.

Usage :
    python add_source.py sources.yaml nitab_asm71507v2.yaml

Securites : refuse si un id existe deja, garde une sauvegarde sources.yaml.bak,
verifie que le fichier final est un YAML valide (sinon restaure la sauvegarde).
Le bloc doit etre une liste YAML indentee de 2 espaces ("  - id: ..."), et `sources:`
doit etre la derniere cle de sources.yaml (c'est le cas du fichier actuel).
"""
import shutil
import sys
from pathlib import Path

import yaml


def main(main_path, snippet_path):
    main_p, snip_p = Path(main_path), Path(snippet_path)
    main_text = main_p.read_text(encoding="utf-8")
    snippet = snip_p.read_text(encoding="utf-8")

    cfg = yaml.safe_load(main_text)
    try:
        new_items = yaml.safe_load(snippet)
    except yaml.YAMLError as exc:
        sys.exit(f"Le bloc {snip_p.name} n'est pas un YAML valide (rien n'est modifie) :\n{exc}")
    if not isinstance(new_items, list):
        sys.exit("Le bloc doit etre une liste YAML de sources (lignes commencant par '- id:').")
    existing = {s["id"] for s in cfg["sources"]}
    dup = [s["id"] for s in new_items if s["id"] in existing]
    if dup:
        sys.exit(f"Deja present dans {main_p.name} : {', '.join(dup)} (rien n'est modifie).")

    backup = main_p.with_suffix(main_p.suffix + ".bak")
    shutil.copy2(main_p, backup)
    if not main_text.endswith("\n"):
        main_text += "\n"
    main_p.write_text(main_text + "\n" + snippet.rstrip("\n") + "\n", encoding="utf-8")

    try:
        check = yaml.safe_load(main_p.read_text(encoding="utf-8"))
        ids = {s["id"] for s in check["sources"]}
        assert all(s["id"] in ids for s in new_items)
    except Exception as exc:
        shutil.copy2(backup, main_p)
        sys.exit(f"YAML invalide apres ajout ({exc}). Fichier restaure depuis {backup.name}.")
    print(f"Ajoute : {', '.join(s['id'] for s in new_items)}  "
          f"(sauvegarde : {backup.name}, {len(check['sources'])} sources au total)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:3])
