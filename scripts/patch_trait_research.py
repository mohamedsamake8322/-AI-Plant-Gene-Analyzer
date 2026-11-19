#!/usr/bin/env python3
"""Ajoute la normalisation './' -> espace côté requête dans trait_research.py (une copie .bak est créée).
Usage : python scripts\\patch_trait_research.py    (depuis C:\\Downloads\\IA)"""
import shutil
from pathlib import Path
p = Path(__file__).resolve().parents[1] / "trait_research.py"
src = p.read_text(encoding="utf-8")
old = "websearch_to_tsquery('simple', %s)"
new = "websearch_to_tsquery('simple', translate(%s, './', '  '))"
n = src.count(old)
if n == 0 and src.count(new) >= 1:
    raise SystemExit("Déjà corrigé.")
if n != 2:
    raise SystemExit(f"Attendu 2 occurrences, trouvé {n} : je ne modifie rien.")
shutil.copy2(p, p.with_suffix(".py.bak"))
p.write_text(src.replace(old, new), encoding="utf-8")
print("trait_research.py corrigé (2 occurrences). Sauvegarde :", p.with_suffix(".py.bak"))
