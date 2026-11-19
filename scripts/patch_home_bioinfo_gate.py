"""Fait reconnaitre USE_BIOINFO=1 par les controles 'base disponible' de views/home.py.
Simulation par defaut ; --apply ecrit (sauvegarde .bak_gate). Idempotent."""
import re, sys, shutil
from pathlib import Path

TARGET = Path("views/home.py")
OLD = b"config.USE_POSTGRES and get_gene_count is not None"
NEW = (b"(config.USE_POSTGRES or __import__('core_engines.bioinfo_reader', "
       b"fromlist=['x']).is_enabled()) and get_gene_count is not None")

def main():
    apply = "--apply" in sys.argv
    if not TARGET.exists():
        print("Introuvable :", TARGET); sys.exit(1)
    data = TARGET.read_bytes()
    n_old = data.count(OLD)
    if NEW in data and n_old == 0:
        print("Deja patche."); return
    if n_old == 0:
        print("Motif non trouve : rien a faire / fichier different."); sys.exit(2)
    print(f"{n_old} occurrence(s) a modifier dans {TARGET}")
    if not apply:
        print("Simulation OK. Relancez avec --apply."); return
    shutil.copy2(TARGET, str(TARGET) + ".bak_gate")
    TARGET.write_bytes(data.replace(OLD, NEW))
    print("Ecrit (sauvegarde views/home.py.bak_gate)")

main()
