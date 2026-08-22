"""
Affiche metadata.total_species et un aperçu de organism_counts SANS charger
tout le fichier master_plant_db.json (qui contient des centaines de milliers
de gènes et peut faire plusieurs centaines de Mo).

Fonctionne parce que "metadata" est le premier champ du JSON, avant la
liste "genes" -- on lit juste assez de caractères pour couvrir metadata,
on s'arrête avant "genes", et on parse seulement ce petit bout.

Usage:
    python check_master_db_metadata.py [chemin_vers_master_plant_db.json]
"""
import json
import sys

path = sys.argv[1] if len(sys.argv) > 1 else r"C:\Downloads\IA\data\clean\master_plant_db.json"

# 5 Mo est largement suffisant pour "metadata" seule (même avec pas mal
# d'espèces) -- si jamais organism_counts contient encore le bug de
# l'ancien regex organism (des centaines de "fausses espèces"), cette
# limite ne suffira pas et on le saura tout de suite via le message
# d'erreur ci-dessous, ce qui est en soi une information utile.
CHUNK_SIZE = 5_000_000

with open(path, encoding="utf-8") as f:
    head = f.read(CHUNK_SIZE)

idx = head.find('"genes"')
if idx == -1:
    print(f"⚠ 'genes' non trouvé dans les {CHUNK_SIZE:,} premiers caractères.")
    print("  Soit le fichier n'a pas ce format, soit metadata est anormalement")
    print("  grosse (par exemple si organism_counts contient encore des")
    print("  centaines d'entrées à cause de l'ancien bug organism).")
    print("  Augmente CHUNK_SIZE dans ce script si besoin.")
    sys.exit(1)

meta_str = head[:idx].rstrip()
if meta_str.endswith(","):
    meta_str = meta_str[:-1]
meta_str += "}"

data = json.loads(meta_str)
metadata = data["metadata"]

print("total_species:", metadata.get("total_species"))
print("total_genes:", metadata.get("total_genes"))
print()
print("Aperçu organism_counts (10 premiers) :")
for name in list(metadata.get("organism_counts", {}).keys())[:10]:
    print(" -", name)
