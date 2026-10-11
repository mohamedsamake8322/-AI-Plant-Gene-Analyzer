"""
Sonde legere la structure de master_plant_db.json SANS le charger entierement
en memoire (fichier attendu ~1.5 Go). Utilise ijson pour lire en streaming.

Etape 1 avant d'ecrire un script d'audit de completude complet : on a besoin
de savoir a quoi ressemble le JSON (liste plate ? dict par espece ? autre ?)
avant de pouvoir parcourir les 300k+ genes correctement.

Usage :
    pip install ijson
    python probe_structure.py "C:\\Downloads\\IA\\Data\\clean\\master_plant_db.json"
"""

import sys
import json
import ijson


def probe(path: str):
    print(f"Analyse de : {path}\n")

    # 1. Regarder juste le tout debut du fichier pour voir si c'est un objet {} ou une liste []
    with open(path, "rb") as f:
        head = f.read(500).decode("utf-8", errors="replace")
    print("--- 500 premiers caracteres du fichier ---")
    print(head)
    print("--- fin extrait ---\n")

    is_list_root = head.lstrip().startswith("[")
    is_dict_root = head.lstrip().startswith("{")

    if is_dict_root:
        print("Racine = objet JSON {}. Clés de premier niveau :")
        with open(path, "rb") as f:
            keys = []
            for prefix, event, value in ijson.parse(f):
                if prefix == "" and event == "map_key":
                    keys.append(value)
                if len(keys) >= 30:  # sécurité, on arrête si trop de clés
                    break
        print(keys)
        print()

        # Pour chaque clé candidate plausible, essayer de voir si c'est une liste de genes
        for key in keys:
            try:
                with open(path, "rb") as f:
                    items = ijson.items(f, f"{key}.item")
                    first = next(items, None)
                if first is not None:
                    kind = type(first).__name__
                    sample_keys = list(first.keys())[:15] if isinstance(first, dict) else None
                    print(f"  '{key}' -> liste, premier element de type {kind}")
                    if sample_keys:
                        print(f"    clés du 1er élément : {sample_keys}")
            except Exception as e:
                print(f"  '{key}' -> pas une liste exploitable ({e})")

    elif is_list_root:
        print("Racine = liste JSON []. Premier élément :")
        with open(path, "rb") as f:
            items = ijson.items(f, "item")
            first = next(items, None)
        if isinstance(first, dict):
            print(f"  type: dict, clés : {list(first.keys())}")
            print(f"  contenu (tronqué) : {json.dumps(first, ensure_ascii=False)[:800]}")
        else:
            print(f"  type: {type(first).__name__}, valeur (tronquée) : {str(first)[:300]}")
    else:
        print("Racine non reconnue (ni [ ni {). Premiers caractères ci-dessus à inspecter à la main.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python probe_structure.py <chemin_master_plant_db.json>")
        sys.exit(1)
    probe(sys.argv[1])
