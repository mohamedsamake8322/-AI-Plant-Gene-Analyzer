"""
inspect_json_structure.py

Affiche la structure du JSON (les clés de haut niveau et un aperçu de ce
qu'elles contiennent) sans jamais charger le fichier entier en mémoire.
Utile pour un gros fichier dont on ne connaît pas la forme exacte avant
d'écrire un script de comparaison.

Usage :
    python inspect_json_structure.py chemin/vers/fichier.json
"""

import sys
from pathlib import Path

import ijson


def main():
    if len(sys.argv) != 2:
        print("Usage: python inspect_json_structure.py <fichier.json>")
        sys.exit(1)

    path = Path(sys.argv[1])
    print(f"Inspection de {path} ...\n")

    with open(path, "rb") as f:
        depth = 0
        top_level_keys = []
        current_key = None
        events_seen = 0
        max_events = 4000  # large marge, mais borné pour rester rapide

        for prefix, event, value in ijson.parse(f):
            events_seen += 1
            if events_seen > max_events:
                print("... (arrêt après", max_events, "événements, structure top-level déjà claire)")
                break

            # On ne s'intéresse qu'au premier niveau de profondeur.
            if event == "start_map":
                depth += 1
            elif event == "end_map":
                depth -= 1
            elif event == "start_array":
                depth += 1
            elif event == "end_array":
                depth -= 1
            elif event == "map_key" and depth == 1:
                current_key = value
                print(f"Clé racine trouvée: {value!r}")
            elif depth == 2 and event in ("string", "number", "boolean", "null"):
                # Aperçu du contenu juste sous une clé racine scalaire.
                print(f"    -> valeur scalaire directe: {value!r}")
            elif depth == 2 and event == "map_key":
                # On est entré dans un sous-objet: montre ses premières clés.
                print(f"    -> sous-clé: {value!r}")

    print("\nSi les vrais enregistrements (gene_id/symbol) apparaissent sous une "
          "des sous-clés listées ci-dessus plutôt qu'à la racine, il faut cibler "
          "cette sous-clé avec ijson (ex: ijson.kvitems(f, 'records') ou "
          "ijson.items(f, 'genes.item')).")


if __name__ == "__main__":
    main()
