"""
Compte combien d'objets 'genes' sont reellement parsables dans
master_plant_db.json avant de tomber sur une eventuelle erreur
(fichier tronque / JSON invalide en fin de fichier).

Usage :
    python count_parsable_genes.py "C:\\Downloads\\IA\\Data\\clean\\master_plant_db.json"
"""

import sys
import ijson


def main(path: str):
    count = 0
    last_gene_id = None
    try:
        with open(path, "rb") as f:
            for gene in ijson.items(f, "genes.item"):
                count += 1
                last_gene_id = gene.get("gene_id")
                if count % 50000 == 0:
                    print(f"  ... {count} gènes lus jusqu'ici (dernier: {last_gene_id})")
    except Exception as e:
        print(f"\nERREUR après {count} gènes lus avec succès.")
        print(f"Dernier gene_id lu correctement : {last_gene_id}")
        print(f"Détail de l'erreur : {e}")
        return

    print(f"\nOK — {count} gènes lus intégralement, sans erreur.")
    print(f"Dernier gene_id : {last_gene_id}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_parsable_genes.py <chemin_master_plant_db.json>")
        sys.exit(1)
    main(sys.argv[1])
