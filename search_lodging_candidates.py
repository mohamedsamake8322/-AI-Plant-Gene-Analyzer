"""
Recherche de gènes candidats pour la verse (lodging) chez le quinoa.

Cible deux axes biologiques :
  1. Lignification / rigidité de la tige (paroi cellulaire secondaire)
  2. Régulation hormonale de l'élongation/rigidité de la tige (gibbérellines, DELLA, brassinostéroïdes)

Usage :
    python search_lodging_candidates.py path/to/master_plant_db.json
    python search_lodging_candidates.py path/to/chenopodium_quinoa_all_sources.json

Sortie : candidats_verse_quinoa.csv dans le dossier courant.
"""

import json
import sys
import csv
from pathlib import Path

KEYWORDS = {
    "lignification": [
        "lignin", "lignification", "cell wall", "secondary cell wall",
        "cellulose synthase", "laccase", "peroxidase", "cinnamyl alcohol",
        "phenylpropanoid", "cad ", "4cl", "ccoaomt", "comt",
    ],
    "rigidite_tige": [
        "stem", "culm", "stalk", "lodging", "verse", "internode",
        "mechanical strength", "stem strength",
    ],
    "hormonal": [
        "gibberellin", "gibberellic", "della", "ga20ox", "ga3ox", "ga2ox",
        "brassinosteroid", "bri1",
    ],
}
ALL_KEYWORDS = [(kw, cat) for cat, kws in KEYWORDS.items() for kw in kws]


def collect_text(gene: dict) -> str:
    """Concatène tous les champs textuels pertinents d'un gène en minuscules."""
    parts = []

    def add(v):
        if isinstance(v, str):
            parts.append(v.lower())
        elif isinstance(v, dict):
            for x in v.values():
                add(x)
        elif isinstance(v, list):
            for x in v:
                add(x)

    for field in ("description", "symbol", "traits"):
        add(gene.get(field))

    annotation = gene.get("annotation") or {}
    for field in ("go_terms", "kegg_pathways", "mapman", "tf_family"):
        add(annotation.get(field))

    add(gene.get("literature"))
    return " | ".join(parts)


def find_matches(text: str) -> list[str]:
    return sorted({f"{cat}:{kw.strip()}" for kw, cat in ALL_KEYWORDS if kw in text})


def load_genes(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    if isinstance(raw, dict) and all(isinstance(v, dict) for v in raw.values()):
        return raw
    if isinstance(raw, dict) and "genes" in raw:
        raw = raw["genes"]
    if isinstance(raw, list):
        return {g.get("gene_id", str(i)): g for i, g in enumerate(raw)}
    raise ValueError("Format de fichier non reconnu")


def main():
    if len(sys.argv) != 2:
        print("Usage: python search_lodging_candidates.py <fichier.json>")
        sys.exit(1)

    path = Path(sys.argv[1])
    genes = load_genes(path)

    results = []
    for gene_id, gene in genes.items():
        organism = str(gene.get("organism", "")).lower()
        if organism and "quinoa" not in organism and "chenopodium" not in organism:
            continue  # ignore si ce n'est manifestement pas du quinoa
        text = collect_text(gene)
        matches = find_matches(text)
        if matches:
            results.append({
                "gene_id": gene_id,
                "symbol": gene.get("symbol", ""),
                "origin": gene.get("origin", ""),
                "categories_matchees": ";".join(sorted({m.split(":")[0] for m in matches})),
                "mots_cles": ";".join(matches),
            })

    out_path = Path("candidats_verse_quinoa.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene_id", "symbol", "origin", "categories_matchees", "mots_cles"])
        writer.writeheader()
        writer.writerows(results)

    print(f"{len(results)} gènes candidats trouvés sur {len(genes)} gènes scannés.")
    print(f"Résultats écrits dans : {out_path.resolve()}")


if __name__ == "__main__":
    main()
