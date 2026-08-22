"""
V2 : corrige le bug "symbol" -> "common_name", et ajoute un diagnostic
sur les gènes origin=sequence_backed dont la séquence est en fait null
(anomalie détectée manuellement sur 3 gènes le 20/08/2026).

Usage :
    python search_lodging_candidates_v2.py path/to/chenopodium_quinoa_all_sources.json
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
NOISY = {"rigidite_tige:stem", "lignification:cell wall"}  # trop génériques


def collect_text(gene: dict) -> str:
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

    for field in ("description", "common_name", "traits"):
        add(gene.get(field))

    annotation = gene.get("annotation") or {}
    for field in ("go_terms", "kegg_pathways", "mapman", "tf_family"):
        add(annotation.get(field))

    add(gene.get("literature"))
    return " | ".join(parts)


def find_matches(text: str) -> list[str]:
    return sorted({f"{cat}:{kw.strip()}" for kw, cat in ALL_KEYWORDS if kw in text})


def has_real_sequence(gene: dict) -> bool:
    seq = gene.get("sequence") or {}
    return bool(seq.get("dna") or seq.get("rna") or seq.get("protein"))


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
        print("Usage: python search_lodging_candidates_v2.py <fichier.json>")
        sys.exit(1)

    genes = load_genes(Path(sys.argv[1]))

    # Diagnostic global sur l'anomalie sequence_backed sans séquence
    sb_genes = [g for g in genes.values() if isinstance(g, dict) and g.get("origin") == "sequence_backed"]
    sb_no_seq = [g for g in sb_genes if not has_real_sequence(g)]
    print(f"Diagnostic : {len(sb_no_seq)}/{len(sb_genes)} gènes 'sequence_backed' ont en fait une séquence null.")

    results = []
    for gene_id, gene in genes.items():
        organism = str(gene.get("organism", "")).lower()
        if organism and "quinoa" not in organism and "chenopodium" not in organism:
            continue
        text = collect_text(gene)
        matches = find_matches(text)
        specific = set(matches) - NOISY
        if specific and gene.get("origin") == "sequence_backed":
            results.append({
                "gene_id": gene_id,
                "common_name": gene.get("common_name", ""),
                "has_real_sequence": has_real_sequence(gene),
                "categories_matchees": ";".join(sorted({m.split(":")[0] for m in specific})),
                "mots_cles": ";".join(sorted(specific)),
            })

    out_path = Path("candidats_verse_quinoa_v2.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene_id", "common_name", "has_real_sequence", "categories_matchees", "mots_cles"])
        writer.writeheader()
        writer.writerows(results)

    print(f"\n{len(results)} gènes candidats (mots-clés spécifiques + sequence_backed).")
    print(f"Dont {sum(1 for r in results if r['has_real_sequence'])} avec une séquence réellement présente.")
    print(f"Résultats écrits dans : {out_path.resolve()}")


if __name__ == "__main__":
    main()
