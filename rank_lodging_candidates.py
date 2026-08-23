"""
Priorise les 213 candidats verse/quinoa vers une short-list de 15-30 gènes.

Score composite :
  - Spécificité des mots-clés (tier A = signal fort, tier B = signal faible)
  - Diversité des catégories matchées (lignification / rigidite_tige / hormonal)
  - Présence de publications (PMID) déjà collectées dans "literature"

Usage :
    python rank_lodging_candidates.py path/to/chenopodium_quinoa_all_sources.json
"""

import json
import sys
import csv
from pathlib import Path

TIER_A = {  # signal biologique fort, spécifique à la voie lignification/hormonal
    "lignin", "lignification", "laccase", "cellulose synthase",
    "secondary cell wall", "phenylpropanoid", "cinnamyl alcohol",
    "4cl", "comt", "ccoaomt", "della", "gibberellin", "gibberellic",
    "ga20ox", "ga3ox", "ga2ox", "brassinosteroid", "bri1",
    "verse", "lodging", "stalk", "mechanical strength", "stem strength",
}
TIER_B = {"peroxidase", "culm", "internode"}  # signal plus faible/générique
NOISY = {"stem", "cell wall"}  # déjà exclus en amont, gardé ici par sécurité

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


def find_matches(text: str) -> set[str]:
    return {f"{cat}:{kw.strip()}" for kw, cat in ALL_KEYWORDS if kw in text}


def has_real_sequence(gene: dict) -> bool:
    seq = gene.get("sequence") or {}
    return bool(seq.get("dna") or seq.get("rna") or seq.get("protein"))


def count_pmids(gene: dict) -> int:
    lit = gene.get("literature")
    if not lit:
        return 0
    if isinstance(lit, list):
        return len([x for x in lit if isinstance(x, dict) and (x.get("pmid") or x.get("id"))])
    return 0


def score_gene(matches: set[str], n_pmids: int) -> tuple[int, int, int]:
    """Retourne (score_total, n_categories, n_tier_a) pour tri."""
    kws = {m.split(":", 1)[1] for m in matches}
    n_tier_a = len(kws & TIER_A)
    n_tier_b = len(kws & TIER_B)
    n_categories = len({m.split(":")[0] for m in matches})

    score = n_tier_a * 3 + n_tier_b * 1
    score += max(0, n_categories - 1) * 2  # bonus diversité (2e/3e catégorie)
    score += min(n_pmids, 3) * 2           # bonus littérature, plafonné
    return score, n_categories, n_tier_a


def load_genes(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    if isinstance(raw, dict) and "genes" in raw:
        raw = raw["genes"]
    if isinstance(raw, dict) and all(isinstance(v, dict) for v in raw.values()):
        return raw
    if isinstance(raw, list):
        return {g.get("gene_id", str(i)): g for i, g in enumerate(raw)}
    raise ValueError("Format de fichier non reconnu")


def main():
    if len(sys.argv) != 2:
        print("Usage: python rank_lodging_candidates.py <fichier.json>")
        sys.exit(1)

    genes = load_genes(Path(sys.argv[1]))

    NOISY_TAGGED = {"rigidite_tige:stem", "lignification:cell wall"}
    rows = []
    for gene_id, gene in genes.items():
        organism = str(gene.get("organism", "")).lower()
        if organism and "quinoa" not in organism and "chenopodium" not in organism:
            continue
        if gene.get("origin") != "sequence_backed":
            continue
        text = collect_text(gene)
        matches = find_matches(text) - NOISY_TAGGED
        if not matches or not has_real_sequence(gene):
            continue

        n_pmids = count_pmids(gene)
        score, n_cat, n_tier_a = score_gene(matches, n_pmids)
        rows.append({
            "gene_id": gene_id,
            "common_name": gene.get("common_name", ""),
            "score": score,
            "n_categories": n_cat,
            "n_tier_a_keywords": n_tier_a,
            "n_pmids": n_pmids,
            "categories": ";".join(sorted({m.split(":")[0] for m in matches})),
            "mots_cles": ";".join(sorted(matches)),
        })

    rows.sort(key=lambda r: r["score"], reverse=True)

    out_path = Path("candidats_verse_quinoa_ranked.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene_id", "common_name", "score", "n_categories",
            "n_tier_a_keywords", "n_pmids", "categories", "mots_cles",
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"{len(rows)} gènes scorés et triés.")
    print(f"Résultats écrits dans : {out_path.resolve()}\n")
    print("Top 30 :")
    for r in rows[:30]:
        print(f"  [{r['score']:2d}] {r['gene_id']} — {r['common_name']} "
              f"(cat: {r['categories']}, pmids: {r['n_pmids']})")


if __name__ == "__main__":
    main()
