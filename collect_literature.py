"""
Recherche PubMed ciblee pour une liste de genes candidats (pas toute la base).

Contrairement a une collecte massive sur 320k genes (trop lent, et trop
ambigu a cause de noms de genes generiques comme PER60/LAC17), ce script
ne cherche que sur les genes que tu lui donnes en entree (ex: les 213
candidats verse/quinoa, ou mieux, les 15-30 finalistes).

Il ne modifie PAS master_plant_db.json. Il ecrit un CSV a part, que tu
relis et valides toi-meme avant de decider quels PMIDs sont pertinents.

IMPORTANT - a lire avant de faire confiance aux resultats :
  Les noms de genes (PER60, LAC17, CESA8...) sont utilises dans plein
  d'especes differentes. Une requete "LAC17 AND quinoa" peut renvoyer 0
  resultat pertinent, ou au contraire un faux positif si "quinoa" apparait
  incidemment dans un article qui ne parle pas du gene. VERIFIE chaque
  PMID manuellement (titre + resume) avant de le citer dans le memoire.
  Ce script sert a degrossir la recherche, pas a la remplacer.

Usage :
    python collect_literature.py candidats_verse_quinoa_ranked.csv --top 30
    python collect_literature.py candidats_verse_quinoa_ranked.csv --top 30 --species "Chenopodium quinoa"
    python collect_literature.py candidats_verse_quinoa_ranked.csv --top 30 --api-key TA_CLE_NCBI

Obtenir une cle API NCBI gratuite (optionnel mais recommande, passe de
3 a 10 requetes/seconde) : https://www.ncbi.nlm.nih.gov/account/settings/
"""

import argparse
import csv
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def clean_gene_symbol(common_name: str) -> str:
    """Extrait le symbole court avant les ':' (ex: 'LAC17: Laccase-17' -> 'LAC17')."""
    if not common_name:
        return ""
    return common_name.split(":")[0].strip()


def http_get_xml(url: str, params: dict, retries: int = 3):
    query = urllib.parse.urlencode(params)
    full_url = f"{url}?{query}"
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(full_url, timeout=15) as resp:
                return ET.fromstring(resp.read())
        except Exception as e:
            if attempt == retries:
                print(f"  [erreur reseau apres {retries} essais] {e}", file=sys.stderr)
                return None
            time.sleep(2 * attempt)


def search_pubmed(query: str, api_key: str | None, max_results: int = 5):
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "xml",
        "retmax": str(max_results),
        "sort": "relevance",
    }
    if api_key:
        params["api_key"] = api_key
    root = http_get_xml(ESEARCH_URL, params)
    if root is None:
        return []
    return [el.text for el in root.findall(".//Id")]


def fetch_summaries(pmids: list[str], api_key: str | None):
    if not pmids:
        return {}
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    if api_key:
        params["api_key"] = api_key
    root = http_get_xml(ESUMMARY_URL, params)
    if root is None:
        return {}
    out = {}
    for doc in root.findall(".//DocSum"):
        pmid = doc.findtext("Id")
        title = ""
        pubdate = ""
        for item in doc.findall("Item"):
            if item.get("Name") == "Title":
                title = item.text or ""
            elif item.get("Name") == "PubDate":
                pubdate = item.text or ""
        out[pmid] = {"title": title, "pubdate": pubdate}
    return out


def load_candidates(path: str, top: int) -> list[dict]:
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows[:top]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates_csv", help="CSV produit par rank_lodging_candidates.py")
    ap.add_argument("--top", type=int, default=30, help="Nombre de gènes à chercher (défaut: 30)")
    ap.add_argument("--species", default="Chenopodium quinoa", help="Nom d'espèce à ajouter à la requête")
    ap.add_argument("--api-key", default=None, help="Clé API NCBI (optionnelle, augmente le débit autorisé)")
    ap.add_argument("--max-per-gene", type=int, default=5, help="Nombre max de PMIDs candidats par gène")
    ap.add_argument("--out", default="literature_candidats_top.csv")
    args = ap.parse_args()

    delay = 0.11 if args.api_key else 0.34  # respecte les limites NCBI (10/s ou 3/s)

    candidates = load_candidates(args.candidates_csv, args.top)
    print(f"{len(candidates)} gènes à chercher sur PubMed (espèce: {args.species})\n")

    out_rows = []
    for i, row in enumerate(candidates, 1):
        gene_id = row.get("gene_id", "")
        symbol = clean_gene_symbol(row.get("common_name", ""))
        if not symbol:
            print(f"[{i}/{len(candidates)}] {gene_id} — pas de common_name exploitable, skip")
            continue

        query = f'({symbol}[Title/Abstract]) AND ("{args.species}"[Title/Abstract] OR quinoa[Title/Abstract])'
        print(f"[{i}/{len(candidates)}] {gene_id} ({symbol}) — requête : {query}")

        pmids = search_pubmed(query, args.api_key, args.max_per_gene)
        time.sleep(delay)

        summaries = fetch_summaries(pmids, args.api_key) if pmids else {}
        time.sleep(delay)

        if not pmids:
            out_rows.append({
                "gene_id": gene_id, "symbol": symbol, "pmid": "",
                "title": "", "pubdate": "", "a_verifier_manuellement": "aucun résultat",
            })
            print("    -> aucun résultat")
        else:
            for pmid in pmids:
                s = summaries.get(pmid, {})
                out_rows.append({
                    "gene_id": gene_id, "symbol": symbol, "pmid": pmid,
                    "title": s.get("title", ""), "pubdate": s.get("pubdate", ""),
                    "a_verifier_manuellement": "OUI - vérifier que l'article parle bien de ce gène chez quinoa",
                })
            print(f"    -> {len(pmids)} PMID(s) candidat(s), À VÉRIFIER manuellement")

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene_id", "symbol", "pmid", "title", "pubdate", "a_verifier_manuellement",
        ])
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"\nTerminé. Résultats bruts (non validés) écrits dans : {args.out}")
    print("Prochaine étape : ouvrir le CSV, lire chaque titre, et cocher/supprimer")
    print("les PMIDs qui ne concernent pas réellement ce gène chez le quinoa.")


if __name__ == "__main__":
    main()
