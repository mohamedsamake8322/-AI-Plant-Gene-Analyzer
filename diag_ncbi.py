#!/usr/bin/env python3
"""
diag_ncbi.py
-------------
Diagnostic INDEPENDANT du pipeline : interroge directement l'API NCBI
Entrez (esearch + efetch) pour verifier que la connexion et les requetes
fonctionnent, sans passer par run_pipeline.py / collect_multi_type.py.

Affiche l'URL exacte utilisee, le code HTTP, et le debut de la reponse
pour chaque etape -- de quoi identifier immediatement : probleme reseau,
cle API absente/invalide, limite de requetes (429), ou requete mal
construite (400).

Usage :
    python diag_ncbi.py
    python diag_ncbi.py --organism "Oryza sativa" --db nucleotide
    python diag_ncbi.py --api-key VOTRE_CLE  # si vous en avez une
"""
from __future__ import annotations

import argparse
import os
import sys
import time

try:
    import requests
except ImportError:
    print("Le module 'requests' n'est pas installe dans cet environnement.")
    print("Lancez : pip install requests")
    sys.exit(1)

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def show_response(label: str, resp: requests.Response) -> None:
    print(f"\n--- {label} ---")
    print(f"  URL appelee : {resp.url}")
    print(f"  Code HTTP   : {resp.status_code}")
    print(f"  Taille reponse : {len(resp.content)} octets")
    body = resp.text
    snippet = body[:800]
    print(f"  Debut de la reponse :\n{snippet}")
    if resp.status_code == 429:
        print("  >>> 429 = limite de requetes NCBI atteinte (rate limit).")
        print("      Sans cle API : 3 requetes/seconde max. Avec cle : 10/s.")
    if resp.status_code == 400:
        print("  >>> 400 = requete mal formee (parametre invalide, terme mal echappe...).")
    if "<ERROR>" in body or "error" in body.lower()[:300]:
        print("  >>> Le corps de la reponse contient le mot 'error' -- a lire attentivement ci-dessus.")


def run_probe(organism: str, db: str, api_key: str | None, email: str | None) -> None:
    print(f"\n{'=' * 70}")
    print(f"TEST : organism='{organism}'  db='{db}'")
    print(f"{'=' * 70}")

    if api_key:
        print(f"Cle API NCBI fournie (longueur {len(api_key)} caracteres).")
    else:
        print("AUCUNE cle API NCBI fournie (limite stricte : 3 requetes/seconde).")
    if email:
        print(f"Email fourni : {email}")
    else:
        print("AUCUN email fourni -- NCBI le recommande fortement, peut bloquer sans lui.")

    term = f"{organism}[Organism]"
    params = {
        "db": db,
        "term": term,
        "retmax": 5,
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key
    if email:
        params["email"] = email

    try:
        t0 = time.time()
        resp = requests.get(ESEARCH_URL, params=params, timeout=20)
        elapsed = time.time() - t0
        print(f"(requete esearch effectuee en {elapsed:.2f}s)")
        show_response("ESEARCH (recherche d'identifiants)", resp)
    except requests.RequestException as e:
        print(f"  ÉCHEC RESEAU sur esearch : {type(e).__name__}: {e}")
        print("  >>> Si vous voyez ceci, le probleme est reseau/proxy/pare-feu,")
        print("      pas le code Python du pipeline.")
        return

    # Essaie d'extraire des IDs pour tester efetch aussi
    ids = []
    try:
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        count = data.get("esearchresult", {}).get("count", "?")
        print(f"\n  Nombre total signale par NCBI (count) : {count}")
        print(f"  Identifiants recuperes (idlist) : {ids}")
    except Exception as e:
        print(f"  Impossible de parser la reponse JSON d'esearch : {e}")

    if not ids:
        print("\n  >>> Aucun identifiant retourne. Causes possibles :")
        print("      - le terme de recherche ne correspond a rien dans cette base")
        print("      - la reponse ci-dessus contient un message d'erreur a lire")
        print("      - throttling silencieux (rare, mais possible sans cle API)")
        return

    time.sleep(0.4)
    fetch_params = {
        "db": db,
        "id": ",".join(ids[:2]),
        "rettype": "fasta",
        "retmode": "text",
    }
    if api_key:
        fetch_params["api_key"] = api_key
    if email:
        fetch_params["email"] = email

    try:
        t0 = time.time()
        resp2 = requests.get(EFETCH_URL, params=fetch_params, timeout=20)
        elapsed = time.time() - t0
        print(f"\n(requete efetch effectuee en {elapsed:.2f}s)")
        show_response("EFETCH (recuperation des sequences)", resp2)
    except requests.RequestException as e:
        print(f"  ÉCHEC RESEAU sur efetch : {type(e).__name__}: {e}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--organism", default="Oryza sativa", help="Nom d'espece a tester")
    p.add_argument("--db", default="nucleotide", choices=["nucleotide", "protein"], help="Base NCBI a interroger")
    p.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"), help="Cle API NCBI (ou variable d'env NCBI_API_KEY)")
    p.add_argument("--email", default=os.environ.get("NCBI_EMAIL"), help="Email NCBI (ou variable d'env NCBI_EMAIL)")
    p.add_argument("--both-species", action="store_true", help="Teste aussi Chenopodium quinoa en plus de --organism")
    args = p.parse_args()

    run_probe(args.organism, args.db, args.api_key, args.email)

    if args.both_species:
        time.sleep(0.5)
        run_probe("Chenopodium quinoa", args.db, args.api_key, args.email)

    print(f"\n{'=' * 70}")
    print("FIN DU DIAGNOSTIC")
    print(f"{'=' * 70}")
    print("Copiez-collez TOUTE cette sortie -- c'est ce qui permettra de")
    print("determiner la cause exacte (reseau, cle API, requete, ou autre).")


if __name__ == "__main__":
    main()
