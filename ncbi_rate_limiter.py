#!/usr/bin/env python3
"""
ncbi_rate_limiter.py
----------------------
Limiteur de debit PARTAGE ENTRE PROCESSUS pour les appels NCBI E-utilities
(esearch, esummary, efetch...). Resout le probleme identifie : avec
plusieurs workers (--workers > 1) ou plusieurs terminaux lances en
parallele (deux especes en meme temps), chaque processus respectait son
propre delai (NCBI_SLEEP) mais ce delai n'etait PAS partage -- le debit
cumule reel envoye a NCBI depassait alors la limite reelle de la cle API,
provoquant des coupures ("Remote end closed connection").

Principe : un compteur de requetes est stocke dans un petit fichier
SQLite partage sur disque. SQLite gere nativement le verrouillage
inter-processus (y compris sous Windows), donc N'IMPORTE QUEL processus
Python sur la machine -- qu'il vienne d'un ProcessPoolExecutor interne ou
d'un terminal PowerShell separe -- voit et respecte le MEME compteur.

La limite NCBI documentee s'applique globalement a TOUTE combinaison
d'appels E-utilities (esearch, esummary, efetch, elink...), pas
endpoint par endpoint -- ce limiteur est donc volontairement agnostique
de l'endpoint : chaque acquire() represente UNE requete E-utilities,
quelle qu'elle soit.

Usage dans le code de collecte (a la place de time.sleep(NCBI_SLEEP)) :

    from ncbi_rate_limiter import acquire
    acquire()  # bloque juste ce qu'il faut pour rester sous la limite,
               # PARTAGE avec tous les autres processus sur la machine
    # ... appel Entrez.esearch / esummary / efetch ...

Configuration (variables d'environnement, aucune valeur figee dans le
code) :
    NCBI_RATE_LIMIT_RPS   : requetes/seconde cible (defaut : 8 -- marge
                            de securite sous la limite NCBI de 10 rps
                            avec cle API, pour absorber la latence entre
                            "on decide d'envoyer" et "NCBI recoit").
    NCBI_RATE_LIMIT_DB    : chemin du fichier SQLite partage (defaut :
                            <racine_du_projet>/.ncbi_rate_limiter.sqlite3)
"""
from __future__ import annotations

import os
import sqlite3
import time
from pathlib import Path

_DEFAULT_DB_PATH = Path(__file__).resolve().parent / ".ncbi_rate_limiter.sqlite3"
_DB_PATH = Path(os.environ.get("NCBI_RATE_LIMIT_DB", str(_DEFAULT_DB_PATH)))
_DEFAULT_RPS = float(os.environ.get("NCBI_RATE_LIMIT_RPS", "8"))
_WINDOW_SECONDS = 1.0
_BUSY_TIMEOUT_MS = 30_000  # attente max sur un verrou SQLite avant erreur


def _connect() -> sqlite3.Connection:
    conn = sqlite3.connect(str(_DB_PATH), timeout=_BUSY_TIMEOUT_MS / 1000)
    conn.execute(f"PRAGMA busy_timeout = {_BUSY_TIMEOUT_MS}")
    conn.execute(
        "CREATE TABLE IF NOT EXISTS requests (ts REAL NOT NULL)"
    )
    return conn


def acquire(max_per_second: float | None = None, window_seconds: float = _WINDOW_SECONDS) -> None:
    """
    Bloque jusqu'a ce qu'il soit sur de pouvoir envoyer UNE requete NCBI
    sans depasser max_per_second, en comptant TOUTES les requetes de TOUS
    les processus utilisant le meme fichier de base (meme cle API,
    partagee entre workers/terminaux).

    A appeler juste avant chaque appel Entrez.esearch/esummary/efetch/etc.
    """
    rps = max_per_second if max_per_second is not None else _DEFAULT_RPS
    if rps <= 0:
        raise ValueError("max_per_second doit etre > 0")

    while True:
        conn = _connect()
        try:
            # BEGIN IMMEDIATE : prend le verrou d'ecriture tout de suite,
            # pour que deux processus ne lisent jamais un etat perime
            # avant d'ecrire (evite la fenetre de "race condition" d'un
            # simple SELECT puis INSERT sans transaction explicite).
            conn.execute("BEGIN IMMEDIATE")
            now = time.time()
            cutoff = now - window_seconds

            conn.execute("DELETE FROM requests WHERE ts < ?", (cutoff,))
            count = conn.execute("SELECT COUNT(*) FROM requests").fetchone()[0]

            if count < rps:
                conn.execute("INSERT INTO requests (ts) VALUES (?)", (now,))
                conn.commit()
                return  # autorise -- la requete peut partir maintenant

            # Quota atteint pour cette fenetre : calcule combien de temps
            # attendre avant que la plus ancienne requete "expire" de la
            # fenetre glissante, puis relache le verrou pour ne pas
            # bloquer les autres processus pendant l'attente.
            oldest = conn.execute("SELECT MIN(ts) FROM requests").fetchone()[0]
            conn.rollback()
            wait = max(0.0, (oldest + window_seconds) - now) + 0.01  # marge
        finally:
            conn.close()

        time.sleep(wait)


def reset() -> None:
    """Vide le compteur -- utile pour les tests, jamais necessaire en usage normal."""
    conn = _connect()
    try:
        conn.execute("DELETE FROM requests")
        conn.commit()
    finally:
        conn.close()


def current_db_path() -> Path:
    """Utile pour verifier/afficher quel fichier est reellement partage."""
    return _DB_PATH


if __name__ == "__main__":
    # Auto-test rapide, mono-processus : verifie que le debit reel
    # n'excede jamais max_per_second sur une fenetre glissante.
    import sys

    n_requests = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rps = float(sys.argv[2]) if len(sys.argv) > 2 else 8

    print(f"Fichier partage : {current_db_path()}")
    print(f"Test : {n_requests} requetes a {rps} req/s max...")
    reset()

    t0 = time.time()
    timestamps = []
    for i in range(n_requests):
        acquire(max_per_second=rps)
        timestamps.append(time.time())
        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{n_requests} requetes envoyees, "
                  f"debit instantane sur la derniere seconde : "
                  f"{sum(1 for t in timestamps if t > timestamps[-1] - 1.0)} req/s")

    elapsed = timestamps[-1] - t0
    # Le vrai critere de correction n'est PAS la moyenne globale (biaisee
    # par la salve initiale autorisee au demarrage d'une fenetre vide) --
    # c'est le nombre maximal de requetes tombant dans N'IMPORTE QUELLE
    # fenetre glissante de 1 seconde sur tout l'historique.
    max_in_any_window = max(
        sum(1 for t in timestamps if t2 - 1.0 < t <= t2) for t2 in timestamps
    )
    print(f"\nTermine en {elapsed:.2f}s | debit moyen sur toute la duree : "
          f"{n_requests / elapsed if elapsed > 0 else float('inf'):.2f} req/s (indicatif seulement)")
    print(f"Debit maximal observe sur une fenetre glissante de 1s : {max_in_any_window} req/s "
          f"(cible : {rps} req/s)")
    if max_in_any_window > rps:
        print("⚠️  Une fenetre glissante depasse la cible -- a examiner.")
    else:
        print("✓ Debit respecte sur toute fenetre glissante de 1s.")
