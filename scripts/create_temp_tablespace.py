#!/usr/bin/env python3
"""Crée un tablespace PostgreSQL 'temp_d' (dossier sur un autre disque) pour les fichiers temporaires.
Prérequis : le dossier existe, est VIDE, et le compte du service PostgreSQL peut y écrire.
Usage : python scripts\\create_temp_tablespace.py D:/pg_temp
Nécessite un rôle superutilisateur ; ne modifie aucune donnée."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect

if len(sys.argv) != 2:
    sys.exit("Usage : python scripts\\create_temp_tablespace.py D:/pg_temp")
loc = sys.argv[1].replace("\\", "/")
c = connect(); c.autocommit = True      # CREATE TABLESPACE refuse d'être dans une transaction
cur = c.cursor()
cur.execute("SELECT current_user, rolsuper FROM pg_roles WHERE rolname = current_user")
user, sup = cur.fetchone()
print("rôle :", user, "| superutilisateur :", sup)
if not sup:
    sys.exit("Ce rôle n'est pas superutilisateur : connectez-vous en 'postgres' (ou demandez-moi l'alternative).")
cur.execute("SELECT spcname, pg_tablespace_location(oid) FROM pg_tablespace WHERE spcname = 'temp_d'")
r = cur.fetchone()
if r:
    print("temp_d existe déjà :", r)
else:
    cur.execute(f"CREATE TABLESPACE temp_d LOCATION '{loc}'")
    print("tablespace temp_d créé sur", loc)
cur.execute(f'GRANT CREATE ON TABLESPACE temp_d TO "{user}"')
cur.execute("GRANT CREATE ON TABLESPACE temp_d TO PUBLIC")
print("droits accordés.")
