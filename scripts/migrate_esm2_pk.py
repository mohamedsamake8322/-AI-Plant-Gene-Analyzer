r"""Aligne la clé primaire de esm2_embeddings sur la conception : (assembly_accession, protein_accession).

Contexte : la clé actuelle est (protein_accession) seule, ce qui interdit de stocker deux variétés (K326, TN90)
qui partagent des protéines identiques. Toutes les requêtes de l'application filtrent déjà par assemblage.

    python scripts\migrate_esm2_pk.py            # simulation (ROLLBACK)
    python scripts\migrate_esm2_pk.py --apply    # applique (une transaction)

Aucune donnée n'est modifiée : seule la contrainte change. Rejouable (déjà à jour -> rien à faire).
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.getcwd())
from core_engines.bioinfo_db import connect  # noqa: E402

WANTED = "PRIMARY KEY (assembly_accession, protein_accession)"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    conn = connect()
    conn.autocommit = False
    cur = conn.cursor()
    try:
        cur.execute("SET LOCAL lock_timeout = '15s'")
        cur.execute("SET LOCAL maintenance_work_mem = '128MB'")
        cur.execute("SELECT conname, pg_get_constraintdef(oid) FROM pg_constraint "
                    "WHERE conrelid = 'esm2_embeddings'::regclass AND contype = 'p'")
        current = cur.fetchone()
        print("Clé primaire actuelle :", current)
        if current and current[1] == WANTED:
            print("Déjà à jour : rien à faire.")
            conn.rollback()
            return 0
        cur.execute("SELECT conname, conrelid::regclass FROM pg_constraint WHERE confrelid = 'esm2_embeddings'::regclass")
        fks = cur.fetchall()
        if fks:
            raise RuntimeError(f"Des clés étrangères dépendent de la table : {fks}. Migration refusée.")
        cur.execute("SELECT count(*), count(DISTINCT (assembly_accession, protein_accession)) FROM esm2_embeddings")
        total, distinct = cur.fetchone()
        print(f"Lignes : {total:,} | couples (assemblage, accession) distincts : {distinct:,}")
        if total != distinct:
            raise RuntimeError("Des doublons (assemblage, accession) existent : migration impossible.")
        start = time.time()
        if current:
            cur.execute(f'ALTER TABLE esm2_embeddings DROP CONSTRAINT "{current[0]}"')
        cur.execute("ALTER TABLE esm2_embeddings ADD CONSTRAINT esm2_embeddings_pkey "
                    "PRIMARY KEY (assembly_accession, protein_accession)")
        print(f"Nouvelle clé créée en {time.time() - start:.1f} s")
        cur.execute("SELECT conname, pg_get_constraintdef(oid) FROM pg_constraint "
                    "WHERE conrelid = 'esm2_embeddings'::regclass AND contype = 'p'")
        print("Clé primaire après    :", cur.fetchone())
        cur.execute("SELECT count(*) FROM esm2_embeddings")
        after = cur.fetchone()[0]
        if after != total:
            raise RuntimeError(f"Nombre de lignes modifié ({total} -> {after}) : annulation.")
        if args.apply:
            conn.commit()
            print("COMMIT effectué.")
        else:
            conn.rollback()
            print("Simulation : ROLLBACK. Relancez avec --apply.")
        return 0
    except Exception as error:  # noqa: BLE001
        conn.rollback()
        print("ERREUR, transaction annulée :", error)
        return 1
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
