r"""Patch de scripts\load_tn90_embeddings.py : accepte les accessions partagées avec une autre variété
(même séquence), et exige la nouvelle clé primaire (assemblage, accession).

    python scripts\patch_tn90_loader.py            # applique (copie .bak créée), idempotent
"""
import ast
import shutil
import sys
from pathlib import Path

TARGET = Path(__file__).resolve().parent / "load_tn90_embeddings.py"

OLD_COLLISION = '''    cursor.execute(
        """
        SELECT COUNT(*)
        FROM esm2_embeddings
        WHERE protein_accession = ANY(%s)
          AND NOT (species = %s AND assembly_accession = %s)
        """,
        (list(index_accessions), SPECIES, ASSEMBLY),
    )
    accession_collisions = cursor.fetchone()[0]
    if accession_collisions:
        raise RuntimeError(
            f"{accession_collisions} TN90 accessions already belong to another embedding row."
        )
'''

NEW_COLLISION = '''    cursor.execute(
        """
        SELECT pg_get_constraintdef(oid)
        FROM pg_constraint
        WHERE conrelid = 'esm2_embeddings'::regclass AND contype = 'p'
        """
    )
    primary_key = cursor.fetchone()
    if not primary_key or primary_key[0] != "PRIMARY KEY (assembly_accession, protein_accession)":
        raise RuntimeError(
            "esm2_embeddings primary key is "
            f"{primary_key[0] if primary_key else 'missing'}; "
            "run scripts/migrate_esm2_pk.py --apply first."
        )
    # The same RefSeq accession can exist in two cultivars (same sequence). That is legitimate as long as
    # the sequence MD5 is identical; a different MD5 under the same accession would be an identity error.
    accession_list = sorted(index_accessions)
    md5_by_accession = {row["protein_accession"]: row["sequence_md5"] for row in index_rows}
    cursor.execute(
        """
        SELECT COUNT(*),
               COUNT(*) FILTER (WHERE embedding.sequence_md5 IS DISTINCT FROM expected.md5)
        FROM esm2_embeddings AS embedding
        JOIN unnest(%s::text[], %s::text[]) AS expected(accession, md5)
          ON expected.accession = embedding.protein_accession
        WHERE NOT (embedding.species = %s AND embedding.assembly_accession = %s)
        """,
        (accession_list, [md5_by_accession[a] for a in accession_list], SPECIES, ASSEMBLY),
    )
    shared_with_other_species, shared_md5_mismatches = cursor.fetchone()
    if shared_md5_mismatches:
        raise RuntimeError(
            f"{shared_md5_mismatches} TN90 accessions exist in another species with a DIFFERENT sequence MD5."
        )
'''

OLD_RETURN = '''        "idempotent_existing_load": existing_tn90_count == EXPECTED_COUNT,
    }
'''
NEW_RETURN = '''        "idempotent_existing_load": existing_tn90_count == EXPECTED_COUNT,
        "accessions_shared_with_other_species": shared_with_other_species,
        "shared_accessions_md5_mismatches": shared_md5_mismatches,
    }
'''


def main() -> int:
    raw = TARGET.read_bytes().decode("utf-8")
    eol = "\r\n" if "\r\n" in raw else "\n"
    text = raw.replace("\r\n", "\n")
    if "migrate_esm2_pk.py --apply first" in text:
        print("Déjà patché : rien à faire.")
        return 0
    for label, old in (("bloc de collision", OLD_COLLISION), ("retour du préflight", OLD_RETURN)):
        if text.count(old) != 1:
            print(f"ERREUR : {label} introuvable exactement une fois ({text.count(old)}). Fichier non modifié.")
            return 1
    patched = text.replace(OLD_COLLISION, NEW_COLLISION, 1).replace(OLD_RETURN, NEW_RETURN, 1)
    try:
        ast.parse(patched)
    except SyntaxError as error:
        print("ERREUR de syntaxe après patch, fichier non modifié :", error)
        return 1
    backup = TARGET.with_suffix(".py.bak")
    shutil.copyfile(TARGET, backup)
    TARGET.write_bytes(patched.replace("\n", eol).encode("utf-8"))
    print(f"Patch appliqué. Sauvegarde : {backup.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
