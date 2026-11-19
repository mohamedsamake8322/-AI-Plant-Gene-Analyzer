import os
import uuid

import pytest

from core_engines.bioinfo_db import connect
from scripts.rebuild_pipeline_sequences import _prune_absent_proteins


def test_interrupted_prune_rolls_back_archive_and_deletes():
    database = os.getenv("BIOINFO_PIPELINE_DBNAME", "")
    if not database.startswith("bioinfo_pipeline_restore_check_"):
        pytest.skip(
            "Set BIOINFO_PIPELINE_DBNAME to a disposable restore-check database."
        )

    assemblage = f"__prune_atomicity_{uuid.uuid4().hex}"
    obsolete = f"OLD_{uuid.uuid4().hex}"
    retained = f"NEW_{uuid.uuid4().hex}"

    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            "INSERT INTO dim_protein (protein_accession, assemblage) VALUES (%s, %s), (%s, %s)",
            (obsolete, assemblage, retained, assemblage),
        )
        cursor.execute(
            """
            INSERT INTO fact_functional_annotation
                (assemblage, protein_accession, source)
            VALUES (%s, %s, 'prune_test'), (%s, %s, 'prune_test')
            """,
            (assemblage, obsolete, assemblage, retained),
        )

    with pytest.raises(RuntimeError, match="simulated interruption"):
        with connect() as connection, connection.cursor() as cursor:
            _prune_absent_proteins(cursor, assemblage, {retained})
            raise RuntimeError("simulated interruption before transaction commit")

    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            "SELECT protein_accession FROM dim_protein WHERE assemblage=%s ORDER BY protein_accession",
            (assemblage,),
        )
        assert cursor.fetchall() == [(item,) for item in sorted((obsolete, retained))]
        cursor.execute(
            """
            SELECT protein_accession FROM fact_functional_annotation
            WHERE assemblage=%s ORDER BY protein_accession
            """,
            (assemblage,),
        )
        assert cursor.fetchall() == [(item,) for item in sorted((obsolete, retained))]
        cursor.execute(
            "SELECT count(*) FROM dim_protein_archive WHERE assemblage=%s",
            (assemblage,),
        )
        assert cursor.fetchone()[0] == 0
        cursor.execute(
            "SELECT count(*) FROM fact_functional_annotation_archive WHERE assemblage=%s",
            (assemblage,),
        )
        assert cursor.fetchone()[0] == 0

    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            "DELETE FROM fact_functional_annotation WHERE assemblage=%s",
            (assemblage,),
        )
        cursor.execute(
            "DELETE FROM dim_protein WHERE assemblage=%s",
            (assemblage,),
        )
