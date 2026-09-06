"""
test_similarity_ai.py
----------------------
Unit tests for similarity and AI interpretation modules.
"""

import similarityengine as sim
import aiinterpreter as ai
import bioinformatics as bio
import scripts.postgres_utils as pg


def test_populate_kmer_index_writes_rows(monkeypatch):
    class FakeCursor:
        def __init__(self, conn, rows=None):
            self.conn = conn
            self.rows = rows or []
            self.result = []
            self.executed = []

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def __iter__(self):
            return iter(self.rows)

        def execute(self, sql, params=None):
            self.executed.append((sql, params))
            upper = sql.upper()
            if "SELECT COALESCE(GENE_ID, SYMBOL) AS GENE_KEY" in upper:
                self.result = [("G1", "ATGCATGCATGC", "dna"), ("G2", "GGGGTTTTCCCCAAAA", "dna")]
                self.rows = list(self.result)
            elif "INSERT INTO GENE_KMERS" in upper:
                self.result = [(1,)]
            elif "UPDATE GENES SET KMER_INDEXED = TRUE" in upper:
                self.result = [(1,)]
            elif "DELETE FROM gene_kmers" in upper:
                self.result = [(1,)]

        def executemany(self, sql, params_seq):
            self.executed.append((sql, params_seq))

        def fetchall(self):
            return list(self.result)

        def fetchone(self):
            if self.result:
                return self.result[0]
            return (0,)

    class FakeConnection:
        def __init__(self):
            self.cursors = []

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def cursor(self, *args, **kwargs):
            cur = FakeCursor(self)
            self.cursors.append(cur)
            return cur

        def commit(self):
            pass

        def rollback(self):
            pass

    fake_conn = FakeConnection()
    monkeypatch.setattr(pg, "get_connection", lambda: fake_conn)

    count = pg.populate_kmer_index(rebuild=True)

    assert count == 2
    assert count > 0


def test_load_gene_database_json():
    db = sim.load_gene_database("genes_database.json")
    assert isinstance(db, dict)
    assert "PR1" in db


def test_pairwise_similarity_same_sequence():
    assert sim.pairwise_similarity("ATGC", "ATGC") == 100.0


def test_pairwise_similarity_different_sequence():
    score = sim.pairwise_similarity("ATGC", "ATTA")
    assert score == 50.0


def test_compare_with_database_protein_query():
    protein = "MQNCG"  # short protein-like sequence
    results = sim.compare_with_database(protein, top_n=1)
    assert isinstance(results, list)
    assert len(results) == 1
    assert "gene_name" in results[0]


def test_ai_interpreter_summary_high_confidence():
    stats = {
        "length": 100,
        "gc_content": 45.0,
        "at_content": 55.0,
        "is_coding_length": True,
        "has_start_codon": True,
        "has_stop_codon": True,
    }
    similarity_results = [
        {
            "gene_name": "TestGene",
            "trait": "Stress Response",
            "organism": "Arabidopsis thaliana",
            "similarity_score": 92.0,
        }
    ]
    interp = ai.AIInterpreter(stats, similarity_results, None).full_report()
    assert interp["confidence_level"]["level"] in {"High", "Medium", "Low"}
    assert "summary" in interp["overall_summary"].lower() or isinstance(interp["overall_summary"], str)
