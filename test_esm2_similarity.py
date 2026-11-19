import numpy as np
import pytest

import core_engines.esm2_similarity as esm2_similarity
from core_engines.esm2_similarity import EMBED_DIM, EmbeddingIndex, _decode_embedding
from core_engines.distance_engine import distance_matrix


def test_embedding_cosine_matrix_preserves_distance_result_shape(monkeypatch):
    vectors = {
        "MKW": np.array([1.0, 0.0], dtype=np.float32),
        "GQE": np.array([0.0, 1.0], dtype=np.float32),
        "MKA": np.array([1.0, 1.0], dtype=np.float32),
    }
    monkeypatch.setattr(esm2_similarity, "embed_sequence", vectors.__getitem__)

    result = distance_matrix(
        [
            {"name": "A", "sequence": " MKW "},
            {"name": "B", "sequence": "GQE"},
            {"name": "C", "sequence": "MKA"},
        ],
        method="embedding_cosine_distance",
        seq_type="protein",
    )

    assert result["method"] == "embedding_cosine_distance"
    assert result["sequence_names"] == ["A", "B", "C"]
    assert result["aligned_sequences"] == ["MKW", "GQE", "MKA"]
    assert result["distance_matrix"][0][0] == 0.0
    assert result["distance_matrix"][0][1] == pytest.approx(1.0)
    assert result["distance_matrix"][0][2] == pytest.approx(1 - 1 / np.sqrt(2))
    assert np.allclose(result["distance_matrix"], np.asarray(result["distance_matrix"]).T)
    assert np.asarray(result["matrix_object"]).shape == (3, 3)


def test_embedding_cosine_matrix_rejects_non_protein_mode():
    with pytest.raises(ValueError, match="requires protein sequences"):
        distance_matrix(
            [{"name": "A", "sequence": "MKW"}, {"name": "B", "sequence": "GGA"}],
            method="embedding_cosine_distance",
            seq_type="dna",
        )


def test_embedding_cosine_matrix_rejects_dna_even_when_caller_sets_protein_mode():
    with pytest.raises(ValueError, match="requires amino-acid sequences"):
        distance_matrix(
            [{"name": "A", "sequence": "ATGCGT"}, {"name": "B", "sequence": "ATGAAA"}],
            method="embedding_cosine_distance",
            seq_type="protein",
        )


def test_embedding_index_ranks_and_filters_species():
    index = EmbeddingIndex(
        ["rice-hit", "maize-hit"],
        np.array(["oryza_sativa", "zea_mays"]),
        np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32),
    )

    hits = index.search(np.array([0.9, 0.1]), top_n=2, species_filter=["zea_mays"])

    assert [hit["protein_accession"] for hit in hits] == ["maize-hit"]
    assert hits[0]["esm2_cosine_similarity"] == pytest.approx(0.1 / np.sqrt(0.82))


def test_embedding_decoder_accepts_bytea_and_legacy_real_array():
    vector = np.arange(EMBED_DIM, dtype=np.float32)

    assert np.array_equal(_decode_embedding(vector.tobytes()), vector)
    assert np.array_equal(_decode_embedding(vector.tolist()), vector)
    with pytest.raises(ValueError, match="Expected an ESM-2 vector"):
        _decode_embedding([1.0, 2.0])


def test_pfam_annotations_are_looked_up_by_assembly_and_accession():
    class FakeCursor:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def execute(self, query, params):
            self.query = query
            self.params = params

        def fetchall(self):
            return [
                ("GCF_RICE", "shared-accession", "PF00001.5", "Domain", "description"),
                ("GCF_MAIZE", "shared-accession", "PF00002.3", "Repeat", None),
            ]

    class FakeConnection:
        def __init__(self):
            self.cursor_instance = FakeCursor()

        def cursor(self):
            return self.cursor_instance

    connection = FakeConnection()
    annotations = esm2_similarity._fetch_pfam_annotations(
        connection,
        [
            {"assembly_accession": "GCF_RICE", "protein_accession": "shared-accession"},
            {"assembly_accession": "GCF_MAIZE", "protein_accession": "shared-accession"},
            {"assembly_accession": "GCF_RICE", "protein_accession": "shared-accession"},
        ],
    )

    assert connection.cursor_instance.params == (
        ["GCF_RICE", "GCF_MAIZE"],
        ["shared-accession", "shared-accession"],
    )
    assert annotations[("GCF_RICE", "shared-accession")] == [
        {
            "term_id": "PF00001.5",
            "term_desc": "Domain description",
            "evidence": None,
        }
    ]
    assert annotations[("GCF_MAIZE", "shared-accession")][0]["term_id"] == "PF00002.3"
    assert "annotation_pfam_domain" in connection.cursor_instance.query


def test_plaza_annotations_require_a_concordant_gene_protein_link():
    class FakeCursor:
        def __init__(self):
            self.queries = []

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def execute(self, query, params=None):
            self.queries.append((query, params))

        def fetchone(self):
            return (False,)

        def fetchall(self):
            return []

    class FakeConnection:
        def __init__(self):
            self.cursor_instance = FakeCursor()

        def cursor(self):
            return self.cursor_instance

    connection = FakeConnection()
    esm2_similarity._fetch_annotations(connection, ["protein-1"])
    query, params = connection.cursor_instance.queries[-1]

    assert params == (["protein-1"], ["go_plaza", "interpro", "mapman"])
    assert "dim_gene_protein AS plaza_link" in query
    assert "plaza_link.gene_source = 'plaza'" in query
    assert "plaza_link.coordinate_concordance = 'concordant'" in query
    assert "annotation.source <> ALL(%s::text[])" in query


def test_embedding_search_loads_database_index_before_expensive_model(monkeypatch):
    def fail_index_load(_species):
        raise RuntimeError("database unavailable")

    def should_not_embed(_sequence):
        raise AssertionError("The model should not load before the database index.")

    monkeypatch.setattr(esm2_similarity, "get_embedding_index", fail_index_load)
    monkeypatch.setattr(esm2_similarity, "get_supported_species", lambda: ["zea_mays"])
    monkeypatch.setattr(esm2_similarity, "embed_sequence", should_not_embed)

    with pytest.raises(RuntimeError, match="database unavailable"):
        esm2_similarity.find_similar_by_embedding("MKW")


def test_embedding_search_reports_truncation_and_annotation_metadata(monkeypatch):
    index = EmbeddingIndex(
        ["rice-protein"],
        np.array(["oryza_sativa"]),
        np.array([[1.0, 0.0]], dtype=np.float32),
    )
    monkeypatch.setattr(esm2_similarity, "get_embedding_index", lambda _species: index)
    monkeypatch.setattr(esm2_similarity, "get_supported_species", lambda: ["oryza_sativa"])
    monkeypatch.setattr(esm2_similarity, "embed_sequence", lambda _sequence: np.array([1.0, 0.0]))
    monkeypatch.setattr(
        esm2_similarity,
        "_fetch_metadata",
        lambda _connection, _accessions: {"rice-protein": {"plaza_gene_id": "Os01g00010"}},
    )
    monkeypatch.setattr(
        esm2_similarity,
        "_fetch_annotations",
        lambda _connection, _accessions: {"rice-protein": {"pfam": [{"term_id": "PF00001"}]}},
    )
    monkeypatch.setattr(
        esm2_similarity,
        "_fetch_pfam_annotations",
        lambda _connection, hits: {
            (hit["assembly_accession"], hit["protein_accession"]):
                [{"term_id": "PF00001"}]
            for hit in hits
        },
    )

    class FakeConnection:
        def close(self):
            pass

    monkeypatch.setattr(esm2_similarity, "_connect", lambda: FakeConnection())
    results = esm2_similarity.find_similar_by_embedding("M" * (esm2_similarity.MAX_RESIDUES + 1))

    assert results.query_length == esm2_similarity.MAX_RESIDUES + 1
    assert results.query_length_used == esm2_similarity.MAX_RESIDUES
    assert results[0]["gene_name"] == "Os01g00010"
    assert results[0]["annotation"]["pfam"][0]["term_id"] == "PF00001"


def test_stored_embedding_neighbor_search_skips_model_inference(monkeypatch):
    query_vector = np.zeros(EMBED_DIM, dtype=np.float32)
    query_vector[0] = 1.0
    neighbor_vector = np.zeros(EMBED_DIM, dtype=np.float32)
    neighbor_vector[0] = 0.8
    neighbor_vector[1] = 0.6
    index = EmbeddingIndex(
        ["candidate", "neighbor"],
        np.array(["zea_mays", "oryza_sativa"]),
        np.stack([query_vector, neighbor_vector]),
        assemblies=["GCF_MAIZE", "GCF_RICE"],
    )

    class FakeCursor:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def execute(self, query, params):
            self.query = query
            self.params = params

        def fetchone(self):
            return ("zea_mays", query_vector.tobytes())

    class FakeConnection:
        def __init__(self):
            self.cursor_instance = FakeCursor()
            self.closed = False

        def cursor(self):
            return self.cursor_instance

        def close(self):
            self.closed = True

    connection = FakeConnection()
    monkeypatch.setattr(esm2_similarity, "_connect", lambda: connection)
    monkeypatch.setattr(esm2_similarity, "get_embedding_index", lambda _species: index)
    monkeypatch.setattr(
        esm2_similarity,
        "get_supported_species",
        lambda: ["zea_mays", "oryza_sativa"],
    )
    monkeypatch.setattr(
        esm2_similarity,
        "_fetch_annotations",
        lambda _connection, accessions: {
            accession: {"pfam": []} for accession in accessions
        },
    )
    monkeypatch.setattr(
        esm2_similarity,
        "_fetch_pfam_annotations",
        lambda _connection, hits: {
            (hit["assembly_accession"], hit["protein_accession"]): []
            for hit in hits
        },
    )
    monkeypatch.setattr(
        esm2_similarity,
        "embed_sequence",
        lambda _sequence: pytest.fail("Stored-vector search must not run ESM-2."),
    )

    results = esm2_similarity.find_similar_by_stored_embedding(
        "GCF_MAIZE",
        "candidate",
        top_n=1,
    )

    assert len(results) == 1
    assert results[0]["accession"] == "neighbor"
    assert results[0]["assembly_accession"] == "GCF_RICE"
    assert results[0]["esm2_cosine_similarity"] == pytest.approx(0.8)
    assert connection.cursor_instance.params == ("GCF_MAIZE", "candidate")
    assert (
        "embedding.assembly_accession = %s" in connection.cursor_instance.query
        and "embedding.protein_accession = %s" in connection.cursor_instance.query
    )
    assert connection.closed


def test_embedding_indexes_are_searched_one_species_at_a_time(monkeypatch):
    query = np.array([1.0, 0.0], dtype=np.float32)
    indexes = {
        "zea_mays": EmbeddingIndex(
            ["maize-hit"],
            np.array(["zea_mays"]),
            np.array([[0.8, 0.6]], dtype=np.float32),
            assemblies=["GCF_MAIZE"],
        ),
        "nicotiana_tabacum_k326": EmbeddingIndex(
            ["k326-hit"],
            np.array(["nicotiana_tabacum_k326"]),
            np.array([[0.6, 0.8]], dtype=np.float32),
            assemblies=["GCF_K326"],
        ),
    }
    loaded_species = []

    def load_species(species):
        loaded_species.append(species)
        return indexes[species]

    monkeypatch.setattr(esm2_similarity, "get_embedding_index", load_species)

    hits = esm2_similarity._search_embedding_indexes(
        query,
        top_n=2,
        species_names=["zea_mays", "nicotiana_tabacum_k326"],
    )

    assert loaded_species == ["zea_mays", "nicotiana_tabacum_k326"]
    assert [hit["species"] for hit in hits] == [
        "zea_mays",
        "nicotiana_tabacum_k326",
    ]