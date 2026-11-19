import trait_research as trait


def test_search_candidates_uses_indexed_postgres_full_text(monkeypatch):
    class FakeCursor:
        description = [
            (name,) for name in (
                "species",
                "assembly_accession",
                "protein_accession",
                "plaza_gene_id",
                "gene_symbol",
                "description",
                "go_text",
                "pfam_text",
                "interpro_text",
                "trait_text",
                "qtl_evidence",
                "search_text",
                "text_rank",
                "total_matches",
            )
        ]

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def execute(self, query, params):
            self.query = query
            self.params = params

        def fetchall(self):
            return [(
                "zea_mays",
                "GCF_902167145.1",
                "NP_001.1",
                "Zm00001eb000170",
                None,
                "lignin biosynthesis",
                "",
                "PF00001",
                "",
                "",
                [],
                "lignin biosynthesis PF00001",
                0.75,
                4,
            )]

    class FakeConnection:
        def __init__(self):
            self.cursor_instance = FakeCursor()
            self.closed = False

        def cursor(self):
            return self.cursor_instance

        def close(self):
            self.closed = True

    connection = FakeConnection()
    monkeypatch.setattr(trait, "_connect_bioinfo_pipeline", lambda: connection)
    template = {
        "keywords": {"lignification": ["lignin"]},
        "noisy": set(),
    }

    result = trait.search_candidates(template, "zea_mays", limit=20)

    assert len(result) == 1
    assert result[0]["gene_id"] == "NP_001.1"
    assert result[0]["matches"] == {"lignification:lignin"}
    assert result[0]["gene"]["assembly_accession"] == "GCF_902167145.1"
    assert "websearch_to_tsquery" in connection.cursor_instance.query
    assert connection.cursor_instance.params == (
        '"lignin"', "zea_mays", '"lignin"', ["%lignin%"], 100
    )
    assert connection.closed


def test_pick_search_term_skips_refseq_and_uncharacterized_descriptions():
    assert trait.pick_search_term(
        {
            "gene_symbol": "XP_012345.1",
            "common_name": "XP_012345.1",
            "description": "uncharacterized protein",
            "annotation": {"go_text": "GO:0008150 biological_process; GO:0003674"},
        },
        "XP_012345.1",
    ) == "biological_process"


def test_candidate_key_contains_assembly_and_accession():
    assert trait._looks_like_systematic_id("NP_001.1")
    assert trait._looks_like_systematic_id("Zm00001eb000170")


def test_species_cli_input_is_normalized_as_a_registry_key():
    assert trait.resolve_species_filter("zea_mays") == "zea_mays"
    assert trait.resolve_species_filter("Oryza sativa") == "oryza_sativa"


def test_candidate_score_uses_specificity_diversity_then_full_text_rank():
    template = {
        "tier_a": {"lignin"},
        "tier_b": {"peroxidase"},
    }
    candidates = [
        {"matches": {"lignification:lignin"}, "text_rank": 0.2},
        {"matches": {"lignification:lignin", "stem_rigidity:peroxidase"}, "text_rank": 0.1},
    ]

    ranked = trait.score_candidates(candidates, template)

    assert ranked[0]["score"] > ranked[1]["score"]
    assert ranked[0]["text_rank"] == 0.1
