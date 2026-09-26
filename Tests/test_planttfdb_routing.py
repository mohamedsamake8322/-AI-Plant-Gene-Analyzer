from collect import collect_planttfdb
from collect.collect_all_sources import _checkpoint_can_be_reused


def test_oryza_sativa_fetches_both_current_planttfdb_datasets(monkeypatch):
    calls = []

    def fake_download(sp_code, species, retmax, is_extended=False):
        calls.append((sp_code, species, retmax, is_extended))
        return [{"gene_id": f"{sp_code}_gene", "source": "planttfdb"}]

    monkeypatch.setattr(collect_planttfdb, "_fetch_via_download", fake_download)

    records = collect_planttfdb.fetch_planttfdb("Oryza sativa", retmax=10)

    assert [call[0] for call in calls] == ["Osi", "Osj"]
    assert all(call[1] == "Oryza sativa" for call in calls)
    assert [record["gene_id"] for record in records] == ["Osi_gene", "Osj_gene"]


def test_oryza_subspecies_map_to_official_planttfdb_codes():
    assert collect_planttfdb.SPECIES_MAP["oryza sativa subsp. indica"] == "Osi"
    assert collect_planttfdb.SPECIES_MAP["oryza sativa subsp. japonica"] == "Osj"
    assert "oryza sativa" not in collect_planttfdb.SPECIES_MAP


def test_recovery_mode_reuses_nonempty_checkpoint_but_not_empty_one(tmp_path):
    complete = tmp_path / "source.json"
    empty = tmp_path / "planttfdb.json"
    complete.write_text("[]", encoding="utf-8")
    empty.write_text("[]", encoding="utf-8")

    assert _checkpoint_can_be_reused(complete, [{"id": 1}], 100, True)
    assert not _checkpoint_can_be_reused(empty, [], 100, True)
    assert not _checkpoint_can_be_reused(complete, [{"id": 1}], 100, False)
