import hashlib
import json
import os
from pathlib import Path

import ijson

from collect.collect_all_sources import _sequence_hash
from scripts.collect_pfam import hmm_index_files
from scripts.protein_json_annotations import annotate_species_file, normalize_protein


def _write_species(path: Path, proteins: list[tuple[str, str]]) -> None:
    genes = [
        {
            "gene_id": gene_id,
            "sequence": {"protein": sequence, "dna_hash": "sha256:dna-hash"},
            "annotation": {"go_terms": [{"id": "GO:0000001"}], "tf_family": None},
        }
        for gene_id, sequence in proteins
    ]
    path.write_text(
        json.dumps({"metadata": {"plant": "Test plant"}, "genes": genes}),
        encoding="utf-8",
    )


def _read_genes(path: Path) -> list[dict]:
    with path.open("rb") as source:
        return list(ijson.items(source, "genes.item"))


def _fake_tools():
    calls: list[list[str]] = []

    def run_batch(query: Path, output: Path, threads: int) -> None:
        identifiers = [
            line[1:]
            for line in query.read_text(encoding="ascii").splitlines()
            if line.startswith(">")
        ]
        calls.append(identifiers)
        output.write_text("\n".join(identifiers), encoding="ascii")

    def parse_output(output: Path) -> dict[str, list[dict]]:
        return {
            identifier: [{"hit": "example"}]
            for identifier in output.read_text(encoding="ascii").splitlines()
            if identifier
        }

    return calls, run_batch, parse_output


def test_protein_hash_uses_existing_sha256_prefix_convention():
    sequence = "MPEPTIDEAC"
    normalized = normalize_protein(sequence.lower())

    assert normalized is not None
    _, digest, schema_hash = normalized
    assert digest == hashlib.sha256(sequence.encode("ascii")).hexdigest()
    assert schema_hash == _sequence_hash(sequence)


def test_hmmer_index_paths_match_hmmpress_suffixes(tmp_path):
    database = tmp_path / "Pfam-A.hmm"

    assert [path.name for path in hmm_index_files(database)] == [
        "Pfam-A.hmm.h3f",
        "Pfam-A.hmm.h3i",
        "Pfam-A.hmm.h3m",
        "Pfam-A.hmm.h3p",
    ]


def test_species_json_deduplicates_sequences_and_preserves_annotations(tmp_path):
    input_path = tmp_path / "test_all_sources.json"
    reference = tmp_path / "reference.dmnd"
    reference.write_text("fake reference", encoding="ascii")
    _write_species(input_path, [
        ("gene-1", "MPEPTIDEAC"),
        ("gene-2", "mpepTideac"),
        ("gene-3", "MPEPTIDEFG"),
    ])
    calls, run_batch, parse_output = _fake_tools()
    output_path = tmp_path / "annotated" / input_path.name

    result = annotate_species_file(
        input_path=input_path,
        output_path=output_path,
        tool_name="diamond-test",
        annotation_key="diamond_hits",
        reference_files=[reference],
        options={},
        batch_size=1,
        threads=1,
        limit=None,
        run_batch=run_batch,
        parse_output=parse_output,
    )

    genes = _read_genes(output_path)
    assert result["protein_records"] == 3
    assert result["unique_sequences"] == 2
    assert result["batches_run"] == 2
    assert sum(map(len, calls)) == 2
    assert genes[0]["sequence"]["protein_hash"] == genes[1]["sequence"]["protein_hash"]
    assert genes[0]["annotation"]["diamond_hits"] == genes[1]["annotation"]["diamond_hits"]
    assert genes[0]["annotation"]["go_terms"] == [{"id": "GO:0000001"}]
    assert "diamond_hits" in genes[2]["annotation"]

    pfam_calls, pfam_run_batch, pfam_parse_output = _fake_tools()
    annotate_species_file(
        input_path=input_path,
        output_path=output_path,
        tool_name="pfam-test",
        annotation_key="pfam_domains",
        reference_files=[reference],
        options={},
        batch_size=1,
        threads=1,
        limit=None,
        run_batch=pfam_run_batch,
        parse_output=pfam_parse_output,
    )

    genes = _read_genes(output_path)
    assert len(pfam_calls) == 2
    assert genes[0]["annotation"]["diamond_hits"] == genes[0]["annotation"]["pfam_domains"]
    assert genes[0]["annotation"]["go_terms"] == [{"id": "GO:0000001"}]


def test_interrupted_species_json_run_resumes_from_committed_batch(tmp_path):
    input_path = tmp_path / "resume_all_sources.json"
    reference = tmp_path / "reference.hmm"
    reference.write_text("fake reference", encoding="ascii")
    _write_species(input_path, [
        ("gene-1", "MPEPTIDEAC"),
        ("gene-2", "MPEPTIDEFG"),
        ("gene-3", "MPEPTIDEHI"),
    ])
    original = input_path.read_bytes()
    failed_calls, run_batch, parse_output = _fake_tools()

    def interrupt_second_batch(query: Path, output: Path, threads: int) -> None:
        if not failed_calls:
            run_batch(query, output, threads)
            return
        raise RuntimeError("simulated interruption")

    try:
        annotate_species_file(
            input_path, input_path, "resume-test", "pfam_domains", [reference], {},
            1, 1, None, interrupt_second_batch, parse_output,
        )
    except RuntimeError as exc:
        assert str(exc) == "simulated interruption"
    else:
        raise AssertionError("the synthetic interruption should have propagated")

    assert input_path.read_bytes() == original
    completed_before_resume = failed_calls[0][0]
    resumed_calls: list[list[str]] = []

    def record_resume(query: Path, output: Path, threads: int) -> None:
        identifiers = [
            line[1:]
            for line in query.read_text(encoding="ascii").splitlines()
            if line.startswith(">")
        ]
        resumed_calls.append(identifiers)
        run_batch(query, output, threads)

    result = annotate_species_file(
        input_path, input_path, "resume-test", "pfam_domains", [reference], {},
        1, 1, None, record_resume, parse_output,
    )

    assert all(completed_before_resume not in batch for batch in resumed_calls)
    assert result["unique_sequences"] == 3
    assert result["annotated_records"] == 3
    assert all("pfam_domains" in gene["annotation"] for gene in _read_genes(input_path))


def test_newer_source_does_not_reuse_stale_annotated_copy(tmp_path):
    input_path = tmp_path / "source_all_sources.json"
    output_path = tmp_path / "annotated" / input_path.name
    output_path.parent.mkdir()
    reference = tmp_path / "reference.dmnd"
    reference.write_text("fake reference", encoding="ascii")
    _write_species(input_path, [("current-gene", "MPEPTIDEAC")])
    _write_species(output_path, [("stale-gene", "MPEPTIDEFG")])
    older_time = input_path.stat().st_mtime_ns - 1_000_000_000
    os.utime(output_path, ns=(older_time, older_time))
    _, run_batch, parse_output = _fake_tools()

    annotate_species_file(
        input_path, output_path, "freshness-test", "diamond_hits", [reference], {},
        1, 1, None, run_batch, parse_output,
    )

    genes = _read_genes(output_path)
    assert [gene["gene_id"] for gene in genes] == ["current-gene"]