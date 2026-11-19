import importlib.util
from pathlib import Path


def _load_registry_module():
    path = Path(__file__).parent / "scripts" / "register_species.py"
    spec = importlib.util.spec_from_file_location("register_species", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_registry_manifest_hash_is_deterministic():
    registry = _load_registry_module()
    entry = {"common_name": "maize", "refseq_accession": "GCF_1"}
    assert registry._manifest_hash(entry) == registry._manifest_hash(dict(entry))


def test_registry_creates_one_key_for_single_reference_assembly():
    registry = _load_registry_module()
    rows = registry._assemblies(
        "zea_mays",
        {
            "reference_assembly": "B73-NAM-5.0",
            "refseq_accession": "GCF_902167145.1",
        },
    )
    assert rows == [
        ("zea_mays", "B73-NAM-5.0", "GCF_902167145.1", "zea_mays")
    ]


def test_registry_rejects_unversioned_assembly_accession():
    registry = _load_registry_module()
    try:
        registry._assemblies(
            "test_species",
            {"reference_assembly": "unversioned", "refseq_accession": "not-an-accession"},
        )
    except ValueError as error:
        assert "valid GCA/GCF accession" in str(error)
    else:
        raise AssertionError("Invalid assembly accession should fail S0.")
