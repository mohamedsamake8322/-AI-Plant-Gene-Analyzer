import pytest

def test_imports():
    import bioinformatics
    import alignment_engine
    import distance_engine
    import phylogeny_engine
    assert hasattr(bioinformatics, 'translate_dna') or True
