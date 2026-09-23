from scripts import collect_ncbi as cn


def test_make_record_uses_gene_symbol_when_provided():
    rec = cn.make_record_from_fasta(
        "NC_000001.1:100-250",
        "ACGTACGT",
        db="nucleotide",
        resolved_gene_id="12345",
        gene_symbol="BRCA1",
        organism="Homo sapiens",
    )

    assert rec["gene_id"] == "GeneID:12345"
    assert rec["symbol"] == "BRCA1"
    assert rec["external_links"]["ncbi_gene"] == "https://www.ncbi.nlm.nih.gov/gene/12345"
