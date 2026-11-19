# Comparaison genes vs bioinfo_pipeline

## Base `plant_gene_analyzer_clean`

### `public.genes` — 168,577 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `id` | integer | 168,577 | 100.0% |
| `gene_id` | text | 168,577 | 100.0% |
| `symbol` | text | 0 | 0.0% |
| `organism` | text | 168,577 | 100.0% |
| `sequence` | text | 105,790 | 62.8% |
| `sequence_type` | text | 105,790 | 62.8% |
| `description` | text | 168,577 | 100.0% |
| `source` | text | 168,577 | 100.0% |
| `source_url` | text | 0 | 0.0% |
| `external_links` | jsonb | 168,577 | 100.0% |
| `expression_profiles` | jsonb | 168,577 | 100.0% |
| `pathways` | jsonb | 168,577 | 100.0% |
| `publications` | jsonb | 168,577 | 100.0% |
| `annotations` | jsonb | 168,577 | 100.0% |
| `traits` | jsonb | 168,577 | 100.0% |
| `length` | integer | 105,790 | 62.8% |
| `date_added` | timestamp with time zone | 0 | 0.0% |
| `sequence_hash` | text | 105,790 | 62.8% |
| `origin` | text | 168,577 | 100.0% |
| `relations` | jsonb | 168,577 | 100.0% |
| `kmer_hashes` | ARRAY | 105,790 | 62.8% |

Répartition par `organism` :

- Zea mays: 96,974
- Chenopodium quinoa: 71,538
- Zea mays subsp. mays: 61
- Zea mays subsp. huehuetenangensis: 2
- Zea mays subsp. mexicana: 1
- Zea mays subsp. parviglumis: 1

## Base `bioinfo_pipeline`

### `public.dim_protein` — 190,607 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `protein_accession` | text | 190,607 | 100.0% |
| `assemblage` | text | 190,607 | 100.0% |
| `plaza_gene_id` | text | 65,204 | 34.2% |
| `chr` | text | 65,204 | 34.2% |
| `start_pos` | bigint | 65,204 | 34.2% |
| `stop_pos` | bigint | 65,204 | 34.2% |
| `strand` | text | 65,204 | 34.2% |

Répartition par `assemblage` :

- nicotiana_tabacum_k326: 64,776
- nicotiana_tabacum_tn90: 60,627
- oryza_sativa: 41,610
- zea_mays: 23,594

### `public.esm2_embeddings` — 191,304 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `protein_accession` | text | 191,304 | 100.0% |
| `species` | text | 191,304 | 100.0% |
| `assembly_accession` | text | 191,304 | 100.0% |
| `embedding` | bytea | 191,304 | 100.0% |

Répartition par `species` :

- nicotiana_tabacum_k326: 91,379
- zea_mays: 57,345
- oryza_sativa: 42,580

### `public.fact_expression` — 833,554 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `assemblage` | text | 833,554 | 100.0% |
| `protein_accession` | text | 833,554 | 100.0% |
| `plaza_gene_id` | text | 833,554 | 100.0% |
| `tissue` | text | 833,554 | 100.0% |
| `tpm_median` | numeric | 833,554 | 100.0% |
| `n_replicates` | integer | 833,554 | 100.0% |
| `source_experiment` | text | 833,554 | 100.0% |

Répartition par `assemblage` :

- zea_mays: 542,662
- oryza_sativa: 290,892

### `public.fact_functional_annotation` — 7,843,548 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `id` | integer | 7,843,548 | 100.0% |
| `assemblage` | text | 7,843,548 | 100.0% |
| `protein_accession` | text | 7,843,548 | 100.0% |
| `source` | text | 7,843,548 | 100.0% |
| `term_id` | text | 7,652,963 | 97.6% |
| `term_desc` | text | 1,880,534 | 24.0% |
| `evidence` | text | 5,963,014 | 76.0% |

Répartition par `assemblage` :

- oryza_sativa: 3,802,094
- zea_mays: 2,753,356
- nicotiana_tabacum_k326: 674,809
- nicotiana_tabacum_tn90: 613,289

### `public.fact_orthology` — 64,383 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `protein_accession` | text | 64,383 | 100.0% |
| `assemblage` | text | 64,383 | 100.0% |
| `plaza_gene_id` | text | 64,383 | 100.0% |
| `orthofam_id` | text | 64,383 | 100.0% |
| `plaza_instance` | text | 64,383 | 100.0% |

Répartition par `assemblage` :

- oryza_sativa: 40,789
- zea_mays: 23,594

### `public.fact_qtl` — 1,051 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `qtl_id` | text | 1,051 | 100.0% |
| `qtl_gene_name` | text | 716 | 68.1% |
| `major_category` | text | 1,051 | 100.0% |
| `trait_category` | text | 1,051 | 100.0% |
| `trait_detail` | text | 1,030 | 98.0% |
| `chr` | text | 1,051 | 100.0% |
| `genome_start` | bigint | 1,051 | 100.0% |
| `genome_end` | bigint | 1,051 | 100.0% |
| `mapping_method` | text | 1,051 | 100.0% |
| `explained_variance` | numeric | 799 | 76.0% |
| `additive_effect` | numeric | 720 | 68.5% |
| `year` | integer | 1,051 | 100.0% |
| `reference` | text | 1,051 | 100.0% |
| `lo` | bigint | 1,051 | 100.0% |
| `hi` | bigint | 1,051 | 100.0% |
| `interval_bp` | bigint | 1,051 | 100.0% |
| `resolution` | text | 1,051 | 100.0% |
| `n_genes_in_interval` | integer | 1,051 | 100.0% |

### `public.fact_qtl_gene_link` — 220,155 lignes

| Colonne | Type | Non-nulles | % rempli |
| --- | --- | --- | --- |
| `qtl_id` | text | 220,155 | 100.0% |
| `gene_id` | text | 220,155 | 100.0% |
| `chr` | text | 220,155 | 100.0% |
| `qtl_start` | bigint | 220,155 | 100.0% |
| `qtl_end` | bigint | 220,155 | 100.0% |
| `gene_start` | bigint | 220,155 | 100.0% |
| `gene_end` | bigint | 220,155 | 100.0% |

## Identifiants en commun

Comparaison des valeurs distinctes (minuscules, suffixe de version retiré) entre colonnes candidates.

| genes (table.colonne) | bioinfo_pipeline (table.colonne) | En commun | Distincts genes | Distincts bio |
| --- | --- | --- | --- | --- |
| `genes.gene_id` | `esm2_embeddings.protein_accession` | 40 | 168,577 | 191,304 |
| `genes.gene_id` | `fact_orthology.protein_accession` | 16 | 168,577 | 61,205 |
| `genes.gene_id` | `fact_functional_annotation.protein_accession` | 16 | 168,577 | 146,502 |
| `genes.gene_id` | `fact_expression.protein_accession` | 16 | 168,577 | 61,528 |
| `genes.gene_id` | `dim_protein.protein_accession` | 16 | 168,577 | 146,502 |
