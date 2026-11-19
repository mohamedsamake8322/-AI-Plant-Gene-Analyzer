# ESM-2 integration

ESM-2 is an additional evidence source. Its cosine similarity is not a sequence-identity percentage, an evolutionary distance, or proof of function. Topic-search scores and alignment scores remain separate.

## Configuration

Install the application requirements and set either a dedicated `BIOINFO_PIPELINE_URL` or the dedicated `BIOINFO_PIPELINE_HOST`, `BIOINFO_PIPELINE_PORT`, `BIOINFO_PIPELINE_USER`, `BIOINFO_PIPELINE_PASSWORD`, and `BIOINFO_PIPELINE_DBNAME` variables. The module loads the project `.env` at import time. Do not point the generic `DATABASE_URL` at this feature: it may refer to the main application database.

The engine lives in `core_engines/esm2_similarity.py` and loads the project-root `.env`. Set `ESM2_MODEL_DIR` to a local ESM-2 model directory. The default is the rice model folder at `C:\Downloads\ESM\ESM_assemblage_ IRGSP-1.0\esm2_rice\model_esm2_650M`. Restart Streamlit after changing these variables.

The database must contain `esm2_embeddings`, `dim_protein`, and `fact_functional_annotation`. Embeddings can be stored as float32 `bytea` (preferred) or as a legacy PostgreSQL `REAL[]`; the loader supports both. The annotation `evidence` column is optional. The first ESM-2 database search loads the full embedding index into process memory, approximately 1 GB for 191,000 1280-dimensional float32 vectors.

## Coverage and interpretation

The reference index covers rice (`oryza_sativa`), maize (`zea_mays`), and
the published tobacco assembly/assemblies registered under
`nicotiana_tabacum`. Query sequences from other plants are compared against
those registry rows only. Sequences longer than 1,022 residues are truncated
for model inference.

The embedding cosine-distance matrix is defined as `1 - cosine_similarity`, is computed without sequence alignment, and is available only for protein sequences. UPGMA or Neighbor-Joining can display its clustering, but the result is phenetic rather than a phylogenetic inference. IQ-TREE continues to use an MSA and substitution model independently.

## Main Sequence Analysis

For protein inputs, the **Similarity** tab offers an on-demand ESM-2 search. It is separate from alignment identity and can still be used when the alignment search is skipped by its cost budget. The UI displays the cosine score as its own metric and includes the returned Pfam, InterPro, GO/PLAZA, and Swiss-Prot annotations. It does not convert the cosine score into percent identity or merge the two rankings.

The button loads the ESM-2 model and the embedding index on first use; both remain cached by the Streamlit process. The search truncates inputs longer than 1,022 residues and reports that truncation. Configure the local model and `BIOINFO_PIPELINE_*` connection variables before using the button; failures are shown in the UI and logged.

## Topic search

The topic-search tab uses the dedicated `bioinfo_pipeline` database, not the
main `genes` database used by `similarityengine.py` and
`scripts/postgres_utils.py`. After importing the normalized tables, provision
its indexed text source once:

```powershell
psql "$env:BIOINFO_PIPELINE_URL" -f scripts/create_trait_search_view.sql
```

The materialized view indexes curated descriptions, GO, Pfam, InterPro, and
linked QTL evidence. KEGG, MapMan, TF-family, and literature fields are not
available in this imported schema. `species_registry` deliberately enables
only rice, maize, and tobacco K326. After updating annotation/QTL imports,
refresh the view with:

```sql
REFRESH MATERIALIZED VIEW CONCURRENTLY v_gene_text;
```

The topic-search ESM-2 button retrieves the candidate vector from
`esm2_embeddings` by `(assembly_accession, protein_accession)` and searches
the cached reference index. It does not recompute the candidate embedding or
load the ESM-2 model. The first neighbor search still loads the reference
matrix into memory.

## Pipeline alignment

The Colab notebooks map to the reproducible pipeline as follows:

| Notebook | Planning step | Scope |
|---|---|---|
| `Phase1.ipynb` | S1/S2 | collect FASTA/GFF3, then load the structural base |
| `Phase2.ipynb` | S3 | Pfam/InterPro/GO functional annotation |
| `Phase3 pour le riz/maïs/tabac.ipynb` | S4 | ESM-2 embeddings per assembly |
| `Phase4.ipynb` | S5 | orthology, expression, QTL and trait enrichment |
| `Phase5.ipynb` | S6 | completeness and quality-control reporting |

S0 is now implemented by `scripts/register_species.py`, which reads
`species.yaml`, validates assembly accessions, and performs an idempotent
upsert into `species_registry`. It never publishes an assembly. S6 is
implemented by `scripts/validate_bioinfo_pipeline.py`; it reports quantitative
coverage for sequences, coordinates, annotations, source-scan evidence,
checksums, sequence fingerprints, and embeddings. S7 is implemented by
`scripts/publish_species.py`: it reruns S6, updates the registry and audit
table in one transaction, refreshes `v_gene_text`, and clears the ESM-2 cache.
An assembly is visible to the app only when its registry row is both
`enabled = true` and `status = 'published'`.

The normalized schema uses `(assembly_accession, protein_accession)` as the
key for `dim_sequence` and `esm2_embeddings`, and `(assemblage,
protein_accession)` for `dim_protein` and functional facts. This is required
for tobacco K326/TN90, where protein accessions can overlap.

### D1 decision: retain CDS

CDS collection remains in scope. `dim_sequence.cds_sequence` is already part
of the normalized schema, and retaining it preserves the GC-content and codon
panels in Statistics without coupling those panels back to the legacy
`genes` JSON. CDS is optional for Trait Search and ESM-2, but mandatory for a
species that declares the Statistics module in its registry metadata.

Example:

```powershell
python scripts/register_species.py --manifest species.yaml
python scripts/validate_bioinfo_pipeline.py --out results/bioinfo_pipeline_qc.json
python scripts/measure_trait_search_recall.py --out results/trait_search_recall.json
python scripts/publish_species.py --species zea_mays
```

The recall command measures full-text `simple` against the historical
case-insensitive substring search for drought/rice and lodging/maize. It must
be run with `BIOINFO_PIPELINE_URL` configured; no counts are fabricated when
the database is unavailable.

The database-backed CLI supports species discovery and targeted searches:

```powershell
python trait_research.py --list-species
python trait_research.py --species zea_mays --topic lodging --top-n 20 --candidate-pool 1000
python trait_research.py --species oryza_sativa --topic drought --top-n 20 --fetch-refs
```

## Checks

Run the focused integration checks with:

```powershell
.\.venv\Scripts\python.exe -m pytest -q test_trait_research.py test_esm2_similarity.py test_phylo_view.py test_i18n.py
```