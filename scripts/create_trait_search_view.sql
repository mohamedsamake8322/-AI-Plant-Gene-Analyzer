-- Run against the dedicated bioinfo_pipeline database after importing
-- dim_protein, fact_functional_annotation, fact_qtl, fact_qtl_gene_link,
-- esm2_embeddings, and initializing annotation_pfam_domain with
-- scripts/load_annotation_evidence.py.

CREATE TABLE IF NOT EXISTS species_registry (
    species TEXT PRIMARY KEY,
    display_name TEXT NOT NULL,
    assemblage TEXT NOT NULL,
    assembly_accession TEXT NOT NULL UNIQUE,
    enabled BOOLEAN NOT NULL DEFAULT FALSE,
    status TEXT NOT NULL DEFAULT 'registered',
    completed_phase SMALLINT NOT NULL DEFAULT 0,
    manifest_hash TEXT,
    modules JSONB NOT NULL DEFAULT '{}'::jsonb,
    source_metadata JSONB NOT NULL DEFAULT '{}'::jsonb,
    tool_versions JSONB NOT NULL DEFAULT '{}'::jsonb,
    source_protein_count BIGINT,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT now(),
    CONSTRAINT species_registry_status_ck
        CHECK (status IN ('registered', 'collecting', 'validating', 'published', 'blocked')),
    CONSTRAINT species_registry_phase_ck
        CHECK (completed_phase BETWEEN 0 AND 7)
);

ALTER TABLE species_registry
    ADD COLUMN IF NOT EXISTS status TEXT NOT NULL DEFAULT 'registered',
    ADD COLUMN IF NOT EXISTS completed_phase SMALLINT NOT NULL DEFAULT 0,
    ADD COLUMN IF NOT EXISTS manifest_hash TEXT,
    ADD COLUMN IF NOT EXISTS modules JSONB NOT NULL DEFAULT '{}'::jsonb,
    ADD COLUMN IF NOT EXISTS source_metadata JSONB NOT NULL DEFAULT '{}'::jsonb,
    ADD COLUMN IF NOT EXISTS tool_versions JSONB NOT NULL DEFAULT '{}'::jsonb,
    ADD COLUMN IF NOT EXISTS source_protein_count BIGINT,
    ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ NOT NULL DEFAULT now();

CREATE INDEX IF NOT EXISTS species_registry_enabled_idx
    ON species_registry (enabled, status, display_name);

DROP MATERIALIZED VIEW IF EXISTS v_gene_text;

CREATE MATERIALIZED VIEW v_gene_text AS
WITH plaza_concordant_proteins AS (
    SELECT DISTINCT
        gene.assembly_name AS assemblage,
        gene.protein_accession
    FROM dim_gene_protein AS gene
    JOIN species_registry AS registry
      ON registry.assemblage = gene.assembly_name
     AND registry.enabled IS TRUE
     AND registry.status = 'published'
    WHERE gene.gene_source = 'plaza'
      AND gene.coordinate_concordance = 'concordant'
),
annotation_values AS (
    SELECT
        annotation.assemblage,
        annotation.protein_accession,
        CASE
            WHEN annotation.source IN ('swissprot', 'description', 'gene_description')
                THEN 'description'
            WHEN annotation.source = 'refseq_product' THEN 'refseq_product'
            WHEN annotation.source IN ('go_plaza', 'go', 'gene_ontology') THEN 'go'
            WHEN annotation.source = 'interpro' THEN 'interpro'
            WHEN annotation.source = 'mapman' THEN 'mapman'
        END AS annotation_kind,
        CASE
            WHEN annotation.source IN ('swissprot', 'description', 'gene_description',
                                       'refseq_product')
                THEN annotation.term_desc
            ELSE concat_ws(' ', annotation.term_id, annotation.term_desc)
        END AS text_value
    FROM fact_functional_annotation AS annotation
    JOIN species_registry AS registry
      ON registry.assemblage = annotation.assemblage
     AND registry.enabled IS TRUE
     AND registry.status = 'published'
    LEFT JOIN plaza_concordant_proteins AS concordant
      ON concordant.assemblage = annotation.assemblage
     AND concordant.protein_accession = annotation.protein_accession
    WHERE annotation.source IN (
        'swissprot', 'description', 'gene_description', 'refseq_product',
        'go_plaza', 'go', 'gene_ontology', 'interpro', 'mapman'
    )
      AND (
          annotation.source NOT IN ('go_plaza', 'interpro', 'mapman')
          OR concordant.protein_accession IS NOT NULL
      )
),
annotation_dedup AS (
    SELECT DISTINCT assemblage, protein_accession, annotation_kind, text_value
    FROM annotation_values
    WHERE text_value IS NOT NULL
),
annotation_rollup AS (
    SELECT
        assemblage,
        protein_accession,
        string_agg(text_value, ' ')
            FILTER (WHERE annotation_kind = 'description') AS description,
        string_agg(text_value, ' ')
            FILTER (WHERE annotation_kind = 'refseq_product') AS refseq_product,
        string_agg(text_value, ' ')
            FILTER (WHERE annotation_kind = 'go') AS go_text,
        string_agg(text_value, ' ')
            FILTER (WHERE annotation_kind = 'interpro') AS interpro_text,
        string_agg(text_value, ' ')
            FILTER (WHERE annotation_kind = 'mapman') AS mapman_text
    FROM annotation_dedup
    GROUP BY assemblage, protein_accession
),
pfam_rollup AS (
    SELECT
        domain.assembly_accession,
        domain.protein_accession,
        string_agg(
            DISTINCT concat_ws(' ', domain.pfam_accession, domain.pfam_name, domain.description),
            ' '
        ) AS pfam_text
    FROM annotation_pfam_domain AS domain
    JOIN species_registry AS registry
      ON registry.assembly_accession = domain.assembly_accession
     AND registry.enabled IS TRUE
     AND registry.status = 'published'
    GROUP BY domain.assembly_accession, domain.protein_accession
),
qtl_rollup AS (
    SELECT
        link.gene_id AS plaza_gene_id,
        string_agg(
            DISTINCT concat_ws(' ', qtl.major_category, qtl.trait_category, qtl.trait_detail),
            ' '
        ) AS trait_text,
        jsonb_agg(DISTINCT jsonb_build_object(
            'qtl_id', qtl.qtl_id,
            'trait_category', qtl.trait_category,
            'trait_detail', qtl.trait_detail,
            'major_category', qtl.major_category,
            'chr', qtl.chr,
            'start', qtl.genome_start,
            'end', qtl.genome_end,
            'explained_variance', qtl.explained_variance,
            'year', qtl.year,
            'resolution', qtl.resolution,
            'reference', qtl.reference
        )) FILTER (WHERE qtl.qtl_id IS NOT NULL) AS qtl_evidence
    FROM fact_qtl_gene_link AS link
    JOIN fact_qtl AS qtl ON qtl.qtl_id = link.qtl_id
    WHERE EXISTS (
        SELECT 1
        FROM dim_gene_protein AS plaza_link
        JOIN species_registry AS registry
          ON registry.assemblage = plaza_link.assembly_name
         AND registry.enabled IS TRUE
         AND registry.status = 'published'
        WHERE plaza_link.gene_source = 'plaza'
          AND plaza_link.coordinate_concordance = 'concordant'
          AND plaza_link.gene_id = link.gene_id
    )
    GROUP BY link.gene_id
),
gene_rollup AS (
    SELECT
        gene.assembly_name AS assemblage,
        gene.protein_accession,
        string_agg(DISTINCT gene.gene_id, ' ' ORDER BY gene.gene_id)
            FILTER (WHERE gene.gene_source = 'plaza') AS plaza_gene_ids,
        string_agg(DISTINCT gene.gene_id, ' ' ORDER BY gene.gene_id) AS gene_ids
    FROM dim_gene_protein AS gene
    JOIN species_registry AS registry
      ON registry.assemblage = gene.assembly_name
     AND registry.enabled IS TRUE
     AND registry.status = 'published'
    WHERE gene.gene_source = 'plaza'
      AND gene.coordinate_concordance = 'concordant'
    GROUP BY gene.assembly_name, gene.protein_accession
),
-- Un gène PLAZA par ligne, puis jointure par égalité (hash join) au lieu de "= ANY(tableau)",
-- qui obligeait PostgreSQL à comparer chaque protéine à chaque gène QTL (boucle imbriquée, très lent).
-- Un seul résultat par protéine : textes et preuves QTL de tous ses gènes concordants.
protein_gene AS (
    SELECT assemblage, protein_accession,
           unnest(string_to_array(plaza_gene_ids, ' ')) AS gid
    FROM gene_rollup
    WHERE plaza_gene_ids IS NOT NULL
),
protein_qtl AS (
    SELECT
        pg.assemblage,
        pg.protein_accession,
        string_agg(DISTINCT q.trait_text, ' ') AS trait_text,
        jsonb_agg(DISTINCT e.value) FILTER (WHERE e.value IS NOT NULL) AS qtl_evidence
    FROM protein_gene AS pg
    JOIN qtl_rollup AS q ON q.plaza_gene_id = pg.gid
    LEFT JOIN LATERAL jsonb_array_elements(q.qtl_evidence) AS e(value) ON TRUE
    GROUP BY pg.assemblage, pg.protein_accession
),
proteins AS (
    SELECT DISTINCT ON (
        registry.species, registry.assembly_accession, protein.protein_accession
    )
        registry.species,
        registry.assembly_accession,
        protein.protein_accession,
        genes.plaza_gene_ids AS plaza_gene_id,
        annotation.description,
        annotation.refseq_product,
        annotation.go_text,
        pfam.pfam_text,
        annotation.interpro_text,
        annotation.mapman_text,
        qtl.trait_text,
        COALESCE(qtl.qtl_evidence, '[]'::jsonb) AS qtl_evidence,
        COALESCE(
            to_jsonb(protein)->>'gene_symbol',
            to_jsonb(protein)->>'gene_id',
            genes.plaza_gene_ids,
            genes.gene_ids
        ) AS gene_symbol
    FROM dim_protein AS protein
    JOIN species_registry AS registry
      ON registry.assemblage = protein.assemblage
     AND registry.enabled IS TRUE
     AND registry.status = 'published'
    JOIN esm2_embeddings AS embedding
      ON embedding.protein_accession = protein.protein_accession
     AND embedding.assembly_accession = registry.assembly_accession
     AND embedding.species = registry.species
    LEFT JOIN annotation_rollup AS annotation
      ON annotation.assemblage = protein.assemblage
     AND annotation.protein_accession = protein.protein_accession
    LEFT JOIN pfam_rollup AS pfam
      ON pfam.assembly_accession = registry.assembly_accession
     AND pfam.protein_accession = protein.protein_accession
    LEFT JOIN gene_rollup AS genes
      ON genes.assemblage = protein.assemblage
     AND genes.protein_accession = protein.protein_accession
    LEFT JOIN protein_qtl AS qtl
      ON qtl.assemblage = protein.assemblage
     AND qtl.protein_accession = protein.protein_accession
    ORDER BY
        registry.species,
        registry.assembly_accession,
        protein.protein_accession,
        (protein.plaza_gene_id IS NULL),
        protein.plaza_gene_id
)
SELECT
    species,
    assembly_accession,
    protein_accession,
    plaza_gene_id,
    gene_symbol,
    description,
    go_text,
    pfam_text,
    interpro_text,
    trait_text,
    qtl_evidence,
    concat_ws(
        ' ', gene_symbol, refseq_product, description, go_text, pfam_text, interpro_text,
        mapman_text, trait_text
    ) AS search_text,
    -- '.' et '/' remplacés par des espaces pour l'index : sinon « salt/drought » ou
    -- « response.drought.stomatal » deviennent un seul mot (chemin / nom d'hôte) et « drought » est introuvable.
    -- La requête de recherche applique la même transformation : translate(%s, './', '  ').
    to_tsvector(
        'simple',
        translate(
            concat_ws(
                ' ', gene_symbol, refseq_product, description, go_text, pfam_text, interpro_text,
                mapman_text, trait_text
            ),
            './', '  '
        )
    ) AS search_vector,
    mapman_text
FROM proteins;

CREATE UNIQUE INDEX v_gene_text_protein_uidx
    ON v_gene_text (species, assembly_accession, protein_accession);
CREATE INDEX v_gene_text_search_idx
    ON v_gene_text USING GIN (search_vector);
CREATE INDEX v_gene_text_species_idx
    ON v_gene_text (species);

-- After importing refreshed annotations/QTLs, run:
-- REFRESH MATERIALIZED VIEW CONCURRENTLY v_gene_text;
