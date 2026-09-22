-- diag_genes_2.sql : 2e diagnostic, LECTURE SEULE (aucune modification, aucune sequence affichee)
-- Lancer depuis PowerShell :
--   psql -h localhost -p 5432 -U postgres -d plant_gene_analyzer_clean -f diag_genes_2.sql -o diag_output_2.txt

\pset pager off
\pset null '(null)'

\qecho === A. Codes de preuve GO (evidence_code) par sequence_type
SELECT g.sequence_type, e->>'evidence_code' AS code, count(*) AS n_annotations
FROM genes g,
     LATERAL jsonb_array_elements(
         CASE WHEN jsonb_typeof(g.annotations->'go_terms') = 'array'
              THEN g.annotations->'go_terms' ELSE '[]'::jsonb END
     ) AS e
GROUP BY 1, 2
ORDER BY 1, 3 DESC;

\qecho === B. Genes ayant au moins une annotation GO non automatique (code different de IEA)
SELECT g.sequence_type, count(*) AS n_genes
FROM genes g
WHERE EXISTS (
    SELECT 1
    FROM jsonb_array_elements(
             CASE WHEN jsonb_typeof(g.annotations->'go_terms') = 'array'
                  THEN g.annotations->'go_terms' ELSE '[]'::jsonb END
         ) AS e
    WHERE e->>'evidence_code' IS DISTINCT FROM 'IEA'
)
GROUP BY 1
ORDER BY 1;

\qecho === C. Origine des annotations (evidence_source), top 15
SELECT e->>'evidence_source' AS evidence_source, count(*) AS n
FROM genes g,
     LATERAL jsonb_array_elements(
         CASE WHEN jsonb_typeof(g.annotations->'go_terms') = 'array'
              THEN g.annotations->'go_terms' ELSE '[]'::jsonb END
     ) AS e
GROUP BY 1
ORDER BY 2 DESC
LIMIT 15;

\qecho === D. Familles PLAZA : couverture parmi les genes avec GO
SELECT sequence_type,
       count(*) AS n_genes_go,
       count(*) FILTER (WHERE relations->>'homologous_family_id' IS NOT NULL) AS avec_famille_homologue,
       count(DISTINCT relations->>'homologous_family_id') AS n_familles_homologues
FROM genes
WHERE jsonb_typeof(annotations->'go_terms') = 'array'
  AND jsonb_array_length(annotations->'go_terms') > 0
GROUP BY 1;

\qecho === E. Plus grandes familles homologues (genes avec GO)
SELECT relations->>'homologous_family_id' AS famille, count(*) AS n
FROM genes
WHERE jsonb_typeof(annotations->'go_terms') = 'array'
  AND jsonb_array_length(annotations->'go_terms') > 0
  AND relations->>'homologous_family_id' IS NOT NULL
GROUP BY 1
ORDER BY 2 DESC
LIMIT 10;

\qecho === F. Formats d identifiants des genes avec GO (A = lettre, 9 = chiffre)
SELECT sequence_type,
       regexp_replace(regexp_replace(gene_id, '[0-9]', '9', 'g'), '[A-Za-z]', 'A', 'g') AS motif,
       count(*) AS n,
       min(gene_id) AS exemple
FROM genes
WHERE jsonb_typeof(annotations->'go_terms') = 'array'
  AND jsonb_array_length(annotations->'go_terms') > 0
GROUP BY 1, 2
ORDER BY 3 DESC
LIMIT 15;

\qecho === G. Plage des id et nombre de date_added renseignes
SELECT min(id) AS id_min, max(id) AS id_max, count(*) AS n,
       count(date_added) AS avec_date_added
FROM genes;

\qecho === H. Dates de recuperation des GO terms (retrieved_at)
SELECT min(e->>'retrieved_at') AS plus_ancien, max(e->>'retrieved_at') AS plus_recent
FROM genes g,
     LATERAL jsonb_array_elements(
         CASE WHEN jsonb_typeof(g.annotations->'go_terms') = 'array'
              THEN g.annotations->'go_terms' ELSE '[]'::jsonb END
     ) AS e;