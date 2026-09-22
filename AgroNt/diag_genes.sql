-- diag_genes.sql : diagnostic en LECTURE SEULE (aucune modification, aucune sequence affichee)
-- Lancer depuis PowerShell :
--   psql -h localhost -p 5432 -U postgres -d plant_gene_analyzer_clean -f diag_genes.sql -o diag_output.txt
-- Si une requete echoue, psql continue avec la suivante : ce n'est pas grave.

\pset pager off
\pset null '(null)'

\qecho === 1. Colonnes et types de la table genes
SELECT column_name, data_type
FROM information_schema.columns
WHERE table_name = 'genes'
ORDER BY ordinal_position;

\qecho === 2. Nombre de genes et taille de la table
SELECT count(*) AS n_genes,
       pg_size_pretty(pg_total_relation_size('genes')) AS taille_totale
FROM genes;

\qecho === 3. Genes par sequence_type
SELECT sequence_type, count(*) AS n
FROM genes
GROUP BY 1
ORDER BY 2 DESC;

\qecho === 4. Genes par organisme et sequence_type (top 20)
SELECT organism, sequence_type, count(*) AS n
FROM genes
GROUP BY 1, 2
ORDER BY 3 DESC
LIMIT 20;

\qecho === 5. Croissance : ajouts par mois et sequence_type (mois les plus recents)
SELECT left(date_added::text, 7) AS mois, sequence_type, count(*) AS n
FROM genes
GROUP BY 1, 2
ORDER BY 1 DESC, 2
LIMIT 40;

\qecho === 6. Provenance : source et origin par sequence_type
SELECT source, origin, sequence_type, count(*) AS n
FROM genes
GROUP BY 1, 2, 3
ORDER BY 4 DESC
LIMIT 30;

\qecho === 7. Formats d identifiants (A = lettre, 9 = chiffre)
SELECT sequence_type,
       regexp_replace(regexp_replace(gene_id, '[0-9]', '9', 'g'), '[A-Za-z]', 'A', 'g') AS motif,
       count(*) AS n,
       min(gene_id) AS exemple
FROM genes
GROUP BY 1, 2
ORDER BY 3 DESC
LIMIT 25;

\qecho === 8. Genes avec GO terms, par sequence_type et organisme
SELECT sequence_type, organism, count(*) AS n
FROM genes
WHERE jsonb_typeof(annotations->'go_terms') = 'array'
  AND jsonb_array_length(annotations->'go_terms') > 0
GROUP BY 1, 2
ORDER BY 3 DESC
LIMIT 20;

\qecho === 9. Cles presentes dans la colonne annotations
SELECT k AS cle, count(*) AS n
FROM genes,
     LATERAL jsonb_object_keys(
         CASE WHEN jsonb_typeof(annotations) = 'object' THEN annotations ELSE '{}'::jsonb END
     ) AS k
GROUP BY 1
ORDER BY 2 DESC
LIMIT 30;

\qecho === 10. Cles a l interieur des elements go_terms
SELECT k AS cle, count(*) AS n
FROM genes g,
     LATERAL jsonb_array_elements(
         CASE WHEN jsonb_typeof(g.annotations->'go_terms') = 'array'
              THEN g.annotations->'go_terms' ELSE '[]'::jsonb END
     ) AS e,
     LATERAL jsonb_object_keys(
         CASE WHEN jsonb_typeof(e) = 'object' THEN e ELSE '{}'::jsonb END
     ) AS k
GROUP BY 1
ORDER BY 2 DESC;

\qecho === 11. Repartition des aspects et des codes de preuve (evidence) GO
SELECT e->>'aspect' AS aspect, e->>'evidence' AS evidence, count(*) AS n
FROM genes g,
     LATERAL jsonb_array_elements(
         CASE WHEN jsonb_typeof(g.annotations->'go_terms') = 'array'
              THEN g.annotations->'go_terms' ELSE '[]'::jsonb END
     ) AS e
GROUP BY 1, 2
ORDER BY 3 DESC
LIMIT 25;

\qecho === 12. Exemple de go_terms (3 genes)
SELECT gene_id, (annotations->'go_terms')->0 AS premier_go_term
FROM genes
WHERE jsonb_typeof(annotations->'go_terms') = 'array'
  AND jsonb_array_length(annotations->'go_terms') > 0
LIMIT 3;

\qecho === 13. Faisabilite du mapping proteine -> ADN : references croisees
SELECT sequence_type,
       count(*) AS n,
       count(*) FILTER (WHERE external_links::text ~* 'embl|ensembl|refseq|genbank|ncbi|phytozome|gramene|maizegdb') AS avec_lien_externe,
       count(*) FILTER (WHERE relations::text ~* 'cds|mrna|transcript|ensembl|refseq|embl') AS avec_relation
FROM genes
GROUP BY 1;

\qecho === 14. Echantillon external_links et relations (5 proteines et 5 ADN avec GO)
(SELECT sequence_type, gene_id,
        left(external_links::text, 300) AS external_links,
        left(relations::text, 300) AS relations
 FROM genes
 WHERE sequence_type = 'protein'
   AND jsonb_typeof(annotations->'go_terms') = 'array'
   AND jsonb_array_length(annotations->'go_terms') > 0
 LIMIT 5)
UNION ALL
(SELECT sequence_type, gene_id,
        left(external_links::text, 300),
        left(relations::text, 300)
 FROM genes
 WHERE sequence_type = 'dna'
   AND jsonb_typeof(annotations->'go_terms') = 'array'
   AND jsonb_array_length(annotations->'go_terms') > 0
 LIMIT 5);

\qecho === 15. Symboles presents a la fois en ADN et en proteine (piste de correspondance)
SELECT count(*) AS symboles_communs
FROM (
    SELECT symbol
    FROM genes
    WHERE symbol IS NOT NULL AND symbol <> ''
    GROUP BY symbol
    HAVING count(DISTINCT sequence_type) = 2
) t;

\qecho === 16. Doublons de sequence (sequence_hash)
SELECT count(*) AS total, count(DISTINCT sequence_hash) AS hash_distincts
FROM genes;

\qecho === 17. Longueur des sequences par type
SELECT sequence_type,
       count(*) AS n,
       min(length) AS min,
       percentile_disc(0.5) WITHIN GROUP (ORDER BY length) AS mediane,
       max(length) AS max
FROM genes
GROUP BY 1;