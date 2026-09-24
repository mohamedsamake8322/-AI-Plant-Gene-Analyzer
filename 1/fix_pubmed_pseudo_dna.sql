-- =============================================================================
-- Nettoyage rétroactif : pseudo-enregistrements PubMed comptés comme ADN
-- =============================================================================
-- Bug corrigé dans collect_pubmed.py : publications_to_gene_record() posait
-- sequence='N' et sequence_type='dna' comme placeholder, ce qui faisait
-- classer ces enregistrements bibliographiques comme de vrais gènes ADN
-- par verify_gene_records.py (et tout autre outil qui se fie à sequence_type).
--
-- Cette requête ne touche QUE les enregistrements correspondant exactement
-- à l'empreinte du bug (source='pubmed' AND sequence='N' AND length=1),
-- pas l'ensemble des enregistrements PubMed -- si un jour un vrai gène
-- avait légitimement 'N' comme séquence d'1 caractère (cas quasi
-- impossible biologiquement), il ne serait pas non plus concerné puisque
-- source doit aussi être 'pubmed'.
-- =============================================================================

-- ─────────────────────────────────────────────────────────────────────────
-- ÉTAPE 1 : PRÉVISUALISATION -- à lancer et vérifier avant l'UPDATE
-- ─────────────────────────────────────────────────────────────────────────
SELECT COUNT(*) AS affected_rows
FROM genes
WHERE source = 'pubmed' AND sequence = 'N' AND length = 1;

SELECT gene_id, symbol, organism, sequence, sequence_type, length, source
FROM genes
WHERE source = 'pubmed' AND sequence = 'N' AND length = 1
LIMIT 10;

-- ─────────────────────────────────────────────────────────────────────────
-- ÉTAPE 2 : CORRECTIF -- ne lancer qu'après avoir vérifié l'étape 1
-- ─────────────────────────────────────────────────────────────────────────
UPDATE genes
SET
    sequence = '',
    sequence_type = '',
    length = 0,
    annotations = COALESCE(annotations, '{}'::jsonb)
        || '{"sequence_available": false, "record_type": "publication_bundle"}'::jsonb,
    record = jsonb_set(
        jsonb_set(
            jsonb_set(
                jsonb_set(
                    COALESCE(record, '{}'::jsonb),
                    '{sequence}', '""'::jsonb
                ),
                '{sequence_type}', '""'::jsonb
            ),
            '{length}', '0'::jsonb
        ),
        '{annotations}',
        (COALESCE(record->'annotations', '{}'::jsonb)
            || '{"sequence_available": false, "record_type": "publication_bundle"}'::jsonb)
    )
WHERE source = 'pubmed' AND sequence = 'N' AND length = 1;

-- ─────────────────────────────────────────────────────────────────────────
-- ÉTAPE 3 : VÉRIFICATION -- doit retourner 0
-- ─────────────────────────────────────────────────────────────────────────
SELECT COUNT(*) AS remaining_buggy_rows
FROM genes
WHERE source = 'pubmed' AND sequence = 'N' AND length = 1;

-- Et pour confirmer que les enregistrements corrigés sont bien repérables :
SELECT COUNT(*) AS fixed_publication_bundles
FROM genes
WHERE annotations->>'record_type' = 'publication_bundle';
