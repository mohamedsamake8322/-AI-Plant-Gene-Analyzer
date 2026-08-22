# Récap projet — Plant Gene Analyzer (mémoire de master)

*Dernière mise à jour : 19/08/2026. À coller en tout premier message d'une nouvelle discussion avec Claude, avant tout fichier `.py`, pour reprendre exactement où on en est.*

## Contexte général

Plateforme bioinformatique Streamlit (`Plant Gene Analyzer`). Backend séparé (`collect/` + `scripts/`) qui alimente PostgreSQL (Neon) via `collect_all_sources.py`. Objectif mémoire : relier un gène à un problème agronomique réel (cas d'usage principal : **verse/lodging chez le quinoa**), avec preuves croisées et sourcées — pas juste comparer des séquences.

## Les 7 cultures finales

**Quinoa** (cas d'usage principal), **riz**, **maïs**, **tomate**, **raisin** (a remplacé le blé, non couvert par PLAZA), **tabac** (intérêt académique, génome hexaploïde-ish), **pomme de terre**.

## Schéma de données — décision structurante (toujours valide, ne pas changer)

Structure en étoile, pas en chaîne : chaque gène a `sequence` (dna/rna/protein), `annotation` (go_terms/kegg_pathways/tf_family/mapman), `traits`, `relations` (orthologs/family IDs), `literature`, `sources_summary`, `origin` (`sequence_backed` vs `plaza_only`), `quality.data_completeness`. Chaque fait individuel porte sa propre `source` + `retrieved_at` — pas de champ `source` global unique.

**Une refonte alternative a été proposée récemment** (`gene_record_schema.md`, organism en objet structuré + provenance séparée) — **rejetée** : le diagnostic de cause racine du document était faux (le bug n'était pas un problème de forme de données, mais de logique d'extraction, déjà corrigé à la source). Ne pas migrer vers cette structure sauf nouvelle raison concrète.

## PLAZA — architecture finale (stable, ne plus retoucher sauf besoin réel)

- Fichiers bruts multi-espèces (`HOMFAM`/`ORTHOFAM`) trop lourds à parser en local → **extraction en 2 passes sur Colab**, résultat mis en cache par espèce (`plaza_<code>_cached.json`)
- **`TARGET_CROP_CODES`** filtre les orthologues gardés aux 7 cultures cibles seulement → cache réduit de ~90% (quelques Go → 100-300 Mo par espèce)
- `CACHE_DIR` séparé de `PLAZA_DIR` (caches filtrés dans un dossier neuf, CSV bruts intacts ailleurs)
- **Codes PLAZA confirmés** : `cqu` (quinoa), `osa` (riz), `zma` (maïs), `sly` (tomate), `vvi` (raisin), `nta` (tabac), `stu` (pomme de terre)
- **Mode `plaza_only`** : gènes PLAZA sans correspondance NCBI/UniProt créés quand même (préfixe `PLAZA:`, `origin: "plaza_only"`, séquence `null`) — couverture 100% assumée, pas cachée
- **Matching fiable** : priorité à l'accession UniProt (`id_conversion.<code>.csv`), repli sur symbole normalisé seulement si échec
- ⚠️ **Piège de nommage à retenir** : le code KEGG du quinoa est `cqi`, **différent** du code PLAZA `cqu` — confusion facile

## Bugs trouvés et corrigés dans cette session (chronologique)

| # | Fichier | Bug | Statut |
|---|---|---|---|
| 1 | `collect_all_sources.py` | DNA/RNA/protein d'un même gène s'écrasaient (`if gid not in all_records`) | ✅ Corrigé |
| 2 | `collect_ncbi.py` | Chromosomes entiers (jusqu'à 100+ Mo) téléchargés puis rejetés après coup | ✅ Corrigé — pré-filtre par taille via `esummary` avant téléchargement |
| 3 | `collect_ncbi.py` | `gene_id` = accession brute → DNA/RNA/protein d'un même gène ne partagent jamais le même ID | 🟡 Tentative de fix via résolution `elink` → Entrez GeneID partagé — **efficacité non confirmée** (dépend de la richesse de curation NCBI par espèce, ex: quinoa peu curaté) |
| 4 | `collect_uniprot.py` | Pagination cassée : `Link` header découpé par virgule, alors que le paramètre `fields=` contient lui-même des virgules → URL de page suivante tronquée | ✅ Corrigé (regex au lieu de split naïf) |
| 5 | `postgres_utils.py` / `load_to_postgres.py` | Toute la couche PostgreSQL supposait l'ancien schéma plat (`sequence` = string) → crash garanti sur le nouveau schéma imbriqué (dict) | ✅ Corrigé — `extract_primary_sequence()`, colonnes `origin`/`relations` ajoutées, gènes `plaza_only` ne sont plus rejetés comme "vides" |
| 6 | `collect_kegg.py` | `KEGG_ORG_MAP` incomplet — tabac et quinoa absents du dictionnaire (pas une vraie absence KEGG) | ✅ Corrigé **par l'utilisateur lui-même** (`nta` pour tabac, `cqi` pour quinoa) |
| 7 | `download_data.py` | Fonction locale `collect_ncbi()` portant le même nom que le module importé `collect_ncbi` → écrasait la référence, plantait à l'exécution | ✅ Corrigé (renommée `collect_ncbi_records`) |
| 8 | `collect_ncbi.py` | `organism` deviné par regex depuis le texte du header FASTA — échouait sur les soumissions GenBank directes dont le header ne commence pas par le nom d'espèce (ex: "Molybdenum cofactor..." pris pour un organisme) | ✅ Corrigé dans le code — utilise maintenant `known_organism` (nom d'espèce déjà connu du contexte d'appel) au lieu de deviner. **⚠️ Fix pas encore vérifié bout en bout — voir "À faire" ci-dessous** |

## Fichiers dans leur état final (tous livrés, tous testés unitairement)

`collect_all_sources.py`, `collect_plaza.py`, `collect_ncbi.py`, `collect_plant_data.py`, `download_data.py`, `collect_uniprot.py`, `postgres_utils.py`, `load_to_postgres.py` — plus deux nouveaux utilitaires : `check_master_db_metadata.py` (lecture rapide de `metadata` sans charger tout le fichier), `repair_organism_field.py` (répare `organism` sur les données déjà collectées, sans recollecter).

**`collect_kegg.py`** : corrigé par l'utilisateur lui-même, pas de version livrée par Claude à re-uploader — juste vérifié.

## État actuel des données

- **319 741 gènes collectés sur 7 espèces** (quinoa, riz, maïs, tomate, raisin, tabac, pomme de terre) dans `data/clean/master_plant_db.json`
- **Pollution `organism` détectée** (`total_species: 109` au lieu de 7) — **pas corrompu**, juste un champ pollué sur une fraction des enregistrements NCBI (protéines surtout). Réparable via `repair_organism_field.py`, **sans recollecter**.
- **Majorité des gènes sont `plaza_only`** (~85-95% selon l'espèce) — c'est voulu (choix de couverture 100% PLAZA), mais à garder en tête pour toute analyse downstream (filtrer sur `origin` avant BLAST/Similarity)

## Neon (PostgreSQL) — point bloquant non résolu

**Quota mensuel de transfert réseau épuisé (5 Go/mois, tier gratuit)** — repéré via un bandeau "Limit reached" dans la console Neon. Cause probable identifiée mais **pas encore corrigée** : `dedupe_by_sequence()` fait une requête `SELECT` vers la base **pour chaque gène, un par un**, avant chaque insertion — avec 319 741 gènes, ça fait beaucoup d'allers-retours réseau. **Prochaine étape technique à faire** : réécrire cette logique pour ne faire qu'une seule requête groupée au début (récupérer tous les hash de séquence existants une fois, comparer localement en Python), plutôt qu'une requête par gène.

## Limitations connues, non bloquantes, pas de plan de fix immédiat

- **NCBI, UniProt, KEGG, PLAZA ne fusionnent pas entre eux** pour un même gène biologique — chacun a son propre espace d'identifiants (accession NCBI, accession UniProt, ID PLAZA, `org:gene_id` KEGG). Seul PLAZA↔UniProt fusionne de façon fiable (via `id_conversion`). Un vrai fix nécessiterait un travail de cross-référencement plus large — mis de côté pour l'instant.
- **`expression_profiles` (GEO/Expression Atlas) disparaît silencieusement** pendant `restructure_to_schema()` — le nouveau schéma imbriqué n'a pas de champ pour ça. Sans conséquence pour l'instant (GEO renvoie toujours 0 dans les runs), mais à corriger si GEO est activé un jour.

## À faire, dans l'ordre

1. **Envoyer `collect_multi_type.py` et `run_pipeline.py`** — nécessaire pour vérifier que le paramètre `organism` (utilisé par le fix #8 ci-dessus) est bien transmis jusqu'à `collect_ncbi.py`, sinon le fix ne prendra pas effet sur les futures collectes
2. **Lancer `repair_organism_field.py`** sur les vraies données, vérifier que `total_species` redescend à 7
3. **Corriger le problème de quota Neon** (regrouper les vérifications de doublons en une seule requête au lieu d'une par gène) avant tout rechargement complet
4. **Recharger PostgreSQL** une fois les données réparées et le quota débloqué
5. **Construire la table de traits manuelle verse/quinoa** (15-30 gènes candidats, sourcés PubMed) — **la vraie pièce maîtresse du mémoire**, reportée depuis le début de tout le travail PLAZA/NCBI/PostgreSQL. Rien ne bloque structurellement ce chantier — c'est juste resté en attente pendant qu'on stabilisait le pipeline de collecte.
