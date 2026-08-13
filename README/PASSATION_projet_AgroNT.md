# État du projet & prochaines étapes — Module IA Interprétable (AgroNT)

*Document de passation — à donner en contexte si tu ouvres une nouvelle conversation.*

## 1. Où on en est, en une phrase

AgroNT est chargé et fonctionnel sur Kaggle (991 261 206 paramètres confirmés),
connecté à un jeu de données de **1 207 gènes complets** (séquence nucléotidique
+ annotation fonctionnelle réelle), construit spécifiquement pour corriger un
défaut de conception découvert dans la collecte initiale. Prochaine étape
immédiate : extraire les embeddings et visualiser leur cohérence avant de
lancer le fine-tuning.

## 2. Contexte du projet (rappel)

- Mémoire de master : *"Développement d'une plateforme bioinformatique modulaire
  pour l'intégration de données génomiques et l'analyse fonctionnelle des
  gènes chez les plantes d'importance agricole."*
- Application : Plant Gene Analyzer (Streamlit), en production
  (`ai-plant-gene-analyzer.streamlit.app`), connectée à une base PostgreSQL
  (Neon, tier gratuit 512 Mo).
- Module actuel en travaux : **AI Interpretation** — objectif : passer d'un
  système à base de règles (actuel) à un système hybride : backbone AgroNT
  (fine-tuné IA3/LoRA) + têtes de classification (GO terms, traits, TF family,
  pathways) + génération de texte ancrée (RAG) sur ces prédictions réelles.
- Stratégie décidée : **qualité plutôt que quantité**. Un petit jeu de données
  propre et équilibré plutôt que les 55 966 gènes bruts collectés initialement.

## 3. Ce qui a été découvert et corrigé (dans l'ordre chronologique)

| # | Problème découvert | Diagnostic | Correction |
|---|---|---|---|
| 1 | Crise de stockage Neon (494/512 Mo) | Table `gene_kmers` (288 Mo) obsolète depuis la migration vers `pg_trgm`, jamais nettoyée | `TRUNCATE gene_kmers` -> 206 Mo. Confirmé sans risque : `pg_trgm` est la voie réellement active (`source: trigram_prefilter`). |
| 2 | Pas de dédoublonnage réel | `postgres_utils.py` dédoublonnait par `gene_id`/`symbol`, jamais par contenu de séquence | Ajout de `sequence_hash` (md5) + `dedupe_by_sequence()`, redirige vers l'enregistrement existant si séquence identique. |
| 3 | Filtre qualité trop strict (2 itérations) | Rejetait les séquences `rna` contenant `T` (convention NCBI standard, pas une erreur) et les codes IUPAC d'ambiguïté (R/Y/S/W/K/M/B/D/H/V) | Alphabet nucléotidique unifié acceptant T+U et tous les codes IUPAC, dans **`postgres_utils.py` ET `build_core_dataset.py`** (les deux endroits, ne pas oublier). |
| 4 | Défaut de conception majeur : 0% de chevauchement entre séquences nucléotidiques (NCBI) et annotations fonctionnelles (UniProt/PlantTFDB) | La collecte "par source" (tout l'ADN d'un côté, toutes les protéines annotées de l'autre) produit deux ensembles de gènes disjoints. Confirmé empiriquement. | **Changement de stratégie : collecte "par gène".** Nouveau script `collect/collect_linked_genes.py` : part d'une protéine déjà annotée (UniProt/PlantTFDB), récupère sa référence croisée NCBI, va chercher CETTE séquence précise. Ne garde le gène que si les deux réussissent. |
| 5 | `collect_uniprot.py` demandait `xref_refseq` à l'API mais ne l'extrayait jamais | Bug de code, pas un manque de données | Extraction de `refseq_nucleotide`/`refseq_protein`/`embl_nucleotide` depuis les références croisées UniProt. |
| 6 | Résolution NCBI par accession échouait pour certains formats (`XM_...`) | Le tag Entrez `[Accession]` échoue silencieusement sur certaines accessions RefSeq "predicted", alors qu'`efetch` direct les trouve sans problème | `fetch_fasta_by_accession()` tente `efetch` direct en premier, `esearch` en repli seulement. |
| 7 | Séquences aberrantes (chromosomes/génomes entiers, jusqu'à 45M pb) | Un gène sans entrée ARNm séparée (souvent mitochondrial/chloroplastique, ex: NAD7) retombe sur le génome entier de l'organite | Seuil de longueur resserré à 50 000 pb (au lieu de 500 000) dans `collect_linked_genes.py`. |
| 8 | Résolution de locus PlantTFDB -> NCBI ne généralise pas | Fonctionne pour Arabidopsis (locus `AT1G01010`, motif reconnu), échoue à 100% (0/8 testés) pour Solanum lycopersicum, Glycine max, Zea mays | Limite assumée et documentée. Pour ces espèces : `--no-tf`, collecte via UniProt seule (suffisant pour go_terms/traits/pathways). |

## 4. Résultat de la collecte (jeu de données final)

Fichier : **`Data/clean/linked_genes_final_clean.json`** (5,95 Mo) — celui à
utiliser partout désormais (passé par `clean_data.py`, cohérent avec ce qui
est déjà chargé en Postgres).

| Métrique | Valeur |
|---|---|
| Total (après dédoublonnage) | 1 207 gènes |
| Arabidopsis thaliana | 474 |
| Zea mays | 256 |
| Oryza sativa | 252 |
| Glycine max | 225 |
| Couverture go_terms | 85,0% |
| Couverture traits | 100% |
| Couverture pathways | 73,8% |
| Couverture tf_family | 15,0% (concentré sur Arabidopsis — limite connue, voir §3.8) |

Déjà chargé dans PostgreSQL (coexiste avec l'ancienne base de 55 986 gènes,
non supprimée — l'appli en production en dépend toujours). Les nouveaux
enregistrements sont identifiables via `source LIKE '%_linked'`.

Document détaillé complet : `README_collecte_liee.md` (fourni séparément,
contient tous les tableaux et le raisonnement complet).

## 5. État technique actuel (Kaggle)

- Notebook : `agroNT_zeroshot_exploration.ipynb` (fourni séparément), importé
  sur Kaggle sous le nom "Embedding AgroNT".
- Dataset Kaggle : contient `linked_genes_final_clean.json`.
- Modèle chargé avec succès : `InstaDeepAI/agro-nucleotide-transformer-1b`,
  991 261 206 paramètres, GPU actif.
- Warning vu au chargement : "unauthenticated requests to the HF Hub" — pas
  bloquant, mais pour éviter les limites de débit sur de futurs
  téléchargements (fine-tuning avec plus d'itérations), envisager de créer un
  token HuggingFace (gratuit) et de l'ajouter comme Kaggle Secret (HF_TOKEN).
- AgroNT : fenêtre de contexte 1024 tokens (~6144 pb, tokenizer 6-mers). Des
  séquences plus longues (ex: MDN1, 16 444 pb) sont tronquées automatiquement
  — compromis accepté pour le test zero-shot, à réévaluer si ça pénalise le
  fine-tuning.

## 6. Prochaines étapes immédiates (dans le notebook)

1. Exécuter les cellules restantes : extraction des embeddings (mean pooling
   sur la dernière couche cachée), sauvegarde (.npy + .csv métadonnées).
2. Visualisation UMAP colorée par organisme et par tf_family.
3. Lire le résultat :
   - Regroupement net par organisme -> normal, pas forcément informatif seul.
   - Un minimum de séparation par tf_family -> signal encourageant pour le
     fine-tuning.
   - Aucune structure -> pas alarmant sur un échantillon aussi petit en
     zero-shot, confirme juste que le fine-tuning est nécessaire (attendu).

## 7. Plan pour la suite (fine-tuning) — pas encore commencé

1. Fine-tuning IA3/LoRA d'AgroNT sur les 1 207 gènes (paramétrique-efficace,
   recommandé par les auteurs d'AgroNT eux-mêmes, adapté à un quota GPU
   Kaggle limité : ~30h/semaine, 12h/session max).
2. Têtes de classification séparées par tâche (un seul backbone partagé) :
   - GO terms (multi-label) — 85% de couverture, meilleur socle
   - Pathways — 73,8%
   - Traits (à consolider — catégories probablement hétérogènes à nettoyer
     avant classification, jamais fait à ce stade)
   - TF family — 15%, entraînable mais concentré sur Arabidopsis ; à
     documenter comme limite dans le mémoire plutôt qu'à cacher
3. Split train/val/test : avec seulement 1 207 gènes, le risque de fuite par
   homologie (deux séquences quasi-identiques dans train ET test) est réduit
   par rapport à l'ancienne base de 55k, mais pas nul — vérifier avant de
   conclure sur les métriques.
4. Couche RAG (génération de texte ancrée) : à construire après avoir des
   prédictions réelles des têtes de classification — utilise les métadonnées
   déjà en base (traçabilité complète via external_links).

## 8. Fichiers clés du projet (inventaire)

Collecte :
- `collect/collect_linked_genes.py` — nouveau collecteur "par gène" (le plus important)
- `collect/collect_uniprot.py` — corrigé (extraction refseq_nucleotide)
- `scripts/collect_ncbi.py` — corrigé (efetch direct en premier)
- `collect/collect_planttfdb.py` — non modifié (fonctionne pour Arabidopsis seulement)

Nettoyage / préparation :
- `scripts/merge_linked_genes.py` — fusion multi-espèces + dédoublonnage
- `scripts/postgres_utils.py` — corrigé (hash, filtre qualité, migration)
- `scripts/load_to_postgres.py` / `scripts/import_species_to_postgres.py` — corrigés
- `scripts/check_db_size.py`, `scripts/export_genes_backup.py` — diagnostics Postgres

Diagnostics ponctuels (jetables, gardés pour référence) :
`diagnose_invalid_chars.py`, `diagnose_label_type_split.py`, `diagnose_join_key.py`,
`diagnose_oryza_resolution.py`, `diagnose_efetch_direct.py`, `probe_species_readiness.py`,
`test_ncbi_resolution_candidates.py`, `inspect_linked_genes.py`

Modélisation :
- `agroNT_zeroshot_exploration.ipynb` — notebook Kaggle actuel

Documentation :
- `README_collecte_liee.md` — historique complet de la collecte (détaillé)
- ce fichier — passation/état général

## 9. Pièges déjà rencontrés — ne pas les refaire

- Ne jamais supposer qu'un filtre qualité "raisonnable sur le papier" est
  correct sans le vérifier sur un échantillon réel (2 itérations ont été
  nécessaires pour le filtre de caractères valides).
- Le clustering d'homologie ou toute boucle appelant find_similar_genes/
  fetch_fasta_by_accession en masse doit être testé sur un petit échantillon
  avant un run complet — le coût par appel réseau (Postgres, NCBI) est
  largement sous-estimable à l'oeil.
- NCBI peut limiter/couper les connexions en cas de sollicitation intensive
  prolongée — augmenter --api-delay et espacer les runs si ça arrive plutôt
  que de re-tenter en boucle.
- Toujours vérifier qu'une correction de bug appliquée dans un script l'est
  aussi dans les autres endroits où la même logique est dupliquée (le filtre
  qualité existait en double, dans postgres_utils.py ET
  build_core_dataset.py — la correction a été oubliée dans l'un des deux la
  première fois).
