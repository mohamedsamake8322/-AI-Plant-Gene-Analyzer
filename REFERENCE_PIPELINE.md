# Référence complète — dossiers `collect/` et `scripts/`

Ce document existe pour une seule raison : que tu n'aies plus jamais à te demander *"quel fichier je lance pour faire X ?"*. Chaque fichier est classé par ce qu'il fait, s'il est à lancer directement ou juste utilisé en interne par un autre script, et ce qui se passe concrètement quand tu le lances.

---

## ⚠️ Le point le plus important : il y a DEUX pipelines de collecte différents

| | `collect/` | `scripts/` |
|---|---|---|
| **Pipeline** | `collect_all_sources.py` — moderne, parallèle, multi-espèces | `run_pipeline.py` + `collect_plant_data.py` — historique, un fichier combiné |
| **Sortie** | Un fichier JSON **par espèce** (`data/clean/species/*.json`) | Un **seul** fichier JSON combiné (`plant_data_clean.json`) |
| **Statut** | ✅ **C'est celui que tu utilises activement** | 🔒 En réserve, ne pas lancer en parallèle du premier |
| **Import Postgres** | `scripts/import_species_to_postgres.py` | `scripts/load_to_postgres.py` |

**Règle d'or : ne jamais mélanger les deux dans le même cycle de collecte.** Le workflow à suivre est détaillé en bas de ce document.

---

## 📁 Dossier `collect/`

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| **`collect_all_sources.py`** | ✅ **Oui — c'est LE point d'entrée** | Orchestre toute la collecte : NCBI, UniProt, KEGG, PlantTFDB, PubMed, Expression Atlas, GEO, en parallèle sur plusieurs espèces. Écrit un fichier par espèce + fusionne en `master_plant_db.json` à la fin. |
| `collect_all_plants.ps1` / `.sh` | Optionnel | Simple wrapper qui appelle `collect_all_sources.py` avec des options prédéfinies. Pratique mais pas obligatoire — tu peux appeler `collect_all_sources.py` directement avec tes propres arguments (ce qu'on fait depuis le début). |
| `run_collect_with_log.py` / `.ps1` | Optionnel | Wrapper qui ajoute une journalisation automatique. Alternative à `Start-Transcript`. |
| `collect_uniprot.py` | ❌ Jamais seul | Utilisé automatiquement par `collect_all_sources.py` pour la source `uniprot`. Contient `fetch_uniprot(species, retmax)`. |
| `collect_kegg.py` | ❌ Jamais seul | Idem pour `kegg`. Contient `fetch_kegg(species, retmax)`. Corrigé récemment (pathways/KO conservés même sans séquence). |
| `collect_planttfdb.py` | ❌ Jamais seul | Idem pour `planttfdb`. Retourne des métadonnées de facteurs de transcription, **sans séquence** dans la plupart des cas (limite connue de PlantTFDB, pas un bug). |
| `collect_pubmed.py` | ❌ Jamais seul | Idem pour `pubmed`. Retourne un paquet de publications par espèce (pas des gènes individuels). Corrigé récemment (ne prétend plus avoir une séquence ADN). |
| `collect_atlas.py` | ❌ Jamais seul | Pont vers `scripts/collect_expression_atlas.py`. Retourne un pseudo-enregistrement par espèce regroupant les expériences Expression Atlas trouvées. |
| `collect_geon.py` | ❌ Jamais seul | Pont vers `scripts/collect_geo.py`. Même principe, pour les datasets NCBI GEO. |
| `collect_ensembl.py` ⚠️ | ❌ Jamais, et **actuellement cassé** | **À renommer en `collect_ensembl_stub.py`** — on l'avait décidé pour éviter une collision de nom avec `scripts/collect_ensembl.py` (deux fichiers différents, même nom). Tant qu'il n'est pas renommé, si tu demandes un jour `--sources ensembl` explicitement, ça plantera. Sans le demander, aucun impact (exclu du défaut). |
| `__pycache__/` | Ignorer | Cache Python auto-généré, jamais à toucher ni committer. |

---

## 📁 Dossier `scripts/`

### Collecte ciblée (usage ponctuel, pas pour la collecte en masse)

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| `collect_ncbi.py` | ✅ Pour une recherche NCBI ponctuelle | Un gène/accession précis, ou un terme de recherche limité. Pas pour collecter toute une espèce. |
| `collect_ensembl.py` | ✅ Pour un gène Ensembl précis | `fetch_gene(species, symbol, feature_id, seq_type)` — un seul gène par son symbole, pas de recherche bulk. |
| `collect_expression_atlas.py` | ✅ Pour chercher des expériences par mot-clé | `search_experiments()`, `gene_profile()`. Utilisé aussi en interne par `collect/collect_atlas.py`. |
| `collect_geo.py` | ✅ Pour chercher des datasets GEO par mot-clé | `search_geo()`, `fetch_geo_summaries()`. Utilisé aussi en interne par `collect/collect_geon.py`. |
| `collect_plant_data.py` | ✅ Wrapper de l'**autre** pipeline | Orchestre GEO+Ensembl+Atlas+NCBI en un seul fichier combiné. Utilisé par `run_pipeline.py`. |
| `download_data.py` | ✅ Optionnel | Collecte brute par source, sans nettoyage immédiat. Rarement nécessaire si tu utilises `collect_all_sources.py`. |

### Le moteur NCBI du pipeline principal

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| `collect_multi_type.py` | ✅ Possible seul, mais surtout utilisé par `collect/collect_all_sources.py` | **C'est lui qui gère la source `ncbi`** dans ton pipeline principal — DNA, RNA, protéines, une espèce à la fois. Contient le correctif `--max-length` qui empêche de récupérer des chromosomes entiers au lieu de gènes. |

### Nettoyage et transformation

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| `clean_data.py` | ✅ Pour l'autre pipeline uniquement | Nettoie/normalise le JSON brut de `run_pipeline.py`. **Pas nécessaire avec `collect_all_sources.py`**, qui nettoie déjà automatiquement. |
| `transform_schema.py` | Rarement | Convertit un ancien format JSON "legacy" vers le schéma professionnel enrichi. |
| `integrate_professional_schema.py` | Rarement | Guide/assistant de migration vers le schéma pro — plutôt un outil de conseil que d'exécution automatique. |
| `migrate_genes_db.py` | Rarement | Migre l'ancienne base locale `genes_database.json` vers un format standard. |
| `rebuild_master.py` | ❌ Pas nécessaire | Reconstruit `master_plant_db.json` depuis les fichiers espèces — mais `collect_all_sources.py` le fait déjà automatiquement à la fin de chaque run. |
| `rebuild_master_stream.py` | ❌ Pas nécessaire | Version streaming du même besoin, pour de très gros volumes. |

### Validation

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| **`verify_gene_records.py`** | ✅ **Oui, systématiquement après chaque collecte/import** | Lit depuis Postgres (`--db`) ou un JSON (`--json-file`), classe les enregistrements (dna/protein/metadata/...), génère CSV/JSON de rapport. Options `--csv-out`, `--stats-out` pour ne pas perdre le résumé dans un terminal trop long. |
| `validate_and_add_gene.py` | ❌ Utilisé en interne | Valide un enregistrement avant insertion dans `genes_database.json`. |
| `validate_bioinformatics.py` | ✅ Ponctuel | Teste les moteurs d'alignement/mutation/phylogénie eux-mêmes (pas tes données) — pertinent quand tu travailles sur `alignment_engine.py` etc., pas sur la collecte. |
| `test_postgres_connection.py` | ✅ Dès qu'il y a un doute de connexion | Diagnostique DNS/TCP/auth vers Postgres, étape par étape. Aurait pu nous faire gagner du temps sur le bug de port 5432/5433 ! À utiliser en premier réflexe la prochaine fois qu'une connexion échoue. |

### Base de données

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| **`import_species_to_postgres.py`** | ✅ **Oui — c'est l'import du pipeline principal** | Lit tous les `data/clean/species/*_all_sources.json` et les upsert dans Postgres. |
| `load_to_postgres.py` | ✅ Pour l'autre pipeline uniquement | Charge le JSON combiné de `run_pipeline.py`/`clean_data.py`. Pas utilisé par ton flux principal. |
| `postgres_utils.py` | ❌ Jamais seul | Bibliothèque partagée (connexion, upsert, création de tables) utilisée par tous les scripts d'import/lecture Postgres ci-dessus. |
| `db_utils.py` | ❌ Jamais seul | Équivalent mais pour la base JSON locale `genes_database.json` (pas Postgres). |

### Analyse (pas de la collecte — pour la suite : `alignment_engine.py` etc.)

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| `run_analysis_suite.py` | ✅ Plus tard | Statistiques, ORF, traduction, alignements, phylogénie sur les données déjà en base. Sujet de notre prochaine session, pas de celle-ci. |

---

## ✅ Le workflow exact à suivre, dans l'ordre

### 1. Collecter
```powershell
python collect/collect_all_sources.py --plant-file toutes_especes.txt --workers 4 --retmax 300 --create-tables
```
**Ce que ça fait** : lance 4 collectes en parallèle (une par espèce à la fois par worker), interroge NCBI/UniProt/KEGG/PlantTFDB/PubMed/Atlas/GEO pour chaque espèce, écrit un fichier JSON par espèce dans `data/clean/species/`, puis fusionne tout dans `data/clean/master_plant_db.json` à la fin. `--create-tables` s'assure que la table Postgres existe (sans importer les données elles-mêmes).

### 2. Nettoyer
**Rien à faire** — déjà inclus dans l'étape 1.

### 3. Valider avant l'import massif
```powershell
.\.venv\Scripts\python.exe .\scripts\verify_gene_records.py --json-file data\clean\master_plant_db.json --stats-out data\clean\stats_avant_import.json
```
**Ce que ça fait** : lit le JSON fusionné, classe chaque enregistrement, écrit un résumé. Sert à repérer un problème (comme les chromosomes entiers de la dernière fois) **avant** de polluer la base Postgres avec.

### 4. Envoyer vers Postgres
```powershell
.\.venv\Scripts\python.exe .\scripts\import_species_to_postgres.py --species-dir data\clean\species
```
**Ce que ça fait** : lit chaque fichier espèce, upsert chaque gène dans la table `genes` (met à jour si le gène existe déjà, insère sinon — sans écraser une séquence existante par du vide, grâce au correctif qu'on a fait).

### 5. Vérifier le résultat final
```powershell
.\.venv\Scripts\python.exe .\scripts\verify_gene_records.py --db --csv-out data\clean\rapport_final.csv --stats-out data\clean\stats_final.json
```
**Ce que ça fait** : relit directement depuis Postgres (source de vérité), confirme que ce qui est en base correspond à ce qui a été collecté.

---

## 🚫 Ce que tu ne dois jamais faire

- Lancer `run_pipeline.py` **en plus** de `collect_all_sources.py` pour la même collecte — ce sont deux systèmes séparés qui ne se complètent pas, ça duplique le travail.
- Lancer `verify_gene_records.py` sans `--json-file` ni `--db` en pensant lire tes données actuelles — il retombe par défaut sur `data/clean/plant_data_clean.json`, un fichier de test périmé s'il existe encore.
- Modifier `collect_uniprot.py`, `collect_kegg.py`, etc. dans `collect/` en pensant que ça affecte l'usage "ponctuel" — ce sont des fichiers différents de leurs équivalents dans `scripts/`, malgré des noms proches.
