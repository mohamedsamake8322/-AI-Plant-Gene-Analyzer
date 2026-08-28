# Référence complète — dossiers `collect/` et `scripts/`

Ce document existe pour une seule raison : que tu n'aies plus jamais à te demander *"quel fichier je lance pour faire X ?"*. Chaque fichier est classé par ce qu'il fait, s'il est à lancer directement ou juste utilisé en interne par un autre script, et ce qui se passe concrètement quand tu le lances.

---

## ⚠️ Le point le plus important : il y a DEUX pipelines de collecte différents

| | `collect/` | `scripts/` |
|---|---|---|
| **Pipeline** | `collect_all_sources.py` — moderne, parallèle, multi-espèces | `run_pipeline.py` + `collect_plant_data.py` — historique, un fichier combiné |
| **Sortie** | Un fichier JSON **par espèce** (`data/clean/species/*.json`) | Un **seul** fichier JSON combiné (`plant_data_clean.json`) |
| **Statut** | ✅ **C'est celui que tu utilises activement** | 🔒 En réserve, ne pas lancer en parallèle du premier |
| **Import Postgres** | `scripts/load_to_postgres.py` sur le master actuel | `scripts/load_to_postgres.py` |

**Règle d'or : ne jamais mélanger les deux dans le même cycle de collecte.** Le workflow à suivre est détaillé en bas de ce document.

---

## 📁 Dossier `collect/`

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| **`collect_all_sources.py`** | ✅ **Oui — c'est LE point d'entrée** | Orchestre toute la collecte : NCBI, UniProt, KEGG, PlantTFDB, PubMed, Expression Atlas et GEO, en parallèle sur plusieurs espèces. Écrit un fichier par espèce puis reconstruit `master_plant_db.json` avec `rebuild_master_safe.py`. |
| `collect_all_plants.ps1` / `.sh` | Optionnel | Simple wrapper qui appelle `collect_all_sources.py` avec des options prédéfinies. Pratique mais pas obligatoire — tu peux appeler `collect_all_sources.py` directement avec tes propres arguments (ce qu'on fait depuis le début). |
| `run_collect_with_log.py` / `.ps1` | Optionnel | Wrapper qui ajoute une journalisation automatique. Alternative à `Start-Transcript`. |
| `collect_uniprot.py` | ❌ Jamais seul | Utilisé automatiquement par `collect_all_sources.py` pour la source `uniprot`. Contient `fetch_uniprot(species, retmax)`. |
| `collect_kegg.py` | ❌ Jamais seul | Idem pour `kegg`. Contient `fetch_kegg(species, retmax)`. Corrigé récemment (pathways/KO conservés même sans séquence). |
| `collect_planttfdb.py` | ❌ Jamais seul | Idem pour `planttfdb`. Retourne des métadonnées de facteurs de transcription, **sans séquence** dans la plupart des cas (limite connue de PlantTFDB, pas un bug). |
| `collect_pubmed.py` | ❌ Jamais seul | Idem pour `pubmed`. Retourne un paquet de publications par espèce (pas des gènes individuels). Corrigé récemment (ne prétend plus avoir une séquence ADN). |
| `collect_atlas.py` | ❌ Jamais seul | Pont vers `scripts/collect_expression_atlas.py`. Retourne un pseudo-enregistrement par espèce regroupant les expériences Expression Atlas trouvées. |
| `collect_geon.py` | ❌ Jamais seul | Pont vers `scripts/collect_geo.py`. Même principe, pour les datasets NCBI GEO. |
| `collect_ensembl.py` | ❌ Pas utilisé par la collecte bulk | Lookup unitaire fonctionnel par symbole ou identifiant. Il ne fournit pas l'interface bulk attendue par `collect_all_sources.py`. Pour `--sources ensembl` dans le pipeline bulk, `collect_ensembl_stub.py` est appelé et retourne zéro résultat. |
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
| **`rebuild_master_safe.py`** | ✅ Ponctuellement | Point unique de vérité pour reconstruire le master. Lit un fichier espèce à la fois, dédoublonne par `(organism, gene_id)`, écrit un temporaire, vérifie le JSON en streaming, puis remplace le master atomiquement. Crée aussi `master_plant_db.integrity.json`. |
| `rebuild_master.py` / `rebuild_master_stream.py` | ❌ À éviter | Anciennes versions de reconstruction qui ne doivent plus être utilisées pour le master actuel. |

### Validation

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| **`verify_gene_records.py`** | ✅ Oui, avec réserve | Lit depuis Postgres (`--db`) ou un JSON (`--json-file`), classe les enregistrements et génère un rapport. Le master actuel utilise `sequence.dna/rna/protein` imbriqué : ce script attend surtout un format plat et ne suffit donc pas à certifier les types de séquences du master. |
| `validate_and_add_gene.py` | ❌ Utilisé en interne | Valide un enregistrement avant insertion dans `genes_database.json`. |
| `validate_bioinformatics.py` | ✅ Ponctuel | Teste les moteurs d'alignement/mutation/phylogénie eux-mêmes (pas tes données) — pertinent quand tu travailles sur `alignment_engine.py` etc., pas sur la collecte. |
| `test_postgres_connection.py` | ✅ Dès qu'il y a un doute de connexion | Diagnostique DNS/TCP/auth vers Postgres, étape par étape. Aurait pu nous faire gagner du temps sur le bug de port 5432/5433 ! À utiliser en premier réflexe la prochaine fois qu'une connexion échoue. |

### Base de données

| Fichier | À lancer directement ? | Rôle |
|---|---|---|
| **`load_to_postgres.py`** | ✅ **Oui — import compatible avec le schéma actuel** | Charge le master JSON, extrait les séquences imbriquées et les upsert dans PostgreSQL. Peut créer les tables avec `--create-tables`. |
| `import_species_to_postgres.py` | ⚠️ À corriger avant usage | Attend encore `sequence` comme chaîne plate et n'est donc pas compatible tel quel avec `sequence.dna/rna/protein`. |
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
python collect/collect_all_sources.py --plant-file toutes_especes.txt --workers 4 --retmax 300
```
**Ce que ça fait** : lance 4 collectes en parallèle, interroge NCBI/UniProt/KEGG/PlantTFDB/PubMed/Atlas/GEO pour chaque espèce, écrit un fichier JSON par espèce dans `data/clean/species/`, puis reconstruit le master avec `rebuild_master_safe.py` à la fin. Les sources Atlas, GEO et PubMed produisent principalement des métadonnées au niveau espèce ou dataset, pas des séquences par gène.

Le pipeline utilise `--skip-existing` par défaut : pour recollecter une espèce déjà présente, ajouter `--force`. Les sources par défaut sont `ncbi,uniprot,kegg,planttfdb,pubmed,atlas,geon`; Ensembl et PLAZA doivent être demandés explicitement et ont des contraintes particulières.

`--create-tables` crée seulement les tables PostgreSQL. Pour importer également les données, utiliser `--load-db --create-tables` après avoir vérifié la connexion.

### 2. Nettoyer
**Rien à faire** — déjà inclus dans l'étape 1.

### 3. Vérifier avant l'import massif
```powershell
.\.venv\Scripts\python.exe .\scripts\verify_gene_records.py --json-file data\clean\master_plant_db.json --stats-out data\clean\stats_avant_import.json
```
**Limite importante** : le master actuel stocke les séquences dans `sequence.dna`, `sequence.rna` et `sequence.protein`. `verify_gene_records.py` attend surtout une séquence plate et son classement des types ne suffit donc pas à valider le master. Contrôler aussi les compteurs d'espèces, `origin`, les séquences imbriquées et les doublons `(organism, gene_id)`.

### 4. Envoyer vers Postgres
```powershell
.\.venv\Scripts\python.exe .\scripts\load_to_postgres.py --json-file data\clean\master_plant_db.json --create-tables
```
**Ce que ça fait** : charge le master, extrait les séquences depuis la structure imbriquée et upsert les enregistrements dans PostgreSQL. `import_species_to_postgres.py` n'est pas compatible tel quel avec le champ `sequence` imbriqué et ne doit pas être utilisé avant adaptation.

### 5. Vérifier le résultat final
```powershell
.\.venv\Scripts\python.exe .\scripts\verify_gene_records.py --db --csv-out data\clean\rapport_final.csv --stats-out data\clean\stats_final.json
```
**Ce que ça fait** : relit directement depuis Postgres (source de vérité), confirme que ce qui est en base correspond à ce qui a été collecté.

### Reconstruction manuelle du master
Si la collecte est terminée mais que le master doit être régénéré :
```powershell
.\.venv\Scripts\python.exe .\rebuild_master_safe.py --species-dir data\clean\species --out data\clean\master_plant_db.json
```
Le script prend tous les fichiers `*_all_sources.json` présents dans le dossier, pas uniquement ceux du dernier run. Il dédoublonne par `(organism, gene_id)`, vérifie le fichier temporaire en streaming, remplace le master atomiquement et écrit `master_plant_db.integrity.json`.

---

## 🚫 Ce que tu ne dois jamais faire

- Lancer `run_pipeline.py` **en plus** de `collect_all_sources.py` pour la même collecte — ce sont deux systèmes séparés qui ne se complètent pas, ça duplique le travail.
- Lancer `verify_gene_records.py` sans `--json-file` ni `--db` en pensant lire tes données actuelles — il retombe par défaut sur `data/clean/plant_data_clean.json`, un fichier de test périmé s'il existe encore.
- Modifier `collect_uniprot.py`, `collect_kegg.py`, etc. dans `collect/` en pensant que ça affecte l'usage "ponctuel" — ce sont des fichiers différents de leurs équivalents dans `scripts/`, malgré des noms proches.
