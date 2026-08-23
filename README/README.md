# ?? Plant Gene Analyzer

Application bioinformatique locale et web pour analyser des s�quences ADN/ARN/prot�ines v�g�tales, comparer les r�sultats � une base de r�f�rence, produire des interpr�tations biologiques et exporter des rapports professionnels.

---

## Vérification des sections d'analyse

L'application Streamlit se lance avec :

```bash
streamlit run app.py
```

Après avoir collé une séquence et cliqué sur **Analyze Sequence**, les onglets
suivants sont disponibles :

| Section | Fonction actuelle |
|---|---|
| **Statistics** | Nettoie et valide la séquence, calcule la longueur, GC%, AT%, ratio GC/AT, fréquences A/T/G/C/N, bases ambiguës, GC3 et les ORF détectés sur les six cadres. Affiche aussi le profil GC par fenêtre et les motifs connus. |
| **Similarity** | Compare la séquence à la base locale avec un alignement global Needleman-Wunsch, affiche les meilleurs scores, la classification, la couverture, les gaps, la carte d'alignement et une jauge de confiance. Le mode Deep Search élargit les candidats. |
| **Mutations** | Compare la requête au meilleur match disponible et liste substitutions, insertions, délétions, identités et positions. Sans meilleur match, l'onglet indique qu'aucun rapport n'est disponible. |
| **Translation** | Traduit une séquence ADN dans le cadre choisi (+1, +2 ou +3), affiche les six cadres, le statut de l'ORF et les séquences complémentaire et reverse-complement. Pour une protéine, la traduction est signalée comme non applicable. |
| **AI Interpretation** | Produit une interprétation déterministe basée sur des règles : profil de séquence, GC, similarité, mutations, potentiel fonctionnel, stress et recommandations agricoles. Il ne s'agit pas d'un appel à une API d'IA externe. |
| **Raw Sequence** | Affiche la séquence nettoyée, sa longueur et ses métriques principales, avec téléchargement FASTA et rapport texte. Une action d'annotation locale est également proposée. |
| **Alignments** | Outil manuel pour deux séquences ou davantage : MSA guidé par étoile, Needleman-Wunsch global et Smith-Waterman local, avec visualisation et scores. |
| **Distance Matrix** | Outil manuel multi-séquences. Aligne les séquences puis calcule une matrice avec la méthode choisie : Hamming, Jukes-Cantor, Kimura ou PAM, exportable en CSV. |
| **Phylogeny** | Outil manuel multi-séquences. Réutilise une matrice Kimura et construit un arbre UPGMA ou Neighbor Joining, affiché en dendrogramme et exportable au format Newick. |
| **Protein Analysis** | Outil manuel pour une protéine : validation, longueur, composition, masse moléculaire, point isoélectrique, hydrophobicité et graphique des acides aminés. |

### Test de référence vérifié

La séquence de contrôle de 300 bp fournie pour la vérification produit avec le
pipeline local :

- type détecté : ADN ;
- GC : 67,67 % ; AT : 32,33 % ;
- 1 ORF détecté et traduction de 22 acides aminés ;
- 2 motifs détectés ;
- aucun meilleur match dans la base locale utilisée au test.

Le dernier point est important : dans ce cas, **Similarity** affiche l'absence
de résultats et **Mutations** ne peut pas comparer la séquence à une référence.
Les onglets **Alignments**, **Distance Matrix** et **Phylogeny** demandent leurs
propres séquences de comparaison dans leurs champs interactifs.

### Limites d'interprétation

Les résultats statistiques, les alignements et les distances sont des calculs
reproductibles. Les phrases de l'onglet **AI Interpretation** sont des règles
heuristiques et ne constituent pas une annotation expérimentale ni une preuve
de fonction ou de résistance agronomique. Les séquences courtes, divergentes ou
sans référence proche doivent donc être interprétées avec prudence.

## ?? D�marrage rapide

### Pr�requis
- Python 3.8+
- pip
- Internet pour installer les d�pendances

### Installation en 2 minutes
```bash
cd c:\Downloads\IA
pip install -r requirements.txt
```

### Lancer l'application
```bash
streamlit run app.py
```

L'application est disponible sur http://localhost:8501.

### Test rapide
1. Ouvrir l'interface Streamlit
2. Charger la d�mo "DREB-like"
3. Cliquer sur "?? Analyze Sequence"
4. Explorer les onglets Statistics, Similarity, AI Interpretation

---

## ? Fonctionnalit�s principales

### Analyse bioinformatique
- Calcul du GC%, AT%, distribution nucl�otidique
- Validation de s�quences et d�tection d'erreurs
- Traduction ADN ? prot�ines sur plusieurs cadres
- D�tection de mutations et comparaison de s�quences
- Recherche de motifs r�gulateurs (TATA-box, CAAT-box, etc.)

### Comparaison � une base de donn�es
- Alignement local et global
- Scores de similarité basés sur pg_trgm (trigram index)
- Meilleurs matchs avec métadonnées (55,979 gènes vegétaux)
- Classification automatique des similarités
- **Visualisations enrichies:**
  - Heatmap de couverture d'alignement
  - Métriques détaillées (gaps%, coverage%, matches/mismatches)
  - Tableau comparatif des 3 meilleurs matches
  - Indicateur de confiance (0-100% composite)
  - Contexte biologique complet

### Interprétation IA
- Analyse règle-based des caractéristiques biologiques
- Prédiction de résistance au stress (sécheresse, chaleur, maladies)
- Suggestions agricoles et fonctionnelles
- **Phase 2 en cours:** AI INTERPRETABLE (classification des similarités, explications NLP)

### Visualisation et export
- Graphiques Plotly interactifs
- Alignements ASCII de type BLAST-like
- Export JSON, CSV et HTML
- Interface Streamlit moderne et sombre

---

## ?? Architecture du projet

```text
Utilisateur
   ?
Streamlit UI (app.py)
   ?
Analyse bioinformatique (bioinformatics.py)
   ?
Comparaison base de donn�es (similarityengine.py)
   ?
Interpr�tation IA (aiinterpreter.py)
   ?
Export / visualisation / logs
```

Le flux principal est le suivant :
1. L'utilisateur fournit une s�quence ADN/ARN/prot�ine
2. L'application valide et nettoie la s�quence
3. Les analyses statistiques et comparatives sont calcul�es
4. Un rapport est g�n�r� avec explications biologiques
5. Les r�sultats peuvent �tre export�s ou sauvegard�s

---

## ?? Modules principaux

### app.py
Interface utilisateur Streamlit. G�re l'entr�e de s�quence, les param�tres, l'ex�cution du pipeline d'analyse, l'affichage des r�sultats et les exports.

### bioinformatics.py
Moteur central de traitement biologique. Contient les fonctions de nettoyage, validation, calcul GC, distribution, traduction, mutations et motifs.

### similarityengine.py
Moteur de comparaison aux g�nes de r�f�rence. Fournit les scores de similarit�, les alignements et la s�lection des meilleurs matchs.

### aiinterpreter.py
Interpr�tation intelligente des r�sultats via logique de r�gles. Produit des commentaires biologiques, des pr�visions fonctionnelles et des recommandations.

### visualization.py
G�n�ration de graphiques interactifs avec Plotly : composition nucl�otidique, GC profile, similarit�, mutations.

### config.py
Configuration centralis�e : chemins, seuils, logging, styles, param�tres d'analyse.

### export_utils.py
Utilitaires d'export multi-format : JSON, CSV, HTML.

### setup_and_run.py
Assistant d'installation et de v�rification de l'environnement.

---
## 📊 Visualisations d'analyse de similarité avancée

Le module `visualization.py` contient 5 nouvelles fonctions dédiées à l'analyse enrichie des résultats de similarité :

### 1. `build_match_context_card(match: dict) → dict`
Extrait le contexte biologique complet d'un gène matché.
- **Données:** gene_name, trait, organism, description, accession, source
- **Utilité:** Comprendre rapidement la biologie du gène matché
- **Affichage:** Markdown formaté dans l'expander de chaque match

### 2. `build_similarity_metrics_table(match: dict, query_len: int) → dict`
Calcule les statistiques d'alignement détaillées.
- **Calculs:**
  - Coverage% = (aligned_length - gaps) / query_len * 100
  - Gap% = gaps / aligned_length * 100
  - Matches/Mismatches = analyse base à base
  - Identity% = similarité du score
- **Utilité:** Valider la qualité de l'alignement (haute couverture + peu de gaps = confiance)
- **Affichage:** Tableau de métriques dans chaque expander

### 3. `plot_alignment_coverage_heatmap(match: dict, query_len: int, window: int = 50) → go.Figure`
Heatmap montrant l'identité en fenêtres glissantes.
- **X-axis:** Position dans l'alignement
- **Y-axis:** Identité% calculée par fenêtre (50bp par défaut)
- **Couleurs:** CORAL (0%) → AMBER (50%) → TEAL (100%)
- **Utilité:** Identifier régions conservées vs variables
- **Affichage:** Graphique Plotly interactif après alignment map

### 4. `build_top3_comparison_table(similarity_results: list[dict], query_len: int) → dict`
Tableau comparatif des 3 meilleurs matches.
- **Colonnes:** Rank | Gene | Similarity | Trait | Organism | Coverage% | Gaps%
- **Utilité:** Détecter rapidement si tous les matches s'accordent (forte confiance) ou divergent (explorer)
- **Affichage:** Table markdown avant les expanders individuels

### 5. `plot_confidence_gauge(metrics: dict) → go.Figure`
Indicateur de confiance composite 0-100%.
- **Formule:** Coverage(40%) + (100-Gaps%)(35%) + Identity(25%)
- **Zones couleur:**
  - TEAL 75-100% = Haute confiance
  - AMBER 50-75% = Confiance moyenne
  - CORAL 0-50% = Basse confiance
- **Annotations:** Détails coverage/gaps/identity en bas du gauge
- **Affichage:** Gauge interactif dans chaque expander

---
## ?? Tests

Le projet contient une suite de tests unitaires pour le moteur bioinformatique.

```bash
pytest test_bioinformatics.py -v
```

Pour une couverture :
```bash
pytest test_bioinformatics.py --cov=bioinformatics
```

---

## ?? Installation d�taill�e

### 1. Ouvrir un terminal
Windows : PowerShell ou CMD

```bash
cd c:\Downloads\IA
```

### 2. Cr�er un environnement virtuel (recommand�)
```bash
python -m venv env
env\Scripts\activate
```

Linux/macOS :
```bash
python3 -m venv env
source env/bin/activate
```

### 3. Installer les d�pendances
```bash
pip install -r requirements.txt
```

### 4. V�rifier l'installation
```bash
python -c "import streamlit, plotly, bioinformatics; print('? OK')"
```

### 5. Lancer l'application
```bash
streamlit run app.py
```

Si le port 8501 est d�j� utilis�, utiliser :
```bash
streamlit run app.py --server.port 8502
```

---

## ?? Utilisation

### Workflow simple
1. Coller une s�quence ADN ou charger un fichier FASTA
2. Choisir les param�tres de l'analyse
3. Cliquer sur "?? Analyze Sequence"
4. Consulter les r�sultats dans les onglets
5. Exporter le rapport si besoin

### Donn�es de d�monstration
Le projet inclut des fichiers de test comme :
- Data/sample_gene.fasta
- Data/sample_protein.fasta

Vous pouvez les charger directement via l'interface pour tester les fonctionnalit�s.

---
## 🗄️ Base de données PostgreSQL

Le projet supporte une base de données PostgreSQL/Neon avec 55,979 gènes végétaux.

### Caractéristiques
- **Taille:** 55,979 gènes avec métadonnées complètes
- **Technologie:** pg_trgm (trigram GIN index) pour recherche rapide
- **Connexion:** psycopg avec pool de connexions (fallback gracieux si psycopg_pool unavailable)
- **Recherche:** Similarité trigram sur colonne `genes.sequence`
- **Performance:** Indexée et optimisée pour requêtes sub-second

### Configuration
Configuration via variables d'environnement (.env):
```
DATABASE_URL=postgresql://user:pass@host:port/dbname
```

### Dégradation gracieuse
Si PostgreSQL est indisponible:
- App affiche `st.error()` explicite
- Pas de fallback silencieux vers base locale 13-gènes
- Force utilisateur à corriger la connexion

### Scripts utilitaires
- `scripts/postgres_utils.py` - Connection pooling, gestion requêtes
- `scripts/load_to_postgres.py` - Import données
- `scripts/test_postgres_connection.py` - Vérification connectivité

---
## ?? Collecte de donn�es v�g�tales

Le projet inclut �galement un pipeline de collecte de donn�es via le dossier `collect` et le dossier `scripts`.

### Collecte multi-sources
Le script principal est :
```bash
python collect/collect_all_sources.py --all-plants --workers 4 --retmax 300
```

### Collecte NCBI multi-type
```bash
python scripts/collect_multi_type.py --plant "Oryza sativa" --retmax 300
```

### Collecte de m�tadonn�es et pipeline complet
Des scripts existent pour :
- GEO / Expression Atlas / Ensembl / NCBI
- UniProt, KEGG, PlantTFDB, PubMed
- Nettoyage des donn�es
- Import dans PostgreSQL
- Reconstruction de bases fusionn�es

Les pipelines sont con�us pour produire des fichiers JSON par esp�ce et un master JSON global.

---

## ?? Sch�ma professionnel (version avanc�e)

Le projet inclut �galement un sch�ma de donn�es professionnel et un chargeur compatible.

### Fichiers associ�s
- professional_schema.py
- professional_loader.py
- scripts/transform_schema.py
- scripts/clean_data.py

### Int�r�t
- Versioning du sch�ma
- M�tadonn�es structur�es (taxonomie, qualit�, analytics)
- Compatibilit� avec l'ancien format JSON
- Pr�paration � une utilisation plus robuste et �volutive

---

## ? Am�liorations majeures impl�ment�es

### Phase 1 : Correctifs critiques (Production stability)
✅ **similarityengine.py** - Correction bug critique
- `find_similar_genes()` appelait des méthodes non-existantes (`pg._sequence_minimizers()`, `pg.KMER_WINDOW`)
- Maintenant utilise `pg.find_candidate_genes_by_kmer()` avec recherche pg_trgm fonctionnelle
- Source indiquée : "trigram_index" au lieu de "kmer_index"

✅ **config.py** - Ajout constantes manquantes
- `CHART_TITLE_COLOR = "#00d9a3"` (causait AttributeError)
- `CHART_LINE_COLOR = "#cccccc"` (causait AttributeError)

✅ **pipeline.py** - Correction erreurs de syntaxe
- Indentation et f-strings malformées aux lignes 120-148
- Validation complète d'import

✅ **postgres_utils.py** - Dégradation gracieuse
- Gestion conditionnelle de `psycopg_pool` (try/except import)
- Classe fallback `_DirectConnectionPool` pour connexions directes
- Fonctionne avec ou sans pool disponible

✅ **requirements.txt** - Dépendance ajoutée
- `psycopg_pool>=1.0.0` (manquait, causait ModuleNotFoundError)

✅ **app.py** - Protection production
- Suppression du fallback silencieux vers 13-gènes JSON
- Affiche `st.error()` si Postgres indisponible
- Garantit utilisation de la vraie base de données (55,979 gènes)

### Phase 2 : Visualisations enrichies (Enhanced similarity analysis)
✅ **visualization.py** - 5 nouvelles fonctions d'analyse :

1. **`build_match_context_card()`** - Contexte biologique du gène
   - Affiche: description, accession, organisme, trait, source
   - Rend claire la biologie du match

2. **`build_similarity_metrics_table()`** - Statistiques d'alignement détaillées
   - Colonnes: Aligned length, matches, mismatches, gaps%, coverage%, identity%
   - Répond: "Est-ce que cet alignement est réel ?"

3. **`plot_alignment_coverage_heatmap()`** - Heatmap de couverture
   - Affiche identité% en fenêtres glissantes (50bp par défaut)
   - Identifie régions conservées vs variables
   - Gradient couleur: CORAL (0%) → AMBER → TEAL (100%)

4. **`build_top3_comparison_table()`** - Tableau comparatif top 3
   - Colonnes: Rank, Gene, Similarity, Trait, Organism, Coverage%, Gaps%
   - Détection rapide de la cohérence biologique entre matches

5. **`plot_confidence_gauge()`** - Indicateur de confiance composite
   - Échelle 0-100 basée sur: Coverage(40%) + Gap penalty(35%) + Identity(25%)
   - Code couleur: TEAL(75-100), AMBER(50-75), CORAL(0-50)
   - Gauge interactif Plotly

✅ **app.py** - Intégration UI des visualisations
- Tab Similarity enrichi avec:
  - Table résumé top 3 avant expanders
  - Heatmap couverture après alignement
  - Gauge confiance + métriques détaillées
  - Contexte biologique du gène
- Chaque expander match affiche toutes les 5 visualisations

### Améliorations antérieures maintenues
- Cache de chargement de la base de données
- Gestion robuste des erreurs et logs
- Export JSON/CSV/HTML
- Configuration centralisée
- Visualisation ASCII des alignements
- Tests unitaires automatisés
- Pipeline de collecte et de nettoyage de données
- Interface Streamlit moderne et sombre

---

## ??? Roadmap

### ✅ Complété (Phase 1-2)
- ✅ Correction bug critique similarityengine.py (pg_trgm functional)
- ✅ Protection production (Postgres failsafe, pas de fallback JSON)
- ✅ Constantes config manquantes (CHART_TITLE_COLOR, CHART_LINE_COLOR)
- ✅ Dependencies robustes (psycopg_pool graceful fallback)
- ✅ 5 nouvelles visualisations (coverage heatmap, metrics, confidence gauge, top3 table, context card)
- ✅ Intégration UI complète (Similarity tab enrichi)

### En cours (Phase 3 - AI INTERPRETABLE)
- 🔄 Classification intelligente des similarités (homolog, domain match, full identity, etc.)
- 🔄 Explications NLP des résultats
- 🔄 Prédictions de stress et recommandations agricoles
- 🔄 Warnings sur résultats low-confidence

### Priorités court terme (Phase 4)
- Support des séquences protéiques avancé
- Support batch / multi-FASTA
- Exports professionnels enrichis (XLSX, rapport HTML)
- API REST pour intégration externe

### Améliorations avancées (Phase 5+)
- Intégration BLAST local ou API NCBI BLAST
- Alignements multiples et phylogénie
- Annotation génomique (GFF/BED)
- Analyse de variantes et d'indels avancée

---

## 📂 Structure du projet

```text
IA/
+-- app.py                          # Interface Streamlit (UI principale)
+-- bioinformatics.py               # Moteur bioinformatique central
+-- similarityengine.py             # Moteur comparaison genes (FIXED: pg_trgm)
+-- aiinterpreter.py                # Interprétation IA (Phase 2 en cours)
+-- visualization.py                # Graphiques Plotly (+ 5 nouvelles fonctions)
+-- config.py                        # Configuration centralisée (+ CHART_* constants)
+-- export_utils.py                 # Export JSON/CSV/HTML
+-- setup_and_run.py                # Assistance installation
+-- professional_schema.py           # Schéma avancé métadonnées
+-- professional_loader.py           # Chargeur schéma avancé
+-- requirements.txt                # Dépendances (+ psycopg_pool)
+-- README.md                        # Documentation (MISE À JOUR)
+-- genes_database.json              # Base 13 gènes (fallback, rarement utilisée)
+-- .env                             # Variables environnement (DATABASE_URL)
+--
+-- collect/                         # Pipeline collecte données
+-- scripts/                         # Utilitaires backend (postgres_utils FIXED)
+-- Tests/                           # Suite de tests
+-- Data/                            # Données de test/demo
+-- Documentation/                   # Documentation avancée
+-- logs/ results/                   # Générés automatiquement
```

---

 ?? D�pannage courant

### Python introuvable
V�rifier :
```bash
python --version
```

### Module introuvable
R�installer les d�pendances :
```bash
pip install -r requirements.txt
```

### Port 8501 d�j� utilis�
```bash
streamlit run app.py --server.port 8502
```

### Fichier genes_database.json absent
V�rifier que le fichier est bien pr�sent � la racine du projet.

### Tests qui �chouent
```bash
pytest test_bioinformatics.py -v -s
```

---

## ?? Documentation unique

Ce README est d�sormais le document principal et unique du projet. Il regroupe les informations de d�marrage, d�installation, d�utilisation, d�architecture, de collecte de donn�es et de roadmap.

---

## ?? R�sum�

Ce projet combine une interface Streamlit simple, un moteur bioinformatique robuste, une base de comparaison de g�nes et des outils de collecte de donn�es v�g�tales. Il est con�u pour �tre � la fois accessible pour un usage rapide et extensible pour des analyses plus avanc�es.

---

Happy analyzing! ??
