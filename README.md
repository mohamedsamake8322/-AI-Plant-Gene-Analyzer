# 🧬 AI-Powered Plant Gene Analyzer

Application bioinformatique locale et web pour analyser des séquences ADN/ARN/protéines végétales, comparer les résultats à une base de référence, produire des interprétations biologiques et exporter des rapports professionnels.

---

## 🚀 Démarrage rapide

### Prérequis
- Python 3.8+
- pip
- Internet pour installer les dépendances

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
2. Charger la démo "DREB-like"
3. Cliquer sur "🔬 Analyze Sequence"
4. Explorer les onglets Statistics, Similarity, AI Interpretation

---

## ✨ Fonctionnalités principales

### Analyse bioinformatique
- Calcul du GC%, AT%, distribution nucléotidique
- Validation de séquences et détection d'erreurs
- Traduction ADN → protéines sur plusieurs cadres
- Détection de mutations et comparaison de séquences
- Recherche de motifs régulateurs (TATA-box, CAAT-box, etc.)

### Comparaison à une base de données
- Alignement local et global
- Scores de similarité
- Meilleurs matchs avec métadonnées
- Classification automatique des similarités

### Interprétation IA
- Analyse règle-based des caractéristiques biologiques
- Prédiction de résistance au stress (sécheresse, chaleur, maladies)
- Suggestions agricoles et fonctionnelles

### Visualisation et export
- Graphiques Plotly interactifs
- Alignements ASCII de type BLAST-like
- Export JSON, CSV et HTML
- Interface Streamlit moderne et sombre

---

## 🧱 Architecture du projet

```text
Utilisateur
   ↓
Streamlit UI (app.py)
   ↓
Analyse bioinformatique (bioinformatics.py)
   ↓
Comparaison base de données (similarityengine.py)
   ↓
Interprétation IA (aiinterpreter.py)
   ↓
Export / visualisation / logs
```

Le flux principal est le suivant :
1. L'utilisateur fournit une séquence ADN/ARN/protéine
2. L'application valide et nettoie la séquence
3. Les analyses statistiques et comparatives sont calculées
4. Un rapport est généré avec explications biologiques
5. Les résultats peuvent être exportés ou sauvegardés

---

## 📁 Modules principaux

### app.py
Interface utilisateur Streamlit. Gère l'entrée de séquence, les paramètres, l'exécution du pipeline d'analyse, l'affichage des résultats et les exports.

### bioinformatics.py
Moteur central de traitement biologique. Contient les fonctions de nettoyage, validation, calcul GC, distribution, traduction, mutations et motifs.

### similarityengine.py
Moteur de comparaison aux gènes de référence. Fournit les scores de similarité, les alignements et la sélection des meilleurs matchs.

### aiinterpreter.py
Interprétation intelligente des résultats via logique de règles. Produit des commentaires biologiques, des prévisions fonctionnelles et des recommandations.

### visualization.py
Génération de graphiques interactifs avec Plotly : composition nucléotidique, GC profile, similarité, mutations.

### config.py
Configuration centralisée : chemins, seuils, logging, styles, paramètres d'analyse.

### export_utils.py
Utilitaires d'export multi-format : JSON, CSV, HTML.

### setup_and_run.py
Assistant d'installation et de vérification de l'environnement.

---

## 🧪 Tests

Le projet contient une suite de tests unitaires pour le moteur bioinformatique.

```bash
pytest test_bioinformatics.py -v
```

Pour une couverture :
```bash
pytest test_bioinformatics.py --cov=bioinformatics
```

---

## 📦 Installation détaillée

### 1. Ouvrir un terminal
Windows : PowerShell ou CMD

```bash
cd c:\Downloads\IA
```

### 2. Créer un environnement virtuel (recommandé)
```bash
python -m venv env
env\Scripts\activate
```

Linux/macOS :
```bash
python3 -m venv env
source env/bin/activate
```

### 3. Installer les dépendances
```bash
pip install -r requirements.txt
```

### 4. Vérifier l'installation
```bash
python -c "import streamlit, plotly, bioinformatics; print('✅ OK')"
```

### 5. Lancer l'application
```bash
streamlit run app.py
```

Si le port 8501 est déjà utilisé, utiliser :
```bash
streamlit run app.py --server.port 8502
```

---

## 🔬 Utilisation

### Workflow simple
1. Coller une séquence ADN ou charger un fichier FASTA
2. Choisir les paramètres de l'analyse
3. Cliquer sur "🔬 Analyze Sequence"
4. Consulter les résultats dans les onglets
5. Exporter le rapport si besoin

### Données de démonstration
Le projet inclut des fichiers de test comme :
- Data/sample_gene.fasta
- Data/sample_protein.fasta

Vous pouvez les charger directement via l'interface pour tester les fonctionnalités.

---

## 🌱 Collecte de données végétales

Le projet inclut également un pipeline de collecte de données via le dossier `collect` et le dossier `scripts`.

### Collecte multi-sources
Le script principal est :
```bash
python collect/collect_all_sources.py --all-plants --workers 4 --retmax 300
```

### Collecte NCBI multi-type
```bash
python scripts/collect_multi_type.py --plant "Oryza sativa" --retmax 300
```

### Collecte de métadonnées et pipeline complet
Des scripts existent pour :
- GEO / Expression Atlas / Ensembl / NCBI
- UniProt, KEGG, PlantTFDB, PubMed
- Nettoyage des données
- Import dans PostgreSQL
- Reconstruction de bases fusionnées

Les pipelines sont conçus pour produire des fichiers JSON par espèce et un master JSON global.

---

## 🧬 Schéma professionnel (version avancée)

Le projet inclut également un schéma de données professionnel et un chargeur compatible.

### Fichiers associés
- professional_schema.py
- professional_loader.py
- scripts/transform_schema.py
- scripts/clean_data.py

### Intérêt
- Versioning du schéma
- Métadonnées structurées (taxonomie, qualité, analytics)
- Compatibilité avec l'ancien format JSON
- Préparation à une utilisation plus robuste et évolutive

---

## ✅ Améliorations majeures implémentées

Le projet a été enrichi avec plusieurs améliorations professionnelles :
- Cache de chargement de la base de données
- Gestion robuste des erreurs et logs
- Export JSON/CSV/HTML
- Configuration centralisée
- Visualisation ASCII des alignements
- Tests unitaires automatisés
- Pipeline de collecte et de nettoyage de données

---

## 🛣️ Roadmap

### Priorités court terme
- Support des séquences protéiques avancé
- Support batch / multi-FASTA
- Base de données plus flexible
- Exports professionnels enrichis (XLSX, rapport HTML plus riche)

### Améliorations avancées
- Intégration BLAST local ou API NCBI BLAST
- Alignements multiples et phylogénie
- Annotation génomique (GFF/BED)
- Analyse de variantes et d'indels

---

## 🧰 Structure du projet

```text
IA/
├── app.py
├── bioinformatics.py
├── similarityengine.py
├── aiinterpreter.py
├── visualization.py
├── config.py
├── export_utils.py
├── setup_and_run.py
├── professional_schema.py
├── professional_loader.py
├── genes_database.json
├── requirements.txt
├── README.md
├── collect/
├── scripts/
├── tests/
└── logs/ / results/ (générés automatiquement)
```

---

## 🐛 Dépannage courant

### Python introuvable
Vérifier :
```bash
python --version
```

### Module introuvable
Réinstaller les dépendances :
```bash
pip install -r requirements.txt
```

### Port 8501 déjà utilisé
```bash
streamlit run app.py --server.port 8502
```

### Fichier genes_database.json absent
Vérifier que le fichier est bien présent à la racine du projet.

### Tests qui échouent
```bash
pytest test_bioinformatics.py -v -s
```

---

## 📚 Documentation unique

Ce README est désormais le document principal et unique du projet. Il regroupe les informations de démarrage, d’installation, d’utilisation, d’architecture, de collecte de données et de roadmap.

---

## 🎯 Résumé

Ce projet combine une interface Streamlit simple, un moteur bioinformatique robuste, une base de comparaison de gènes et des outils de collecte de données végétales. Il est conçu pour être à la fois accessible pour un usage rapide et extensible pour des analyses plus avancées.

---

Happy analyzing! 🧬
