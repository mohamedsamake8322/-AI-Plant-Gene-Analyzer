# AI-Powered Plant Gene Analyzer — Overview

## Objectif du projet

Cette application analyse des séquences ADN végétales et génère des résultats bioinformatiques complets, des comparaisons de similarité, des traductions protéiques et des interprétations biologiques intelligentes.

Elle est destinée à la fois aux utilisateurs débutants et aux professionnels qui veulent un outil local, sécurisé et simple pour tester des gènes et fragments de séquence.

---

## Ce que l'app peut analyser aujourd'hui

### 1. Séquences nucléotidiques (ADN)
- Nettoyage et validation des séquences ADN brutes ou FASTA
- Calcul du contenu GC et AT
- Distribution des nucléotides A/T/G/C/N
- Statistiques de séquence (longueur, ratio GC/AT, codons de départ/fin)
- Recherche de motifs réglementaires connus
- Traduction en protéine dans les 3 cadres de lecture (+1, +2, +3)

### 2. Traduction protéique
- Génération de séquences protéiques à partir d'un ADN codant
- Indication de l'existence d'une ORF complète
- Affichage du résultat pour chaque cadre de lecture

### 3. Analyse de similarité contre base locale
- Comparaison du query contre `genes_database.json`
- Similarité globale directe et similitude par fenêtre glissante
- Classement des meilleurs matches
- Affichage de la meilleure correspondance
- Visualisation de l'alignement en ASCII

### 4. Interprétation basée sur des règles
- Évaluation du profil de séquence (court, moyen, long)
- Analyse du contenu GC pour estimer tolérances au stress
- Interprétation de la similarité et du match de gènes
- Rapport sur les mutations et leur gravité
- Recommandations agricoles générales

### 5. Export
- Generation de rapports JSON
- Downloads CSV
- Downloads HTML

---

## Limites actuelles

### Séquences prises en charge
- Séquences ADN uniquement
- Pas de support natif pour séquences protéiques en entrée
- Pas d'analyse directe de génomes entiers ou d'assemblages complets
- Base de références limitée à `genes_database.json` (quelques gènes de démonstration)

### Capacités bioinformatiques
- Pas de BLAST ou alignement multiple sophistiqué
- Pas d'analyse de domaine protéique ou annotation fonctionnelle avancée
- Pas de prédiction de structure ou de phylogénie
- Pas de traitement de grands fichiers FASTA multi-séquences

---

## Fichiers principaux et rôles

### `app.py`
- Interface Streamlit
- Chargement du database
- Entrée utilisateur + UI
- Orchestration de l'analyse
- Présentation des résultats
- Export des rapports

### `bioinformatics.py`
- Nettoyage et validation de séquence
- Statistiques de base et distribution modernes
- Calcul GC/AT
- Traduction codon/protéine
- Détection de mutations
- Recherche de motifs et fonctions utilitaires

### `similarityengine.py`
- Chargement de la base de gènes
- Similarité paire-à-paire
- Fenêtre glissante pour alignements locaux
- Classement des meilleures correspondances
- Génération d'alignement ASCII

### `aiinterpreter.py`
- Interprétation basée sur règles
- Génération de rapports textuels
- Évaluation de la similarité, GC, mutations
- Recommandations agricoles

### `visualization.py`
- Graphiques Plotly pour l'interface
- Pies, barres, jauges, profils GC, scores de similarité, mutations

### `export_utils.py`
- Export JSON/CSV/HTML
- Sauvegarde dans `results/`
- Génération de rapport structurés

### `config.py`
- Paramètres constants et seuils
- Chemins des fichiers
- Configuration du logging
- Paramètres UI

### `genes_database.json`
- Base de gènes de référence
- Utilisée pour les matches de similarité

### `requirements.txt`
- Bibliothèques Python nécessaires

### `test_bioinformatics.py`
- Tests unitaires pour le moteur bioinformatique

---

## Grande lignes d'architecture

1. **Entrée utilisateur**
   - Texte / fichier FASTA / démo
2. **Nettoyage + validation**
   - `bio.clean_sequence`
   - `bio.validate_sequence`
3. **Analyse bioinformatique principale**
   - Statistiques, GC, traduction, motifs
4. **Comparaison de similarité**
   - `sim.compare_with_database`
   - `sim.pairwise_similarity`
   - `sim.sliding_window_similarity`
5. **Interprétation AI**
   - `ai_interp.interpret`
6. **Affichage & export**
   - Graphiques Streamlit
   - JSON/CSV/HTML

---

## Améliorations recommandées pour dépasser les limites

### Niveau débutant
- Ajouter un mode "pas à pas" avec explications simples
- Ajouter des guides contextuels dans l'UI pour chaque onglet
- Fournir une section "Exemple de séquence" pour chaque gène
- Ajouter un guide de lecture des résultats (`README.md` plus pédagogique)

### Niveau professionnel
- Ajouter une prise en charge des entrées protéiques (FASTA AA)
- Ajouter une option "Analyse génomique" pour séquences très longues
- Intégrer un module BLAST local / appel NCBI BLAST API
- Ajouter alignement multiple et clustering phylogénétique
- Permettre l'import de FASTA multi-séquences
- Ajouter gestion de métadonnées d'échantillon et export XLSX
- Ajouter des rapports d'annotation génomique (+ GFF/BED)

---

## Priorités d'amélioration

1. **Support des séquences protéiques**
2. **Support FASTA multi-séquences**
3. **Mise en place d'un backend d'annotation (BLAST / HMMER)**
4. **Amélioration du database format**
   - `genes_database.json` → `genes_database.fasta` + métadonnées
5. **Ajout d'une page d'aide/guide intégré**
6. **Ajout de tests sur `similarityengine.py` et `aiinterpreter.py`**

---

## Conclusion

Ton application est déjà solide pour l'analyse d'ADN court et la comparaison de gènes. Pour dépasser les limites, il faut maintenant ajouter des capacités protéiques, génomiques et d'annotation, tout en gardant l'expérience utilisateur accessible.
