# Roadmap — AI-Powered Plant Gene Analyzer

## Vision

Faire de cette application un outil local professionnel pour l'analyse de séquences végétales, tout en conservant une interface simple pour les débutants.

---

## Priorités court terme

1. **Support des entrées protéiques**
   - Reconnaissance automatique des séquences AA
   - Analyse de composition d'acides aminés
   - Comparaison entre protéines et ADN par traduction

2. **Support multi-FASTA / batch**
   - Chargement de fichiers FASTA multi-séquences
   - Sélection de la séquence à analyser
   - Gestion de lots et rapports groupés

3. **Base de données plus flexible**
   - Prise en charge de JSON et FASTA pour la base de référence
   - Métadonnées enrichies (`sequence_type`, `trait`, `organism`, `accession`)

4. **Exports professionnels**
   - Ajout de l'export XLSX
   - Rapport HTML plus riche
   - Package d'export automatique JSON/CSV/HTML/XLSX

---

## Améliorations UX

- Mode pas à pas pour débutants
- Aide intégrée dans l'interface Streamlit
- Explications sur les metrics biologiques
- Pré-chargement de séquences de démonstration par trait

---

## Améliorations avancées

- Intégrer BLAST local ou appel NCBI BLAST API
- Alignement multiple et construction d'arbres phylogénétiques
- Annotation génomique (GFF/BED)
- Analyse de variantes et d'indels
- Intégration de données contextualisées (climat, culture, épigénétique)

---

## Tests et qualité

- Ajouter tests unitaires pour tous les modules
- Ajouter tests d'intégration pour app.py / pipeline
- Mettre en place un pipeline CI avec `pytest` et `pytest-cov`
- Documenter les règles de développement et le style de code
