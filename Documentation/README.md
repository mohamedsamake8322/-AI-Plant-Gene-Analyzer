🧬 AI-Powered Bioinformatics Analysis Platform — Évolution du Projet

Le projet a évolué d’une simple application d’analyse ADN vers une plateforme bioinformatique computationnelle modulaire.

Objectifs principaux:
- Alignement de séquences (pairwise & MSA)
- Calcul de distances évolutives
- Construction d’arbres phylogénétiques
- Analyse protéique
- Annotation biologique
- Interprétation intelligente des résultats

Architecture proposée:

Core Engines/
- `alignment_engine.py` — Needleman–Wunsch, Smith–Waterman, MSA
- `distance_engine.py` — Hamming, Jukes-Cantor, Kimura, PAM, distance matrix
- `phylogeny_engine.py` — UPGMA, Neighbor-Joining, Newick export

Legacy Modules/
- Modules existants migrés progressivement: `bioinformatics.py`, `similarityengine.py`, `visualization.py`, `export_utils.py`, `app.py`.

Étapes immédiates recommandées:
1. Factoriser le code des modules existants vers `Core Engines`.
2. Ajouter tests unitaires pour chaque moteur.
3. Améliorer les visualisations interactives (Plotly dendrogrammes, alignement coloré).
4. Préparer une API pour consommation programmatique (FastAPI / Flask) si besoin.

Futurs ajouts possibles:
- Engine de protéomique (PAM, BLOSUM, pI, poids moléculaire)
- Annotation via bases publiques (GenBank, UniProt)
- Export Newick et formats standards
- Intégration CI (GitHub Actions) et packaging
