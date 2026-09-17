# Guide utilisateur et technique de l'analyse de sequence

Ce dossier explique le parcours complet d'une sequence dans **AI Plant Gene Analyzer**.
Il s'adresse aux debutants, aux utilisateurs intermediaires et aux utilisateurs avances.

## Parcours general

```text
Sequence collee ou fichier FASTA
        |
        v
Nettoyage + validation + detection ADN/proteine
        |
        v
Pipeline d'analyse
        |
        +--> Statistics
        +--> Similarity
        +--> Mutations
        +--> Translation
        +--> AI Interpretation
        +--> Raw Sequence / exports
```

Le point d'entree de l'application est `app.py`, mais l'interface principale est implementee dans `views/home.py`.
L'orchestration scientifique est dans `pipeline.py`.

## Demarrage rapide

```powershell
cd C:\Downloads\IA
.\.venv\Scripts\python.exe -m streamlit run app.py
```

1. Ouvrir `http://localhost:8501`.
2. Choisir le type d'entree dans la barre laterale : `Auto detect`, `DNA` ou `Protein`.
3. Coller une sequence ou charger un fichier `.fasta`, `.fa` ou `.txt`.
4. Choisir le cadre de lecture si l'entree est de l'ADN : `+1`, `+2`, `+3`, `-1`, `-2`, `-3`.
5. Cliquer sur **Analyze Sequence**.
6. Lire les onglets dans l'ordre recommande : Statistics, Similarity, Mutations, Translation.

## Documents de ce dossier

- [01 - Entree et pipeline](01_entree_et_pipeline.md)
- [02 - Statistics](02_statistics.md)
- [03 - Similarity](03_similarity.md)
- [04 - Mutations](04_mutations.md)
- [05 - Translation](05_translation.md)
- [06 - Limites et depannage](06_limites_et_depannage.md)

## Regle importante d'interpretation

Une statistique calculee sur la sequence est une mesure directe. Une similarite est une comparaison avec une reference. Une mutation n'est interpretable comme telle que si la reference est suffisamment proche ou explicitement fournie par l'utilisateur.

Le pipeline refuse donc de produire un rapport de mutations lorsqu'un meilleur match est sous `MIN_MUTATION_REFERENCE_IDENTITY` (85 % par defaut). Cette protection evite de presenter une divergence evolutive normale comme une liste de mutations.

## Modules principaux

| Module | Responsabilite |
|---|---|
| `sequence_loader.py` | Lecture FASTA, metadonnees d'en-tete et detection ADN/proteine |
| `pipeline.py` | Orchestration de l'analyse |
| `bioinformatics.py` | Statistiques, traduction, ORF et mutations de base |
| `similarityengine.py` | Recherche et classement des references |
| `variant_analysis.py` | Consequences biologiques, codons, acides amines, BLOSUM62 et frameshifts |
| `alignment_engine.py` | Needleman-Wunsch, Smith-Waterman et matrices de scores |
| `visualization.py` | Graphiques Plotly |
| `export_utils.py` | Exports JSON, CSV, VCF, FASTA, GFF3, HTML et XLSX |
| `i18n.py` | Traductions francais, anglais et turc |

## Niveau de confiance

Les calculs d'alignement, de GC, de traduction et de variants sont deterministes. Les textes de l'onglet **AI Interpretation** sont produits par des regles locales et ne remplacent pas une annotation experimentale, une base clinique ou une validation biologique.
