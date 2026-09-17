# 1. Entree et pipeline

## 1.1 Types d'entree

L'application accepte :

- une sequence brute collee dans la zone de texte ;
- un fichier FASTA avec une ou plusieurs sequences ;
- un fichier `.fa` ou `.txt` contenant une sequence ;
- une sequence d'ADN ou une sequence proteique.

Le choix `Auto detect` utilise `sequence_loader.detect_sequence_type()`.
Une sequence contenant uniquement des caracteres ADN est traitee comme ADN. Une sequence contenant des caracteres d'acides amines est traitee comme proteine.

## 1.2 Metadonnees FASTA

Un en-tete comme celui-ci peut fournir des metadonnees :

```text
>gene_01 | organism=Chenopodium quinoa | contig=chr3 | GC=62%
ATG...
```

`sequence_loader.parse_header_metadata()` transforme les champs en dictionnaire. Ces metadonnees peuvent etre utilisees pour l'organisme, le GC, la longueur et le contig.

Une sequence collee sans en-tete possede `metadata={}`. Elle n'a donc pas de coordonnee genomique reelle. Les exports VCF utilisent alors un contig local `sequence` et le signalent comme coordonnee non genomique.

## 1.3 Nettoyage et validation

`pipeline.analyze_sequence_record()` appelle `bioinformatics.clean_sequence()` puis `bioinformatics.validate_sequence()`.

Le nettoyage :

- met les lettres en majuscules ;
- retire les espaces et retours a la ligne ;
- retire l'en-tete FASTA deja interprete ;
- conserve les caracteres ADN ambigus autorises ou les residus proteiques valides.

Les longueurs minimales par defaut sont :

- ADN : 10 bases pour l'analyse applicative ;
- proteine : 5 acides amines.

La limite dure de longueur generale est `MAX_SEQUENCE_LENGTH`.
La limite d'alignement est plus stricte : `MAX_ALIGNMENT_SEQUENCE_LENGTH`.

## 1.4 Parametres de la barre laterale

- `Top database matches` : nombre de meilleurs matchs affiches, de 1 a 8.
- `Deep search` : elargit la recherche de candidats, mais peut etre plus lente.
- `Sliding window (GC profile)` : taille de fenetre GC, de 5 a 60 bases.
- `Reading frame` : `+1`, `+2`, `+3`, `-1`, `-2`, `-3`.
- `Input type` : detection automatique, ADN ou proteine.
- `Explicit reference` : sequence de reference facultative pour l'analyse des mutations.

## 1.5 Execution

Pour chaque sequence, `views/home.py` appelle le wrapper cache `_cached_analyze()`, qui appelle ensuite `pipeline.analyze_sequence_record()`.

Le cache evite de recalculer les alignements lorsqu'un widget Streamlit declenche seulement un rerun d'interface.

Le resultat contient notamment :

```python
{
    "stats": {...},
    "similarity_results": [...],
    "mutation_report": {...} ou None,
    "variant_report": {...} ou None,
    "translation": {...} ou None,
    "orfs": [...],
    "metadata_warnings": [...]
}
```
