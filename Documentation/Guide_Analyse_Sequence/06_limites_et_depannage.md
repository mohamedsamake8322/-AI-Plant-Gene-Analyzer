# 6. Limites et depannage

## Erreur de cadre de lecture

Erreur typique :

```text
frame must be one of 1, 2, 3, -1, -2, or -3
```

Verifier :

- que l'interface utilise `config.READING_FRAMES` ;
- qu'aucun appel recent ne transmet `0`, `4` ou `-4` ;
- que `translate_dna()` recoit un cadre signe.

`frame=0` reste accepte comme alias historique de `+1` dans les appels programmatiques.

## Mutations absentes

Ce n'est pas toujours une erreur.

Cas possibles :

1. aucun meilleur match n'a ete trouve ;
2. le meilleur match est sous 85 % d'identite ;
3. la recherche a ete ignoree car la sequence est trop longue ;
4. le budget d'alignement a ete depasse ;
5. fournir une reference explicite pour comparer directement les sequences.

## Sequence trop longue

Les statistiques de base peuvent rester possibles, mais les alignements couteux sont interrompus au-dela de `MAX_ALIGNMENT_SEQUENCE_LENGTH` ou du budget `MAX_ALIGNMENT_CELL_BUDGET`.

Cette protection evite une consommation excessive de memoire dans l'alignement dynamique.

## VCF sans contig

Un VCF genere depuis une sequence collee sans en-tete FASTA ne peut pas inventer une coordonnee genomique.

Dans ce cas :

- le contig est `sequence` ;
- les positions sont locales a la reference fournie ;
- l'en-tete du VCF le precise.

Pour obtenir des coordonnees genomiques, fournir un FASTA avec un en-tete contenant un champ de contig reconnu, par exemple `contig=chr3`, puis utiliser une reference genome compatible.

## Similarite faible

Un score faible peut signifier :

- gene nouveau ou tres divergent ;
- sequence non codante ;
- fragment court ;
- mauvais type ADN/proteine ;
- reference absente de la base ;
- candidat masque par un prefiltre.

Le mode Deep Search peut aider, mais il ne transforme pas une reference distante en reference de mutation valide.

## Interpretation AI

L'onglet AI Interpretation est un moteur de regles local. Il ne consulte pas une API externe et ne fournit pas une preuve experimentale.
Les recommandations agricoles doivent etre traitees comme des pistes de travail.

## Contrôles techniques recommandes

```powershell
cd C:\Downloads\IA
.\.venv\Scripts\python.exe -m py_compile pipeline.py bioinformatics.py variant_analysis.py
.\.venv\Scripts\python.exe -m pytest -q test_bioinformatics.py test_pipeline_similarity.py test_i18n.py
```

Pour verifier uniquement la traduction :

```powershell
.\.venv\Scripts\python.exe -c "import bioinformatics as b; print(b.translate_all_frames('ATGAAATAGCCC'))"
```
