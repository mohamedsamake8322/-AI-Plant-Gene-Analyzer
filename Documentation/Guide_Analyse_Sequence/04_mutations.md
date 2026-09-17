# 4. Mutations

## Objectif

L'onglet **Mutations** compare une requete a une reference et decrit les differences observees.
Une difference est techniquement reelle dans l'alignement, mais son interpretation comme mutation depend de la qualite et de la proximite de la reference.

## 4.1 Quelle reference est utilisee ?

Priorite :

1. la sequence saisie dans **Explicit reference** ;
2. sinon le meilleur match de Similarity si son identite est au moins de 85 % ;
3. sinon aucun rapport mutationnel n'est produit.

Une sequence collee sans en-tete FASTA ne fournit pas de contig genomique. Les positions sont locales a la sequence comparee.

## 4.2 Resume et indicateurs

L'interface affiche :

- substitutions ;
- indels regroupes ;
- transitions ;
- transversions ;
- taux de substitution ;
- identite sur bases alignees uniquement ;
- identite sur l'alignement complet, gaps compris ;
- positions comparees sans gaps.

Difference entre les identites :

- `aligned bases only` = matches / colonnes non-gap ;
- `full alignment` = matches / longueur totale de l'alignement.

Avec un indel, ces deux pourcentages peuvent diverger fortement.

## 4.3 Substitutions ADN

Une substitution est une colonne ou la base requete est differente de la base reference.

- Transition : A↔G ou C↔T.
- Transversion : changement entre une purine et une pyrimidine.

Le tableau indique :

- position reference ;
- position requete ;
- base reference ;
- base requete ;
- type transition/transversion ;
- codon avant/apres lorsque le cadre est exploitable ;
- acide amine avant/apres ;
- consequence ;
- impact ;
- score BLOSUM62.

## 4.4 Consequences ADN

Pour un cadre de lecture donne, `variant_analysis.analyze_variants()` classe les substitutions en :

- `silent` : le codon change mais l'acide amine reste identique ;
- `missense` : l'acide amine change ;
- `nonsense` : un codon stop est gagne ;
- `readthrough` : un codon stop est perdu ;
- `incomplete_codon` : la substitution ne permet pas un codon complet ;
- `downstream_of_frameshift` : la lecture est decalee par un indel en amont.

Exemple :

```text
Codon : ACT -> ACG
AA    : T   -> T
Effet : silent
```

## 4.5 BLOSUM62 et impact

Pour les acides amines standards, BLOSUM62 fournit un score de substitution.

Dans l'interface :

- score positif : impact classe conservative ;
- score negatif : impact classe radical ;
- score nul : impact classe neutral ;
- residu stop ou ambigu : not_scored lorsque la paire n'est pas dans la matrice standard.

Cette classification est un indicateur de compatibilite physicochimique, pas une preuve de conservation fonctionnelle.

## 4.6 Indels et frameshifts

Les colonnes consecutives d'une insertion ou d'une deletion sont regroupees en un evenement biologique.

- longueur multiple de 3 en ADN : indel in-frame ;
- longueur non multiple de 3 : frameshift probable.

La carte affiche les substitutions comme losanges et les indels comme triangles.
Le clic sur un element affiche sa position, ses bases, sa longueur et sa consequence.

## 4.7 Filtre important

Le filtre **Show important mutations only** conserve principalement :

- missense ;
- nonsense ;
- readthrough ;
- radical ;
- variants downstream d'un frameshift ;
- indels frameshift.

Il sert a reduire la table aux variants qui meritent une revue prioritaire.

## 4.8 Frequence par region

La fenetre de frequence regroupe les variants par intervalles de la sequence.
Le tableau indique :

- region en bp ;
- substitutions ;
- indels ;
- nombre total de variants ;
- consequences observees dans la region.

Cette vue aide a repérer des zones concentrees en differences, mais elle ne remplace pas une analyse de couverture ou de qualite de variant calling.

## 4.9 Exports

### CSV detaille

L'export contient notamment :

- positions reference/requete ;
- alleles ;
- type ;
- consequence ;
- impact class ;
- BLOSUM62 ;
- codons ;
- acides amines ;
- longueur et frameshift des indels.

### VCF

Le VCF exporte les substitutions compatibles avec un format VCF.
Si un contig existe dans les metadonnees FASTA, il est utilise.
Sinon le fichier porte une mention indiquant que le contig `sequence` correspond a des coordonnees locales et non genomiques.

## Fonctions principales

- `bioinformatics.detect_mutations()`
- `variant_analysis.analyze_variants()`
- `variant_analysis.classify_dna_substitution()`
- `variant_analysis.classify_protein_substitution()`
- `variant_analysis.group_indels()`
- `visualization.plot_mutation_map()`
- `export_utils.export_mutations_csv()`
- `export_utils.export_mutations_vcf()`
