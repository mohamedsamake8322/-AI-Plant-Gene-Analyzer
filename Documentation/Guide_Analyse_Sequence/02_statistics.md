# 2. Statistics

## Objectif

L'onglet **Statistics** decrit la sequence elle-meme, sans avoir besoin d'une reference externe.
Il est donc normalement disponible meme si aucune similarite n'est trouvee.

## 2.1 Sequence ADN

L'application affiche ou calcule :

- longueur en paires de bases ;
- contenu GC ;
- contenu AT ;
- ratio GC/AT ;
- comptage A, T, G, C et N ;
- pourcentage de bases ambiguës ;
- GC de troisieme position (`GC3`) ;
- longueur et nombre d'ORF ;
- presence de ATG quelque part dans les six cadres ;
- presence d'un codon stop quelque part dans les six cadres ;
- presence d'un ORF complet start-to-stop dans un meme cadre.

Attention : `Contains ATG in any frame` et `Contains stop codon in any frame` sont deux observations independantes. Elles ne prouvent pas qu'un gene complet existe. Pour cette conclusion, il faut regarder `Complete ORF found`.

## 2.2 Profil GC

Le profil GC utilise une fenetre glissante selectionnee dans la barre laterale.
Le graphique montre le pourcentage GC local le long de la sequence.

Un profil local ne dit pas qu'une region est bonne ou mauvaise. Il decrit seulement la composition observee.

Le profil GC skew et AT skew montrent aussi l'equilibre entre G/C et A/T dans chaque fenetre.

## 2.3 Methylation et qualite

Pour l'ADN vegetal, l'application peut afficher les contextes de cytosine :

- CG ;
- CHG ;
- CHH, avec H = A, T ou C.

Le rapport de qualite indique notamment le pourcentage de N et si la sequence passe le filtre courant.

## 2.4 Codons, motifs et restrictions

Pour l'ADN, l'onglet peut aussi montrer :

- usage des codons ;
- comparaison d'usage des codons avec une reference d'organisme si disponible ;
- Codon Adaptation Index, quand une reference est disponible ;
- motifs regulateurs connus, par exemple TATA-box et CAAT-box ;
- sites de restriction connus ;
- suggestions de primers avec estimation de Tm.

Ces resultats sont des calculs ou des detections de motifs. Ils ne constituent pas une preuve d'expression ou de fonction.

## 2.5 Sequence proteique

Pour une proteine, Statistics affiche :

- longueur en acides amines ;
- nombre de residus uniques ;
- distribution des acides amines ;
- residu le plus abondant ;
- masse moleculaire estimee ;
- point isoelectrique estime ;
- hydrophobicite/GRAVY ;
- indices d'instabilite et d'aliphaticite.

Le GC n'est pas applicable directement a une sequence proteique.

## Fonctions principales

- `bioinformatics.sequence_statistics()`
- `bioinformatics.generate_protein_statistics()`
- `bioinformatics.nucleotide_distribution()`
- `bioinformatics.gc_skew_profile()`
- `bioinformatics.cytosine_methylation_context()`
- `bioinformatics.find_orfs()`
