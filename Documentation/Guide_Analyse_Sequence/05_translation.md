# 5. Translation

## Objectif

L'onglet **Translation** traduit une sequence ADN en acides amines selon un cadre de lecture signe.
Il ne traduit pas une entree deja proteique.

## 5.1 Cadres disponibles

- `+1` : brin fourni, decalage 0 ;
- `+2` : brin fourni, decalage 1 ;
- `+3` : brin fourni, decalage 2 ;
- `-1` : reverse-complement, decalage 0 ;
- `-2` : reverse-complement, decalage 1 ;
- `-3` : reverse-complement, decalage 2.

Le signe est important. Un cadre negatif n'est pas simplement un cadre positif avec un nombre different : il est calcule sur le reverse-complement.

## 5.2 Resultats du cadre choisi

L'interface affiche :

- cadre selectionne ;
- longueur de la proteine produite ;
- nombre de codons complets ;
- statut `Complete` ou `Open` ;
- bases restantes hors codon complet ;
- sequence proteique, avec stop si present ;
- carte des codons ;
- position de debut et de fin de chaque codon ;
- acide amine ;
- indication du codon stop.

Un cadre `Open` signifie seulement qu'aucun stop n'a ete rencontre dans le parcours traduit. Cela ne prouve pas que la sequence code une proteine fonctionnelle.

## 5.3 Comparaison des six cadres

La table Six-frame comparison resume les six traductions :

- cadre ;
- brin ;
- longueur proteique ;
- presence d'un stop ;
- nombre de codons complets ;
- bases restantes.

Le cadre recommande est une suggestion informatique basee sur l'existence d'un stop et la longueur traduite. Ce n'est pas une preuve de l'ORF biologique reel.

## 5.4 ORF predits

La section Predicted ORFs affiche les regions ouvertes detectees :

- brin/cadre ;
- debut et fin en nucleotides ;
- longueur ;
- longueur proteique ;
- caractere complet ou tronque.

Un export GFF3 des ORF est disponible lorsque des ORF sont detectes.

## 5.5 Sequences complementaires

L'onglet affiche :

- le complement du segment ;
- le reverse-complement du segment.

Ces vues sont utiles pour controler manuellement l'orientation et comprendre les cadres negatifs.

## 5.6 Export proteique

La proteine du cadre choisi peut etre telechargee en FASTA.
Le nom du fichier conserve le cadre, par exemple `translated_frame_+1.fasta` ou `translated_frame_-1.fasta`.

## Fonctions principales

- `bioinformatics.translate_dna()`
- `bioinformatics.translate_all_frames()`
- `bioinformatics.translation_codon_rows()`
- `bioinformatics.reverse_complement()`
- `pipeline.analyze_sequence_record()`

## Compatibilite des anciens appels

Les anciens appels programmatiques avec `frame=0` sont conserves comme alias de `+1`.
L'interface et le nouveau contrat public utilisent cependant exclusivement les cadres signes `+1..+3/-1..-3`.
