# 3. Similarity

## Objectif

L'onglet **Similarity** compare la sequence utilisateur avec des sequences candidates de la base de genes vegetaux.
La similarite sert a trouver une reference ou un homologue. Elle ne constitue pas automatiquement une identification certaine.

## 3.1 Selection des candidats

Le pipeline peut utiliser :

- la base locale JSON ;
- PostgreSQL lorsque la connexion est disponible ;
- un prefiltre k-mer ;
- un prefiltre par longueur ;
- le mode Deep Search pour elargir l'exploration.

Les filtres sont des optimisations de performance. Ils peuvent reduire le nombre de candidats evalues.
L'interface indique la source, le nombre de candidats et le temps de recherche.

## 3.2 Score principal

Le score principal est l'identite globale issue d'un alignement Needleman-Wunsch.
Il compare la sequence complete a la reference candidate.

La recherche peut aussi calculer :

- identite locale Smith-Waterman ;
- couverture locale de la requete ;
- nombre de gaps ;
- longueur alignee ;
- correspondances et mismatches ;
- score d'alignement ;
- approximation de significativite locale.

L'identite globale reste le score principal du classement. Une region locale parfaite ne remplace pas une faible couverture de la sequence complete.

## 3.3 Lecture des resultats

Chaque resultat contient notamment :

- nom ou gene ID ;
- trait ;
- organisme ;
- accession ;
- score de similarite ;
- methode d'alignement ;
- description ;
- couverture locale ;
- carte d'alignement.

La section de conservation affiche une comparaison multiple independante entre la requete et les meilleurs matchs. Elle peut placer les gaps differemment de l'alignement utilise pour le classement.

## 3.4 Limites et protection Mutation

Un meilleur match distant peut etre biologiquement utile comme homologue, mais pas comme reference de mutation.

Par defaut :

```text
MIN_MUTATION_REFERENCE_IDENTITY = 85.0
```

Si le meilleur match est sous ce seuil :

- Similarity affiche toujours le resultat ;
- Statistics reste disponible ;
- Mutations ne produit pas de faux rapport de mutations ;
- un message indique qu'aucune reference suffisamment proche n'a ete trouvee.

Pour comparer malgre tout une sequence distante, utiliser le champ **Explicit reference** et fournir une vraie reference adaptee.

## Fonctions principales

- `similarityengine.find_similar_genes()`
- `similarityengine.find_similar_genes_deep()`
- `similarityengine.compare_with_database()`
- `similarityengine.aligned_similarity()`
- `alignment_engine.needleman_wunsch()`
- `alignment_engine.smith_waterman()`
- `visualization.plot_similarity_scores()`
