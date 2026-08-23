# Plant Gene Analyzer — Vérification des analyses

## Lancer l'application

```bash
streamlit run app.py
```

Ouvrir ensuite http://localhost:8501, coller une séquence et cliquer sur **Analyze Sequence**.

## Sections disponibles

| Section | Fonction |
|---|---|
| **Statistics** | Nettoyage, validation, longueur, GC%, AT%, ratio GC/AT, fréquences A/T/G/C/N, bases ambiguës, GC3, ORF et motifs. Affiche aussi le profil GC par fenêtre. |
| **Similarity** | Recherche dans la base locale avec alignement global Needleman-Wunsch. Affiche les meilleurs matchs, scores, couverture, gaps, carte d'alignement et confiance. |
| **Mutations** | Compare la séquence au meilleur match et affiche substitutions, insertions, délétions, identités et positions. Nécessite un match de référence. |
| **Translation** | Traduit l'ADN dans le cadre choisi, affiche les six cadres, le statut ORF, le complément et le reverse-complement. |
| **AI Interpretation** | Génère une interprétation déterministe par règles : profil, GC, similarité, mutations, fonction, stress et recommandations agricoles. Ce n'est pas une API d'IA externe. |
| **Raw Sequence** | Affiche la séquence nettoyée et permet de télécharger un FASTA ou un rapport texte. Une annotation locale est également disponible. |
| **Alignments** | Compare manuellement plusieurs séquences avec MSA, Needleman-Wunsch global et Smith-Waterman local. |
| **Distance Matrix** | Aligne plusieurs séquences et calcule une matrice avec Hamming, Jukes-Cantor, Kimura ou PAM. Export CSV disponible. |
| **Phylogeny** | Construit un arbre UPGMA ou Neighbor Joining à partir d'une matrice Kimura. Export Newick disponible. |
| **Protein Analysis** | Analyse une protéine : validation, longueur, composition, masse moléculaire, point isoélectrique et hydrophobicité. |

## Outils indépendants

Les outils **Alignments**, **Distance Matrix**, **Phylogeny** et **Protein
Analysis** ne dépendent pas du résultat de l'analyse principale. Ils sont
également affichés sur l'écran d'accueil, avant toute analyse, et possèdent
leurs propres champs de séquence et boutons d'exécution.

Les autres sections (**Statistics**, **Similarity**, **Mutations**,
**Translation**, **AI Interpretation** et **Raw Sequence**) restent liées à la
séquence principale : elles apparaissent après l'action **Analyze Sequence**.
## Test avec la séquence de référence

Pour la séquence ADN de 300 bp testée :

- GC : **67,67 %**
- AT : **32,33 %**
- 1 ORF détecté
- Traduction de 22 acides aminés
- 2 motifs détectés
- Aucun meilleur match trouvé dans la base locale du test

L'absence de match explique pourquoi l'onglet **Mutations** ne produit pas de rapport automatique. Les sections **Alignments**, **Distance Matrix** et **Phylogeny** utilisent leurs propres champs de saisie pour les séquences de comparaison.

## Interprétation

Les statistiques, distances et alignements sont des calculs reproductibles. Les textes de **AI Interpretation** sont heuristiques et ne constituent pas une preuve expérimentale de fonction, de résistance ou de trait agronomique.

## Vérifications effectuées

- Serveur Streamlit : endpoint de santé HTTP 200
- Compilation Python : réussie
- Pipeline ADN : réussie
- Alignement global/local et MSA : réussis
- Matrice de distances : réussie
- UPGMA et Neighbor Joining : réussis
- Traduction et analyse protéique : réussies
