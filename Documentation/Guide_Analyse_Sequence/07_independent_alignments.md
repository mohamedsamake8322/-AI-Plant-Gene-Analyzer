# 7. Independent Tools - Alignments

## Architecture

L'onglet **Alignments** est une interface independante du pipeline principal. Il ne depend pas de la base PostgreSQL ni du meilleur match de Similarity.

Il reutilise le moteur partage `alignment_engine.py` :

- `needleman_wunsch()` pour le global ;
- `smith_waterman()` pour le local ;
- `star_alignment()` pour le MSA ;
- `alignment_statistics()` pour les metriques ;
- `get_score()` et les matrices de substitution.

Aucun second moteur d'alignement n'est cree dans l'interface.

## Fonctionnalites actuellement disponibles

### MSA

- Entree FASTA multiple ou sequences separees ;
- type `Auto`, `DNA` ou `Protein` ;
- choix de la sequence de reference du star alignment ;
- avertissement explicite sur l'effet du choix de reference ;
- nombre de sequences ;
- longueur alignee ;
- score de conservation ;
- affichage Plotly ;
- exemple pret a l'emploi ;
- avertissement reutilisant `bio.detect_low_complexity_regions()` ;
- consensus et colonnes variables ;
- navigation par fenetre ;
- exports FASTA aligne, Clustal, CSV, PHYLIP et NEXUS ;
- metadonnees de reproductibilite et hash SHA-256 des entrees.

Le MSA actuel est un **star alignment reference-guided**. Ce n'est pas encore un MSA progressif ou iteratif de type MAFFT/MUSCLE.

### Pairwise

- type Auto, ADN ou proteine ;
- Needleman-Wunsch global ;
- Smith-Waterman local ;
- score ;
- identite ;
- matches ;
- mismatches ;
- gaps ;
- couverture des deux sequences pour le global ;
- longueur de l'alignement local ;
- comparaison separee des scores global et local ;
- matrices proteiques BLOSUM45, BLOSUM62, BLOSUM80, PAM30, PAM70 et PAM250 ;
- penalites gap-open et gap-extension configurables.

## Garde-fous

Les memes limites que le reste de l'application sont reutilisees :

- `config.MAX_ALIGNMENT_SEQUENCE_LENGTH` ;
- `config.MAX_ALIGNMENT_CELL_BUDGET`.

Le MSA estime le cout cumule des alignements reference-vers-sequence. Une execution est refusee si la longueur ou le budget est depasse.

## Interpretation debutant

- **Global** : aligne les sequences du debut a la fin.
- **Local** : cherche la meilleure region commune.
- **MSA** : aligne plusieurs sequences dans une grille commune.
- **Match** : caracteres identiques dans une colonne.
- **Mismatch** : caracteres differents dans une colonne.
- **Gap** : insertion d'un tiret pour maintenir la correspondance.

## Phase 2 restante

- coloration ADN/proteine plus detaillee.
- conservation visuelle avancee et scores par colonne.

## Phase 3 restante

- sites variables et informatifs pour la phylogenie ;
- rapport de reproductibilite enrichi avec version applicative ;
- MAFFT, MUSCLE ou Clustal Omega en integration optionnelle.

MAFFT et MUSCLE sont des binaires externes. Leur integration devra gerer la disponibilite du binaire, les sous-processus, les erreurs, les versions et les differences Windows/Linux. Ils ne doivent pas etre ajoutes comme une simple dependance Python sans strategie de deploiement.

## Limites actuelles

- Le star alignment depend du choix de reference.
- Il ne remplace pas un MSA progressif professionnel pour de nombreux sequences divergentes.
- Les coordonnees de debut/fin du traceback local ne sont pas encore exposees par le contrat du moteur.
- Aucune E-value n'est calculee : une mesure de significativite ne doit pas etre ajoutee avec une formule non calibree.
