# Section Statistics

## Objectif

La section **Statistics** fournit une analyse descriptive et reproductible d'une séquence végétale. Elle fonctionne avec des séquences d'ADN et de protéines. Le contenu affiché dépend automatiquement du type de séquence détecté.

Les statistiques sont calculées localement par le pipeline bioinformatique. Les comparaisons à une espèce utilisent les agrégations PostgreSQL disponibles lorsque l'organisme est identifié dans les métadonnées de la séquence.

## Accès

1. Lancer l'application :

   ```bash
   streamlit run app.py
   ```

2. Charger un fichier FASTA/FA/TXT ou coller une séquence.
3. Cliquer sur **Analyze Sequence**.
4. Ouvrir l'onglet **Statistics**.

---

## Analyse ADN

### Composition nucléotidique

La section affiche :

- un graphique de composition des bases ;
- le nombre de bases `A`, `T`, `G` et `C` ;
- la longueur en paires de bases ;
- le pourcentage de GC ;
- le pourcentage de AT ;
- le ratio GC/AT ;
- la distribution de composition par fenêtre glissante.

Le profil GC local aide à repérer les variations de composition le long de la séquence. La taille de fenêtre est contrôlée par la configuration de l'application.

### Référence GC par espèce

Lorsque l'organisme est connu, le GC de la séquence peut être comparé à la moyenne des séquences du même organisme présentes dans la base.

La comparaison est affichée seulement si le nombre de séquences de référence atteint le minimum requis, actuellement `n >= 10`. Dans le cas contraire :

- le gauge générique reste affiché ;
- un message explicite indique que la référence espèce est indisponible ;
- aucune moyenne artificielle à `0 %` n'est affichée.

Le même mécanisme de référence est utilisé pour le codon usage et la longueur.

### GC skew et AT skew

Un profil par fenêtre affiche :

- le GC skew :

  ```text
  (G - C) / (G + C)
  ```

- le AT skew :

  ```text
  (A - T) / (A + T)
  ```

Une valeur n'est pas calculée lorsque le dénominateur est nul. Cela évite de transformer une fenêtre sans `G/C` ou sans `A/T` en valeur trompeuse.

### Contextes de méthylation de la cytosine

Pour les séquences d'ADN, les cytosines sont classées dans les contextes végétaux suivants :

- **CG** ;
- **CHG**, où `H` est `A`, `T` ou `C` ;
- **CHH**, où `H` est `A`, `T` ou `C`.

L'interface affiche pour chaque contexte :

- le nombre de cytosines détectées ;
- le pourcentage parmi les cytosines classifiables.

Il s'agit d'une estimation basée sur la séquence fournie. Elle ne remplace pas une mesure expérimentale de méthylation.

### Contrôle qualité

Le rapport qualité indique si la séquence respecte les règles utilisées par la collecte :

- longueur minimale ;
- proportion maximale de `N` ou de bases ambiguës ;
- raison du rejet lorsqu'une règle n'est pas respectée.

Le seuil de référence actuel est centralisé dans `quality_rules.py`, afin que l'ingestion et l'interface appliquent les mêmes règles.

### Complexité faible

Le scanner recherche des régions répétitives ou de faible complexité. La section affiche un avertissement lorsque la séquence contient de telles régions et indique leur couverture approximative en pourcentage.

Cette information est importante pour interpréter la similarité : une correspondance située principalement dans une région répétitive peut ne pas représenter une homologie biologique réelle.

### Codon usage

Pour les séquences d'ADN, les codons complets sont comptés dans le cadre de lecture utilisé par le calcul.

Le tableau affiche jusqu'aux dix codons les plus divergents lorsque la référence espèce est disponible :

- fréquence dans la séquence analysée ;
- fréquence moyenne dans l'organisme ;
- différence entre les deux fréquences.

Si la séquence n'est pas un multiple de trois, les bases terminales incomplètes sont exclues du calcul et un avertissement est affiché.

#### Codon Adaptation Index

Le **CAI** est affiché comme un score unique lorsque la référence codon de l'organisme est disponible. Il est calculé comme une moyenne géométrique des fréquences relatives des codons synonymes.

Dans l'implémentation actuelle, la référence est une distribution globale par organisme. Elle ne constitue pas encore un jeu dédié de gènes fortement exprimés, tels que des gènes ribosomiques ou de ménage. Le CAI doit donc être interprété comme une approximation comparative.

### Statistiques détaillées et ORF

Le résumé des six cadres de lecture, affiché avant l'analyse puis réutilisable
dans les résultats, distingue explicitement :

- **Complete ORFs** : ORF `ATG -> codon stop` trouvé dans le même cadre ;
- **Truncated ORFs** : ORF `ATG -> fin de séquence` lorsqu'aucun codon stop
  n'est rencontré avant la fin.

La colonne **Longest ORF (bp, stop codon included)** inclut le codon stop dans
la longueur lorsqu'il est présent. Cette convention explique un écart de trois
bases avec un outil qui mesure uniquement jusqu'au dernier codon avant l'arrêt.

La zone détaillée indique notamment :

- si la longueur est divisible par trois ;
- la présence d'un codon de démarrage `ATG` dans un cadre quelconque ;
- la présence d'un codon stop dans un cadre quelconque ;
- la présence d'un ORF complet avec démarrage et arrêt dans le même cadre.

La présence indépendante d'un `ATG` et d'un codon stop ne prouve pas qu'ils appartiennent au même ORF. L'application affiche une précision dans ce cas et le résultat complet est disponible dans les informations de traduction/ORF.

### Motifs réglementaires

Les motifs connus détectés dans la séquence sont listés avec :

- leur nom ;
- leur séquence ;
- leurs positions de début et de fin ;
- la séquence effectivement trouvée.

L'absence de motif connu ne signifie pas que la séquence ne possède aucune fonction réglementaire.

### Sites de restriction

Un panneau dédié scanne les sites de plusieurs enzymes courantes, notamment :

- EcoRI ;
- BamHI ;
- HindIII ;
- NdeI ;
- XhoI ;
- KpnI ;
- NotI ;
- SmaI ;
- PstI.

Chaque résultat indique l'enzyme, le motif, la position et la séquence trouvée. Les motifs chevauchants sont pris en compte.

Les sites de restriction ne sont pas calculés pour les séquences protéiques.

### Primer design hints

Pour une séquence d'au moins 25 bases, l'application propose des candidats naturels aux extrémités :

- les 20 premières bases pour l'amorce avant ;
- les 20 dernières bases, complémentées et inversées, pour l'amorce arrière ;
- une estimation de Tm par la règle de Wallace : `2 °C` par base `A/T` et `4 °C` par base `G/C` ;
- la présence ou non d'un GC clamp en extrémité 3'.

Ces valeurs sont des indications préliminaires et ne remplacent pas une vérification complète des amorces : dimères, hairpins, spécificité, concentration en sels et conditions de PCR ne sont pas modélisés ici.

Aucune estimation n'est affichée pour une séquence trop courte ou une protéine.

---

## Analyse protéique

Lorsque l'entrée est détectée comme une protéine, les éléments ADN ne sont pas affichés. La section présente :

- la composition en acides aminés ;
- la longueur en acides aminés ;
- le nombre de résidus uniques ;
- la diversité des résidus ;
- le résidu le plus abondant ;
- un graphique de composition ;
- le score GRAVY, basé sur l'échelle d'hydrophobicité de Kyte-Doolittle ;
- l'indice d'instabilité ;
- l'indice aliphatique ;
- les motifs protéiques connus lorsqu'ils sont détectés.

Les caractères ambigus tels que `X` sont exclus des calculs biochimiques qui nécessitent un résidu connu. Ils ne doivent pas être interprétés comme des acides aminés hydrophobes ou chargés.

Les sites de restriction, la méthylation, le GC skew, le codon usage et les amorces PCR sont non applicables aux protéines.

---

## Références par organisme

Les comparaisons espèce utilisent un helper commun :

```python
get_organism_reference(organism, metric, min_n=10)
```

Le résultat normalisé contient :

```python
{
    "available": bool,
    "value": value_or_none,
    "n": sample_count,
    "fallback_reason": message_or_none,
}
```

Ce contrat est utilisé pour :

- le GC content ;
- le codon usage ;
- la longueur.

Lorsque `available` vaut `False`, l'interface affiche un message de fallback standard et omet la comparaison plutôt que d'afficher une valeur inventée.

---

## Export et paragraphe méthodes

Les résultats peuvent être exportés dans les formats déjà proposés par l'application :

- JSON ;
- CSV ;
- HTML ;
- XLSX ;
- FASTA ;
- GFF3 pour les séquences ADN.

Le CSV est un export plat `Metric,Value`, pratique pour un script ou une importation rapide. Le XLSX est un classeur structuré pour le tri, le filtrage et la comparaison manuelle. Ils ne sont donc pas strictement redondants.

Le FASTA contient la séquence analysée, avec un en-tête minimal indiquant l'identifiant, la longueur et le GC. Il peut être réutilisé directement dans BLAST, un aligneur ou un autre pipeline bioinformatique.

Le GFF3 décrit les features détectées avec des coordonnées 1-based inclusives :

- ORF, avec le brin et le cadre lorsqu'ils sont connus ;
- motifs réglementaires ;
- sites de restriction.

Les motifs et sites de restriction utilisent le brin `.` car le scanner actuel signale les occurrences dans la séquence fournie sans attribuer de brin. Le fichier peut être chargé dans IGV, Artemis ou un visualiseur GFF3 compatible.

Les exports comprennent les nouveaux résultats lorsqu'ils sont disponibles, notamment les sites de restriction, les amorces et les propriétés protéiques.

Le bouton **Copy methods paragraph** génère également un fichier texte contenant un paragraphe prêt à copier dans un rapport ou un article. Le paragraphe inclut uniquement les données disponibles : une référence espèce absente est omise proprement.

Exemple de formulation :

```text
The 7107 bp sequence exhibited a GC content of 41.93% (species average: 42.10%, n=27) with a length +118 bp from the species mean (6989 bp, n=27), consistent with a complete open reading frame (start-to-stop, same frame).
```

Le paragraphe est une synthèse automatique des calculs. Il doit être relu et adapté au protocole expérimental avant publication.

---

## Limites générales

- Les statistiques dépendent de la qualité et de la représentativité de la séquence fournie.
- Une référence espèce n'est disponible que si l'organisme est identifié et suffisamment représenté dans la base.
- Les moyennes espèce ne constituent pas nécessairement une population taxonomiquement homogène.
- Les prédictions d'ORF, de CAI, de motifs et de sites de restriction sont in silico.
- Les estimations Tm et les indices protéiques sont des approximations, pas des résultats expérimentaux.
- L'absence de résultat ne prouve pas l'absence biologique du phénomène recherché.

## Modules principaux

- `app.py` : affichage Streamlit et formatage des résultats ;
- `pipeline.py` : orchestration de l'analyse ;
- `bioinformatics.py` : calculs de séquence ;
- `organism_reference.py` : gestion commune des références par organisme ;
- `scripts/postgres_utils.py` : agrégations PostgreSQL ;
- `export_utils.py` : exports JSON, CSV, HTML et XLSX.
