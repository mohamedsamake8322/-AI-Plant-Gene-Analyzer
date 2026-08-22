# Plant Gene Analyzer — Structure canonique proposée pour un enregistrement de gène

## Principes directeurs

1. **`organism` devient un objet structuré, jamais une chaîne libre.** C'est la cause racine des deux bugs qu'on vient de corriger (code court PLAZA au lieu du nom complet, description de protéine devinée depuis un header FASTA). Une chaîne libre n'a aucune garde-fou ; un objet avec des clés typées, si.
2. **Chaque champ garde sa provenance.** Tu fusionnes 6+ sources (NCBI, UniProt, KEGG, PLAZA, PlantTFDB, Ensembl). Sans traçabilité par champ, un futur bug de fusion redevient indétectable — exactement le problème qu'on a eu aujourd'hui.
3. **`pathways` et `traits` sont des listes d'objets typés, jamais des chaînes ou des IDs bruts** — pour que tes têtes de classification (GO terms, trait, TF family, pathway) puissent consommer directement sans reparser.
4. **Un seul hash de déduplication (`sequence_hash`)**, déjà en place côté infra — à garder comme clé de dédup canonique.
5. **Séparer "ce qui sert à l'entraînement" de "ce qui sert à la traçabilité/qualité"** pour que le code de training puisse ignorer proprement les métadonnées non pertinentes.

## Structure proposée

```json
{
  "gene_id": "XM_026013608.1",
  "symbol": "LOX2",

  "organism": {
    "scientific_name": "Chenopodium quinoa",
    "common_name": "quinoa",
    "taxon_id": 63459,
    "codes": {
      "kegg": "cqi",
      "plaza": "cqu",
      "ncbi": "chenopodium_quinoa"
    }
  },

  "sequence": {
    "value": "ATGGCTTCG...",
    "type": "dna",
    "length": 1824,
    "sequence_hash": "a3f9e1c2..."
  },

  "description": "Lipoxygenase 2, involved in jasmonic acid biosynthesis",

  "traits": [
    {"trait": "stress response", "category": "biological", "evidence": "GO annotation", "source": "uniprot", "retrieved_at": "2026-08-19T10:22:00Z"},
    {"trait": "lodging resistance candidate", "category": "biological", "evidence": "manual curation", "source": "user", "retrieved_at": "2026-08-19T10:22:00Z"}
  ],

  "pathways": [
    {"id": "map00591", "name": "Linoleic acid metabolism", "source": "kegg"},
    {"id": "map04075", "name": "Plant hormone signal transduction", "source": "kegg"}
  ],

  "go_terms": [
    {"id": "GO:0009611", "name": "response to wounding", "aspect": "P", "source": "uniprot"}
  ],

  "tf_family": {
    "family": "bHLH",
    "source": "planttfdb"
  },

  "expression_profiles": [
    {"tissue": "stem", "condition": "panicle load stress", "source": "expression_atlas", "url": "..."}
  ],

  "publications": [
    {"pmid": "31234567", "title_ref": "internal note, not full title", "source": "ncbi"}
  ],

  "external_links": {
    "ncbi": "https://www.ncbi.nlm.nih.gov/nuccore/XM_026013608.1",
    "uniprot": "https://www.uniprot.org/uniprotkb/A0A...",
    "kegg": "https://www.genome.jp/dbget-bin/www_bget?cqi:...",
    "plaza": "https://bioinformatics.psb.ugent.be/plaza/..."
  },

  "quality": {
    "reviewed": false,
    "protein_existence": 2,
    "annotation_score": 3,
    "sources_seen": ["ncbi", "uniprot", "kegg", "plaza"]
  },

  "provenance": {
    "organism": {"source": "ncbi", "retrieved_at": "2026-08-19T10:22:00Z"},
    "sequence": {"source": "ncbi", "retrieved_at": "2026-08-19T10:22:00Z"},
    "pathways": {"source": "kegg", "retrieved_at": "2026-08-19T10:24:00Z"}
  },

  "date_added": "2026-08-19T10:22:00Z",
  "schema_version": "2.1"
}
```

## Pourquoi ces choix précis

- **`organism.codes`** centralise tous les codes courts par base (KEGG `cqi`, PLAZA `cqu`, etc.) plutôt que de les laisser fuiter dans le champ `organism` principal comme on l'a vu aujourd'hui. Un futur bug similaire deviendrait immédiatement visible (un code dans `scientific_name` sauterait aux yeux en revue).
- **`traits[].category`** règle directement ton problème ouvert de filtrage des mots-clés UniProt : chaque trait est tagué `biological` vs `technical` dès la collecte, au lieu d'un filtrage a posteriori fragile.
- **`pathways`** est déjà dans le format `{"id", "name", "source"}` que tu voulais vérifier — donc cette structure valide directement ta question en attente.
- **`provenance`** par champ (pas juste un `source` global sur tout l'enregistrement) — utile parce que tes enregistrements sont fusionnés depuis plusieurs sources à des moments différents ; ça permet un audit champ par champ si un bug de fusion réapparaît.
- **`quality`** séparé de `traits`/`pathways` — pour que ton code d'entraînement (têtes de classification) puisse ignorer ce sous-objet sans risque de fuite de métadonnées non pertinentes dans les features.
- **`sequence` en sous-objet** plutôt que des champs à plat (`sequence`, `sequence_type`, `length`, `sequence_hash` séparés) — regroupe tout ce qui concerne la séquence elle-même, plus facile à valider comme unité.

## Migration depuis ta structure actuelle

Pas besoin de tout réécrire à la main : un script de migration peut lire `master_plant_db.json` existant et remapper vers cette structure, avec une règle spéciale pour `organism` — valider chaque valeur contre une liste blanche de tes ~10 espèces réelles, et marquer (`quality.needs_review: true`) tout enregistrement dont l'`organism` ne matche pas, plutôt que de deviner. Ça te donnerait un rapport direct des enregistrements pollués à corriger (probablement concentrés dans les runs NCBI protéine, vu le bug qu'on vient de trouver).

Tu veux que je t'écrive ce script de migration ?
