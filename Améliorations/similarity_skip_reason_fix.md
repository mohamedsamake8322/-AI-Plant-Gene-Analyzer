# Correctif — Onglet Similarity : distinguer "non lancée" de "lancée, vide"

## Constat

`app.py`, onglet Similarity (~ligne 1258) :
```python
if not similarity_results:
    st.warning("No similarity results available.")
```
Ce message s'affiche identiquement dans deux cas très différents :
1. La recherche a été **sautée** (séquence > `config.MAX_ALIGNMENT_SEQUENCE_LENGTH`,
   voir `pipeline.py`) — actuellement seule trace : une phrase libre dans
   `metadata_warnings`, jamais vérifiée par cet endroit précis du code.
2. La recherche a **réellement eu lieu** et n'a trouvé aucun match.

Un chercheur qui ne remarque pas l'avertissement affiché plus haut sur la
page peut conclure à tort "ce gène n'a pas d'homologue en base", alors que
la vraie réponse est "la recherche n'a jamais eu lieu".

## Correctif

- **`pipeline.py`** : ajouter un champ structuré explicite au dict retourné,
  ex. `"similarity_skipped_reason": "sequence_too_long"` (ou `None` si la
  recherche a bien eu lieu). Le peupler au même endroit où
  `pipeline_warnings.append(...)` est déjà appelé pour ce cas — ne pas
  dupliquer la condition de seuil, la poser une seule fois.
- **`app.py`**, onglet Similarity : remplacer le message générique par une
  vérification de ce champ :
  ```python
  skipped_reason = result.get("similarity_skipped_reason")
  if skipped_reason == "sequence_too_long":
      st.warning("Similarity search was not run (sequence exceeds the alignment length threshold) — this is not the same as 'no matches found'.")
  elif not similarity_results:
      st.warning("No similarity results available.")
  else:
      ...
  ```
- **Bénéfice structurel** : en passant par un champ dédié plutôt qu'un
  texte libre, toute future raison de skip (ex: base de données non
  chargée, erreur technique déjà catchée ailleurs dans pipeline.py) peut
  s'ajouter au même endroit sans réintroduire ce même problème d'ambiguïté.
