import hashlib
import json
from collections import Counter

d = json.load(open("Data/clean/linked_genes.json", encoding="utf-8"))
genes = d["genes"]
meta = d.get("metadata", {})

print(f"Nombre total de gènes : {len(genes)}")
print(f"Stratégie : {meta.get('strategy')}\n")

# --- 1. Échantillons complets à lire à l'œil ---------------------------------
print("=" * 70)
print("ÉCHANTILLON 1 : un gène lié via UniProt")
print("=" * 70)
uniprot_examples = [g for g in genes if g.get("source") == "uniprot+ncbi_linked"]
if uniprot_examples:
    g = uniprot_examples[0]
    print(f"gene_id (accession NCBI) : {g.get('gene_id')}")
    print(f"symbol                   : {g.get('symbol')}")
    print(f"organism                 : {g.get('organism')}")
    print(f"sequence_type            : {g.get('sequence_type')}")
    print(f"sequence (longueur {len(g.get('sequence',''))})  : {g.get('sequence','')[:80]}...")
    print(f"description              : {g.get('description')}")
    print(f"traits                   : {g.get('traits')}")
    print(f"pathways                 : {g.get('pathways')}")
    ann = g.get("annotations") or {}
    print(f"annotations.go_terms     : {ann.get('go_terms', [])[:3]} ... ({len(ann.get('go_terms', []))} total)")
    print(f"annotations.tf_family    : {ann.get('tf_family')}")
    print(f"protein_sequence présente: {bool(g.get('protein_sequence'))} (len={len(g.get('protein_sequence') or '')})")
    print(f"external_links           : {g.get('external_links')}")
else:
    print("Aucun exemple trouvé (source=uniprot+ncbi_linked) !")

print("\n" + "=" * 70)
print("ÉCHANTILLON 2 : un gène lié via PlantTFDB")
print("=" * 70)
tf_examples = [g for g in genes if g.get("source") == "planttfdb+ncbi_linked"]
if tf_examples:
    g = tf_examples[0]
    print(f"gene_id (accession NCBI) : {g.get('gene_id')}")
    print(f"symbol                   : {g.get('symbol')}")
    print(f"organism                 : {g.get('organism')}")
    print(f"sequence_type            : {g.get('sequence_type')}")
    print(f"sequence (longueur {len(g.get('sequence',''))})  : {g.get('sequence','')[:80]}...")
    print(f"description              : {g.get('description')}")
    print(f"traits                   : {g.get('traits')}")
    ann = g.get("annotations") or {}
    print(f"annotations.tf_family    : {ann.get('tf_family')}")
    print(f"annotations.is_tf        : {ann.get('is_transcription_factor')}")
    print(f"external_links           : {g.get('external_links')}")
else:
    print("Aucun exemple trouvé (source=planttfdb+ncbi_linked) !")

# --- 2. Statistiques globales -------------------------------------------------
print("\n" + "=" * 70)
print("STATISTIQUES GLOBALES")
print("=" * 70)

seq_types = Counter(g.get("sequence_type") for g in genes)
print(f"Répartition sequence_type : {dict(seq_types)}")

lengths = [len(g.get("sequence") or "") for g in genes]
lengths_nonzero = [l for l in lengths if l > 0]
print(f"Longueur séquence : min={min(lengths_nonzero) if lengths_nonzero else 0} "
      f"max={max(lengths_nonzero) if lengths_nonzero else 0} "
      f"moyenne={sum(lengths_nonzero)/len(lengths_nonzero) if lengths_nonzero else 0:.0f}")
empty_seq = sum(1 for l in lengths if l == 0)
print(f"Séquences vides : {empty_seq} / {len(genes)}")

has_go = sum(1 for g in genes if (g.get("annotations") or {}).get("go_terms"))
has_tf = sum(1 for g in genes if (g.get("annotations") or {}).get("tf_family"))
has_traits = sum(1 for g in genes if g.get("traits"))
has_pathways = sum(1 for g in genes if g.get("pathways"))
print(f"\nCouverture des labels sur ce jeu de données :")
print(f"  go_terms  : {has_go}/{len(genes)} ({has_go/len(genes)*100:.1f}%)")
print(f"  tf_family : {has_tf}/{len(genes)} ({has_tf/len(genes)*100:.1f}%)")
print(f"  traits    : {has_traits}/{len(genes)} ({has_traits/len(genes)*100:.1f}%)")
print(f"  pathways  : {has_pathways}/{len(genes)} ({has_pathways/len(genes)*100:.1f}%)")

# --- 3. Doublons de séquence (attendu, sera géré par postgres_utils au chargement) --
hashes = Counter()
for g in genes:
    seq = str(g.get("sequence") or "").upper().strip()
    h = hashlib.sha256(seq.encode("utf-8")).hexdigest()
    hashes[h] += 1
dup_groups = sum(1 for c in hashes.values() if c > 1)
dup_total = sum(c - 1 for c in hashes.values() if c > 1)
print(f"\nDoublons de séquence exacts : {dup_groups} groupes, {dup_total} enregistrements redondants "
      f"(seront fusionnés automatiquement au chargement Postgres)")

# --- 4. Anomalies à surveiller ------------------------------------------------
print("\n" + "=" * 70)
print("ANOMALIES POTENTIELLES")
print("=" * 70)
no_organism = sum(1 for g in genes if not g.get("organism") or g.get("organism") == "Unknown")
print(f"Gènes sans organisme renseigné : {no_organism}")
no_symbol = sum(1 for g in genes if not g.get("symbol"))
print(f"Gènes sans symbol : {no_symbol}")
very_short = sum(1 for l in lengths_nonzero if l < 50)
print(f"Séquences très courtes (<50pb, suspectes) : {very_short}")
