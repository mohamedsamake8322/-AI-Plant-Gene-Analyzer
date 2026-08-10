import json
from collections import Counter

d = json.load(open("Data/clean/master_plant_db.json", encoding="utf-8"))
genes = d["genes"]

VALID = {
    "dna": set("ACGTRYSWKMBDHVN"),
    "rna": set("ACGURYSWKMBDHVN"),
}

offending_chars = Counter()
examples = {}
mismatch_type_but_valid_combined = 0
n_checked = 0

for g in genes:
    seq_type = str(g.get("sequence_type") or "").lower()
    if seq_type not in ("dna", "rna"):
        continue
    seq = str(g.get("sequence") or "").upper().strip()
    if not seq:
        continue
    n_checked += 1

    invalid = set(seq) - VALID[seq_type]
    if invalid:
        for c in invalid:
            offending_chars[c] += 1
            if c not in examples:
                examples[c] = {"gene_id": g.get("gene_id"), "source": g.get("source"),
                                "sequence_type": seq_type, "snippet": seq[:60]}
        # est-ce que ça passerait si on acceptait T et U indifféremment des deux côtés ?
        combined = VALID["dna"] | VALID["rna"]
        if not (set(seq) - combined):
            mismatch_type_but_valid_combined += 1

print(f"Séquences dna/rna vérifiées : {n_checked}")
print(f"\nCaractères problématiques rencontrés (nombre de séquences affectées) :")
for c, count in offending_chars.most_common(20):
    print(f"  '{c}' : {count} séquences")

print(f"\nSéquences invalides pour LEUR type déclaré mais valides si dna/rna étaient fusionnés (probable erreur de labellisation T/U) : {mismatch_type_but_valid_combined}")

print("\nExemples (un par caractère problématique) :")
for c, info in list(examples.items())[:10]:
    print(f"  '{c}' -> gene_id={info['gene_id']} source={info['source']} type={info['sequence_type']}")
    print(f"        extrait: {info['snippet']}")
