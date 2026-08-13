import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "collect"), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import collect_uniprot
import collect_planttfdb

# Especes presentes a la fois dans TAXON_MAP (uniprot) et SPECIES_MAP (planttfdb),
# pas encore testees (Arabidopsis/Oryza/Solanum lycopersicum deja faites).
CANDIDATES = [
    "Glycine max",
    "Zea mays",
    "Vitis vinifera",
    "Medicago sativa",
    "Hordeum vulgare",
    "Sorghum bicolor",
]

SAMPLE_SIZE = 20

print(f"Sonde rapide sur {len(CANDIDATES)} especes candidates (echantillon de {SAMPLE_SIZE} chacune)\n")

results = []
for species in CANDIDATES:
    print(f"=== {species} ===")

    # Partie UniProt : combien ont une reference croisee exploitable ?
    try:
        uniprot_recs = collect_uniprot.fetch_uniprot(species, retmax=SAMPLE_SIZE, reviewed_only=False)
        with_go = [r for r in uniprot_recs if (r.get("annotations") or {}).get("go_terms")]
        with_xref = [r for r in with_go
                     if (r.get("external_links") or {}).get("refseq_nucleotide")
                     or (r.get("external_links") or {}).get("embl_nucleotide")]
        uniprot_rate = len(with_xref) / len(with_go) * 100 if with_go else 0
        print(f"  UniProt : {len(uniprot_recs)} recuperes, {len(with_go)} avec go_terms, "
              f"{len(with_xref)} avec reference croisee NCBI ({uniprot_rate:.0f}%)")
    except Exception as e:
        uniprot_rate = 0
        print(f"  UniProt : ECHEC ({e})")

    # Partie PlantTFDB : le telechargement fonctionne-t-il, et a quoi ressemblent les gene_id ?
    try:
        tf_recs = collect_planttfdb.fetch_planttfdb(species, retmax=SAMPLE_SIZE)
        sample_ids = [r.get("gene_id") for r in tf_recs[:5]]
        print(f"  PlantTFDB : {len(tf_recs)} recuperes. Exemples de gene_id : {sample_ids}")
    except Exception as e:
        tf_recs = []
        print(f"  PlantTFDB : ECHEC ({e})")

    results.append({
        "species": species,
        "uniprot_xref_rate": uniprot_rate,
        "planttfdb_count": len(tf_recs),
    })
    print()

print("=" * 70)
print("RESUME (trie par taux de reference croisee UniProt, decroissant)")
print("=" * 70)
for r in sorted(results, key=lambda x: -x["uniprot_xref_rate"]):
    verdict = "BON CANDIDAT" if r["uniprot_xref_rate"] > 70 and r["planttfdb_count"] > 0 else "risque"
    print(f"  {r['species']:20} uniprot_xref={r['uniprot_xref_rate']:5.0f}%  "
          f"planttfdb_recuperes={r['planttfdb_count']:3}  -> {verdict}")
