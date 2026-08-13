import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "collect"), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import collect_planttfdb
import collect_ncbi

CANDIDATES = ["Glycine max", "Zea mays"]
SAMPLE_SIZE = 8

for species in CANDIDATES:
    print(f"=== {species} : test de resolution NCBI reelle sur {SAMPLE_SIZE} gene_id PlantTFDB ===")
    tf_recs = collect_planttfdb.fetch_planttfdb(species, retmax=SAMPLE_SIZE)

    ok = 0
    for r in tf_recs:
        gene_locus = r.get("gene_id")
        result = collect_ncbi.fetch_fasta_by_accession(
            gene_locus, db="nucleotide", plants_only=True, organism=species, max_length=50_000,
        )
        status = "OK" if result else "ECHEC"
        if result:
            ok += 1
        print(f"  {gene_locus:25} -> {status}")

    print(f"  Taux de reussite : {ok}/{len(tf_recs)} ({ok/len(tf_recs)*100 if tf_recs else 0:.0f}%)\n")
