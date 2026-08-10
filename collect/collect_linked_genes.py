#!/usr/bin/env python3
"""
collect_linked_genes.py
-------------------------
Collecte "par gène", pas "par source". Corrige le défaut de conception
identifié : collecter l'ADN/ARN (NCBI) et les annotations (UniProt) de
façon indépendante produit deux ensembles de gènes disjoints (0% de
recouvrement mesuré sur la base actuelle). Ici, pour CHAQUE gène, on va
chercher sa séquence ET son annotation ENSEMBLE, et on ne le garde que si
les deux sont trouvées.

Logique :
  1. Interroger UniProt pour les protéines déjà bien annotées d'une espèce
     (go_terms non vide -- c'est notre signal de qualité).
  2. Pour chacune, lire sa référence croisée RefSeq/EMBL (ajoutée dans
     collect_uniprot.py -- external_links.refseq_nucleotide /
     embl_nucleotide) : c'est l'accession NCBI de LA séquence nucléotidique
     de CE MÊME gène.
  3. Aller chercher cette séquence précise chez NCBI via
     collect_ncbi.fetch_fasta_by_accession (pas une recherche large par
     espèce -- une recherche ciblée par accession).
  4. Ne garder le gène que si l'étape 3 réussit. Sinon, il est compté
     comme "sans référence croisée exploitable" et exclu -- pas de
     séquence orpheline, pas de label orphelin.

Résultat : un jeu de données plus petit, mais où CHAQUE gène a une
séquence nucléotidique ET une annotation fonctionnelle réelle -- ce que
la base actuelle n'a pas.

Usage :
    python collect/collect_linked_genes.py --species "Arabidopsis thaliana" --retmax 300
    python collect/collect_linked_genes.py --species "Arabidopsis thaliana" "Oryza sativa" "Solanum lycopersicum" --retmax 500
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "collect"), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import collect_uniprot  # noqa: E402
import collect_ncbi  # noqa: E402
import collect_planttfdb  # noqa: E402


def _sequence_type_from_accession(accession: str, embl_molecule_type: str | None = None) -> str:
    """RefSeq mRNA accessions start with NM_/XM_ ; genomic/other nucleotide
    accessions default to dna. EMBL/GenBank/DDBJ accessions (e.g. from the
    embl_nucleotide fallback) don't follow that NCBI prefix convention at
    all, so without embl_molecule_type they'd silently default to "dna"
    even when UniProt already told us it's an mRNA -- use that hint when
    the accession-prefix heuristic doesn't apply.
    """
    prefix = accession.split("_")[0].upper() if "_" in accession else ""
    if prefix in ("NM", "XM"):
        return "rna"
    if embl_molecule_type and "rna" in embl_molecule_type.lower():
        return "rna"
    return "dna"


def build_linked_gene(uniprot_rec: dict, organism: str, api_delay: float, max_gene_length: int) -> tuple[dict | None, str]:
    """Try to attach a real nucleotide sequence to one already-annotated
    UniProt record. Returns (merged_record_or_None, reason).
    """
    links = uniprot_rec.get("external_links") or {}
    refseq_accession = links.get("refseq_nucleotide")
    nuc_accession = refseq_accession or links.get("embl_nucleotide")
    if not nuc_accession:
        return None, "no_nucleotide_xref"
    # Only meaningful when we actually fell back to the EMBL xref -- RefSeq
    # accessions are already classified correctly by their NM_/XM_ prefix.
    embl_molecule_type = None if refseq_accession else links.get("embl_molecule_type")

    time.sleep(api_delay)
    ncbi_records = collect_ncbi.fetch_fasta_by_accession(
        nuc_accession, db="nucleotide", plants_only=True, organism=organism, max_length=max_gene_length,
    )
    if not ncbi_records:
        print(
            f"    [DEBUG] NCBI fetch failed for UniProt {uniprot_rec.get('gene_id')} -> "
            f"nucleotide accession {nuc_accession} (organism={organism})"
        )
        return None, "ncbi_fetch_failed"

    header, seq = ncbi_records[0]
    if not seq:
        return None, "empty_ncbi_sequence"

    merged = {
        "gene_id": nuc_accession,
        "symbol": uniprot_rec.get("symbol"),
        "organism": organism,
        "sequence": seq.upper().replace(" ", ""),
        "sequence_type": _sequence_type_from_accession(nuc_accession, embl_molecule_type),
        "description": uniprot_rec.get("description"),
        "traits": uniprot_rec.get("traits") or [],
        "pathways": uniprot_rec.get("pathways") or [],
        "publications": [],
        "annotations": uniprot_rec.get("annotations") or {},
        "external_links": {
            **links,
            "ncbi_nuccore": f"https://www.ncbi.nlm.nih.gov/nuccore/{nuc_accession.split('.')[0]}",
            "linked_from": "uniprot_go_terms",
        },
        "expression_profiles": [],
        "protein_sequence": uniprot_rec.get("sequence"),  # conservee pour d'autres usages de la plateforme
        "protein_accession": uniprot_rec.get("gene_id"),
        "source": "uniprot+ncbi_linked",
        "date_added": datetime.now(timezone.utc).isoformat(),
    }
    return merged, "ok"


def build_linked_gene_from_tf(tf_rec: dict, organism: str, api_delay: float, max_gene_length: int) -> tuple[dict | None, str]:
    """Meme logique que build_linked_gene, mais pour un enregistrement
    PlantTFDB. Pas de refseq_nucleotide ici (PlantTFDB n'a pas cette
    reference croisee) -- on tente de resoudre gene_id directement, car
    c'est souvent deja un identifiant de locus (ex: AT1G01010, ou
    Solyc04g007000.1 pour la tomate) que resolve_accession_id() sait
    deja traiter comme tel (voir le pattern generalise
    [A-Za-z]{2,}\\d+[A-Za-z]\\d+(\\.\\d+)? dans collect_ncbi.py).
    """
    gene_locus = tf_rec.get("gene_id")
    if not gene_locus:
        return None, "no_gene_locus"

    time.sleep(api_delay)
    ncbi_records = collect_ncbi.fetch_fasta_by_accession(
        gene_locus, db="nucleotide", plants_only=True, organism=organism, max_length=max_gene_length,
    )
    if not ncbi_records:
        print(
            f"    [DEBUG] NCBI fetch failed for PlantTFDB {gene_locus} -> "
            f"could not resolve nucleotide sequence for {organism}"
        )
        return None, "ncbi_fetch_failed"

    header, seq = ncbi_records[0]
    if not seq:
        return None, "empty_ncbi_sequence"

    accession = header.split()[0]
    merged = {
        "gene_id": accession,
        "symbol": tf_rec.get("symbol"),
        "organism": organism,
        "sequence": seq.upper().replace(" ", ""),
        "sequence_type": _sequence_type_from_accession(accession),
        "description": tf_rec.get("description"),
        "traits": tf_rec.get("traits") or [],
        "pathways": tf_rec.get("pathways") or [],
        "publications": tf_rec.get("publications") or [],
        "annotations": tf_rec.get("annotations") or {},  # contient deja tf_family
        "external_links": {
            **(tf_rec.get("external_links") or {}),
            "ncbi_nuccore": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession.split('.')[0]}",
            "linked_from": "planttfdb_gene_locus",
        },
        "expression_profiles": [],
        "source": "planttfdb+ncbi_linked",
        "date_added": datetime.now(timezone.utc).isoformat(),
    }
    return merged, "ok"


def collect_species(species: str, retmax: int, api_delay: float, max_gene_length: int, include_tf: bool = True) -> tuple[list[dict], dict]:
    print(f"\n=== {species} ===")
    print(f"  Interrogation UniProt (proteines annotees, cible {retmax})...")
    uniprot_recs = collect_uniprot.fetch_uniprot(species, retmax=retmax, reviewed_only=False)

    annotated = [r for r in uniprot_recs if (r.get("annotations") or {}).get("go_terms")]
    print(f"  {len(uniprot_recs)} proteines recuperees, {len(annotated)} avec go_terms non vide")

    linked = []
    reasons = {"ok": 0, "no_nucleotide_xref": 0, "ncbi_fetch_failed": 0, "empty_ncbi_sequence": 0}

    for i, rec in enumerate(annotated, 1):
        merged, reason = build_linked_gene(rec, species, api_delay, max_gene_length)
        reasons[reason] += 1
        if merged:
            linked.append(merged)
        if i % 50 == 0:
            print(f"    [{i}/{len(annotated)}] traites -- {reasons['ok']} lies avec succes jusqu'ici")

    stats = {
        "species": species,
        "uniprot_fetched": len(uniprot_recs),
        "uniprot_with_go_terms": len(annotated),
        **reasons,
        "final_linked_genes_uniprot": len(linked),
    }

    if include_tf:
        print(f"  Interrogation PlantTFDB (facteurs de transcription, cible {retmax})...")
        tf_recs = collect_planttfdb.fetch_planttfdb(species, retmax=retmax)
        print(f"  {len(tf_recs)} facteurs de transcription recuperes")

        tf_reasons = {"ok": 0, "no_gene_locus": 0, "ncbi_fetch_failed": 0, "empty_ncbi_sequence": 0}
        tf_linked = []
        for i, tf_rec in enumerate(tf_recs, 1):
            merged, reason = build_linked_gene_from_tf(tf_rec, species, api_delay, max_gene_length)
            tf_reasons[reason] += 1
            if merged:
                tf_linked.append(merged)
            if i % 50 == 0:
                print(f"    [{i}/{len(tf_recs)}] TF traites -- {tf_reasons['ok']} lies avec succes jusqu'ici")

        linked.extend(tf_linked)
        stats["planttfdb_fetched"] = len(tf_recs)
        stats["planttfdb_tf_reasons"] = tf_reasons
        stats["final_linked_genes_planttfdb"] = len(tf_linked)

    stats["final_linked_genes"] = len(linked)
    print(f"  -> {stats}")
    return linked, stats


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", nargs="+", required=True,
                    help='Une ou plusieurs especes, ex: --species "Arabidopsis thaliana" "Oryza sativa"')
    p.add_argument("--retmax", type=int, default=300, help="Proteines UniProt a examiner par espece (defaut: 300)")
    p.add_argument("--api-delay", type=float, default=0.4,
                    help="Pause entre chaque appel NCBI, secondes (defaut: 0.4 -- respecte la limite NCBI "
                         "de 3 requetes/s sans cle API ; reduire si NCBI_API_KEY est configure dans .env)")
    p.add_argument("--no-tf", action="store_true", help="Desactiver la collecte PlantTFDB (UniProt seul)")
    p.add_argument("--max-gene-length", type=int, default=50_000,
                    help="Longueur max acceptee pour une sequence 'gene' (defaut: 50000pb -- resserre le "
                         "defaut de collect_ncbi.py, 500000, qui laisse passer les genomes d'organites "
                         "mitochondriaux/chloroplastiques mal resolus en genes individuels, ex: NC_037304.1 "
                         "= genome mitochondrial complet d'Arabidopsis, 367808pb, retourne pour le gene NAD7)")
    p.add_argument("--out", default=str(ROOT / "Data" / "clean" / "linked_genes.json"))
    args = p.parse_args()

    all_linked = []
    all_stats = []
    for species in args.species:
        linked, stats = collect_species(species, args.retmax, args.api_delay, args.max_gene_length, include_tf=not args.no_tf)
        all_linked.extend(linked)
        all_stats.append(stats)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps({
            "metadata": {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "strategy": "per-gene linked collection (uniprot annotation -> ncbi sequence by accession)",
                "species_stats": all_stats,
                "total_genes": len(all_linked),
            },
            "genes": all_linked,
        }, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(f"\n{'=' * 60}")
    print(f"TERMINE : {len(all_linked)} genes complets (sequence + annotation) -> {out_path}")
    for s in all_stats:
        total_candidates = s['uniprot_with_go_terms'] + s.get('planttfdb_fetched', 0)
        print(f"  {s['species']:30} {s['final_linked_genes']:5} genes lies "
              f"(sur {total_candidates} candidats : {s['uniprot_with_go_terms']} UniProt + "
              f"{s.get('planttfdb_fetched', 0)} PlantTFDB)")


if __name__ == "__main__":
    main()