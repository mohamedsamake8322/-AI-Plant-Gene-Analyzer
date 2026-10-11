#!/usr/bin/env python3
"""
diagnose_uniprot_ncbi_bridge.py — Compare les vraies valeurs des deux côtés
du pont UniProt<->NCBI pour comprendre pourquoi 0 fusion a eu lieu malgré
un test synthétique concluant.

Usage:
    python diagnose_uniprot_ncbi_bridge.py "data\\clean\\species\\chenopodium_quinoa_all_sources.json"
"""

import sys
import json


def flatten(obj, out):
    if isinstance(obj, dict):
        if "gene_id" in obj:
            out.append(obj)
        else:
            for v in obj.values():
                flatten(v, out)
    elif isinstance(obj, list):
        for item in obj:
            flatten(item, out)


def main():
    json_path = sys.argv[1]
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)
    records = []
    flatten(data, records)

    print(f"{len(records)} enregistrements au total\n")

    # 1) Un enregistrement NCBI typique -- on veut voir TOUTES ses clés
    #    top-level pour trouver où l'accession vit réellement.
    ncbi_like = [r for r in records if str(r.get("gene_id", "")).startswith("GeneID:")]
    print(f"Enregistrements 'GeneID:' trouvés : {len(ncbi_like)}")
    if ncbi_like:
        sample = ncbi_like[0]
        print("Clés top-level d'un enregistrement NCBI fusionné :", list(sample.keys()))
        print("Exemple complet :")
        print(json.dumps(sample, indent=2, ensure_ascii=False)[:1500])
    print()

    # 2) Un enregistrement UniProt typique -- on veut voir external_links
    #    tel qu'il existe VRAIMENT (nom des clés, format des valeurs).
    uniprot_like = [
        r for r in records
        if r.get("external_links", {}).get("refseq_nucleotide")
        or r.get("external_links", {}).get("refseq_protein")
    ]
    print(f"Enregistrements UniProt avec un refseq_nucleotide/protein renseigné : {len(uniprot_like)}")
    if uniprot_like:
        sample = uniprot_like[0]
        print("gene_id UniProt :", sample.get("gene_id"))
        print("external_links  :", sample.get("external_links"))
    print()

    # 3) Comparaison directe : est-ce que la valeur EXACTE de
    #    refseq_nucleotide/protein correspond à un "accession" existant
    #    quelque part dans les vrais enregistrements (peu importe la clé) ?
    all_accessions_anywhere = set()
    for r in records:
        for key in ("accession", "ncbi_accession", "source_accession"):
            v = r.get(key)
            if v:
                all_accessions_anywhere.add(v)
        ext = r.get("external_links", {}) or {}
        for key in ("accession", "ncbi_accession"):
            v = ext.get(key)
            if v:
                all_accessions_anywhere.add(v)

    print(f"Valeurs d'accession trouvées (toutes clés confondues) : {len(all_accessions_anywhere)}")
    print("Exemples :", list(all_accessions_anywhere)[:5])

    if uniprot_like:
        for r in uniprot_like[:5]:
            ext = r.get("external_links", {})
            for candidate in (ext.get("refseq_nucleotide"), ext.get("refseq_protein")):
                if candidate:
                    exact_match = candidate in all_accessions_anywhere
                    versionless = candidate.split(".")[0]
                    versionless_match = any(
                        a.split(".")[0] == versionless for a in all_accessions_anywhere
                    )
                    print(f"  UniProt dit '{candidate}' -> match exact: {exact_match} | "
                          f"match sans version: {versionless_match}")


if __name__ == "__main__":
    main()
