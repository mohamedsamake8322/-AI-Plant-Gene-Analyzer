#!/usr/bin/env python3
"""
reconcile_uniprot_ncbi_kegg.py — Réconciliation post-hoc, SANS re-télécharger
NCBI ni KEGG (déjà sur disque, corrects). Ne re-télécharge QUE les
external_links UniProt (rapide), qui avaient été perdus par
restructure_to_schema() avant qu'on ne les utilise pour fusionner.

Principe :
  1. Charge le fichier espèce déjà collecté (toutes les séquences/annotations
     sont déjà bonnes -- rien à refaire côté NCBI/KEGG/PLAZA).
  2. Repère les enregistrements NCBI dont gene_id EST déjà l'accession brute
     (la majorité -- seuls ceux fusionnés via GeneID ont perdu ce lien, et
     ceux-là on ne peut de toute façon plus les récupérer sans re-fetch NCBI).
  3. Repère les enregistrements KEGG (gene_id du type "cqi:xxxxx").
  4. Re-télécharge UNIQUEMENT les external_links des enregistrements UniProt
     déjà présents (requête légère, par lot, via l'API REST UniProt).
  5. Fusionne chaque UniProt dans l'enregistrement NCBI/KEGG correspondant
     s'il y a match, supprime le doublon UniProt autonome.
  6. Écrit le fichier réconcilié (l'original n'est jamais modifié).

Usage:
    python reconcile_uniprot_ncbi_kegg.py "data\\clean\\species\\chenopodium_quinoa_all_sources.json"
"""

import sys
import json
import re
import time
import argparse
from pathlib import Path

sys.path.insert(0, "collect")
sys.path.insert(0, "scripts")

import requests
import collect_uniprot as cu

NCBI_ACCESSION_PATTERNS = [
    re.compile(r"^[NX]M_\d+\.\d+$"),
    re.compile(r"^[NXY]P_\d+\.\d+$"),
    re.compile(r"^[A-Z]{2,4}\d{4,8}\.\d+$"),
]
KEGG_PATTERN = re.compile(r"^cqi:\d+$")
UNIPROT_PATTERN = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$")

UNIPROT_BATCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def flatten(obj, out, path=None):
    """Aplatit en gardant, pour chaque record, la liste des clés parentes
    (nécessaire pour pouvoir le supprimer/remplacer au bon endroit ensuite)."""
    if isinstance(obj, dict):
        if "gene_id" in obj:
            out.append((path or [], obj))
        else:
            for k, v in obj.items():
                flatten(v, out, (path or []) + [k])
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            flatten(item, out, (path or []) + [i])


def fetch_uniprot_external_links(accessions: list[str], species: str) -> dict[str, dict]:
    """Re-télécharge UNIQUEMENT les external_links pour une liste précise
    d'accessions UniProt déjà connues -- pas de recherche par espèce, pas de
    nouvelle collecte de séquence. Retourne {accession: external_links}."""
    result: dict[str, dict] = {}
    batch_size = 100
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        query = "(" + " OR ".join(f"accession:{a}" for a in batch) + ")"
        params = {"query": query, "format": "json", "size": batch_size}
        try:
            resp = requests.get(UNIPROT_BATCH_URL, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            print(f"  ⚠ Lot {i//batch_size + 1} : échec ({e}), ignoré")
            continue
        for entry in data.get("results", []):
            parsed = cu._parse_entry(entry, species)
            if parsed and parsed.get("gene_id"):
                result[parsed["gene_id"]] = parsed.get("external_links", {})
        done = min(i + batch_size, len(accessions))
        print(f"  → {done}/{len(accessions)} accessions UniProt re-vérifiées...")
        time.sleep(0.3)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_path")
    parser.add_argument("--out", default=None, help="Fichier de sortie (défaut: <input>.reconciled.json)")
    parser.add_argument("--species", default="Chenopodium quinoa")
    args = parser.parse_args()

    in_path = Path(args.json_path)
    out_path = Path(args.out) if args.out else in_path.with_suffix(".reconciled.json")

    with open(in_path, encoding="utf-8") as f:
        data = json.load(f)

    flat_records = []
    flatten(data, flat_records)
    print(f"{len(flat_records)} enregistrements chargés depuis {in_path}\n")

    # Index construits SANS aucun réseau, directement depuis le fichier.
    ncbi_index: dict[str, dict] = {}
    kegg_index: dict[str, dict] = {}
    uniprot_records: list[dict] = []

    for path, rec in flat_records:
        gid = str(rec.get("gene_id", ""))
        if any(p.match(gid) for p in NCBI_ACCESSION_PATTERNS):
            ncbi_index[gid] = rec
        elif KEGG_PATTERN.match(gid):
            kegg_index[gid] = rec
        elif UNIPROT_PATTERN.match(gid):
            uniprot_records.append(rec)

    print(f"NCBI (accession brute, réconciliables sans re-fetch) : {len(ncbi_index)}")
    print(f"KEGG (gene_id natif)                                  : {len(kegg_index)}")
    print(f"UniProt (à re-vérifier)                                : {len(uniprot_records)}\n")

    if not uniprot_records:
        print("Rien à réconcilier (aucun enregistrement UniProt trouvé). Arrêt.")
        return

    accessions = [r["gene_id"] for r in uniprot_records]
    print("Re-téléchargement des external_links UniProt (léger, pas de séquence)...")
    ext_links_by_accession = fetch_uniprot_external_links(accessions, args.species)

    merged_into_ncbi = 0
    merged_into_kegg = 0
    to_remove_paths = []

    for rec_path, rec in [(p, r) for p, r in flat_records if r in uniprot_records]:
        gid = rec["gene_id"]
        ext = ext_links_by_accession.get(gid, {})
        target = None

        for candidate in (ext.get("refseq_nucleotide"), ext.get("refseq_protein")):
            if candidate and candidate in ncbi_index:
                target = ncbi_index[candidate]
                merged_into_ncbi += 1
                break

        if target is None:
            for kegg_ref in ext.get("kegg_gene_refs") or []:
                if kegg_ref in kegg_index:
                    target = kegg_index[kegg_ref]
                    merged_into_kegg += 1
                    break

        if target is not None:
            # Fusionne la séquence protéine + annotations dans l'enregistrement
            # cible, sans écraser ce qui existe déjà.
            if rec.get("sequence", {}).get("protein") and not target.get("sequence", {}).get("protein"):
                target.setdefault("sequence", {})["protein"] = rec["sequence"]["protein"]
            for go in rec.get("annotation", {}).get("go_terms", []):
                target.setdefault("annotation", {}).setdefault("go_terms", []).append(go)
            target.setdefault("sources_summary", [])
            if "uniprot" not in target["sources_summary"]:
                target["sources_summary"] = sorted(set(target["sources_summary"]) | {"uniprot"})
            target.setdefault("uniprot_accession", gid)
            to_remove_paths.append(rec_path)

    print(f"\nFusionné dans NCBI : {merged_into_ncbi}")
    print(f"Fusionné dans KEGG  : {merged_into_kegg}")
    print(f"Total réconcilié    : {merged_into_ncbi + merged_into_kegg} / {len(uniprot_records)}")

    # Retire les enregistrements UniProt maintenant fusionnés ailleurs, pour
    # ne pas les garder en double dans le fichier de sortie.
    removed = 0
    for rec_path, rec in flat_records:
        if rec_path in to_remove_paths:
            continue
    # Suppression réelle : reconstruit récursivement en filtrant les records
    # dont le chemin est dans to_remove_paths.
    def prune(obj, path=None):
        nonlocal removed
        path = path or []
        if isinstance(obj, dict):
            if "gene_id" in obj:
                return None if path in to_remove_paths else obj
            out = {}
            for k, v in obj.items():
                pruned = prune(v, path + [k])
                if pruned is None and isinstance(v, dict) and "gene_id" in v:
                    removed += 1
                    continue
                out[k] = pruned
            return out
        if isinstance(obj, list):
            out = []
            for i, item in enumerate(obj):
                pruned = prune(item, path + [i])
                if pruned is None and isinstance(item, dict) and "gene_id" in item:
                    removed += 1
                    continue
                out.append(pruned)
            return out
        return obj

    reconciled = prune(data)
    print(f"\nEnregistrements UniProt autonomes supprimés (fusionnés ailleurs) : {removed}")

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(reconciled, f, ensure_ascii=False, indent=2)
    print(f"\n✓ Fichier réconcilié écrit : {out_path}")
    print("  (le fichier original n'a pas été modifié)")


if __name__ == "__main__":
    main()
