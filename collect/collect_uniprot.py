#!/usr/bin/env python3
"""
UniProt collector — reviewed (Swiss-Prot) + unreviewed (TrEMBL) plant proteomes.
Returns a list of gene records compatible with the pipeline's canonical schema.
"""

from __future__ import annotations

import json
import os
import time
import requests
import request_utils as rq
from typing import Optional

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
KEGG_LINK_API = "https://rest.kegg.jp/link/pathway"

# GO aspect codes as embedded in UniProt's "GoTerm" property prefix
# (e.g. "C:chloroplast envelope" -> aspect "C"). Kept human-readable.
GO_ASPECT_LABELS = {
    "C": "cellular_component",
    "F": "molecular_function",
    "P": "biological_process",
}

# Map common plant names to UniProt taxonomy IDs
TAXON_MAP: dict[str, str] = {
    "arabidopsis thaliana": "3702",
    "oryza sativa": "39947",
    "zea mays": "4577",
    "triticum aestivum": "4565",
    "glycine max": "3847",
    "solanum lycopersicum": "4081",
    "solanum tuberosum": "4113",
    "vitis vinifera": "29760",
    "hordeum vulgare": "4513",
    "sorghum bicolor": "4558",
    "medicago sativa": "3879",
    "phaseolus vulgaris": "3885",
    "helianthus annuus": "4232",
    "daucus carota": "4039",
    "lactuca sativa": "4236",
    "allium cepa": "35883",
    "brassica oleracea": "3712",
    "cucumis sativus": "3659",
    "malus domestica": "3750",
    "prunus persica": "3760",
    "citrus sinensis": "2711",
    "fragaria ananassa": "3747",
    "olea europaea": "4146",
}


def resolve_kegg_pathways(
    kegg_gene_refs: list[str],
    delay: float = 0.4,
    cache_path: Optional[str] = None,
) -> dict[str, list[str]]:
    """
    Résout une liste de références de gène KEGG (ex: 'ath:AT1G67120')
    vers leurs vraies voies métaboliques (ex: 'path:ath04075').

    IMPORTANT: une référence croisée KEGG côté UniProt (xref_kegg) est un
    identifiant de GÈNE, pas un identifiant de voie métabolique -- il faut
    un second appel à l'API KEGG (endpoint link/pathway) pour obtenir les
    voies réelles associées à ce gène.

    Respecte la limite KEGG de 3 requêtes/seconde (usage académique).
    Batch de 10 identifiants par requête (max autorisé par l'API).
    Tolérant aux erreurs : un batch en échec est marqué comme résolu-vide
    plutôt que de faire planter toute la collecte.

    Si cache_path est fourni, les résolutions sont chargées/sauvegardées
    sur disque -- évite de re-interroger KEGG pour des références déjà
    résolues lors d'un run précédent (collecte interrompue, debug, ajout
    d'espèces supplémentaires).
    """
    resolved: dict[str, list[str]] = {}
    if cache_path and os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8") as f:
            resolved = json.load(f)

    # NOTE: "resolved" ne contient QUE des résultats confirmés par KEGG
    # (y compris une vraie liste vide = confirmé sans voie associée).
    # Un échec réseau va dans "failed_refs" et n'est jamais persisté dans
    # le cache -- sans cette distinction, une erreur transitoire lors d'un
    # run serait mémorisée comme "aucune voie" pour toujours, et ne serait
    # plus jamais réessayée lors des runs suivants. Critique à l'échelle
    # d'une collecte massive, où des erreurs réseau ponctuelles sont
    # statistiquement inévitables.
    unique_refs = sorted(set(kegg_gene_refs) - set(resolved.keys()))
    failed_refs: set[str] = set()

    for i in range(0, len(unique_refs), 10):
        chunk = unique_refs[i:i + 10]
        query = "+".join(chunk)
        url = f"{KEGG_LINK_API}/{query}"
        try:
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
            for ref in chunk:
                resolved.setdefault(ref, [])
            if resp.text.strip():
                for line in resp.text.strip().split("\n"):
                    gene_id, pathway_id = line.split("\t")
                    resolved.setdefault(gene_id, []).append(pathway_id)
        except (requests.RequestException, ValueError) as e:
            print(f"  [KEGG] Erreur résolution pathway pour {chunk}: {e} -- sera réessayé au prochain run")
            failed_refs.update(chunk)
        time.sleep(delay)

        if cache_path and (i // 10) % 20 == 0:
            with open(cache_path, "w", encoding="utf-8") as f:
                json.dump(resolved, f)

    if cache_path:
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(resolved, f)

    if failed_refs:
        print(f"  [KEGG] {len(failed_refs)} référence(s) en échec réseau -- "
              f"non mises en cache, seront retentées au prochain run")

    # Pour l'usage immédiat (le run en cours), on retourne quand même une
    # entrée (vide) pour les refs en échec -- mais uniquement en mémoire,
    # jamais persistée sur disque.
    result = dict(resolved)
    for ref in failed_refs:
        result.setdefault(ref, [])

    return result


def fetch_uniprot(
    species: str,
    retmax: int = 300,
    reviewed_only: bool = False,
    kegg_cache_path: Optional[str] = None,
) -> list[dict]:
    """
    Fetch protein records from UniProt for a given plant species.

    Args:
        species: Scientific name (e.g. "Arabidopsis thaliana")
        retmax: Maximum number of records to retrieve
        reviewed_only: If True, only fetch Swiss-Prot (reviewed) entries
        kegg_cache_path: Optional path to a JSON cache for KEGG pathway
            resolution, reused across species/runs to avoid redundant
            API calls.

    Returns:
        List of normalized gene records
    """
    taxon_id = TAXON_MAP.get(species.lower())
    if taxon_id:
        query = f"taxonomy_id:{taxon_id}"
    else:
        query = f'organism_name:"{species}"'

    if reviewed_only:
        query += " AND reviewed:true"

    params = {
        "query": query,
        "format": "json",
        "size": min(retmax, 500),
        "fields": (
            "accession,gene_names,organism_name,protein_name,"
            "sequence,go_id,go,ft_domain,cc_function,"
            "xref_kegg,xref_ensembl,xref_refseq,xref_embl,reviewed,"
            "protein_existence,annotation_score,length,"
            "cc_subcellular_location,keyword,feature_count"
        ),
    }

    records: list[dict] = []
    next_link: Optional[str] = None
    fetched = 0

    while fetched < retmax:
        try:
            if next_link:
                resp = rq.get(next_link, timeout=30)
            else:
                resp = rq.get(UNIPROT_API, params=params, timeout=30)
            data = resp.json()
            results = data.get("results", [])

            for entry in results:
                if fetched >= retmax:
                    break
                rec = _parse_entry(entry, species)
                if rec:
                    records.append(rec)
                    fetched += 1

            # Pagination via Link header
            link_header = resp.headers.get("Link", "")
            next_link = _parse_next_link(link_header)
            if not next_link or not results:
                break

            time.sleep(0.3)  # Be polite to UniProt API

        except (requests.RequestException, ValueError) as e:
            print(f"  [UniProt] Request error for {species}: {e}")
            break

    # Résolution des vraies voies métaboliques KEGG -- post-traitement
    # batché sur l'ensemble des enregistrements collectés, plus efficace
    # qu'un appel API par gène pendant la boucle de collecte ci-dessus.
    all_kegg_refs = [
        ref for rec in records
        for ref in rec["external_links"].get("kegg_gene_refs", [])
    ]
    if all_kegg_refs:
        pathway_map = resolve_kegg_pathways(all_kegg_refs, cache_path=kegg_cache_path)
        n_resolved = sum(1 for v in pathway_map.values() if v)
        print(f"  [KEGG] {n_resolved} / {len(pathway_map)} références de gène "
              f"résolues vers au moins une vraie voie métabolique")
        for rec in records:
            refs = rec["external_links"].get("kegg_gene_refs", [])
            resolved_pathways = sorted(set(
                pid for ref in refs for pid in pathway_map.get(ref, [])
            ))
            rec["pathways"] = [{"id": pid, "source": "kegg"} for pid in resolved_pathways]

    return records


def _parse_entry(entry: dict, species: str) -> dict | None:
    """Parse a single UniProt JSON entry into pipeline canonical format."""
    accession = entry.get("primaryAccession", "")
    if not accession:
        return None

    # Gene names
    gene_names = entry.get("genes", [])
    symbol = ""
    if gene_names:
        gn = gene_names[0]
        symbol = (
            gn.get("geneName", {}).get("value", "")
            or (gn.get("synonyms", [{}])[0].get("value", "") if gn.get("synonyms") else "")
            or accession
        )

    # Protein name.
    # NOTE: la version précédente utilisait un ternaire dont la précédence
    # d'opérateur ("A or B if C else D" == "(A or B) if C else D" en
    # Python) effaçait silencieusement recommendedName.fullName dès que
    # submittedName était absent -- le cas le plus courant pour une entrée
    # Swiss-Prot bien annotée. Réécrit explicitement pour éviter ce piège.
    pn = entry.get("proteinDescription", {})
    rec_name = pn.get("recommendedName", {})
    recommended_full_name = rec_name.get("fullName", {}).get("value", "")
    if recommended_full_name:
        description = recommended_full_name
    elif pn.get("submittedName"):
        description = pn["submittedName"][0].get("fullName", {}).get("value", "")
    else:
        description = ""

    # Sequence
    seq_obj = entry.get("sequence", {})
    seq = seq_obj.get("value", "")
    if not seq:
        return None

    # Organism
    organism = entry.get("organism", {}).get("scientificName", species)

    # GO terms.
    # NOTE: la version précédente stockait "aspect": props["GoEvidenceType"]
    # (ex: "IDA:TAIR") sous une clé nommée "aspect", alors que le VRAI
    # aspect GO (Cellular Component / Molecular Function / Biological
    # Process) est en fait le préfixe caché dans "GoTerm" (ex: "C:nucleus").
    # Corrigé ici : term est nettoyé du préfixe, aspect contient le vrai
    # C/F/P (libellé lisible), et le code d'évidence + sa source sont
    # exposés séparément plutôt que conflatés dans "aspect".
    go_terms = []
    uniProtKB_cross_refs = entry.get("uniProtKBCrossReferences", [])
    for ref in uniProtKB_cross_refs:
        if ref.get("database") == "GO":
            go_id = ref.get("id", "")
            props = {p["key"]: p["value"] for p in ref.get("properties", [])}
            raw_term = props.get("GoTerm", "")  # e.g. "C:chloroplast envelope"
            aspect_code, sep, clean_term = raw_term.partition(":")
            if not sep:
                # Pas de préfixe reconnu -- garde la chaîne brute telle quelle
                aspect_code, clean_term = "", raw_term
            evidence_raw = props.get("GoEvidenceType", "")  # e.g. "IDA:TAIR"
            evidence_code, _, evidence_source = evidence_raw.partition(":")
            go_terms.append({
                "id": go_id,
                "term": clean_term.strip() or raw_term,
                "aspect": GO_ASPECT_LABELS.get(aspect_code, aspect_code),
                "evidence_code": evidence_code,
                "evidence_source": evidence_source,
            })

    # KEGG cross-refs -- ce sont des identifiants de GÈNE (ex: "ath:AT1G67120"),
    # PAS des identifiants de voie métabolique. Résolus vers les vraies
    # voies après coup par resolve_kegg_pathways(), appelée depuis
    # fetch_uniprot() une fois toute la collecte terminée.
    kegg_ids = [
        ref.get("id") for ref in uniProtKB_cross_refs
        if ref.get("database") == "KEGG"
    ]

    # Ensembl cross-refs
    ensembl_ids = [
        ref.get("id") for ref in uniProtKB_cross_refs
        if ref.get("database") == "Ensembl"
    ]

    def _valid_refseq_nucleotide_id(nuc_id: str) -> bool:
        return nuc_id.startswith(("NM_", "XM_", "NR_", "XR_"))

    # RefSeq cross-refs -- THE join key back to NCBI nucleotide records.
    # UniProt's RefSeq cross-reference "id" is the PROTEIN accession
    # (NP_/XP_...), but it carries a "NucleotideSequenceId" property that
    # is usually the corresponding transcript accession (NM_/XM_/NR_/XR_).
    # Some plant entries still expose genomic genome references (NC_/CM_),
    # which are too large for our per-gene pipeline and must be ignored.
    refseq_protein_ids = []
    refseq_nucleotide_ids = []
    for ref in uniProtKB_cross_refs:
        if ref.get("database") == "RefSeq":
            if ref.get("id"):
                refseq_protein_ids.append(ref["id"])
            props = {p["key"]: p["value"] for p in ref.get("properties", [])}
            nuc_id = props.get("NucleotideSequenceId")
            if nuc_id and _valid_refseq_nucleotide_id(nuc_id):
                refseq_nucleotide_ids.append(nuc_id)

    # EMBL/GenBank cross-refs -- secondary fallback join key when a gene
    # has no RefSeq entry. Prefer transcript/RNA accessions and ignore
    # genomic mappings like Genomic_DNA, which often point to whole
    # chromosomes/contigs rather than a single-gene transcript.
    embl_candidates = []
    for ref in uniProtKB_cross_refs:
        if ref.get("database") == "EMBL" and ref.get("id"):
            props = {p["key"]: p["value"] for p in ref.get("properties", [])}
            molecule_type = props.get("MoleculeType", "")
            normalized = molecule_type.lower()
            if "genomic" in normalized and "rna" not in normalized:
                continue
            if "rna" in normalized or "transcript" in normalized:
                priority = 0
            else:
                priority = 1
            embl_candidates.append((ref["id"], molecule_type, priority))
    embl_candidates.sort(key=lambda c: c[2])
    embl_nucleotide_ids = [c[0] for c in embl_candidates]
    embl_molecule_types = [c[1] for c in embl_candidates]

    # Keywords -- garde maintenant la catégorie UniProt (ex: "Biological
    # process", "Domain", "PTM", "Molecular function") en plus du nom.
    # Utile pour un futur nettoyage de "traits", qui mélange actuellement
    # des voies de signalisation, propriétés structurelles et fonctions
    # dans la même liste plate (limite documentée séparément).
    keywords_detailed = [
        {"name": kw.get("name", ""), "category": kw.get("category", "")}
        for kw in entry.get("keywords", [])
    ]
    keywords = [kw["name"] for kw in keywords_detailed]

    # Subcellular location
    comments = entry.get("comments", [])
    subcell = []
    func_desc = ""
    for comment in comments:
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                loc_val = loc.get("location", {}).get("value", "")
                if loc_val:
                    subcell.append(loc_val)
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                func_desc = texts[0].get("value", "")

    reviewed = entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"

    return {
        "gene_id": accession,
        "symbol": symbol or accession,
        "organism": organism,
        "sequence": seq.upper(),
        "sequence_type": "protein",
        "description": description or func_desc,
        "length": len(seq),
        "source": "uniprot",
        "annotations": {
            "go_terms": go_terms,
            "keywords": keywords,
            "keywords_detailed": keywords_detailed,
            "subcellular_location": subcell,
            "function": func_desc,
            "annotation_score": entry.get("annotationScore", 0),
            "protein_existence": entry.get("proteinExistence", ""),
            "reviewed": reviewed,
        },
        "external_links": {
            "uniprot": f"https://www.uniprot.org/uniprotkb/{accession}",
            "accession": accession,
            "kegg": kegg_ids[0] if kegg_ids else None,
            # Toutes les références de gène KEGG (pas juste la première),
            # utilisées par resolve_kegg_pathways() dans fetch_uniprot()
            # pour retrouver les vraies voies métaboliques après coup.
            "kegg_gene_refs": kegg_ids,
            "ensembl": ensembl_ids[0] if ensembl_ids else None,
            # Clé de jointure vers NCBI -- consommée par le futur collecteur
            # NCBI "fetch-by-accession" pour récupérer la séquence ADN/ARN
            # DU MÊME GÈNE que cette entrée protéine annotée.
            "refseq_nucleotide": refseq_nucleotide_ids[0] if refseq_nucleotide_ids else None,
            "refseq_protein": refseq_protein_ids[0] if refseq_protein_ids else None,
            "embl_nucleotide": embl_nucleotide_ids[0] if embl_nucleotide_ids else None,
            "embl_molecule_type": embl_molecule_types[0] if embl_molecule_types else None,
        },
        "traits": keywords[:10],  # top keywords as traits
        "expression_profiles": [],
        # Rempli après coup par resolve_kegg_pathways(), appelée depuis
        # fetch_uniprot() une fois toute la collecte de l'espèce terminée.
        "pathways": [],
        "publications": [],
    }


def _parse_next_link(link_header: str) -> str | None:
    """Extract 'next' URL from Link header."""
    if not link_header:
        return None
    for part in link_header.split(","):
        part = part.strip()
        if 'rel="next"' in part:
            url = part.split(";")[0].strip().strip("<>")
            return url
    return None


if __name__ == "__main__":
    results = fetch_uniprot(
        "Arabidopsis thaliana",
        retmax=5,
        kegg_cache_path="kegg_pathway_cache.json",
    )
    print(json.dumps(results[:2], indent=2))
    print(f"Total fetched: {len(results)}")