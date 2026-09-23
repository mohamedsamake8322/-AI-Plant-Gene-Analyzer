#!/usr/bin/env python3
"""
analyze_collection.py
----------------------
Analyse en LECTURE SEULE d'un ou plusieurs fichiers *_all_sources.json
produits par collect_all_sources.py. Ne modifie rien, n'affiche aucune
sequence complete (seulement des longueurs et des comptes).

Usage :
    python analyze_collection.py oryza_sativa_all_sources.json chenopodium_quinoa_all_sources.json
    python analyze_collection.py oryza_sativa_all_sources.json --out rapport_oryza.txt

Si --out n'est pas fourni, le rapport est ecrit dans <nom_du_fichier>.report.txt
a cote du fichier source, et aussi affiche a l'ecran.
"""
from __future__ import annotations

import argparse
import collections
import json
import statistics
import sys
from pathlib import Path


def pct(n, total):
    return f"{n} ({100 * n / total:.1f}%)" if total else f"{n} (0%)"


def analyze_file(path: Path, out) -> None:
    def w(line=""):
        print(line, file=out)

    w(f"{'=' * 70}")
    w(f"FICHIER : {path.name}")
    w(f"{'=' * 70}")

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        w(f"ERREUR : impossible de lire/parser le fichier : {e}")
        return

    meta = data.get("metadata", {})
    genes = data.get("genes", [])
    n = len(genes)

    w("\n--- Metadonnees du fichier ---")
    for k in ("plant", "common_name", "category", "generated_at", "count"):
        if k in meta:
            w(f"  {k}: {meta[k]}")
    if "sources" in meta:
        w(f"  sources demandees: {meta['sources']}")
    if "source_counts" in meta:
        w("  source_counts (bruts, avant fusion finale):")
        for k, v in meta["source_counts"].items():
            w(f"    {k}: {v}")
    if meta.get("errors"):
        w(f"  ERREURS DE COLLECTE ({len(meta['errors'])}):")
        for e in meta["errors"]:
            w(f"    - {e}")

    w(f"\n--- Volume ---\n  Genes dans le fichier: {n}")
    if n == 0:
        w("  Fichier vide, arret de l'analyse.")
        return

    # ------------------------------------------------------------------
    # 1. origin (sequence_backed vs plaza_only)
    # ------------------------------------------------------------------
    w("\n--- 1. Origine des enregistrements ---")
    origin_counts = collections.Counter(g.get("origin", "(absent)") for g in genes)
    for k, v in origin_counts.most_common():
        w(f"  {k}: {pct(v, n)}")

    seq_backed = [g for g in genes if g.get("origin") == "sequence_backed"]

    # ------------------------------------------------------------------
    # 2. Sequences presentes par type (dna/rna/protein), sur l'ensemble
    #    ET seulement parmi les sequence_backed
    # ------------------------------------------------------------------
    w("\n--- 2. Sequences presentes par type ---")
    for label, pool in (("Tous les genes", genes), ("Genes sequence_backed", seq_backed)):
        w(f"  [{label}, n={len(pool)}]")
        for st in ("dna", "rna", "protein"):
            has = sum(1 for g in pool if (g.get("sequence") or {}).get(st))
            w(f"    a une sequence {st}: {pct(has, len(pool)) if pool else '0 (0%)'}")
        both_dna_prot = sum(
            1 for g in pool
            if (g.get("sequence") or {}).get("dna") and (g.get("sequence") or {}).get("protein")
        )
        w(f"    a la fois dna ET protein: {pct(both_dna_prot, len(pool)) if pool else '0 (0%)'}")
        no_seq_at_all = sum(
            1 for g in pool
            if not any((g.get("sequence") or {}).get(st) for st in ("dna", "rna", "protein"))
        )
        w(f"    AUCUNE sequence (dna/rna/protein toutes vides): {pct(no_seq_at_all, len(pool)) if pool else '0 (0%)'}")

    # ------------------------------------------------------------------
    # 3. external_links -- LE point critique identifie precedemment
    # ------------------------------------------------------------------
    w("\n--- 3. external_links (cle de correction n0 1) ---")
    has_ext = sum(1 for g in genes if g.get("external_links"))
    w(f"  Genes avec external_links non vide: {pct(has_ext, n)}")
    if has_ext == 0:
        w("  >>> external_links est VIDE PARTOUT : confirme le bug (non recopie")
        w("      dans restructure_to_schema). A corriger AVANT de relancer une collecte,")
        w("      sinon les identifiants refseq_nucleotide/embl_nucleotide seront")
        w("      calcules puis perdus, comme dans la base actuelle.")
    else:
        keys_seen = collections.Counter()
        for g in genes:
            ext = g.get("external_links") or {}
            for k, v in ext.items():
                if v:
                    keys_seen[k] += 1
        w("  Cles presentes dans external_links (genes concernes):")
        for k, v in keys_seen.most_common():
            w(f"    {k}: {pct(v, n)}")

    # ------------------------------------------------------------------
    # 4. GO terms : couverture, par presence de sequence dna/protein
    # ------------------------------------------------------------------
    w("\n--- 4. Couverture GO terms ---")
    def go_terms_of(g):
        return ((g.get("annotation") or {}).get("go_terms")) or []

    has_go = [g for g in genes if go_terms_of(g)]
    w(f"  Genes avec au moins un GO term: {pct(len(has_go), n)}")

    has_go_and_dna = sum(1 for g in has_go if (g.get("sequence") or {}).get("dna"))
    has_go_and_protein = sum(1 for g in has_go if (g.get("sequence") or {}).get("protein"))
    has_go_no_seq = sum(
        1 for g in has_go
        if not any((g.get("sequence") or {}).get(st) for st in ("dna", "rna", "protein"))
    )
    w(f"    parmi eux, avec sequence dna:      {pct(has_go_and_dna, len(has_go)) if has_go else '0 (0%)'}")
    w(f"    parmi eux, avec sequence protein:  {pct(has_go_and_protein, len(has_go)) if has_go else '0 (0%)'}")
    w(f"    parmi eux, SANS AUCUNE sequence:   {pct(has_go_no_seq, len(has_go)) if has_go else '0 (0%)'}")
    if has_go_no_seq:
        w("    >>> Ces genes ont un GO term mais aucune sequence exploitable pour")
        w("        un modele (AgroNT ou proteique) : inutilisables tels quels pour")
        w("        l'entrainement, seulement pour la similarite/l'annotation.")

    # Aspect + evidence code
    aspect_counter = collections.Counter()
    evidence_counter = collections.Counter()
    evidence_source_counter = collections.Counter()
    for g in genes:
        for go in go_terms_of(g):
            aspect_counter[go.get("aspect") or "(absent)"] += 1
            evidence_counter[go.get("evidence_code") or "(absent)"] += 1
            evidence_source_counter[go.get("evidence_source") or "(absent)"] += 1
    w("  Repartition par aspect (comptage d'annotations, pas de genes):")
    for k, v in aspect_counter.most_common():
        w(f"    {k}: {v}")
    w("  Repartition par evidence_code:")
    total_ev = sum(evidence_counter.values())
    for k, v in evidence_counter.most_common():
        w(f"    {k}: {pct(v, total_ev)}")
    non_iea = total_ev - evidence_counter.get("IEA", 0)
    w(f"  Annotations non-IEA (plus fiables): {pct(non_iea, total_ev)}")

    # ------------------------------------------------------------------
    # 5. sources_summary : verifier l'attribution ncbi/uniprot (correction n0 2)
    # ------------------------------------------------------------------
    w("\n--- 5. Attribution des sources (sources_summary) ---")
    src_counter = collections.Counter()
    for g in genes:
        for s in g.get("sources_summary") or []:
            src_counter[s] += 1
    for k, v in src_counter.most_common():
        w(f"  {k}: {pct(v, n)}")
    # Cas suspect : source ncbi mais aucune sequence dna/rna
    suspect_ncbi = [
        g for g in genes
        if "ncbi" in (g.get("sources_summary") or [])
        and not (g.get("sequence") or {}).get("dna")
        and not (g.get("sequence") or {}).get("rna")
    ]
    w(f"  Genes marques 'ncbi' mais SANS sequence dna/rna: {pct(len(suspect_ncbi), n)}")
    if suspect_ncbi:
        w("    >>> Signe du bug d'attribution (correction n0 2) : la source 'ncbi' est")
        w("        ajoutee des qu'une sequence proteine existe, meme venue d'UniProt.")
        exemples = [g.get("gene_id") for g in suspect_ncbi[:5]]
        w(f"    Exemples de gene_id concernes: {exemples}")

    # ------------------------------------------------------------------
    # 6. Doublons de symbole (risque de fusion croisee, correction n0 3)
    # ------------------------------------------------------------------
    w("\n--- 6. Doublons de symbole (risque de fusion croisee) ---")
    def norm_symbol(s):
        return "".join(ch for ch in (s or "").lower() if ch.isalnum())

    sym_to_ids = collections.defaultdict(set)
    for g in genes:
        sym = norm_symbol(g.get("symbol"))
        if sym:
            sym_to_ids[sym].add(g.get("gene_id"))
    collisions = {s: ids for s, ids in sym_to_ids.items() if len(ids) > 1}
    total_genes_in_collision = sum(len(ids) for ids in collisions.values())
    w(f"  Symboles normalises partages par plusieurs gene_id: {len(collisions)}")
    w(f"  Genes concernes par une collision de symbole: {pct(total_genes_in_collision, n)}")
    if collisions:
        exemples = list(collisions.items())[:5]
        for sym, ids in exemples:
            w(f"    symbole '{sym}': {sorted(ids)[:6]}")

    # ------------------------------------------------------------------
    # 7. Familles PLAZA (orthologues / homologues) -- utile pour le split
    # ------------------------------------------------------------------
    w("\n--- 7. Familles PLAZA (relations) ---")
    has_ortho_fam = sum(1 for g in genes if (g.get("relations") or {}).get("orthologous_family_id"))
    has_hom_fam = sum(1 for g in genes if (g.get("relations") or {}).get("homologous_family_id"))
    w(f"  Genes avec orthologous_family_id: {pct(has_ortho_fam, n)}")
    w(f"  Genes avec homologous_family_id:  {pct(has_hom_fam, n)}")
    fam_counter = collections.Counter(
        g["relations"]["homologous_family_id"]
        for g in genes
        if (g.get("relations") or {}).get("homologous_family_id")
    )
    if fam_counter:
        top = fam_counter.most_common(5)
        w(f"  Plus grandes familles homologues: {top}")

    # ------------------------------------------------------------------
    # 8. TF family / traits / kegg / mapman -- couverture rapide
    # ------------------------------------------------------------------
    w("\n--- 8. Autres annotations ---")
    has_tf = sum(1 for g in genes if (g.get("annotation") or {}).get("tf_family"))
    has_kegg = sum(1 for g in genes if (g.get("annotation") or {}).get("kegg_pathways"))
    has_mapman = sum(1 for g in genes if (g.get("annotation") or {}).get("mapman"))
    has_traits = sum(1 for g in genes if g.get("traits"))
    w(f"  Genes avec tf_family:      {pct(has_tf, n)}")
    w(f"  Genes avec kegg_pathways:  {pct(has_kegg, n)}")
    w(f"  Genes avec mapman:         {pct(has_mapman, n)}")
    w(f"  Genes avec traits:         {pct(has_traits, n)}")

    # ------------------------------------------------------------------
    # 9. Longueurs de sequence (pas les sequences elles-memes)
    # ------------------------------------------------------------------
    w("\n--- 9. Longueurs de sequence ---")
    for st in ("dna", "protein"):
        lens = [len(g["sequence"][st]) for g in genes if (g.get("sequence") or {}).get(st)]
        if lens:
            w(f"  {st}: n={len(lens)} min={min(lens)} mediane={statistics.median(lens):.0f} "
              f"moyenne={statistics.mean(lens):.0f} max={max(lens)}")
            over_6144 = sum(1 for l in lens if l > 6144)
            if st == "dna":
                w(f"    > 6144 nt (troncature AgroNT): {pct(over_6144, len(lens))}")
        else:
            w(f"  {st}: aucune sequence")

    # ------------------------------------------------------------------
    # 10. Doublons exacts de sequence (via dna_hash)
    # ------------------------------------------------------------------
    w("\n--- 10. Doublons de sequence ADN (dna_hash) ---")
    hashes = [g["sequence"].get("dna_hash") for g in genes if (g.get("sequence") or {}).get("dna_hash")]
    if hashes:
        w(f"  Sequences ADN avec hash: {len(hashes)}, hashes distincts: {len(set(hashes))}")
        dupes = len(hashes) - len(set(hashes))
        w(f"  Doublons exacts: {dupes}")
    else:
        w("  Aucun dna_hash present.")

    # ------------------------------------------------------------------
    # 11. Score de completude (si present)
    # ------------------------------------------------------------------
    w("\n--- 11. Score de completude (quality.data_completeness) ---")
    scores = [g["quality"]["data_completeness"] for g in genes if (g.get("quality") or {}).get("data_completeness") is not None]
    if scores:
        w(f"  n={len(scores)} min={min(scores):.2f} mediane={statistics.median(scores):.2f} "
          f"moyenne={statistics.mean(scores):.2f} max={max(scores):.2f}")
        low = sum(1 for s in scores if s < 0.3)
        w(f"  Genes avec score < 0.3: {pct(low, len(scores))}")
    else:
        w("  Champ absent.")

    # ------------------------------------------------------------------
    # 12. Exemples concrets pour verification manuelle
    # ------------------------------------------------------------------
    w("\n--- 12. Exemples pour verification manuelle (identifiants seulement) ---")
    ex_go_no_seq = [g.get("gene_id") for g in has_go if not any((g.get("sequence") or {}).get(st) for st in ("dna", "rna", "protein"))][:5]
    w(f"  GO term mais aucune sequence: {ex_go_no_seq}")
    ex_ext = [g.get("gene_id") for g in genes if g.get("external_links")][:5]
    w(f"  Avec external_links non vide: {ex_ext}")
    ex_plaza_only = [g.get("gene_id") for g in genes if g.get("origin") == "plaza_only"][:5]
    w(f"  origin=plaza_only: {ex_plaza_only}")

    w("\n(Fin de l'analyse pour ce fichier)\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("files", nargs="+", help="Un ou plusieurs fichiers *_all_sources.json")
    p.add_argument("--out", help="Fichier de sortie unique pour tous les rapports (sinon: un .report.txt par fichier + affichage ecran)")
    args = p.parse_args()

    if args.out:
        out_path = Path(args.out)
        with out_path.open("w", encoding="utf-8") as out:
            tee = _Tee(out, sys.stdout)
            for f in args.files:
                analyze_file(Path(f), tee)
        print(f"\n-> Rapport combine ecrit dans {out_path}")
    else:
        for f in args.files:
            in_path = Path(f)
            out_path = in_path.with_suffix(in_path.suffix + ".report.txt")
            with out_path.open("w", encoding="utf-8") as out:
                tee = _Tee(out, sys.stdout)
                analyze_file(in_path, tee)
            print(f"-> Rapport ecrit dans {out_path}")


class _Tee:
    """Ecrit simultanement dans plusieurs flux (fichier + ecran)."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)

    def flush(self):
        for st in self.streams:
            st.flush()


if __name__ == "__main__":
    main()
