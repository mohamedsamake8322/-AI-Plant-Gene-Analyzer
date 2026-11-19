"""Lecture des donnees de bioinfo_pipeline pour la barre laterale, la recherche et les statistiques.

Remplace, derriere l'interrupteur USE_BIOINFO, les fonctions de scripts/postgres_utils.py qui lisaient
l'ancienne table `genes`. Memes noms, memes formats de retour.

Interrupteur : variable d'environnement USE_BIOINFO=1 (actif) / 0 (ancien chemin).
Pour activer par defaut apres validation : mettre USE_BIOINFO_DEFAULT = "1" ci-dessous.
Les statistiques (GC, codons, longueurs) viennent de la table pre-calculee `species_cds_stats`
(construite par scripts/build_species_stats.py), jamais calculees a la volee.
"""
from __future__ import annotations

import math
import os

USE_BIOINFO_DEFAULT = "0"

SPECIES_LABELS = {
    "oryza_sativa": "Oryza sativa",
    "zea_mays": "Zea mays",
    "nicotiana_tabacum_k326": "Nicotiana tabacum K326",
    "nicotiana_tabacum_tn90": "Nicotiana tabacum TN90",
}


def is_enabled() -> bool:
    return os.getenv("USE_BIOINFO", USE_BIOINFO_DEFAULT).strip() == "1"


def species_label(species: str, display_name: str | None = None) -> str:
    if species in SPECIES_LABELS:
        return SPECIES_LABELS[species]
    return (species or display_name or "").replace("_", " ").capitalize()


def _connect():
    from core_engines.bioinfo_db import connect
    return connect()


def published_species(cur) -> list[dict]:
    cur.execute("SELECT species, display_name, assembly_accession FROM species_registry "
                "WHERE enabled AND status = 'published' ORDER BY species")
    return [{"species": s, "display_name": d, "assembly_accession": a, "label": species_label(s, d)}
            for s, d, a in cur.fetchall()]


def resolve_species(rows: list[dict], organism: str | None) -> list[dict]:
    """Retourne les especes publiees correspondant a `organism` (cle, etiquette, nom d'affichage ou genre-espece)."""
    s = (organism or "").strip().lower()
    if not s:
        return []
    for pred in (
        lambda r: r["species"].lower() == s,
        lambda r: r["label"].lower() == s,
        lambda r: (r["display_name"] or "").lower() == s,
    ):
        hit = [r for r in rows if pred(r)]
        if hit:
            return hit
    hit = [r for r in rows if r["label"].lower().startswith(s)]
    if hit:
        return hit
    return [r for r in rows if s.startswith(" ".join(r["label"].lower().split()[:2]))]


# ------------------------------------------------------------------ barre laterale
_COND = ("(v.protein_accession ILIKE %(p)s OR v.gene_symbol ILIKE %(p)s OR v.description ILIKE %(p)s "
         "OR v.search_vector @@ websearch_to_tsquery('simple', translate(%(q)s, './', '  ')))")


def get_gene_count() -> int:
    with _connect() as conn, conn.cursor() as cur:
        cur.execute("SELECT count(*) FROM v_gene_text")
        return int(cur.fetchone()[0])


def search_gene_metadata(query: str | None, limit: int = 20, offset: int = 0) -> list[dict]:
    q = (query or "").strip()
    where = f"WHERE {_COND}" if q else ""
    params = {"p": f"%{q}%", "q": q, "limit": int(limit), "offset": int(offset)}
    with _connect() as conn, conn.cursor() as cur:
        cur.execute(
            f"""SELECT v.protein_accession, v.gene_symbol, v.species, r.display_name, v.description,
                       v.assembly_accession, v.trait_text, length(s.protein_sequence)
                FROM v_gene_text v
                JOIN species_registry r ON r.species = v.species
                LEFT JOIN dim_sequence s ON s.assembly_accession = v.assembly_accession
                                        AND s.protein_accession = v.protein_accession
                {where}
                ORDER BY v.species, v.protein_accession LIMIT %(limit)s OFFSET %(offset)s""", params)
        rows = cur.fetchall()
    return [{
        "gene_id": acc, "symbol": symbol, "organism": species_label(species, display),
        "description": description or "", "source": f"RefSeq {assembly}",
        "traits": [trait] if trait else [], "length": length,
    } for acc, symbol, species, display, description, assembly, trait, length in rows]


def count_gene_metadata_matches(query: str | None) -> int:
    q = (query or "").strip()
    if not q:
        return get_gene_count()
    with _connect() as conn, conn.cursor() as cur:
        cur.execute(f"SELECT count(*) FROM v_gene_text v WHERE {_COND}", {"p": f"%{q}%", "q": q})
        return int(cur.fetchone()[0])


# ------------------------------------------------------------------ statistiques par organisme
class StatsUnavailable(RuntimeError):
    pass


def _stat_rows(organism: str) -> list[dict]:
    with _connect() as conn, conn.cursor() as cur:
        species = resolve_species(published_species(cur), organism)
        if not species:
            return []
        try:
            cur.execute("SELECT assembly_accession, n_sequences, gc_mean, gc_m2, len_mean, len_m2, codon_counts "
                        "FROM species_cds_stats WHERE assembly_accession = ANY(%s)",
                        ([s["assembly_accession"] for s in species],))
        except Exception as exc:  # noqa: BLE001
            conn.rollback()
            raise StatsUnavailable(
                "Table species_cds_stats absente : lancez scripts\\build_species_stats.py --apply") from exc
        return [dict(zip(("assembly", "n", "gc_mean", "gc_m2", "len_mean", "len_m2", "codons"), r))
                for r in cur.fetchall()]


def _pool(rows: list[dict], mean_key: str, m2_key: str) -> tuple[float, float, int]:
    n = mean = m2 = 0.0
    for r in rows:
        nb, mb, m2b = r["n"], r[mean_key], r[m2_key]
        if not nb:
            continue
        delta = mb - mean
        tot = n + nb
        mean += delta * nb / tot
        m2 += m2b + delta * delta * n * nb / tot
        n = tot
    sd = math.sqrt(m2 / (n - 1)) if n > 1 else 0.0
    return mean, sd, int(n)


def get_gc_content_stats_for_organism(organism: str) -> dict:
    rows = _stat_rows(organism) if organism else []
    mean, sd, n = _pool(rows, "gc_mean", "gc_m2")
    return {"mean_gc": round(mean, 2), "stdev_gc": round(sd, 2), "n_sequences": n}


def get_length_stats_for_organism(organism: str) -> dict:
    rows = _stat_rows(organism) if organism else []
    mean, sd, n = _pool(rows, "len_mean", "len_m2")
    return {"mean_length": round(mean, 2), "stdev_length": round(sd, 2), "n_sequences": n}


def get_codon_usage_for_organism(organism: str) -> dict:
    rows = _stat_rows(organism) if organism else []
    total: dict[str, float] = {}
    for r in rows:
        for codon, count in (r["codons"] or {}).items():
            total[codon] = total.get(codon, 0.0) + float(count)
    s = sum(total.values())
    return {c: v / s for c, v in total.items()} if s else {}


def get_codon_reference_for_organism(organism: str) -> dict:
    rows = _stat_rows(organism) if organism else []
    return {"value": get_codon_usage_for_organism(organism) if rows else {}, "n": int(sum(r["n"] for r in rows))}
