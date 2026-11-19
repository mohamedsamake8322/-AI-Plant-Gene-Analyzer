#!/usr/bin/env python
"""Pre-calcule les statistiques de CDS par espece (GC, longueur, codons) dans species_cds_stats.

Source : dim_sequence.cds_sequence, hors CDS divergentes de leur proteine (cds_matches_protein non
'exact' ni 'near'), comme prevu pour les statistiques de codons. Lecture en flux : memoire faible.
Simulation par defaut (affiche les statistiques, n'ecrit rien).

Usage (depuis C:\\Downloads\\IA, venv active) :
    python scripts\\build_species_stats.py
    python scripts\\build_species_stats.py --apply        # cree/remplit la table species_cds_stats
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time

sys.path.insert(0, ".")

DDL = """CREATE TABLE IF NOT EXISTS species_cds_stats (
    assembly_accession text PRIMARY KEY,
    species text NOT NULL,
    n_sequences integer NOT NULL,
    n_excluded_divergent integer NOT NULL,
    gc_mean double precision, gc_m2 double precision,
    len_mean double precision, len_m2 double precision,
    codon_counts jsonb,
    method text,
    computed_at timestamptz DEFAULT now()
)"""
CODON = re.compile(r"^[ATGC]{3}$")
METHOD = ("per-CDS GC%% mean (Welford); lengths in nt; codon counts in frame 1, ACGT-only; "
          "CDS with cds_matches_protein NOT IN ('exact','near') excluded")


def compute(conn, assembly: str) -> dict:
    n = n_div = 0
    gc_mean = gc_m2 = len_mean = len_m2 = 0.0
    codons: dict[str, int] = {}
    with conn.cursor(name=f"cds_{assembly.replace('.', '_')}") as cur:
        cur.itersize = 2000
        cur.execute("SELECT cds_sequence, cds_matches_protein FROM dim_sequence "
                    "WHERE assembly_accession = %s AND cds_sequence IS NOT NULL AND cds_sequence <> ''", (assembly,))
        for seq, match in cur:
            if match not in ("exact", "near"):
                n_div += 1
                continue
            seq = seq.upper()
            L = len(seq)
            gc = 100.0 * (seq.count("G") + seq.count("C")) / L
            n += 1
            d = gc - gc_mean
            gc_mean += d / n
            gc_m2 += d * (gc - gc_mean)
            d = L - len_mean
            len_mean += d / n
            len_m2 += d * (L - len_mean)
            for i in range(0, L - 2, 3):
                c = seq[i:i + 3]
                if CODON.match(c):
                    codons[c] = codons.get(c, 0) + 1
    return {"n": n, "n_div": n_div, "gc_mean": gc_mean, "gc_m2": gc_m2, "len_mean": len_mean,
            "len_m2": len_m2, "codons": codons}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    from core_engines.bioinfo_db import connect
    with connect() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT species, assembly_accession FROM species_registry "
                        "WHERE enabled AND status = 'published' ORDER BY 1")
            species = cur.fetchall()
        conn.commit()
        results = []
        for sp, asm in species:
            t0 = time.time()
            r = compute(conn, asm)
            conn.commit()
            sd_gc = (r["gc_m2"] / (r["n"] - 1)) ** 0.5 if r["n"] > 1 else 0
            sd_len = (r["len_m2"] / (r["n"] - 1)) ** 0.5 if r["n"] > 1 else 0
            print(f"{sp:<26} CDS {r['n']:>7,} (exclues divergentes {r['n_div']:,})  GC {r['gc_mean']:.2f} +/- {sd_gc:.2f} %"
                  f"  longueur {r['len_mean']:.0f} +/- {sd_len:.0f} nt  codons {sum(r['codons'].values()):,}  ({time.time() - t0:.0f} s)")
            results.append((sp, asm, r))
        if not a.apply:
            print("\nSimulation terminee. Relancez avec --apply pour ecrire species_cds_stats.")
            return
        with conn.cursor() as cur:
            cur.execute(DDL)
            for sp, asm, r in results:
                cur.execute(
                    """INSERT INTO species_cds_stats (assembly_accession, species, n_sequences, n_excluded_divergent,
                           gc_mean, gc_m2, len_mean, len_m2, codon_counts, method, computed_at)
                       VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, now())
                       ON CONFLICT (assembly_accession) DO UPDATE SET species=EXCLUDED.species,
                           n_sequences=EXCLUDED.n_sequences, n_excluded_divergent=EXCLUDED.n_excluded_divergent,
                           gc_mean=EXCLUDED.gc_mean, gc_m2=EXCLUDED.gc_m2, len_mean=EXCLUDED.len_mean,
                           len_m2=EXCLUDED.len_m2, codon_counts=EXCLUDED.codon_counts, method=EXCLUDED.method,
                           computed_at=now()""",
                    (asm, sp, r["n"], r["n_div"], r["gc_mean"], r["gc_m2"], r["len_mean"], r["len_m2"],
                     json.dumps(r["codons"]), METHOD % ()))
        conn.commit()
        print(f"\nEcrit : species_cds_stats ({len(results)} lignes).")


if __name__ == "__main__":
    main()
