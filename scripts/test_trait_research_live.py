"""Test automatique de trait_research.py sur la VRAIE base (lecture seule).

Lancer depuis C:\\Downloads\\IA, venv actif :
    python scripts\\test_trait_research_live.py

Sortie : une ligne PASS / FAIL / WARN par contrôle, puis un résumé.
Aucune écriture en base. Les fichiers d'export de test vont dans un dossier temporaire.
"""
from __future__ import annotations

import os
import re
import sys
import tempfile
import time
import traceback
import zipfile

ROOT = os.getcwd()
sys.path.insert(0, ROOT)

RESULTS = {"PASS": 0, "FAIL": 0, "WARN": 0}


def report(status: str, name: str, detail: str = "") -> None:
    RESULTS[status] += 1
    print(f"{status:4}  {name}" + (f"   -> {detail}" if detail else ""), flush=True)


def check(name: str, condition: bool, detail: str = "", warn_only: bool = False) -> bool:
    if condition:
        report("PASS", name, detail)
    else:
        report("WARN" if warn_only else "FAIL", name, detail)
    return bool(condition)


def in_range(value: float, low: float, high: float) -> bool:
    return low <= value <= high


def tiers(candidates: list[dict]) -> dict[str, int]:
    out: dict[str, int] = {}
    for c in candidates:
        out[c["tier"]] = out.get(c["tier"], 0) + 1
    return out


def timed(label: str, fn, *args, **kwargs):
    start = time.time()
    result = fn(*args, **kwargs)
    elapsed = time.time() - start
    print(f"      ({label}: {elapsed:.1f} s)", flush=True)
    return result, elapsed


def main() -> int:
    try:
        import trait_research as tr
    except Exception:
        traceback.print_exc()
        print("Impossible d'importer trait_research (lancez le script depuis C:\\Downloads\\IA).")
        return 2

    # ── 0. Espèces présentes dans la vue ────────────────────────────────────
    print("\n== 0. Base et vue ==")
    conn = tr._connect_bioinfo_pipeline()
    cur = conn.cursor()
    cur.execute("SELECT species, count(*) FROM v_gene_text GROUP BY species ORDER BY species")
    species_rows = cur.fetchall()
    conn.close()
    species_counts = {s: n for s, n in species_rows}
    print("      espèces dans v_gene_text :", species_counts)
    rice = next((s for s in species_counts if "oryza" in s), None)
    maize = next((s for s in species_counts if "zea" in s), None)
    check("riz présent dans la vue", rice is not None, str(rice))
    check("maïs présent dans la vue", maize is not None, str(maize))
    conn = tr._connect_bioinfo_pipeline()
    cur = conn.cursor()
    cur.execute("SELECT to_regclass('public.v_gene_func') IS NOT NULL")
    has_func = cur.fetchone()[0]
    if has_func:
        cur.execute("SELECT (SELECT count(*) FROM v_gene_func), (SELECT count(*) FROM v_gene_text)")
        n_f, n_t = cur.fetchone()
        check("v_gene_func présente et synchronisée avec v_gene_text", n_f == n_t, f"{n_f} / {n_t}")
    else:
        check("v_gene_func présente (recherche rapide)", False,
              "absente : python scripts\\create_gene_func_view.py", warn_only=True)
    conn.close()
    if rice:
        check("riz : 42 580 protéines", species_counts[rice] == 42580, str(species_counts[rice]))
    if maize:
        check("maïs : 57 345 protéines", species_counts[maize] == 57345, str(species_counts[maize]))

    # ── 1. Riz : salinité (cas vu à l'écran) ────────────────────────────────
    if rice:
        print("\n== 1. Riz — Salinity (thème) ==")
        bundle, _ = timed("salinité riz", tr.run_advanced_search, "topic", "salinity", [rice])
        check("thème reconnu", bundle.get("error") is None, str(bundle.get("error")))
        res = bundle["species"][rice]
        cands = res["candidates"]
        t = tiers(cands)
        n_valid = t.get("convergent", 0) + t.get("functional", 0)
        print("      niveaux :", t, "| correspondances SQL :", res["total_matches"],
              "| fonctionnelles SQL :", res["functional_matches"])
        check("fonctionnels validés ≈ 554 (500–620) [748 avant retrait du mot-clé « proline »]", in_range(n_valid, 500, 620), str(n_valid))
        check("convergents ≈ 42 (30–60)", in_range(t.get("convergent", 0), 30, 60), str(t.get("convergent", 0)))
        check("aucun 'positional' sans l'option", t.get("positional", 0) == 0)
        check("fonctionnel SQL ≥ validés", res["functional_matches"] >= n_valid)
        top10 = [c["gene"].get("gene_symbol") or c["gene"].get("common_name") for c in cands[:10]]
        skc1 = any(c["gene"].get("protein_accession") == "XP_015631953.1" for c in cands[:10])
        check("SKC1/OsHKT1;5 (XP_015631953.1) dans le top 10", skc1, str(top10))
        check("aucune balise HTML dans les descriptions",
              not any(re.search(r"</?(sup|sub|i|b)>", c["gene"].get("description", "")) for c in cands))
        sc = [c["score"] for c in cands]
        check("scores triés par niveau (les convergents d'abord)",
              all(tr.TIER_RANK[a["tier"]] >= tr.TIER_RANK[b["tier"]] for a, b in zip(cands, cands[1:])))

        # Règle PLAZA : texte InterPro / MapMan / trait seulement si lien concordant (plaza_gene_id rempli)
        bad = [c["gene"]["protein_accession"] for c in cands
               if (c["gene"]["annotation"].get("interpro_text") or c["gene"]["annotation"].get("mapman_text")
                   or c["gene"].get("trait_text")) and not c["gene"].get("plaza_gene_id")]
        check("règle PLAZA : aucun texte InterPro/MapMan/QTL sans lien concordant", not bad, str(bad[:5]))

        # Chaque mot-clé annoncé figure bien dans l'annotation du candidat
        miss = 0
        for c in [x for x in cands if x["tier"] != "positional"][:200]:
            g = c["gene"]
            blob = " ".join([g.get("description", ""), g["annotation"].get("go_text", ""),
                             g["annotation"].get("pfam_text", ""), g["annotation"].get("interpro_text", ""),
                             g["annotation"].get("mapman_text", ""), g.get("common_name", "") or ""]).lower()
            kws = {m.split(":", 1)[1] for m in c["matches"]}
            if kws and not any(k.lower() in blob or k.lower() in (g.get("trait_text") or "").lower() for k in kws):
                miss += 1
        check("mots-clés annoncés retrouvés dans le texte (200 premiers)", miss == 0, f"{miss} incohérents")

        # explications
        lines = tr.explain_candidate(cands[0])
        check("explication non vide pour le n°1", len(lines) >= 2, f"{len(lines)} lignes")

        # Exports
        csv_bytes = tr.candidates_to_csv(cands[:50], {"search": bundle["label"]})
        check("CSV exporté", len(csv_bytes) > 500 and b"XP_" in csv_bytes, f"{len(csv_bytes)} octets")
        out = os.path.join(tempfile.gettempdir(), "trait_test.docx")
        try:
            tr.generate_advanced_docx(bundle["label"], {"search": bundle["label"]}, cands, out, top_n=20)
            ok = os.path.exists(out) and zipfile.is_zipfile(out) and os.path.getsize(out) > 5000
            check("DOCX exporté et valide", ok, f"{os.path.getsize(out)} octets" if os.path.exists(out) else "absent")
        except Exception as e:  # noqa: BLE001
            check("DOCX exporté et valide", False, repr(e))

    # ── 2. Riz : sécheresse, avec/sans positionnel ──────────────────────────
    if rice:
        print("\n== 2. Riz — Drought (avec et sans QTL seuls) ==")
        b1, _ = timed("sécheresse riz sans positionnel", tr.run_advanced_search, "topic", "drought", [rice])
        r1 = b1["species"][rice]
        n1 = tiers(r1["candidates"])
        valid1 = n1.get("convergent", 0) + n1.get("functional", 0)
        print("      niveaux :", n1, "| total SQL :", r1["total_matches"])
        check("sécheresse riz : fonctionnels validés (thème complet) 600–1 100", in_range(valid1, 600, 1100), str(valid1))
        b2, _ = timed("sécheresse riz avec positionnel", tr.run_advanced_search, "topic", "drought", [rice], True, 2000)
        r2 = b2["species"][rice]
        n2 = tiers(r2["candidates"])
        print("      niveaux :", n2, "| total SQL :", r2["total_matches"])
        check("total SQL du thème complet 10 000–17 500", in_range(r2["total_matches"], 10000, 17500),
              str(r2["total_matches"]))
        check("QTL seuls ≈ 97 % du total", r2["total_matches"] > 0 and
              (r2["total_matches"] - r2["functional_matches"]) / r2["total_matches"] > 0.9,
              f"{r2['functional_matches']} fonctionnels / {r2['total_matches']}")
        check("avec positionnel : des 'positional' apparaissent", n2.get("positional", 0) > 0, str(n2))
        check("les fonctionnels passent avant les QTL seuls",
              all(tr.TIER_RANK[a["tier"]] >= tr.TIER_RANK[b["tier"]]
                  for a, b in zip(r2["candidates"], r2["candidates"][1:])))

    # ── 3. Maïs : sécheresse et verse ───────────────────────────────────────
    if maize:
        print("\n== 3. Maïs ==")
        bm, _ = timed("sécheresse maïs", tr.run_advanced_search, "topic", "drought", [maize])
        rm = bm["species"][maize]
        nm = tiers(rm["candidates"])
        print("      niveaux :", nm, "| total SQL :", rm["total_matches"])
        check("maïs sécheresse (thème complet) : 400–800 protéines fonctionnelles", in_range(nm.get("functional", 0), 400, 800),
              str(nm.get("functional", 0)))
        check("maïs : aucun convergent (pas de couche QTL)", nm.get("convergent", 0) == 0, str(nm))
        bl, _ = timed("verse maïs", tr.run_advanced_search, "topic", "lodging", [maize])
        check("maïs verse : exécution sans erreur", bl.get("error") is None)
        nl = tiers(bl["species"][maize]["candidates"])
        print("      niveaux :", nl)

    # ── 4. Mots-clés libres ────────────────────────────────────────────────
    if rice:
        print("\n== 4. Mots-clés libres ==")
        for q, low in (("stomatal closure", 1), ("salt/drought", 1), ("aquaporin", 5)):
            b, _ = timed(f"mots-clés « {q} »", tr.run_advanced_search, "keywords", q, [rice])
            if b.get("error"):
                check(f"« {q} »", False, str(b["error"]))
                continue
            cs = b["species"][rice]["candidates"]
            check(f"« {q} » : ≥ {low} résultat", len(cs) >= low, f"{len(cs)} candidats, SQL {b['species'][rice]['total_matches']}")
        b, _ = timed("match_all", tr.run_advanced_search, "keywords", "aquaporin, water", [rice], False, 1000, True)
        ba, _ = timed("match_any", tr.run_advanced_search, "keywords", "aquaporin, water", [rice], False, 1000, False)
        check("'tous les mots' ≤ 'au moins un mot'",
              len(b["species"][rice]["candidates"]) <= len(ba["species"][rice]["candidates"]),
              f"{len(b['species'][rice]['candidates'])} ≤ {len(ba['species'][rice]['candidates'])}")

    # ── 5. Identifiants ─────────────────────────────────────────────────────
    if rice:
        print("\n== 5. Recherche par identifiant ==")
        for q, key, low in (("PF00069", "pfam_text", 100), ("GO:0006629", "go_text", 10)):
            b, _ = timed(f"identifiant {q}", tr.run_advanced_search, "identifier", q, [rice])
            if b.get("error"):
                check(q, False, str(b["error"]))
                continue
            cs = b["species"][rice]["candidates"]
            check(f"{q} : ≥ {low} protéines", len(cs) >= low, str(len(cs)))
            if cs and key == "pfam_text":
                ok = all("PF00069" in c["gene"]["annotation"].get("pfam_text", "") for c in cs[:100])
                check("PF00069 : chaque protéine porte bien le domaine (100 premières)", ok)
        b, _ = timed("accession", tr.run_advanced_search, "identifier", "XP_015631953.1", [rice])
        cs = b["species"][rice]["candidates"] if not b.get("error") else []
        check("accession XP_015631953.1 retrouvée", len(cs) == 1, str(len(cs)))

    # ── 6. Cas limites ──────────────────────────────────────────────────────
    print("\n== 6. Cas limites ==")
    b = tr.run_advanced_search("topic", "zzzz inconnu", [rice or "x"])
    check("thème inconnu → erreur propre", b.get("error") == "topic_not_recognized", str(b.get("error")))
    b = tr.run_advanced_search("keywords", "   ", [rice or "x"])
    check("mots-clés vides → erreur propre", b.get("error") == "no_keywords", str(b.get("error")))
    if rice and maize:
        b, el = timed("deux espèces", tr.run_advanced_search, "topic", "drought", [rice, maize])
        merged = tr.merge_ranked(b["species"])
        check("fusion riz + maïs", {c["gene"]["organism"] for c in merged} == {rice, maize},
              str({c["gene"]["organism"] for c in merged}))
        check("recherche deux espèces < 20 s (objectif après optimisation)", el < 20, f"{el:.1f} s", warn_only=True)

    print("\n== Résumé ==")
    print(f"PASS {RESULTS['PASS']}   WARN {RESULTS['WARN']}   FAIL {RESULTS['FAIL']}")
    return 1 if RESULTS["FAIL"] else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(3)
