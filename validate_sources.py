"""Valide sources.yaml et applique la porte de licence.

Usage :
    python validate_sources.py sources.yaml            # structure + etat des licences
    python validate_sources.py sources.yaml --check-urls   # teste aussi les URLs (reseau)

Sortie : code 0 si la structure est valide (meme si des licences restent a verifier).
Les sources `ready` sont celles dont license_checked est true et status != blocked.
"""
import sys
from pathlib import Path
import urllib.error
import urllib.request

import yaml

REQUIRED = ["id", "species", "phase", "role", "url", "license_expected",
            "license_checked", "redistribution", "status"]
STATUS = {"to_verify", "verified", "blocked"}
REDIST = {"open", "private_only", "unknown"}
CAPABILITIES = {"sequence", "annotation", "go", "ontology", "domains", "id_map", "ortho",
                "pathway", "resistance", "expression", "variation", "qtl", "phenotype", "training"}


def check_url(url, timeout=30, retries=2):
    """HEAD d'abord ; si refuse (403/405) ou erreur, GET avec Range 0-0. Reessaie en cas de timeout."""
    import time
    last = None
    for attempt in range(retries + 1):
        for method in ("HEAD", "GET"):
            headers = {"User-Agent": "Mozilla/5.0 (sources-validator)"}
            if method == "GET":
                headers["Range"] = "bytes=0-0"
            req = urllib.request.Request(url, method=method, headers=headers)
            try:
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    return f"OK {r.status} ({method})"
            except urllib.error.HTTPError as exc:
                last = f"HTTP {exc.code} ({method})"
                if exc.code in (403, 405, 501):
                    continue          # essayer GET
                return last
            except Exception as exc:
                reason = str(exc)
                if "CERTIFICATE_VERIFY_FAILED" in reason:
                    return "SSL: certificat non valide (voir le message ; ne pas desactiver la verification globalement)"
                last = f"ERREUR ({reason[:80]})"
                break                 # timeout/reseau : reessayer
        time.sleep(2 * (attempt + 1))
    return last


def main(path, check_urls=False, species_path=None):
    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    species_path = species_path or str(Path(path).with_name("species.yaml"))
    with open(species_path, encoding="utf-8") as f:
        known_species = set(yaml.safe_load(f)["species"])

    errors, seen = [], set()
    for s in cfg.get("sources", []):
        sid = s.get("id", "<sans id>")
        for key in REQUIRED:
            if key not in s:
                errors.append(f"{sid}: champ manquant `{key}`")
        if sid in seen:
            errors.append(f"{sid}: id en double")
        seen.add(sid)
        if s.get("status") not in STATUS:
            errors.append(f"{sid}: status invalide ({s.get('status')})")
        if s.get("redistribution") not in REDIST:
            errors.append(f"{sid}: redistribution invalide ({s.get('redistribution')})")
        sp = s.get("species", [])
        for one in ([sp] if isinstance(sp, str) else sp):
            if one not in known_species:
                errors.append(f"{sid}: espece inconnue dans species.yaml ({one})")
        caps = s.get("capability")
        if not caps:
            errors.append(f"{sid}: champ `capability` manquant")
        else:
            for c in caps:
                if c not in CAPABILITIES:
                    errors.append(f"{sid}: capability inconnue ({c})")
        if s.get("status") == "verified" and not s.get("license_checked"):
            errors.append(f"{sid}: status verified mais licence non verifiee")

    ready = [s["id"] for s in cfg["sources"]
             if s.get("license_checked") and s.get("status") != "blocked"]
    pending = [s["id"] for s in cfg["sources"] if not s.get("license_checked")]

    print(f"{len(cfg['sources'])} sources, {len(ready)} autorisees au telechargement, "
          f"{len(pending)} en attente de verification de licence.")
    if pending:
        print("  En attente :", ", ".join(pending))

    publishable = [s["id"] for s in cfg["sources"]
                   if s.get("license_checked") and s.get("redistribution") == "open"
                   and s.get("status") != "blocked"]
    print(f"Publiables dans le dataset public : {len(publishable)}"
          + (" (" + ", ".join(publishable) + ")" if publishable else ""))

    by_cap = {}
    for s in cfg["sources"]:
        for c in s.get("capability", []):
            by_cap.setdefault(c, []).append(s["id"])
    print("\nSources par capacite :")
    for c in sorted(by_cap):
        print(f"  {c}: {len(by_cap[c])}")

    if check_urls:
        print("\nTest des URLs :")
        for s in cfg["sources"]:
            print(f"  {s['id']}: {check_url(s['url'])}")

    if errors:
        print("\nERREURS :")
        for e in errors:
            print("  -", e)
        return 1
    print("\nStructure valide.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], "--check-urls" in sys.argv))
