"""Remplace le chargement des index ESM-2 (relecture PostgreSQL a chaque recherche) par un cache disque memmap.

  - une matrice float32 DEJA normalisee par espece : pipeline_data/esm2_cache/<assembly>.f32.npy
  - la liste des accessions + une signature de la base (nombre, date de verification, revision) : <assembly>.meta.json
  - le cache est reconstruit tout seul si la base change ; les scores sont ceux de la base (memes vecteurs).
Simulation par defaut ; --apply ecrit (sauvegarde .bak_cache). Idempotent.
Usage : python scripts\\patch_esm2_cache.py [--apply] [--file core_engines\\esm2_similarity.py]
"""
import re
import shutil
import sys
from pathlib import Path

INIT_OLD = "def __init__(self, accessions, species, matrix: np.ndarray, assemblies=None):"
INIT_NEW = "def __init__(self, accessions, species, matrix: np.ndarray, assemblies=None, prenormalized: bool = False):"

NORM_OLD = """        norms = np.linalg.norm(matrix, axis=1, keepdims=True)
        np.divide(matrix, np.clip(norms, 1e-12, None), out=matrix)
        self._normalized = matrix
"""
NORM_NEW = """        if not prenormalized:
            norms = np.linalg.norm(matrix, axis=1, keepdims=True)
            np.divide(matrix, np.clip(norms, 1e-12, None), out=matrix)
        self._normalized = matrix
"""

START = '@st.cache_resource(show_spinner=False, max_entries=1)\ndef get_embedding_index('
END = "\n\n\ndef _search_species("

NEW_BLOCK = '''def _esm2_cache_dir():
    from pathlib import Path

    path = Path(__file__).resolve().parent.parent / "pipeline_data" / "esm2_cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _species_signature(cursor, species: str, assembly_accession: str) -> dict:
    cursor.execute(
        """
        SELECT COUNT(*),
               COALESCE(MAX(md5_verified_at)::text, ''),
               COALESCE(MIN(esm2_model_revision), ''),
               COALESCE(MAX(esm2_model_revision), '')
        FROM esm2_embeddings
        WHERE species = %s AND assembly_accession = %s
        """,
        (species, assembly_accession),
    )
    count, verified_max, rev_min, rev_max = cursor.fetchone()
    return {"count": int(count), "verified_max": verified_max,
            "revision": [rev_min, rev_max], "dim": int(EMBED_DIM), "format": 1}


def _build_disk_cache(connection, species: str, assembly_accession: str, signature: dict, npy_path, meta_path) -> None:
    """Lit la base UNE fois et ecrit la matrice normalisee (float32) + la liste des accessions."""
    import json
    import os

    row_count = signature["count"]
    if meta_path.exists():
        meta_path.unlink()
    tmp_npy = npy_path.with_name(npy_path.name + ".tmp")
    tmp_meta = meta_path.with_name(meta_path.name + ".tmp")
    matrix = np.lib.format.open_memmap(str(tmp_npy), mode="w+", dtype=np.float32, shape=(row_count, EMBED_DIM))
    accessions = []
    row_index = 0
    with connection.cursor(name="esm2_cache_cursor") as cursor:
        cursor.itersize = 2000
        cursor.execute(
            """
            SELECT protein_accession, species, embedding
            FROM esm2_embeddings
            WHERE species = %s AND assembly_accession = %s
            ORDER BY protein_accession
            """,
            (species, assembly_accession),
        )
        while rows := cursor.fetchmany(2000):
            for accession, species_name, embedding in rows:
                if row_index >= row_count:
                    raise RuntimeError("The esm2_embeddings table changed while the cache was being built.")
                if species_name != species:
                    raise RuntimeError("Embedding species changed while the cache was being built.")
                accessions.append(accession)
                matrix[row_index] = _decode_embedding(embedding)
                row_index += 1
    if row_index != row_count:
        raise RuntimeError(f"Embedding rows for {species} changed while caching: expected {row_count}, got {row_index}.")
    for start in range(0, row_count, 20000):
        block = matrix[start:start + 20000]
        norms = np.linalg.norm(block, axis=1, keepdims=True)
        np.divide(block, np.clip(norms, 1e-12, None), out=block)
    matrix.flush()
    del matrix
    with open(tmp_meta, "w", encoding="utf-8") as handle:
        json.dump({"signature": signature, "species": species, "assembly_accession": assembly_accession,
                   "accessions": accessions}, handle)
    os.replace(tmp_npy, npy_path)
    os.replace(tmp_meta, meta_path)


def _load_index_uncached(species: str) -> EmbeddingIndex:
    """Charge l'index d'une espece depuis le cache disque (le reconstruit si la base a change)."""
    import json

    connection = _connect()
    try:
        with connection.cursor() as cursor:
            cursor.execute(
                """
                SELECT assembly_accession
                FROM species_registry
                WHERE species = %s AND enabled IS TRUE AND status = 'published'
                """,
                (species,),
            )
            registry_row = cursor.fetchone()
            if not registry_row:
                raise LookupError(f"Species is not published: {species}")
            assembly_accession = registry_row[0]
            signature = _species_signature(cursor, species, assembly_accession)
        if not signature["count"]:
            raise RuntimeError(f"No ESM-2 embeddings found for published species {species}.")

        cache_dir = _esm2_cache_dir()
        npy_path = cache_dir / f"{assembly_accession}.f32.npy"
        meta_path = cache_dir / f"{assembly_accession}.meta.json"

        def _try_open():
            if not (npy_path.exists() and meta_path.exists()):
                return None
            try:
                with open(meta_path, encoding="utf-8") as handle:
                    info = json.load(handle)
                if info.get("signature") != signature:
                    return None
                matrix = np.load(str(npy_path), mmap_mode="r")
                if matrix.shape != (signature["count"], EMBED_DIM) or len(info["accessions"]) != signature["count"]:
                    return None
                return matrix, info["accessions"]
            except Exception:
                return None

        opened = _try_open()
        if opened is None:
            _build_disk_cache(connection, species, assembly_accession, signature, npy_path, meta_path)
            opened = _try_open()
            if opened is None:
                raise RuntimeError(f"ESM-2 disk cache for {species} could not be built.")
    finally:
        connection.close()

    matrix, accessions = opened
    row_count = len(accessions)
    return EmbeddingIndex(accessions, np.full(row_count, species), matrix,
                          assemblies=[assembly_accession] * row_count, prenormalized=True)


@st.cache_resource(show_spinner=False, max_entries=4)
def get_embedding_index(species: str) -> EmbeddingIndex:
    """Index d'une espece publiee : matrice memmap (cache disque), partagee entre les sessions."""
    return _load_index_uncached(species)'''


def main():
    apply = "--apply" in sys.argv
    target = Path("core_engines/esm2_similarity.py")
    if "--file" in sys.argv:
        target = Path(sys.argv[sys.argv.index("--file") + 1])
    if not target.exists():
        print("Introuvable :", target)
        sys.exit(1)
    text = target.read_bytes().decode("utf-8")
    crlf = "\r\n" in text
    work = text.replace("\r\n", "\n")
    if "_load_index_uncached" in work:
        print("Deja patche.")
        return
    problems = []
    if INIT_OLD not in work:
        problems.append("constructeur EmbeddingIndex")
    if NORM_OLD not in work:
        problems.append("bloc de normalisation")
    i, j = work.find(START), work.find(END)
    if i < 0 or j < i:
        problems.append("fonction get_embedding_index")
    if problems:
        print("Ancres introuvables : " + ", ".join(problems) + " -> fichier different de celui attendu, rien ecrit.")
        sys.exit(2)
    out = work.replace(INIT_OLD, INIT_NEW, 1).replace(NORM_OLD, NORM_NEW, 1)
    i, j = out.find(START), out.find(END)
    out = out[:i] + NEW_BLOCK + out[j:]
    compile(out, str(target), "exec")
    print(f"Patch pret pour {target} (ancres OK, syntaxe OK).")
    if not apply:
        print("Simulation OK. Relancez avec --apply pour ecrire.")
        return
    shutil.copy2(target, str(target) + ".bak_cache")
    target.write_bytes((out.replace("\n", "\r\n") if crlf else out).encode("utf-8"))
    print(f"Ecrit (sauvegarde {target}.bak_cache)")


main()
