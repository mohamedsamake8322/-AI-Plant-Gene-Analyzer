"""Corrige _build_disk_cache dans core_engines/esm2_similarity.py (deja patche) : plus de memmap pendant la construction
(erreur Windows 'fichier utilise par un autre processus'). Simulation par defaut ; --apply ecrit (.bak_cache2). Idempotent."""
import shutil
import sys
from pathlib import Path

NEW = r'''def _build_disk_cache(connection, species: str, assembly_accession: str, signature: dict, npy_path, meta_path) -> None:
    """Lit la base UNE fois et ecrit la matrice normalisee (float32) + la liste des accessions.

    Ecriture par blocs dans un fichier ordinaire (pas de memmap pendant la construction : sous Windows
    un memmap encore ouvert empeche le renommage final).
    """
    import json
    import os

    row_count = signature["count"]
    if meta_path.exists():
        meta_path.unlink()
    tmp_npy = npy_path.with_name(npy_path.name + ".tmp")
    tmp_meta = meta_path.with_name(meta_path.name + ".tmp")
    accessions = []
    row_index = 0
    block_rows = 2000
    buffer = np.empty((block_rows, EMBED_DIM), dtype=np.float32)
    filled = 0

    def _flush(handle):
        nonlocal filled
        if not filled:
            return
        block = buffer[:filled]
        norms = np.linalg.norm(block, axis=1, keepdims=True)
        np.divide(block, np.clip(norms, 1e-12, None), out=block)
        handle.write(block.tobytes())
        filled = 0

    with open(tmp_npy, "wb") as handle:
        np.lib.format.write_array_header_1_0(
            handle, {"descr": "<f4", "fortran_order": False, "shape": (row_count, EMBED_DIM)}
        )
        with connection.cursor(name="esm2_cache_cursor") as cursor:
            cursor.itersize = block_rows
            cursor.execute(
                """
                SELECT protein_accession, species, embedding
                FROM esm2_embeddings
                WHERE species = %s AND assembly_accession = %s
                ORDER BY protein_accession
                """,
                (species, assembly_accession),
            )
            while rows := cursor.fetchmany(block_rows):
                for accession, species_name, embedding in rows:
                    if row_index >= row_count:
                        raise RuntimeError("The esm2_embeddings table changed while the cache was being built.")
                    if species_name != species:
                        raise RuntimeError("Embedding species changed while the cache was being built.")
                    accessions.append(accession)
                    buffer[filled] = _decode_embedding(embedding)
                    filled += 1
                    row_index += 1
                    if filled == block_rows:
                        _flush(handle)
        _flush(handle)
    if row_index != row_count:
        raise RuntimeError(f"Embedding rows for {species} changed while caching: expected {row_count}, got {row_index}.")
    with open(tmp_meta, "w", encoding="utf-8") as handle:
        json.dump({"signature": signature, "species": species, "assembly_accession": assembly_accession,
                   "accessions": accessions}, handle)
    os.replace(tmp_npy, npy_path)
    os.replace(tmp_meta, meta_path)


'''


def main():
    apply = "--apply" in sys.argv
    target = Path("core_engines/esm2_similarity.py")
    text = target.read_bytes().decode("utf-8")
    crlf = "\r\n" in text
    work = text.replace("\r\n", "\n")
    if "write_array_header_1_0" in work:
        print("Deja corrige.")
        return
    i, j = work.find("def _build_disk_cache("), work.find("def _load_index_uncached(")
    if i < 0 or j < i:
        print("Ancres introuvables (appliquez d'abord patch_esm2_cache.py) : rien ecrit.")
        sys.exit(2)
    out = work[:i] + NEW + work[j:]
    compile(out, str(target), "exec")
    print("Correction prete (syntaxe OK).")
    if not apply:
        print("Simulation OK. Relancez avec --apply.")
        return
    shutil.copy2(target, str(target) + ".bak_cache2")
    target.write_bytes((out.replace("\n", "\r\n") if crlf else out).encode("utf-8"))
    print("Ecrit (sauvegarde " + str(target) + ".bak_cache2)")


main()
