import hashlib
import json


def build_analysis_signature(
    sequence: str,
    reference_sequence: str = "",
    input_type: str = "DNA",
    reading_frame: int = 1,
    top_n_matches: int = 10,
    similarity_deep_search: bool = False,
) -> str:
    """Stable hash for the current analysis request.

    The UI must not reuse stale results when the pasted sequence or analysis
    settings change. This signature is intentionally based only on the inputs
    that materially affect the analysis, not on the computed output.
    """
    payload = {
        "sequence": (sequence or "").strip(),
        "reference_sequence": (reference_sequence or "").strip(),
        "input_type": input_type,
        "reading_frame": reading_frame,
        "top_n_matches": top_n_matches,
        "similarity_deep_search": similarity_deep_search,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def should_reuse_cached_results(current_signature: str | None, last_signature: str | None, current_sequence: str | None) -> bool:
    """True when the cached result still matches the current input."""
    return bool(current_sequence) and bool(last_signature) and bool(current_signature) and current_signature == last_signature
