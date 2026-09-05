"""Shared reference-result handling for organism-level comparisons."""


def get_organism_reference(
    organism: str | None,
    metric: str,
    min_n: int = 10,
    fetcher=None,
) -> dict:
    """Normalize organism reference data and its standard fallback message."""
    if not organism or fetcher is None:
        return {"available": False, "value": None, "n": 0, "fallback_reason": f"No organism reference available for {metric}."}
    try:
        raw = fetcher(organism) or {}
    except Exception:
        raw = {}
    value = raw.get("value", raw.get(metric, raw.get(f"mean_{metric}")))
    n = int(raw.get("n", raw.get("n_sequences", 0)) or 0)
    available = value is not None and n >= min_n
    return {
        "available": available,
        "value": value if available else None,
        "n": n,
        "fallback_reason": None if available else f"Insufficient organism-specific data for {metric} (n={n}, minimum={min_n}); using no species comparison.",
    }