def _alignment_window_slider_bounds(alignment_length):
    """Compute valid Streamlit slider bounds for a visible alignment window.

    The previous implementation used a 20-column minimum with a max derived from
    the actual alignment length. When the alignment length was exactly 20, both
    bounds collapsed to the same value and Streamlit raised:
    "Slider `min_value` must be less than the `max_value`."
    """
    length = max(1, int(alignment_length or 1))

    min_value = 1 if length <= 20 else 20
    max_value = max(min_value + 1, min(200, length))
    if max_value <= min_value:
        max_value = min_value + 1

    value = min(max_value, max(min_value, min(60, length)))
    if value < min_value:
        value = min_value

    return min_value, max_value, value
