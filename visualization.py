"""
visualization.py
----------------
Visualization module for the Plant Gene Analyzer.
Produces Plotly and Matplotlib figures from analysis results.
All functions return Plotly figures (compatible with st.plotly_chart).

Layout colors (background, grid, font) and the nucleotide palette are read
from config.py, which is the single source of truth for the app's dark
bio-tech theme (also used by style.css). Semantic accent colors used only
for chart-specific meaning (e.g. "high similarity" vs "low similarity") are
defined locally below.
"""

import re
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

import config


# ─── Semantic accents (chart-only meaning, not part of the base UI theme) ─────
TEAL = config.CHART_TITLE_COLOR   # "#00d9a3" — primary accent / good score
CYAN = "#4fc3f7"     # secondary accent
AMBER = "#ffd166"    # caution / mid-range
CORAL = "#ff6b6b"    # mutation / mismatch / low-score
MINT = "#69f0ae"     # success / match
SLATE = "#5f7a86"    # muted / neutral (e.g. ambiguous bases)

THEME = dict(
    paper=config.CHART_PAPER,
    plot_bg=config.CHART_BG,
    font_color=config.CHART_FONT_COLOR,
    grid_color=config.CHART_GRID_COLOR,
    line_color=config.CHART_LINE_COLOR,
    title_color=config.CHART_TITLE_COLOR,
)

# ─── Color palette ─────────────────────────────────────────────────────────────
NUCLEOTIDE_COLORS: dict[str, str] = config.NUCLEOTIDE_COLORS


def _base_layout(title: str = "") -> dict:
    """Shared Plotly layout for the app's dark theme."""
    return dict(
        title=dict(text=title, font=dict(color=THEME["title_color"], size=16, family="Space Grotesk, Arial")),
        paper_bgcolor=THEME["paper"],
        plot_bgcolor=THEME["plot_bg"],
        font=dict(color=THEME["font_color"], family="Inter, Arial"),
        margin=dict(l=40, r=40, t=60, b=40),
        # tickfont is set explicitly here (not left to inherit from the
        # figure-level `font` above) because it wasn't: axis tick labels
        # were rendering at low contrast in production (barely-visible
        # gray gene names on the Similarity bar chart) even though
        # THEME["font_color"] is a light color meant for a dark background.
        # Plotly's font inheritance for tick labels can be overridden by
        # page-level CSS targeting the SVG text elements it renders into
        # (a likely candidate here, alongside style.css also breaking the
        # Material Symbols icon font elsewhere in the app) — setting it
        # explicitly per-axis is a defensive fix that doesn't depend on
        # sorting out that CSS interaction first.
        xaxis=dict(
            gridcolor=THEME["grid_color"], zerolinecolor=THEME["grid_color"], linecolor=THEME["line_color"],
            tickfont=dict(color=THEME["font_color"]),
        ),
        yaxis=dict(
            gridcolor=THEME["grid_color"], zerolinecolor=THEME["grid_color"], linecolor=THEME["line_color"],
            tickfont=dict(color=THEME["font_color"]),
        ),
    )


def _normalize_plotly_color(color_value: str) -> str:
    """Normalize matplotlib-style color names for Plotly."""
    named_colors = {
        'C0': TEAL,
        'C1': CYAN,
        'C2': AMBER,
        'C3': CORAL,
        'C4': '#9467bd',
        'C5': '#8c564b',
        'C6': '#e377c2',
        'C7': SLATE,
        'C8': '#bcbd22',
        'C9': '#17becf',
    }
    if color_value in named_colors:
        return named_colors[color_value]
    if isinstance(color_value, str) and color_value.startswith('C') and color_value[1:].isdigit():
        return named_colors.get(color_value, TEAL)
    return color_value


# ─── Nucleotide distribution ───────────────────────────────────────────────────

def plot_nucleotide_pie(dist: dict) -> go.Figure:
    """
    Pie chart of nucleotide composition.

    Args:
        dist: output of bioinformatics.nucleotide_distribution()
    """
    counts = dist["counts"]
    labels = [k for k, v in counts.items() if v > 0]
    values = [counts[k] for k in labels]
    colors = [NUCLEOTIDE_COLORS.get(k, SLATE) for k in labels]

    fig = go.Figure(
        go.Pie(
            labels=labels,
            values=values,
            marker=dict(colors=colors, line=dict(color="#0d1b2a", width=2)),
            textinfo="label+percent",
            hovertemplate="<b>%{label}</b><br>Count: %{value}<br>%{percent}<extra></extra>",
            hole=0.4,
        )
    )
    fig.update_layout(
        **_base_layout("Nucleotide Composition"),
        showlegend=True,
        legend=dict(font=dict(color=THEME["font_color"])),
        annotations=[
            dict(
                text="DNA<br>Base",
                x=0.5, y=0.5,
                font=dict(size=13, color=TEAL),
                showarrow=False,
            )
        ],
    )
    return fig


def plot_nucleotide_bar(dist: dict) -> go.Figure:
    """
    Bar chart of nucleotide counts. Categories with a zero count are omitted
    (e.g. "N" on a clean sequence) so they don't take up visual space with no
    information to show.
    """
    counts = dist["counts"]
    nucleotides = [n for n, v in counts.items() if v > 0]
    values = [counts[n] for n in nucleotides]
    colors = [NUCLEOTIDE_COLORS.get(n, SLATE) for n in nucleotides]

    fig = go.Figure(
        go.Bar(
            x=nucleotides,
            y=values,
            marker=dict(color=colors, line=dict(color="#0d1b2a", width=1)),
            text=values,
            textposition="outside",
            textfont=dict(color=THEME["font_color"]),
            hovertemplate="<b>%{x}</b><br>Count: %{y}<extra></extra>",
        )
    )
    layout = _base_layout("Nucleotide Counts")
    layout["yaxis"]["title"] = "Count"
    layout["xaxis"]["title"] = "Nucleotide"
    fig.update_layout(**layout)
    return fig


def plot_amino_acid_bar(dist: dict) -> go.Figure:
    """Bar chart of amino acid composition for protein sequences."""
    counts = dist["counts"]
    residues = [aa for aa, count in counts.items() if count > 0]
    values = [counts[aa] for aa in residues]
    colors = [TEAL if aa in {"A", "G", "V", "L", "I", "M"} else CYAN for aa in residues]

    fig = go.Figure(
        go.Bar(
            x=residues,
            y=values,
            marker=dict(color=colors, line=dict(color="#0d1b2a", width=1)),
            text=values,
            textposition="outside",
            textfont=dict(color=THEME["font_color"]),
            hovertemplate="<b>%{x}</b><br>Count: %{y}<extra></extra>",
        )
    )
    layout = _base_layout("Amino Acid Composition")
    layout["yaxis"]["title"] = "Count"
    layout["xaxis"]["title"] = "Amino Acid"
    fig.update_layout(**layout)
    return fig


# ─── GC content gauge ─────────────────────────────────────────────────────────

def plot_gc_gauge(
    gc_percent: float,
    reference_low: float = 35.0,
    reference_high: float = 65.0,
) -> go.Figure:
    """
    Gauge chart for GC content.

    The reference bands (low / balanced / high) are descriptive only — GC
    content varies naturally by species and gene region, so a "high" reading
    is not inherently good or bad. The caption below the gauge makes this
    explicit instead of relying on a red/green traffic-light color scheme.
    """
    fig = go.Figure(
        go.Indicator(
            mode="gauge+number+delta",
            value=gc_percent,
            number=dict(suffix="%", font=dict(color=TEAL, size=36)),
            delta=dict(
                reference=50,
                increasing=dict(color=AMBER),
                decreasing=dict(color=CYAN),
            ),
            gauge=dict(
                axis=dict(
                    range=[0, 100],
                    tickcolor=THEME["font_color"],
                    tickfont=dict(color=THEME["font_color"]),
                ),
                bar=dict(color=TEAL, thickness=0.25),
                bgcolor=THEME["plot_bg"],
                bordercolor=TEAL,
                steps=[
                    dict(range=[0, reference_low], color="rgba(79,195,247,0.15)"),
                    dict(range=[reference_low, reference_high], color="rgba(0,217,163,0.15)"),
                    dict(range=[reference_high, 100], color="rgba(255,209,102,0.15)"),
                ],
                threshold=dict(
                    line=dict(color=AMBER, width=3),
                    thickness=0.8,
                    value=(reference_low + reference_high) / 2,
                ),
            ),
            title=dict(
                text="GC Content",
                font=dict(color=TEAL, size=16),
            ),
        )
    )
    fig.update_layout(
        paper_bgcolor=THEME["paper"],
        font=dict(color=THEME["font_color"], family="JetBrains Mono, Consolas, monospace"),
        margin=dict(l=20, r=20, t=60, b=50),
        height=300,
        annotations=[
            dict(
                text="Reference bands only — a high or low GC% reflects species/gene "
                     "characteristics, not sequence quality.",
                x=0.5, y=-0.08,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=10, color=SLATE),
            )
        ],
    )
    return fig


# ─── Similarity scores ─────────────────────────────────────────────────────────

def plot_similarity_scores(similarity_results: list[dict]) -> go.Figure:
    """
    Horizontal bar chart comparing similarity scores against database genes.

    Args:
        similarity_results: list from similarityengine.compare_with_database()
    """
    if not similarity_results:
        return go.Figure()

    def clean_gene_name(raw_name: str) -> str:
        return re.sub(r"^_[a-z]{2,20}[_-]", "", raw_name or "", flags=re.IGNORECASE)

    genes = [clean_gene_name(r["gene_name"]) for r in similarity_results]
    scores = [r["similarity_score"] for r in similarity_results]
    traits = [r["trait"] for r in similarity_results]
    raw_genes = [r["gene_name"] for r in similarity_results]

    colors = [
        TEAL if s >= 75 else AMBER if s >= 50 else CORAL
        for s in scores
    ]

    fig = go.Figure(
        go.Bar(
            x=scores,
            y=genes,
            orientation="h",
            marker=dict(
                color=colors,
                line=dict(color="#0d1b2a", width=1),
            ),
            text=[f"{s:.1f}%" for s in scores],
            textposition="outside",
            textfont=dict(color=THEME["font_color"]),
            customdata=list(zip(raw_genes, traits)),
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "Gene: %{y}<br>"
                "Trait: %{customdata[1]}<br>"
                "Similarity: %{x:.1f}%<extra></extra>"
            ),
        )
    )
    layout = _base_layout("Database Similarity Scores")
    layout["xaxis"]["title"] = "Similarity (%)"
    layout["xaxis"]["range"] = [0, max(100, max(scores, default=0) + 8)]
    layout["xaxis"]["ticksuffix"] = "%"
    layout["yaxis"]["title"] = "Gene"
    layout["yaxis"]["tickfont"] = dict(size=12)
    layout["margin"] = dict(l=170, r=75, t=55, b=55)
    fig.update_layout(
        **layout,
        height=max(320, len(genes) * 78),
        bargap=0.28,
    )
    return fig


# ─── Alignment visualization ───────────────────────────────────────────────────

def plot_alignment(alignment_map: dict, max_chars: int = 60) -> go.Figure:
    """
    Display query-reference alignment as an annotated heatmap.

    Args:
        alignment_map: output of similarityengine._build_alignment_map()
        max_chars:     max characters to show
    """
    query = alignment_map.get("query", "")[:max_chars]
    ref = alignment_map.get("reference", "")[:max_chars]
    match_line = alignment_map.get("match_line", "")[:max_chars]

    if not query or not ref:
        return go.Figure()

    n = len(match_line)
    positions = list(range(1, n + 1))

    match_values = [1 if c == "|" else 0 for c in match_line]
    q_bases = list(query[:n])
    r_bases = list(ref[:n])

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=positions,
            y=match_values,
            marker=dict(
                color=[TEAL if v else CORAL for v in match_values],
                line=dict(width=0),
            ),
            name="Match / Mismatch",
            hovertemplate=(
                "Position: %{x}<br>"
                "Query: %{customdata[0]}<br>"
                "Reference: %{customdata[1]}<extra></extra>"
            ),
            customdata=list(zip(q_bases, r_bases)),
            showlegend=True,
        )
    )

    layout = _base_layout("Sequence Alignment Map")
    layout["xaxis"]["title"] = "Position (bp)"
    layout["yaxis"]["title"] = "Match"
    layout["yaxis"]["tickvals"] = [0, 1]
    layout["yaxis"]["ticktext"] = ["Mismatch", "Match"]
    layout["showlegend"] = False
    fig.update_layout(**layout, height=250)
    return fig


def plot_alignment_overview(alignment_map: dict) -> go.Figure:
    """
    Full-length companion to plot_alignment().

    plot_alignment() only ever renders the first `max_chars` (default 60)
    positions of the alignment -- on a multi-thousand-bp query, any
    mismatch/gap beyond that window is invisible even though the metrics
    below (Mismatches: N, Gaps: N) correctly count it. Confirmed in
    practice: an 8-substitution test sequence (mutations at positions
    861-7066) rendered as a perfect, all-green 0-60bp bar chart, silently
    hiding every one of the 8 real differences.

    This renders every mismatch/gap position as a sparse marker across the
    FULL alignment length instead of one bar per base -- scales to any
    sequence length (a 60,000bp alignment with 20 mismatches draws 20
    points, not 60,000 bars), and is meant to sit alongside
    plot_alignment()'s zoomed view, not replace it: this answers "where in
    the whole sequence are the differences", the zoomed bar chart answers
    "what exactly changed at this one position".
    """
    query = alignment_map.get("query", "")
    ref = alignment_map.get("reference", "")
    match_line = alignment_map.get("match_line", "")

    if not match_line:
        return go.Figure()

    total_len = len(match_line)
    mismatch_x, mismatch_hover = [], []
    gap_x, gap_hover = [], []

    for i, (m, q, r) in enumerate(zip(match_line, query, ref)):
        pos = i + 1
        if m == " ":  # gap column (query or reference has "-")
            gap_x.append(pos)
            gap_hover.append(f"Position: {pos}<br>Query: {q}<br>Reference: {r}<br>Gap")
        elif m == "X":  # mismatch, no gap
            mismatch_x.append(pos)
            mismatch_hover.append(f"Position: {pos}<br>Query: {q}<br>Reference: {r}<br>Mismatch")

    fig = go.Figure()

    # Baseline showing the full alignment span, so an empty result (perfect
    # match everywhere) still renders an informative "nothing to see here"
    # line rather than a blank chart.
    fig.add_trace(
        go.Scatter(
            x=[1, total_len],
            y=[0, 0],
            mode="lines",
            line=dict(color=THEME["grid_color"], width=2),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    if mismatch_x:
        fig.add_trace(
            go.Scatter(
                x=mismatch_x,
                y=[0] * len(mismatch_x),
                mode="markers",
                marker=dict(color=CORAL, size=9, symbol="line-ns", line=dict(width=2, color=CORAL)),
                name="Mismatch",
                hovertext=mismatch_hover,
                hoverinfo="text",
            )
        )
    if gap_x:
        fig.add_trace(
            go.Scatter(
                x=gap_x,
                y=[0] * len(gap_x),
                mode="markers",
                marker=dict(color=AMBER, size=9, symbol="line-ns", line=dict(width=2, color=AMBER)),
                name="Gap",
                hovertext=gap_hover,
                hoverinfo="text",
            )
        )

    layout = _base_layout(f"Full-Length Difference Map ({total_len:,} bp/aa)")
    layout["xaxis"]["title"] = "Position (bp)"
    layout["xaxis"]["range"] = [0, total_len + 1]
    layout["yaxis"]["visible"] = False
    layout["yaxis"]["range"] = [-1, 1]
    layout["showlegend"] = bool(mismatch_x or gap_x)
    fig.update_layout(**layout, height=150)
    if not mismatch_x and not gap_x:
        fig.add_annotation(
            text="No mismatches or gaps across the full alignment",
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False, font=dict(color=TEAL, size=12),
        )
    return fig


# ─── Mutation map ──────────────────────────────────────────────────────────────

def plot_mutation_map(
    mutation_report: dict,
    seq_length: int,
    labels: dict[str, str] | None = None,
) -> go.Figure:
    """
    Scatter plot of mutation positions along the sequence.

    Args:
        mutation_report: output of bioinformatics.detect_mutations()
        seq_length:      length of the compared region
    """
    labels = labels or {
        "title": "Mutation Map",
        "transition": "Transitions (yellow)",
        "transversion": "Transversions (red)",
        "indel": "Indels (gray)",
        "legend_title": "Substitution type",
        "xaxis": "Position (bp)",
        "position": "Position",
        "reference_position": "Reference",
        "query_position": "Query",
        "change": "Change",
        "type": "Type",
        "event": "Event",
        "bases": "Bases",
        "consequence": "Consequence",
        "frameshift": "Frameshift",
        "in_frame": "In-frame",
    }
    mutations = mutation_report.get("mutations", [])
    indels = mutation_report.get("indel_blocks", mutation_report.get("indels", []))

    if not mutations and not indels:
        fig = go.Figure()
        fig.update_layout(
            **_base_layout("Mutation Map"),
            annotations=[
                dict(
                    text="No mutations detected",
                    x=0.5, y=0.5,
                    xref="paper", yref="paper",
                    showarrow=False,
                    font=dict(color=MINT, size=16),
                )
            ],
        )
        return fig

    grouped_mutations = {
        "transition": {"positions": [], "hover_texts": []},
        "transversion": {"positions": [], "hover_texts": []},
    }
    for m in mutations:
        query_position = m.get("position_query")
        reference_position = m.get("position_reference")
        position = query_position or reference_position
        if position is None:
            continue

        mutation_type = str(m.get("type", "transversion")).lower()
        group = "transition" if mutation_type == "transition" else "transversion"
        position_label = (
            f"{labels['reference_position']} {reference_position} / {labels['query_position']} {query_position}"
            if reference_position is not None and query_position is not None
            else f"{labels['position']} {position}"
        )
        grouped_mutations[group]["positions"].append(position)
        grouped_mutations[group]["hover_texts"].append(
            f"<b>{position_label}</b><br>"
            f"{labels['change']} : {m.get('reference', '?')} → {m.get('query', '?')}<br>"
            f"{labels['type']} : {group.capitalize()}"
        )

    indel_positions = []
    indel_hover_texts = []
    for indel in indels:
        position = indel.get("start_position_query") or indel.get("start_position_reference")
        if position is None:
            continue
        indel_positions.append(position)
        indel_type = str(indel.get("type", "indel")).capitalize()
        status = labels["frameshift"] if indel.get("frameshift") else labels["in_frame"]
        indel_hover_texts.append(
            f"<b>{labels['position']} {position}</b><br>"
            f"{labels['event']} : {indel_type}<br>"
            f"{labels['bases']} : {indel.get('bases', '?')} ({indel.get('length', 1)} bp)<br>"
            f"{labels['consequence']} : {status}"
        )

    if not any(item["positions"] for item in grouped_mutations.values()) and not indel_positions:
        fig = go.Figure()
        fig.update_layout(
            **_base_layout("Mutation Map"),
            annotations=[
                dict(
                    text="No mutation positions available",
                    x=0.5, y=0.5,
                    xref="paper", yref="paper",
                    showarrow=False,
                    font=dict(color=AMBER, size=16),
                )
            ],
        )
        return fig

    fig = go.Figure()

    for group, color, label in (
        ("transition", AMBER, labels["transition"]),
        ("transversion", CORAL, labels["transversion"]),
    ):
        positions = grouped_mutations[group]["positions"]
        if not positions:
            continue
        fig.add_trace(
            go.Scatter(
                x=positions,
                y=[1] * len(positions),
                mode="markers",
                marker=dict(
                    size=14,
                    color=color,
                    symbol="diamond",
                    line=dict(color="#0d1b2a", width=1),
                ),
                text=grouped_mutations[group]["hover_texts"],
                hovertemplate="%{text}<extra></extra>",
                name=label,
            )
        )

    if indel_positions:
        fig.add_trace(
            go.Scatter(
                x=indel_positions,
                y=[0.72] * len(indel_positions),
                mode="markers",
                marker=dict(
                    size=16,
                    color=SLATE,
                    symbol="triangle-up",
                    line=dict(color="#0d1b2a", width=1),
                ),
                text=indel_hover_texts,
                hovertemplate="%{text}<extra></extra>",
                name=labels["indel"],
            )
        )

    fig.add_trace(
        go.Scatter(
            x=[0, seq_length],
            y=[1, 1],
            mode="lines",
            line=dict(color=TEAL, width=2),
            name="Sequence",
            hoverinfo="skip",
        )
    )

    layout = _base_layout(
        f"{labels['title']} — {len(mutations)} substitution(s), {len(indels)} indel(s)"
    )
    layout["xaxis"]["title"] = labels["xaxis"]
    layout["yaxis"]["visible"] = False
    layout["yaxis"]["range"] = [0.45, 1.2]
    layout["margin"] = dict(l=55, r=190, t=55, b=45)
    layout["showlegend"] = True
    layout["legend"] = dict(
        title=labels["legend_title"],
        orientation="v",
        y=1,
        x=1.02,
        xanchor="left",
        yanchor="top",
    )
    fig.update_layout(**layout, height=220)
    return fig


# ─── GC sliding window ────────────────────────────────────────────────────────

def plot_gc_sliding_window(sequence: str, window: int = 20) -> go.Figure:
    """
    Plot GC content along the sequence using a sliding window.

    Args:
        sequence: cleaned DNA string
        window:   window size in bp
    """
    if len(sequence) < window:
        return go.Figure()

    positions: list[int] = []
    gc_values: list[float] = []

    for i in range(0, len(sequence) - window + 1, max(1, window // 4)):
        chunk = sequence[i : i + window]
        gc = (chunk.count("G") + chunk.count("C")) / len(chunk) * 100
        positions.append(i + window // 2)
        gc_values.append(round(gc, 2))

    fig = go.Figure()

    fig.add_hrect(y0=65, y1=100, fillcolor="rgba(255,209,102,0.08)", line_width=0)
    fig.add_hrect(y0=0, y1=35, fillcolor="rgba(79,195,247,0.08)", line_width=0)
    fig.add_hrect(y0=35, y1=65, fillcolor="rgba(0,217,163,0.06)", line_width=0)

    fig.add_trace(
        go.Scatter(
            x=positions,
            y=gc_values,
            mode="lines",
            line=dict(color=TEAL, width=2),
            fill="tozeroy",
            fillcolor="rgba(0,217,163,0.12)",
            name=f"GC% (w={window}bp)",
            hovertemplate="Position %{x}<br>GC: %{y:.1f}%<extra></extra>",
        )
    )

    fig.add_hline(y=50, line=dict(color=AMBER, dash="dash", width=1))

    layout = _base_layout(f"GC Content Profile (window = {window} bp)")
    layout["xaxis"]["title"] = "Position (bp)"
    layout["yaxis"]["title"] = "GC (%)"
    layout["yaxis"]["range"] = [0, 100]
    layout["yaxis"]["ticksuffix"] = "%"
    fig.update_layout(
        **layout,
        height=320,
        annotations=[
            dict(
                text="Shaded bands are reference ranges, not a pass/fail score.",
                x=0.5, y=1.12,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=10, color=SLATE),
            )
        ],
    )
    return fig


def plot_gc_skew_profile(profile: list[dict]) -> go.Figure:
    """Plot GC and AT skew values around the zero reference line."""
    fig = go.Figure()
    positions = [item["position"] for item in profile]
    fig.add_trace(go.Scatter(
        x=positions,
        y=[item["gc_skew"] for item in profile],
        mode="lines",
        name="GC skew",
        line=dict(color=TEAL, width=2),
        connectgaps=False,
    ))
    fig.add_trace(go.Scatter(
        x=positions,
        y=[item["at_skew"] for item in profile],
        mode="lines",
        name="AT skew",
        line=dict(color=AMBER, width=2),
        connectgaps=False,
    ))
    fig.add_hline(y=0, line=dict(color=SLATE, dash="dash", width=1))
    layout = _base_layout("GC / AT Skew Profile")
    layout["xaxis"]["title"] = "Position (bp)"
    layout["yaxis"]["title"] = "Skew"
    layout["yaxis"]["range"] = [-1, 1]
    fig.update_layout(**layout, height=300)
    return fig


# ─── Multiple Sequence Alignment table (colored) ───────────────────────────
def plot_msa_table(aligned_sequences: list, labels: list | None = None) -> go.Figure:
    """
    Render an MSA as a colorized table using Plotly `go.Table`.

    Args:
        aligned_sequences: list of aligned sequence strings (same length)
        labels: optional list of row labels (sequence names)
    """
    if not aligned_sequences:
        return go.Figure()

    align_len = len(aligned_sequences[0])
    # Normalize sequences to same length
    rows = [list(s.ljust(align_len, '-')) for s in aligned_sequences]

    # Build color map per cell
    fill_colors = []
    for row in rows:
        colors = [NUCLEOTIDE_COLORS.get(base.upper(), SLATE) for base in row]
        fill_colors.append(colors)

    # Build header and cell values: show columns as positions
    header_values = [f"Pos {i+1}" for i in range(align_len)]
    # transpose rows -> columns for go.Table cells expects list of columns
    cell_values = list(map(list, zip(*rows)))
    cell_colors = list(map(list, zip(*fill_colors)))

    header_row_height = 28
    cell_row_height = 24

    fig = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=["Sequence"] + header_values,
                    fill_color="#0d1b2a",
                    align="center",
                    font=dict(color=TEAL, size=12),
                    line_color="rgba(0,217,163,0.25)",
                    height=header_row_height,
                ),
                cells=dict(
                    values=[labels or [f"Seq {i+1}" for i in range(len(rows))]] + cell_values,
                    fill_color=[["#0d1b2a"] * len(cell_values[0])] * 1 + cell_colors,
                    align="center",
                    font=dict(color="#061019", size=11),
                    line_color="rgba(0,217,163,0.15)",
                    height=cell_row_height,
                ),
            )
        ]
    )
    # Explicit row heights above make this deterministic instead of leaving
    # Plotly's own default table-row sizing to disagree with whatever height
    # we set here -- with many columns (even after windowing to a
    # reasonable width upstream), a mismatch between the two previously
    # left the actual colored cells squeezed into a thin strip inside a
    # much taller, mostly-blank figure.
    total_height = 40 + header_row_height + cell_row_height * len(rows)
    fig.update_layout(paper_bgcolor=THEME["paper"], height=max(200, total_height))
    return fig


# ─── Dendrogram from scipy dendrogram-data ──────────────────────────────────
def plot_dendrogram(dendro: dict, labels: list | None = None) -> go.Figure:
    """
    Build a simple interactive tree diagram from scipy-style dendrogram data.

    Args:
        dendro: dict containing 'icoord' and 'dcoord' lists (as returned by scipy.dendrogram)
        labels: optional list of leaf labels in order
    """
    if not dendro or 'icoord' not in dendro:
        return go.Figure()

    icoord = dendro.get('icoord', [])
    dcoord = dendro.get('dcoord', [])
    color_list = dendro.get('color_list', [])

    fig = go.Figure()
    node_points = set()
    leaf_x = []

    for xs, ys, col in zip(icoord, dcoord, color_list if color_list else [TEAL] * len(icoord)):
        xs = list(xs)
        ys = list(ys)
        for x, y in zip(xs, ys):
            node_points.add((x, y))
            if y == 0:
                leaf_x.append(x)

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode='lines',
                line=dict(color=_normalize_plotly_color(col), width=3),
                hoverinfo='none',
                showlegend=False,
            )
        )

    internal_nodes = sorted({(x, y) for x, y in node_points if y != 0}, key=lambda p: (p[1], p[0]))
    leaves = sorted({(x, y) for x, y in node_points if y == 0}, key=lambda p: p[0])

    if internal_nodes:
        fig.add_trace(
            go.Scatter(
                x=[x for x, _ in internal_nodes],
                y=[y for _, y in internal_nodes],
                mode='markers',
                marker=dict(size=8, color=CYAN, symbol='circle'),
                hoverinfo='none',
                showlegend=False,
            )
        )

    if leaves:
        fig.add_trace(
            go.Scatter(
                x=[x for x, _ in leaves],
                y=[y for _, y in leaves],
                mode='markers',
                marker=dict(size=10, color=MINT, symbol='circle'),
                hoverinfo='none',
                showlegend=False,
            )

        )

    leaf_x_unique = sorted(set(leaf_x))
    if labels and len(labels) == len(leaf_x_unique):
        fig.add_trace(
            go.Scatter(
                x=leaf_x_unique,
                y=[0] * len(leaf_x_unique),
                mode='text',
                text=labels,
                textposition='bottom center',
                textfont=dict(color=THEME["font_color"], size=11),
                hoverinfo='skip',
                showlegend=False,
            )
        )

    fig.update_layout(
        **_base_layout("Phylogenetic Tree"),
        xaxis=dict(showticklabels=False, zeroline=False, showgrid=False),
        yaxis=dict(title="Distance", zeroline=False, showgrid=True),
        height=420,
    )
    return fig


# ─── Enhanced Similarity Analysis (Phases 1-5) ────────────────────────────────

def build_match_context_card(match: dict) -> dict:
    """
    (1) Biological context: extract gene metadata for display.
    
    Returns a dict with gene name, trait, organism, description, and source
    so the UI can show "what is this gene and why should I care?"
    """
    return {
        "gene_name": match.get("gene_name", "N/A"),
        "trait": match.get("trait", "Unknown trait"),
        "organism": match.get("organism", "Unknown organism"),
        "description": match.get("description", "(No description available)"),
        "accession": match.get("accession", "N/A"),
        "source": match.get("source", "Unknown source"),
    }


def build_similarity_metrics_table(match: dict, query_len: int) -> dict:
    """
    (2) Detailed alignment metrics: gaps%, coverage%, matches/mismatches.
    
    Extracts alignment statistics to answer: "Is this alignment real or noise?"
    High coverage + low gaps = trustworthy match.
    """
    alignment = match.get("alignment", {})
    if not isinstance(alignment, dict):
        return {}

    # similarityengine.compare_with_database nests the aligned strings (and
    # precomputed match/mismatch/gap counts) inside alignment_map -- there
    # is no top-level "seq1_aligned"/"seq2_aligned" key on `alignment`
    # itself. Reading alignment_map also lets us reuse the counts
    # alignment_engine.alignment_statistics() already computed instead of
    # re-deriving them from the strings a second time.
    alignment_map = alignment.get("alignment_map", {})
    if not isinstance(alignment_map, dict):
        return {}

    aligned_query = alignment_map.get("query", "")
    aligned_reference = alignment_map.get("reference", "")

    if not aligned_query or not aligned_reference:
        return {}

    aln_len = len(aligned_query)
    matches = alignment_map.get("identity_count", 0)
    mismatches = alignment_map.get("mismatch_count", 0)
    gaps = alignment_map.get("gap_count", 0)
    coverage = (aln_len - gaps) / query_len * 100 if query_len > 0 else 0

    return {
        "alignment_length": aln_len,
        "matches": matches,
        "mismatches": mismatches,
        "total_gaps": gaps,
        "gap_percent": (gaps / aln_len * 100) if aln_len > 0 else 0,
        "coverage_percent": coverage,
        "identity_percent": match.get("similarity_score", 0),
    }


def plot_alignment_coverage_heatmap(match: dict, query_len: int, window: int = 50) -> go.Figure:
    """
    (3) Heatmap: identity percent in sliding windows across alignment.
    
    Shows which regions are highly conserved vs variable.
    Helps identify if similarity is uniform or concentrated in domains.
    """
    alignment = match.get("alignment", {})
    alignment_map = alignment.get("alignment_map", {}) if isinstance(alignment, dict) else {}
    seq1 = alignment_map.get("query", "") if isinstance(alignment_map, dict) else ""
    seq2 = alignment_map.get("reference", "") if isinstance(alignment_map, dict) else ""
    
    if not seq1 or not seq2 or len(seq1) < window:
        fig = go.Figure()
        fig.update_layout(**_base_layout("Alignment Coverage"))
        return fig
    
    positions = []
    identities = []
    
    for i in range(0, len(seq1) - window + 1, max(1, window // 2)):
        chunk1 = seq1[i:i+window]
        chunk2 = seq2[i:i+window]
        matches = sum(1 for a, b in zip(chunk1, chunk2) if a == b and a != "-")
        identity = matches / len(chunk1) * 100 if chunk1 else 0
        positions.append(i + window // 2)
        identities.append(identity)
    
    if not positions:
        fig = go.Figure()
        fig.update_layout(**_base_layout("Alignment Coverage"))
        return fig
    
    fig = go.Figure(
        go.Bar(
            x=positions,
            y=identities,
            marker=dict(
                color=identities,
                colorscale=[[0, CORAL], [0.5, AMBER], [1, TEAL]],
                showscale=True,
                colorbar=dict(title=dict(text="Identity %", side="right")),
            ),
            text=[f"{v:.0f}%" for v in identities],
            textposition="outside",
            hovertemplate="Position: %{x}<br>Identity: %{y:.1f}%<extra></extra>",
        )
    )
    layout = _base_layout(f"Alignment Coverage — Identity by Window ({window}bp)")
    layout["xaxis"]["title"] = "Position in Alignment"
    layout["yaxis"]["title"] = "Local Identity (%)"
    layout["yaxis"]["range"] = [0, 105]
    fig.update_layout(**layout, height=300)
    return fig


def plot_confidence_gauge(metrics: dict) -> go.Figure:
    """
    Primary number is raw Identity, NOT a composite.

    Previously the gauge's headline number was a weighted composite
    (coverage*0.4 + gap_score*0.35 + identity*0.25) labeled "Match
    Confidence", with the actual identity percent relegated to 10pt
    caption text below the gauge. On a real 8-substitution test case
    (coverage 100%, gaps 0%, identity 99.89%) the composite rounded to a
    dominant "100%" -- visually indistinguishable from a truly identical
    sequence, even though 8 real differences exist. This inverts that:
    identity (the metric the rest of the app already treats as the
    primary, ranking-relevant number -- see similarity_score/"Similarity
    (global)" elsewhere) is now the number the gauge actually displays.
    The composite is kept as a secondary, explicitly-labeled figure for
    users who want a single "should I trust this alignment" heuristic,
    but it can never again be mistaken for the identity score itself.
    """
    coverage = metrics.get("coverage_percent", 0)
    gap_percent = metrics.get("gap_percent", 0)
    identity = metrics.get("identity_percent", 0)

    coverage_score = min(100, coverage * 1.1)
    gap_score = max(0, 100 - gap_percent * 5)
    identity_score = identity
    composite = min(100, max(0, coverage_score * 0.4 + gap_score * 0.35 + identity_score * 0.25))

    color = TEAL if identity >= 90 else AMBER if identity >= 70 else CORAL

    fig = go.Figure(
        go.Indicator(
            mode="gauge+number",
            value=identity,
            number=dict(suffix="%", font=dict(color=color, size=32)),
            gauge=dict(
                axis=dict(range=[0, 100], tickcolor=THEME["font_color"]),
                bar=dict(color=color, thickness=0.25),
                steps=[
                    dict(range=[0, 70], color="rgba(255,107,107,0.1)"),
                    dict(range=[70, 90], color="rgba(255,209,102,0.1)"),
                    dict(range=[90, 100], color="rgba(0,217,163,0.1)"),
                ],
            ),
            title=dict(text="Match Identity", font=dict(color=color, size=14)),
        )
    )
    fig.update_layout(
        paper_bgcolor=THEME["paper"],
        font=dict(color=THEME["font_color"]),
        margin=dict(l=20, r=20, t=60, b=50),
        height=280,
        annotations=[
            dict(
                text=(
                    f"Coverage: {coverage:.0f}% | Gaps: {gap_percent:.1f}% | "
                    f"Composite confidence*: {composite:.0f}%"
                ),
                x=0.5, y=-0.2,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=10, color=SLATE),
            ),
            dict(
                text="*weighted blend of coverage, gaps and identity — not a substitute for Identity above",
                x=0.5, y=-0.32,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=8, color=SLATE),
            ),
        ],
    )
    return fig


def build_top3_comparison_table(similarity_results: list[dict], query_len: int) -> dict:
    """
    (4) Side-by-side comparison of top 3 matches.
    
    Allows user to see if all matches agree on gene function (high confidence)
    or if results are scattered (low confidence / need deep search).
    
    Returns a dict with 'header' and 'rows' suitable for Streamlit table display.
    """
    table_data = {
        "header": ["Rank", "Gene", "Similarity", "Trait", "Organism", "Global coverage", "Gaps"],
        "rows": []
    }
    
    for rank, match in enumerate(similarity_results[:3], 1):
        metrics = build_similarity_metrics_table(match, query_len)
        coverage = metrics.get("coverage_percent", 0)
        gaps = metrics.get("gap_percent", 0)
        
        table_data["rows"].append({
            "rank": rank,
            "gene": match.get("gene_name", "N/A"),
            "similarity": f"{match.get('similarity_score', 0):.1f}%",
            "trait": match.get("trait", "Unknown")[:40],
            "organism": match.get("organism", "Unknown")[:20],
            "coverage": f"{coverage:.0f}%",
            "gaps": f"{gaps:.1f}%",
        })
    
    return table_data