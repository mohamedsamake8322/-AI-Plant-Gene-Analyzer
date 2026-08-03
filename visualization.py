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
        xaxis=dict(gridcolor=THEME["grid_color"], zerolinecolor=THEME["grid_color"], linecolor=THEME["line_color"]),
        yaxis=dict(gridcolor=THEME["grid_color"], zerolinecolor=THEME["grid_color"], linecolor=THEME["line_color"]),
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

def plot_gc_gauge(gc_percent: float) -> go.Figure:
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
                    dict(range=[0, 35], color="rgba(79,195,247,0.15)"),
                    dict(range=[35, 65], color="rgba(0,217,163,0.15)"),
                    dict(range=[65, 100], color="rgba(255,209,102,0.15)"),
                ],
                threshold=dict(
                    line=dict(color=AMBER, width=3),
                    thickness=0.8,
                    value=50,
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

    genes = [r["gene_name"] for r in similarity_results]
    scores = [r["similarity_score"] for r in similarity_results]
    traits = [r["trait"] for r in similarity_results]

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
            customdata=traits,
            hovertemplate=(
                "<b>%{y}</b><br>"
                "Trait: %{customdata}<br>"
                "Similarity: %{x:.1f}%<extra></extra>"
            ),
        )
    )
    layout = _base_layout("Database Similarity Scores")
    layout["xaxis"]["title"] = "Similarity (%)"
    layout["xaxis"]["range"] = [0, 105]
    layout["xaxis"]["ticksuffix"] = "%"
    layout["yaxis"]["title"] = "Gene"
    fig.update_layout(**layout, height=max(300, len(genes) * 60))
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


# ─── Mutation map ──────────────────────────────────────────────────────────────

def plot_mutation_map(mutation_report: dict, seq_length: int) -> go.Figure:
    """
    Scatter plot of mutation positions along the sequence.

    Args:
        mutation_report: output of bioinformatics.detect_mutations()
        seq_length:      length of the compared region
    """
    mutations = mutation_report.get("mutations", [])

    if not mutations:
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

    positions = []
    colors = []
    hover_texts = []
    for m in mutations:
        position = m.get("position")
        if position is None:
            position = m.get("position_query", m.get("position_reference"))
        if position is None:
            continue

        positions.append(position)
        mutation_type = m.get("type", "unknown")
        colors.append(AMBER if mutation_type == "transition" else CORAL)
        hover_texts.append(
            f"Pos {position}: {m.get('reference', '?')} → {m.get('query', '?')} ({mutation_type})"
        )

    if not positions:
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

    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1] * len(positions),
            mode="markers",
            marker=dict(
                size=14,
                color=colors,
                symbol="diamond",
                line=dict(color="#0d1b2a", width=1),
            ),
            text=hover_texts,
            hovertemplate="%{text}<extra></extra>",
            name="Mutations",
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
        f"Mutation Map — {len(mutations)} mutation(s) detected"
    )
    layout["xaxis"]["title"] = "Position (bp)"
    layout["yaxis"]["visible"] = False
    layout["showlegend"] = True
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

    fig = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=["Sequence"] + header_values,
                    fill_color="#0d1b2a",
                    align="center",
                    font=dict(color=TEAL, size=12),
                    line_color="rgba(0,217,163,0.25)",
                ),
                cells=dict(
                    values=[labels or [f"Seq {i+1}" for i in range(len(rows))]] + cell_values,
                    fill_color=[["#0d1b2a"] * len(cell_values[0])] * 1 + cell_colors,
                    align="center",
                    font=dict(color="#061019", size=11),
                    line_color="rgba(0,217,163,0.15)",
                ),
            )
        ]
    )
    fig.update_layout(paper_bgcolor=THEME["paper"], height=max(200, 40 * len(rows)))
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