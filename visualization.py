"""
visualization.py
----------------
Visualization module for the AI-Powered Plant Gene Analyzer.
Produces Plotly and Matplotlib figures from analysis results.
All functions return Plotly figures (compatible with st.plotly_chart).
"""

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

import config


# ─── Color palette ─────────────────────────────────────────────────────────────
NUCLEOTIDE_COLORS: dict[str, str] = config.NUCLEOTIDE_COLORS


def _base_layout(title: str = "") -> dict:
    """Shared Plotly layout for the app theme."""
    return dict(
        title=dict(text=title, font=dict(color="#0f3d0f", size=16)),
        paper_bgcolor=config.CHART_PAPER,
        plot_bgcolor=config.CHART_BG,
        font=dict(color=config.CHART_FONT_COLOR, family="Arial"),
        margin=dict(l=40, r=40, t=60, b=40),
        xaxis=dict(gridcolor=config.CHART_GRID_COLOR, zerolinecolor=config.CHART_GRID_COLOR, linecolor="#cccccc"),
        yaxis=dict(gridcolor=config.CHART_GRID_COLOR, zerolinecolor=config.CHART_GRID_COLOR, linecolor="#cccccc"),
    )


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
    colors = [NUCLEOTIDE_COLORS.get(k, "#9e9e9e") for k in labels]

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
        legend=dict(font=dict(color=config.CHART_FONT_COLOR)),
        annotations=[
            dict(
                text="DNA<br>Base",
                x=0.5, y=0.5,
                font=dict(size=13, color="#69f0ae"),
                showarrow=False,
            )
        ],
    )
    return fig


def plot_nucleotide_bar(dist: dict) -> go.Figure:
    """
    Bar chart of nucleotide counts.
    """
    counts = dist["counts"]
    nucleotides = list(counts.keys())
    values = [counts[n] for n in nucleotides]
    colors = [NUCLEOTIDE_COLORS.get(n, "#9e9e9e") for n in nucleotides]

    fig = go.Figure(
        go.Bar(
            x=nucleotides,
            y=values,
            marker=dict(color=colors, line=dict(color="#0d1b2a", width=1)),
            text=values,
            textposition="outside",
            textfont=dict(color=config.CHART_FONT_COLOR),
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
    colors = ["#00c853" if aa in {"A", "G", "V", "L", "I", "M"} else "#2979ff" for aa in residues]

    fig = go.Figure(
        go.Bar(
            x=residues,
            y=values,
            marker=dict(color=colors, line=dict(color="#0d1b2a", width=1)),
            text=values,
            textposition="outside",
            textfont=dict(color=config.CHART_FONT_COLOR),
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
    """
    if gc_percent >= 60:
        bar_color = "#ff3d00"
    elif gc_percent <= 35:
        bar_color = "#2979ff"
    else:
        bar_color = "#00c853"

    fig = go.Figure(
        go.Indicator(
            mode="gauge+number+delta",
            value=gc_percent,
            number=dict(suffix="%", font=dict(color="#69f0ae", size=36)),
            delta=dict(
                reference=50,
                increasing=dict(color="#ff3d00"),
                decreasing=dict(color="#2979ff"),
            ),
            gauge=dict(
                axis=dict(
                    range=[0, 100],
                    tickcolor=config.CHART_FONT_COLOR,
                    tickfont=dict(color=config.CHART_FONT_COLOR),
                ),
                bar=dict(color=bar_color, thickness=0.25),
                bgcolor=config.CHART_BG,
                bordercolor="#69f0ae",
                steps=[
                    dict(range=[0, 35], color="rgba(41,121,255,0.15)"),
                    dict(range=[35, 65], color="rgba(0,200,83,0.15)"),
                    dict(range=[65, 100], color="rgba(255,61,0,0.15)"),
                ],
                threshold=dict(
                    line=dict(color="#ffd600", width=3),
                    thickness=0.8,
                    value=50,
                ),
            ),
            title=dict(
                text="GC Content",
                font=dict(color="#69f0ae", size=16),
            ),
        )
    )
    fig.update_layout(
        paper_bgcolor=config.CHART_PAPER,
        font=dict(color=config.CHART_FONT_COLOR, family="Courier New"),
        margin=dict(l=20, r=20, t=60, b=20),
        height=280,
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
        "#00c853" if s >= 75 else "#ffd600" if s >= 50 else "#ff3d00"
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
            textfont=dict(color=config.CHART_FONT_COLOR),
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
                color=["#00c853" if v else "#ff3d00" for v in match_values],
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
                    font=dict(color="#69f0ae", size=16),
                )
            ],
        )
        return fig

    positions = [m["position"] for m in mutations]
    colors = [
        "#ffd600" if m["type"] == "transition" else "#ff3d00"
        for m in mutations
    ]
    hover_texts = [
        f"Pos {m['position']}: {m['reference']} → {m['query']} ({m['type']})"
        for m in mutations
    ]

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
            line=dict(color="#69f0ae", width=2),
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

    fig.add_hrect(y0=60, y1=100, fillcolor="rgba(255,61,0,0.1)", line_width=0)
    fig.add_hrect(y0=0, y1=35, fillcolor="rgba(41,121,255,0.1)", line_width=0)
    fig.add_hrect(y0=35, y1=60, fillcolor="rgba(0,200,83,0.05)", line_width=0)

    fig.add_trace(
        go.Scatter(
            x=positions,
            y=gc_values,
            mode="lines",
            line=dict(color="#00c853", width=2),
            fill="tozeroy",
            fillcolor="rgba(0,200,83,0.1)",
            name=f"GC% (w={window}bp)",
            hovertemplate="Position %{x}<br>GC: %{y:.1f}%<extra></extra>",
        )
    )

    fig.add_hline(y=50, line=dict(color="#ffd600", dash="dash", width=1))

    layout = _base_layout(f"GC Content Profile (window = {window} bp)")
    layout["xaxis"]["title"] = "Position (bp)"
    layout["yaxis"]["title"] = "GC (%)"
    layout["yaxis"]["range"] = [0, 100]
    layout["yaxis"]["ticksuffix"] = "%"
    fig.update_layout(**layout, height=300)
    return fig
