"""phylo_view.py - professional, self-explanatory phylogeny output (Plotly + Streamlit).

Usage in the app (replaces the current tree plots for IQ-TREE / UPGMA / NJ):

    from phylo_view import render_phylo_result
    render_phylo_result(newick, method="iqtree",            # "iqtree" | "upgma" | "nj"
                        meta={"model": "F81+F", "bootstrap": 1000, "n_sites": 619,
                              "aligner": "MAFFT v7.505 (--auto)", "engine": "IQ-TREE v2.0.7"},
                        dist_names=names, dist_matrix=matrix)   # last two are optional

Only dependencies: plotly (+ streamlit for render_phylo_result). Newick is parsed here,
so quoted names ('cqi:1 | organism=...'), internal support labels and negative lengths work.
"""
from __future__ import annotations

import html
import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import plotly.graph_objects as go

# --------------------------------------------------------------------------- Newick
@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    support: Optional[float] = None
    children: List["Node"] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0
    color: str = ""

    @property
    def is_leaf(self) -> bool:
        return not self.children


_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def parse_newick(text: str) -> Node:
    s = re.sub(r"\s*[\r\n]+\s*", "", text.strip()).rstrip(";")
    if not s:
        raise ValueError("Newick tree is empty.")
    pos = 0

    def label() -> str:
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            pos += 1
            out = []
            while pos < len(s):
                if s[pos] == "'":
                    if pos + 1 < len(s) and s[pos + 1] == "'":
                        out.append("'")
                        pos += 2
                        continue
                    pos += 1
                    break
                out.append(s[pos])
                pos += 1
            return "".join(out)
        start = pos
        while pos < len(s) and s[pos] not in ",():;":
            pos += 1
        return s[start:pos].strip()

    def node() -> Node:
        nonlocal pos
        if pos >= len(s):
            raise ValueError("Invalid Newick: unexpected end of input.")
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                n.children.append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                break
            if pos >= len(s) or s[pos] != ")":
                raise ValueError("Invalid Newick: missing ')'")
            pos += 1
            lab = label()
            if lab:
                if re.fullmatch(r"[\d.eE+\-/ ]+", lab) and _NUM.findall(lab):
                    n.support = float(_NUM.findall(lab)[-1])  # '95.2/100' -> UFBoot 100
                else:
                    n.name = lab
        else:
            n.name = label()
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
            n.length = float(s[start:pos])
        return n

    root = node()
    if pos != len(s):
        raise ValueError("Invalid Newick: unexpected content after the tree.")
    return root


def _leaves(n: Node) -> List[Node]:
    return [n] if n.is_leaf else [l for c in n.children for l in _leaves(c)]


def _internal(n: Node) -> List[Node]:
    return [] if n.is_leaf else [n] + [i for c in n.children for i in _internal(c)]


# --------------------------------------------------------------- midpoint rooting
def midpoint_root(root: Node) -> Node:
    """Root an unrooted tree (root with >= 3 children) at the middle of its longest leaf-leaf path."""
    nodes: List[Node] = []
    adj: Dict[int, List[Tuple[int, float, Optional[float]]]] = {}

    def walk(n: Node, parent: Optional[int], plen: float, psup: Optional[float]) -> int:
        nid = len(nodes)
        nodes.append(n)
        adj[nid] = []
        if parent is not None:
            adj[nid].append((parent, plen, psup))
            adj[parent].append((nid, plen, psup))
        for c in n.children:
            walk(c, nid, max(c.length, 0.0), c.support)
        return nid

    walk(root, None, 0.0, None)
    leaf_ids = [i for i, n in enumerate(nodes) if n.is_leaf]

    def dists(src: int):
        dist, prev, stack = {src: 0.0}, {src: None}, [src]
        while stack:
            u = stack.pop()
            for v, l, _ in adj[u]:
                if v not in dist:
                    dist[v], prev[v] = dist[u] + l, u
                    stack.append(v)
        return dist, prev

    best = (-1.0, leaf_ids[0], leaf_ids[0])
    for a in leaf_ids:
        d, _ = dists(a)
        for b in leaf_ids:
            if d[b] > best[0]:
                best = (d[b], a, b)
    total, a, b = best
    _, prev = dists(a)
    path = [b]
    while path[-1] != a:
        path.append(prev[path[-1]])
    path.reverse()

    half, acc = total / 2.0, 0.0
    u, v, elen, esup = path[0], path[1], 0.0, None
    for u, v in zip(path, path[1:]):
        elen, esup = next((l, s) for w, l, s in adj[u] if w == v)
        if acc + elen >= half - 1e-15:
            break
        acc += elen
    d1 = max(half - acc, 0.0)

    def build(nid: int, parent: int, length: float, support: Optional[float]) -> Node:
        src = nodes[nid]
        out = Node(name=src.name, length=length, support=support)
        out.children = [build(w, nid, l, s) for w, l, s in adj[nid] if w != parent]
        if len(out.children) == 1:
            child = out.children[0]
            child.length += out.length
            if out.support is not None:
                child.support = out.support
            return child
        return out

    new_root = Node()
    new_root.children = [build(u, v, d1, esup), build(v, u, max(elen - d1, 0.0), esup)]
    return new_root


def _ladderize(n: Node) -> int:
    size = 1 if n.is_leaf else sum(_ladderize(c) for c in n.children)
    n.children.sort(key=lambda c: len(_leaves(c)))
    return size


def _layout(root: Node) -> None:
    cursor = [0]

    def rec(n: Node, x: float) -> None:
        n.x = x
        if n.is_leaf:
            n.y = float(cursor[0])
            cursor[0] += 1
        else:
            for c in n.children:
                rec(c, x + max(c.length, 0.0))
            n.y = (min(c.y for c in n.children) + max(c.y for c in n.children)) / 2

    rec(root, 0.0)


def prepare_tree(newick: str, *, midpoint: bool = False) -> Tuple[Node, bool]:
    root = parse_newick(newick)
    rerooted = midpoint or len(root.children) >= 3
    if rerooted:
        root = midpoint_root(root)
    _ladderize(root)
    _layout(root)
    return root, rerooted


def leaf_order(newick: str, *, midpoint: bool = False) -> List[str]:
    """Leaf names from top to bottom of the drawn tree (use it to order a distance matrix)."""
    root, _ = prepare_tree(newick, midpoint=midpoint)
    return [l.name for l in _leaves(root)]


# ------------------------------------------------------------------------- insights
def split_label(name: str) -> Tuple[str, str]:
    """'cqi:110688118 | organism=Chenopodium quinoa' -> ('cqi:110688118', 'Chenopodium quinoa')"""
    m = re.search(r"organism=([^|;]+)", name)
    species = m.group(1).strip() if m else ""
    head = name.split("|")[0].strip() or name
    return head, species


def _patristic(root: Node):
    chains: Dict[int, List[Node]] = {}

    def rec(n: Node, chain: List[Node]) -> None:
        chain = chain + [n]
        if n.is_leaf:
            chains[id(n)] = chain
        for c in n.children:
            rec(c, chain)

    rec(root, [])
    leaves = _leaves(root)
    out = []
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):
            ca, cb = chains[id(leaves[i])], chains[id(leaves[j])]
            k = 0
            while k < min(len(ca), len(cb)) and ca[k] is cb[k]:
                k += 1
            lca = ca[k - 1]
            out.append((leaves[i].x + leaves[j].x - 2 * lca.x, leaves[i].name, leaves[j].name))
    return out


def _pretty_distance(name: str) -> str:
    return {"kimura": "Kimura 2-parameter (K2P)", "k2p": "Kimura 2-parameter (K2P)", "p": "p-distance",
            "p-distance": "p-distance", "jc": "Jukes-Cantor", "jukes-cantor": "Jukes-Cantor"}.get((name or "").strip().lower(), name or "K2P")


def _tied(rows: list, smallest: bool, tol: float = 0.02) -> list:
    ref = (min if smallest else max)(r[0] for r in rows)
    return [r for r in rows if abs(r[0] - ref) <= tol * ref]


def tree_insights(newick: str, support_cut: float = 95.0, dist_names: Optional[List[str]] = None,
                  dist_matrix: Optional[List[List[float]]] = None, method: Optional[str] = None,
                  fit_warn: float = 0.15, midpoint: bool = False) -> dict:
    """Summary of the tree. If the observed distance matrix is given, 'closest/farthest' use it (the tree only
    approximates the data) and the gap between tree and data is reported."""
    parsed = parse_newick(newick)
    has_negative = any(node.length < 0 for node in _leaves(parsed) + _internal(parsed))
    root, rerooted = prepare_tree(newick, midpoint=midpoint)
    leaves, internals = _leaves(root), _internal(root)
    n = len(leaves)
    ins: dict = {"n": n, "warnings": [], "closest": [], "farthest": [], "split": None, "weak": 0, "fit": None,
                 "source": "tree"}
    obs: Dict[frozenset, float] = {}
    if dist_names and dist_matrix and set(dist_names) == {l.name for l in leaves}:
        obs = {frozenset((dist_names[i], dist_names[j])): dist_matrix[i][j]
               for i in range(len(dist_names)) for j in range(i + 1, len(dist_names))}
        ins["source"] = "observed"
    if n <= 300:
        rows = [(obs.get(frozenset((a, b)), dt), a, b, dt, obs.get(frozenset((a, b)))) for dt, a, b in _patristic(root)]
        if rows:
            ins["closest"], ins["farthest"] = _tied(rows, True), _tied(rows, False)
        fits = [(abs(dt - do) / do, a, b, dt, do) for _, a, b, dt, do in rows if do and do > 1e-6]
        if fits:
            ins["fit"] = max(fits)
    if len(root.children) == 2:
        sup = next((c.support for c in root.children if c.support is not None), None)
        groups = [[l.name for l in _leaves(c)] for c in root.children]
        group_heights = {
            tuple(group): max(leaf.y for leaf in leaves if leaf.name in group)
            for group in groups
        }
        groups.sort(key=lambda group: group_heights[tuple(group)], reverse=True)
        ins["split"] = (groups, sup)
    skip = {id(c) for c in root.children} if rerooted else set()
    ins["weak"] = sum(1 for i in internals[1:] if i.support is not None and i.support < support_cut and id(i) not in skip)
    if rerooted and ins["split"] and ins["split"][1] is not None and ins["split"][1] < support_cut:
        ins["weak"] += 1
    has_support = any(i.support is not None for i in internals)

    if n < 4:
        ins["warnings"].append(f"Only {n} sequences: at least 4 are needed to say anything about branching order.")
    elif n == 4:
        ins["warnings"].append("Exactly 4 sequences: the unrooted tree has a single internal split (3 possible groupings). Read that split with the pairwise distances; add sequences and an outgroup for a robust inference.")
    if method in ("upgma", "nj") and not has_support:
        ins["warnings"].append("No branch support was computed for this method, so the reliability of each grouping is unknown. Use IQ-TREE (UFBoot) to quantify it.")
    if ins["weak"]:
        ins["warnings"].append(f"{ins['weak']} internal branch(es) have support below {support_cut:g}: treat their grouping as uncertain.")
    if ins["fit"] and ins["fit"][0] > fit_warn:
        r, a, b, dt, do = ins["fit"]
        hint = " UPGMA's constant-rate assumption is the usual cause: compare with IQ-TREE or Neighbor-Joining." if method == "upgma" else ""
        ins["warnings"].append(f"The tree does not reproduce the data well: up to {r*100:.0f}% gap for {split_label(a)[0]} - {split_label(b)[0]} "
                               f"({dt:.4f} in the tree vs {do:.4f} observed).{hint}")
    near_zero = [split_label(l.name)[0] for l in leaves if l.length < 1e-5]
    if near_zero:
        ins["warnings"].append("Terminal branch length is about 0 for " + ", ".join(near_zero)
                               + ": no private changes relative to the inferred ancestor of its group (identical or ancestral-like). Check for duplicates.")
    if has_negative:
        ins["warnings"].append("Negative branch length(s) in the input (typical NJ artefact on noisy distances); drawn as 0.")
    if ins["farthest"] and ins["farthest"][0][0] > 0.5:
        ins["warnings"].append("Some distances exceed 0.5 substitutions/site: sites may be saturated and the topology less reliable.")
    return ins


def _pairs_txt(rows: list, k: int = 3) -> str:
    out = [f"{html.escape(split_label(a)[0])} & {html.escape(split_label(b)[0])}" for _, a, b, *_ in rows[:k]]
    return "; ".join(out) + (f"; +{len(rows) - k} more" if len(rows) > k else "")


def findings(ins: dict, support_label: str = "support", dist_label: str = "distance",
             support_cut: float = 95.0) -> List[str]:
    src = f"observed {dist_label}" if ins.get("source") == "observed" else "distance along the tree"
    out = []
    if ins["closest"]:
        d = ins["closest"][0][0]
        many = len(ins["closest"]) > 1
        out.append(f"**Closest pair{'s (tied)' if many else ''}:** {_pairs_txt(ins['closest'])} - {d:.4f} substitutions/site "
                   f"(about {d*100:.1f} differences per 100 sites; {src}).")
    if ins["farthest"] and ins["n"] > 2:
        d = ins["farthest"][0][0]
        many = len(ins["farthest"]) > 1
        out.append(f"**Most divergent pair{'s (tied)' if many else ''}:** {_pairs_txt(ins['farthest'])} - {d:.4f} substitutions/site ({src}).")
    if ins["split"]:
        groups, sup = ins["split"]
        txt = " | ".join("{" + ", ".join(html.escape(split_label(x)[0]) for x in g) + "}" for g in groups)
        moderate_cut = 80.0 if support_cut >= 95.0 else 50.0
        tail = f" - {support_label} {sup:g} ({'reliable' if sup >= support_cut else 'moderate' if sup >= moderate_cut else 'weak'})" if sup is not None else ""
        out.append(f"**Main split:** {txt}{tail}.")
    if ins.get("fit"):
        r, a, b, dt, do = ins["fit"]
        out.append(f"**Fit to the data:** the worst gap between a distance read on the tree and the observed one is {r*100:.0f}% "
                   f"({split_label(a)[0]} & {split_label(b)[0]}); the smaller, the more faithful the tree.")
    return out


# ------------------------------------------------------------------------- figure
THEMES = {
    "dark": dict(bg="rgba(0,0,0,0)", solid="#171a14", fg="#E6EDF3", muted="#8B98A5", grid="rgba(255,255,255,0.07)",
                 branch="#9FB3C8", clades=["#4FC3F7", "#FFB74D"]),
    "light": dict(bg="#FFFFFF", solid="#FFFFFF", fg="#1F2933", muted="#6B7785", grid="#E5E9EF",
                  branch="#52606D", clades=["#1F77B4", "#D9822B"]),
}
SUPPORT_BANDS = [(95.0, "#2E9E6B", "&#8805; 95  reliable"), (80.0, "#E0A100", "80-94  moderate"), (-1e9, "#D64545", "&lt; 80  weak")]
FONT = "Inter, 'Segoe UI', Helvetica, Arial, sans-serif"
PLOTLY_CONFIG = dict(displaylogo=False, modeBarButtonsToRemove=["select2d", "lasso2d"],
                     toImageButtonOptions=dict(format="svg", filename="phylogeny", scale=2))


def _support_bands(support_label: str) -> list[tuple[float, str, str]]:
    if support_label.casefold() in {"ufboot", "ufboot2", "ultrafast bootstrap"}:
        return SUPPORT_BANDS
    return [(70.0, "#2E9E6B", "&#8805; 70  reliable"), (50.0, "#E0A100", "50-69  moderate"), (-1e9, "#D64545", "&lt; 50  weak")]


def _band(v: float, support_label: str = "UFBoot") -> str:
    return next(c for t, c, _ in _support_bands(support_label) if v >= t)


def _nice(x: float) -> float:
    if x <= 0:
        return 0.01
    e = 10 ** math.floor(math.log10(x))
    return next(m * e for m in (1, 2, 5, 10) if m * e >= x * 0.7)


def make_tree_figure(newick: str, *, theme: str = "dark", title: str = "Phylogenetic tree", subtitle: str = "",
                     axis_mode: str = "bar", support_label: str = "UFBoot", group_colors: bool = True,
                     midpoint: bool = False) -> go.Figure:
    """axis_mode: 'bar' = scale bar (additive trees)  |  'height' = node-height axis (ultrametric / UPGMA)."""
    P = THEMES[theme]
    root, rerooted = prepare_tree(newick, midpoint=midpoint)
    leaves, internals = _leaves(root), _internal(root)
    n = len(leaves)
    xmax = max(l.x for l in leaves) or 1e-9

    def paint(node: Node, color: str) -> None:
        node.color = color
        for c in node.children:
            paint(c, color)

    root.color = P["branch"]
    for i, c in enumerate(root.children):
        paint(c, P["clades"][i % 2] if (group_colors and len(root.children) == 2) else P["branch"])

    heads = [split_label(l.name) for l in leaves]
    species = {s for _, s in heads}
    common = species.pop() if len(species) == 1 and heads[0][1] else ""
    if common and not subtitle.endswith(common):
        subtitle = (subtitle + " &#183; " if subtitle else "") + f"<i>{common}</i>"

    fig = go.Figure()
    # ---- branches (grouped by colour -> few traces)
    segs: Dict[str, Tuple[list, list]] = {}
    for p in internals:
        for c in p.children:
            xs, ys = segs.setdefault(c.color, ([], []))
            xs += [p.x, c.x, None, p.x, p.x, None]          # horizontal branch + connector from parent centre
            ys += [c.y, c.y, None, p.y, c.y, None]
    for color, (xs, ys) in segs.items():
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=2.4), hoverinfo="skip", showlegend=False))

    # ---- tips
    def tip_text(l: Node) -> str:
        head, sp = split_label(l.name)
        head = html.escape(head)
        sp = html.escape(sp)
        return f"<b>{head}</b>" + (f"  <i style='color:{P['muted']}'>{sp}</i>" if sp and not common else "")

    fig.add_trace(go.Scatter(
        x=[l.x for l in leaves], y=[l.y for l in leaves], mode="markers+text",
        text=[tip_text(l) for l in leaves], textposition="middle right", textfont=dict(size=13, color=P["fg"], family=FONT),
        marker=dict(size=9, color=[l.color for l in leaves], line=dict(color=P["solid"], width=1.5)),
        hovertext=[f"<b>{html.escape(l.name)}</b><br>Terminal branch: {max(l.length, 0):.6f}<br>Distance from root: {l.x:.6f}" for l in leaves],
        hovertemplate="%{hovertext}<extra></extra>", showlegend=False, cliponaxis=False))

    # ---- support markers (one marker for the root split when the tree was midpoint-rooted)
    sup_pts: List[Tuple[Node, float]] = []
    skip: set = set()
    if rerooted and len(root.children) == 2 and all(c.support is not None for c in root.children) \
            and root.children[0].support == root.children[1].support:
        sup_pts.append((root, root.children[0].support))
        skip = {id(c) for c in root.children}
    sup_pts += [(i, i.support) for i in internals[1:] if i.support is not None and id(i) not in skip]
    if sup_pts:
        fig.add_trace(go.Scatter(
            x=[p.x for p, _ in sup_pts], y=[p.y for p, _ in sup_pts], mode="markers+text",
            text=[f"{v:g}" for _, v in sup_pts], textposition=["middle left" if p is root else "top left" for p, _ in sup_pts],
            textfont=dict(size=12, color=P["fg"], family=FONT),
            marker=dict(size=13, color=[_band(v, support_label) for _, v in sup_pts], line=dict(color=P["solid"], width=2)),
            hovertext=[f"{support_label} support: <b>{v:g}</b>" + ("<br>(root split, same edge on both sides)" if p is root else "") for p, v in sup_pts],
            hovertemplate="%{hovertext}<extra></extra>", showlegend=False))
        for _, color, name in _support_bands(support_label):
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=name, marker=dict(size=11, color=color),
                                     legendgroup="sup", legendgrouptitle_text=f"{support_label} support"))

    # ---- scale bar / axis
    max_chars = max(len(re.sub(r"<[^>]+>", "", tip_text(l))) for l in leaves)
    frac = min(0.55, max_chars * 7.2 / 950)
    x_hi = xmax / (1 - frac)
    y_lo = -2.4 if axis_mode == "bar" else -1.0
    xaxis = dict(visible=False, range=[-0.045 * xmax, x_hi], fixedrange=False)
    if axis_mode == "bar":
        bar = _nice(xmax / 4)
        fig.add_trace(go.Scatter(x=[0, bar, None, 0, 0, None, bar, bar], y=[-1.2, -1.2, None, -1.05, -1.35, None, -1.05, -1.35],
                                 mode="lines", line=dict(color=P["muted"], width=2), hoverinfo="skip", showlegend=False))
        fig.add_annotation(x=bar / 2, y=-1.85, text=f"{bar:g} substitutions / site", showarrow=False,
                           font=dict(size=12, color=P["muted"], family=FONT), xref="x", yref="y")
    else:
        step = _nice(xmax / 5)
        hs = [i * step for i in range(int(xmax / step) + 1)]
        xaxis = dict(visible=True, range=[-0.045 * xmax, x_hi], tickmode="array", tickvals=[xmax - h for h in hs],
                     ticktext=[f"{h:g}" for h in hs], showgrid=True, gridcolor=P["grid"], zeroline=False, showline=True,
                     linecolor=P["muted"], tickfont=dict(size=12, color=P["muted"], family=FONT),
                     title=dict(text="Node height (substitutions / site)", font=dict(size=12, color=P["muted"], family=FONT)))

    fig.update_layout(
        height=max(380, 150 + 52 * n), margin=dict(l=24, r=24, t=104, b=24 if axis_mode == "bar" else 56),
        paper_bgcolor=P["bg"], plot_bgcolor=P["bg"], font=dict(family=FONT, color=P["fg"]),
        title=dict(text=f"<b>{title}</b>" + (f"<br><sup><span style='color:{P['muted']}'>{subtitle}</span></sup>" if subtitle else ""),
                   x=0.01, xanchor="left", y=0.96, yanchor="top", font=dict(size=18, family=FONT)),
        xaxis=xaxis, yaxis=dict(visible=False, range=[y_lo, n - 0.4], fixedrange=True),
        legend=dict(orientation="h", x=0.99, xanchor="right", y=1.12, yanchor="bottom", font=dict(size=12, color=P["muted"]),
                    groupclick="toggleitem"), hoverlabel=dict(font_family=FONT), showlegend=bool(sup_pts))
    return fig


def make_distance_heatmap(names: List[str], matrix: List[List[float]], order: Optional[List[str]] = None,
                          theme: str = "dark", title: str = "Observed pairwise distances (substitutions / site)") -> go.Figure:
    P = THEMES[theme]
    idx = list(range(len(names)))
    if order and set(order) == set(names):
        idx = [names.index(o) for o in order]
    labs = [split_label(names[i])[0] for i in idx]
    z = [[None if a == b else matrix[i][j] for b, j in enumerate(idx)] for a, i in enumerate(idx)]
    text = [["" if v is None else f"{v:.4f}" for v in r] for r in z]
    scale = [[0, "#24504F"], [1, "#1E9486"]] if theme == "dark" else [[0, "#EAF4F4"], [1, "#59BFB1"]]
    fig = go.Figure(go.Heatmap(z=z, x=labs, y=labs, text=text, texttemplate="%{text}", zmin=0, colorscale=scale, showscale=False,
                               textfont=dict(size=13, color="#FFFFFF" if theme == "dark" else "#1F2933"), xgap=3, ygap=3,
                               hovertemplate="%{y} - %{x}<br>%{z:.5f}<extra></extra>"))
    fig.update_layout(title=dict(text=f"<b>{title}</b>", x=0.01, y=0.96, yanchor="top", font=dict(size=16, family=FONT)),
                      height=190 + 62 * len(idx), margin=dict(l=24, r=24, t=64, b=64), paper_bgcolor=P["bg"], plot_bgcolor=P["bg"],
                      font=dict(family=FONT, color=P["fg"]), yaxis=dict(autorange="reversed", showgrid=False, zeroline=False, showline=False, ticks=""),
                      xaxis=dict(side="bottom", showgrid=False, zeroline=False, showline=False, ticks=""))
    return fig


# ------------------------------------------------------------------------- texts
METHODS = {
    "iqtree": dict(label="IQ-TREE", axis="bar"),
    "nj": dict(label="Neighbor-Joining", axis="bar"),
    "upgma": dict(label="UPGMA", axis="height"),
}

HOW_TO_READ = {
    "common": [
        "Each **tip** is one of your sequences. Tips that join close to the tips are more similar.",
        "Only **horizontal distances** carry information: compare them with the scale bar (substitutions per site; 0.01 = 1 difference per 100 sites).",
        "Vertical order does not matter: branches can be rotated around any node without changing the tree.",
        "Colours only mark the two main groups on each side of the display root - they have no taxonomic meaning.",
    ],
    "iqtree": [
        "The circle shows **UFBoot support** for the split it sits on: green >= 95 is reliable, amber 80-94 is moderate, red < 80 is weak.",
        "IQ-TREE estimates an **unrooted** tree. The root is placed at the midpoint of the longest path for display only; it does not say which sequence is ancestral.",
    ],
    "nj": [
        "Neighbor-Joining is a **distance method**: branch lengths come from the pairwise distance matrix and are additive.",
        "The tree is **unrooted**; the root is a display anchor. No support values are shown unless a bootstrap was run.",
    ],
    "upgma": [
        "UPGMA **assumes a constant substitution rate (molecular clock)**: all tips end at the same level and a node's height is half the distance between the groups it joins.",
        "If sequences evolve at different rates, UPGMA can give a wrong topology; compare with IQ-TREE or Neighbor-Joining.",
    ],
}

REFERENCES = {
    "mafft": "Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7. Mol Biol Evol 30:772-780.",
    "iqtree": "Minh BQ et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol 37:1530-1534.",
    "modelfinder": "Kalyaanamoorthy S et al. (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nat Methods 14:587-589.",
    "ufboot": "Hoang DT et al. (2018) UFBoot2: improving the ultrafast bootstrap approximation. Mol Biol Evol 35:518-522.",
    "k2p": "Kimura M (1980) A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. J Mol Evol 16:111-120.",
    "nj": "Saitou N, Nei M (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. Mol Biol Evol 4:406-425.",
    "upgma": "Sokal RR, Michener CD (1958) A statistical method for evaluating systematic relationships. Univ Kansas Sci Bull 38:1409-1438.",
}


def methods_paragraph(method: str, meta: dict) -> str:
    aligner = meta.get("aligner")
    pre = f"Sequences were aligned with {aligner}. " if aligner else ""
    if method == "iqtree":
        return (f"{pre}The best-fit substitution model ({meta.get('model', 'n/a')}) was selected with ModelFinder (BIC) and the maximum-likelihood tree was "
                f"inferred with {meta.get('engine', 'IQ-TREE')}; branch support was assessed with {meta.get('bootstrap', 1000)} ultrafast bootstrap replicates (UFBoot2). "
                "The tree is unrooted and drawn midpoint-rooted for display only.")
    algo = "UPGMA" if method == "upgma" else "Neighbor-Joining"
    return f"{pre}Pairwise {_pretty_distance(meta.get('distance_method', 'K2P'))} distances were computed and a {algo} tree was built."


def references_for(method: str, meta: dict) -> List[str]:
    keys = {"iqtree": ["iqtree", "modelfinder", "ufboot"], "nj": ["nj"], "upgma": ["upgma"]}[method]
    if method != "iqtree" and "kimura" in str(meta.get("distance_method", "kimura")).lower():
        keys = ["k2p"] + keys
    if "mafft" in str(meta.get("aligner", "")).lower():
        keys = ["mafft"] + keys
    return [REFERENCES[k] for k in keys]


# ------------------------------------------------------------------ Streamlit panel
def _plot(st, fig: go.Figure, key: str) -> None:
    """theme=None keeps this module's own fonts/colours; width='stretch' with a fallback for older Streamlit."""
    try:
        st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG, key=key, theme=None)
    except TypeError:
        st.plotly_chart(fig, use_container_width=True, config=PLOTLY_CONFIG, key=key, theme=None)


def render_phylo_result(newick: str, *, method: str = "iqtree", meta: Optional[dict] = None, theme: str = "dark",
                        dist_names: Optional[List[str]] = None, dist_matrix: Optional[List[List[float]]] = None,
                        key: str = "phylo") -> None:
    import streamlit as st

    meta = meta or {}
    spec = METHODS[method]
    dist_label = _pretty_distance(meta.get("distance_method", "K2P"))
    sup_label = "UFBoot" if method == "iqtree" else "bootstrap"
    support_cut = 95.0 if method == "iqtree" else 70.0
    midpoint = method in {"iqtree", "nj"}
    ins = tree_insights(newick, support_cut=support_cut, dist_names=dist_names,
                        dist_matrix=dist_matrix, method=method, midpoint=midpoint)
    root_note = {"iqtree": "unrooted estimate, midpoint root for display",
                 "nj": "unrooted estimate, midpoint root for display", "upgma": "constant-rate (clock) assumption"}[method]
    mid = f"{meta.get('model', '')} &#183; UFBoot x{meta.get('bootstrap', 1000)}" if method == "iqtree" else f"{dist_label} distances"
    subtitle = f"{spec['label']} &#183; {mid} &#183; {root_note}"

    c = st.columns(4)
    c[0].metric("Method", spec["label"])
    c[1].metric("Model" if method == "iqtree" else "Distance", str(meta.get("model") or "n/a") if method == "iqtree" else dist_label)
    c[2].metric("Sequences", ins["n"])
    c[3].metric("Alignment", f"{meta['n_sites']:,} bp" if meta.get("n_sites") else "-")

    _plot(st, make_tree_figure(newick, theme=theme, subtitle=subtitle, axis_mode=spec["axis"],
                               support_label=sup_label, midpoint=midpoint), f"{key}_tree")
    st.caption("Camera icon = download as SVG (vector, ready for a paper or slide).")

    left, right = st.columns(2)
    with left:
        st.markdown("##### Key findings")
        for line in findings(ins, sup_label, dist_label, support_cut) or ["Not enough sequences to summarise."]:
            st.markdown(f"- {line}")
    with right:
        st.markdown("##### How to read this tree")
        for line in HOW_TO_READ["common"] + HOW_TO_READ[method]:
            st.markdown(f"- {line}")
    for w in ins["warnings"]:
        st.warning(w)

    if dist_names and dist_matrix:
        with st.expander("Observed pairwise distances (same order as the tree)"):
            _plot(st, make_distance_heatmap(dist_names, dist_matrix, order=leaf_order(newick, midpoint=midpoint), theme=theme), f"{key}_heat")
    with st.expander("Methods paragraph and references (copy into your report)"):
        st.code(methods_paragraph(method, meta) + "\n\nReferences\n" + "\n".join(f"- {r}" for r in references_for(method, meta)), language=None)
