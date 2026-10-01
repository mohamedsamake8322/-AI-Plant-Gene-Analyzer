"""
phylogeny_engine.py
-------------------
Phylogenetic tree construction and visualization.
Implements UPGMA and Neighbor-Joining algorithms.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform
import json


class PhyloNode:
    """Simple tree node for phylogenetic trees."""
    
    def __init__(self, name: str = "", distance: float = 0.0, children: List = None):
        self.name = name
        self.distance = distance
        self.children = children if children else []
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "distance": self.distance,
            "children": [child.to_dict() for child in self.children],
        }


# ─── UPGMA (Unweighted Pair Group Method with Arithmetic Mean) ────────────────

def upgma(distance_matrix: np.ndarray, names: List[str]) -> Dict:
    """
    UPGMA hierarchical clustering algorithm.
    Constructs an ultrametric tree (all leaves equidistant from root).
    
    Args:
        distance_matrix: Pairwise distance matrix (numpy array)
        names: Sequence/species names
    
    Returns:
        Dict with tree structure, linkage matrix, and dendrogram data
    """
    n = len(names)
    
# Ensure symmetric distance matrix and numeric values
    dm = np.array(distance_matrix, dtype=float)
    if dm.ndim != 2 or dm.shape[0] != dm.shape[1]:
        raise ValueError("Distance matrix must be square.")
    dm = (dm + dm.T) / 2.0  # Make symmetric
    np.fill_diagonal(dm, 0.0)
    
    # Convert distance matrix to condensed form for scipy
    condensed_dist = squareform(dm)
    
    # Perform hierarchical clustering
    Z = linkage(condensed_dist, method='average')
    
    # Convert to tree structure
    tree = to_tree(Z, rd=False)
    
    # Build dendrogram data
    dendro = dendrogram(Z, labels=names, no_plot=True)
    
    newick = linkage_to_newick(Z, names)
    return {
        "algorithm": "UPGMA",
        "sequence_names": names,
        "linkage_matrix": Z.tolist(),
        "dendrogram_data": {
            "icoord": dendro["icoord"],
            "dcoord": dendro["dcoord"],
            "leaves": dendro["leaves"],
            "color_list": dendro.get("color_list", []),
        },
        "newick": newick,
        "tree_type": "Ultrametric (clock-like)",
    }


# ─── NEIGHBOR-JOINING ─────────────────────────────────────────────────────────

def neighbor_joining(distance_matrix: np.ndarray, names: List[str]) -> Dict:
    """
    Neighbor-Joining algorithm (Saitou & Nei, 1987).
    Constructs an additive tree allowing different evolutionary rates.
    
    Args:
        distance_matrix: Pairwise distance matrix (numpy array)
        names: Sequence/species names
    
    Returns:
        Dict with tree structure and branch lengths
    """
    dm = np.asarray(distance_matrix, dtype=float)
    if dm.ndim != 2 or dm.shape[0] != dm.shape[1] or dm.shape[0] != len(names):
        raise ValueError("Distance matrix must be square and match the number of names.")
    if len(names) < 2:
        raise ValueError("At least two sequences are required for Neighbor-Joining.")

    dm = (dm + dm.T) / 2
    np.fill_diagonal(dm, 0.0)
    active_nodes = list(range(len(names)))
    node_names = list(names)
    edges = []

    while len(active_nodes) > 2:
        active_count = len(active_nodes)
        row_sums = {
            node: sum(dm[node, other] for other in active_nodes)
            for node in active_nodes
        }
        min_pair = min(
            ((i, j) for index, i in enumerate(active_nodes) for j in active_nodes[index + 1:]),
            key=lambda pair: (active_count - 2) * dm[pair[0], pair[1]]
            - row_sums[pair[0]] - row_sums[pair[1]],
        )
        node_i, node_j = min_pair
        pair_distance = dm[node_i, node_j]
        branch_i = 0.5 * pair_distance + (
            row_sums[node_i] - row_sums[node_j]
        ) / (2 * (active_count - 2))
        branch_j = pair_distance - branch_i

        parent_name = f"Node_{len(node_names)}"
        edges.append({
            "parent": parent_name,
            "child_1": node_names[node_i],
            "child_2": node_names[node_j],
            "branch_1": float(branch_i),
            "branch_2": float(branch_j),
        })

        remaining = [node for node in active_nodes if node not in (node_i, node_j)]
        new_idx = len(dm)
        expanded = np.zeros((new_idx + 1, new_idx + 1), dtype=float)
        expanded[:new_idx, :new_idx] = dm
        for node in remaining:
            new_distance = 0.5 * (dm[node_i, node] + dm[node_j, node] - pair_distance)
            expanded[new_idx, node] = new_distance
            expanded[node, new_idx] = new_distance

        dm = expanded
        active_nodes = remaining + [new_idx]
        node_names.append(parent_name)

    final_distance = dm[active_nodes[0], active_nodes[1]]
    final_edge = {
        "parent": "Root",
        "child_1": node_names[active_nodes[0]],
        "child_2": node_names[active_nodes[1]],
        "branch_1": float(final_distance / 2),
        "branch_2": float(final_distance / 2),
    }
    edges.append(final_edge)
    
    newick = nj_edges_to_newick(edges)
    return {
        "algorithm": "Neighbor-Joining",
        "sequence_names": names,
        "edges": edges,
        "newick": newick,
        "tree_type": "Additive (non-clock)",
    }


def linkage_to_newick(linkage_matrix: np.ndarray, labels: List[str]) -> str:
    """Convert a SciPy linkage matrix to Newick with correct branch lengths.

    SciPy stores each node's height in ``node.dist``. Newick needs the
    distance from a node to its parent. SciPy's average-linkage height is
    the full pairwise distance at a UPGMA merge, so each edge uses half of
    the height difference. The root has no parent and must not receive a
    branch length.
    """
    tree = to_tree(linkage_matrix, rd=False)

    def _node_to_newick(node, parent_height: float | None = None) -> str:
        branch_length = None if parent_height is None else max(parent_height - node.dist, 0.0) / 2
        if node.is_leaf():
            suffix = "" if branch_length is None else f":{branch_length:.6f}"
            return f"{_format_newick_label(labels[node.id])}{suffix}"
        left = _node_to_newick(node.get_left(), node.dist)
        right = _node_to_newick(node.get_right(), node.dist)
        subtree = f"({left},{right})"
        return subtree if branch_length is None else f"{subtree}:{branch_length:.6f}"

    return _node_to_newick(tree) + ";"


def _format_newick_label(label: str) -> str:
    """Quote labels that Newick parsers could split or normalize."""
    if not label or any(char.isspace() or char in "_():,;[]'" for char in label):
        return "'" + label.replace("'", "''") + "'"
    return label


def nj_edges_to_newick(edges: List[Dict]) -> str:
    """Build a Newick string from Neighbor-Joining edge list."""
    if not edges:
        return ";"

    children: Dict[str, List[tuple[str, float]]] = {}
    nodes = set()
    for edge in edges:
        parent = edge["parent"]
        nodes.add(parent)
        nodes.add(edge["child_1"])
        nodes.add(edge["child_2"])
        children.setdefault(parent, []).append((edge["child_1"], float(edge["branch_1"])))
        children.setdefault(parent, []).append((edge["child_2"], float(edge["branch_2"])))

    child_nodes = {c for kids in children.values() for c, _ in kids}
    roots = [n for n in nodes if n not in child_nodes] or ["Root"]

    def _render(node: str) -> str:
        if node not in children:
            return _format_newick_label(node)
        parts = [_render(child) + f":{length:.6f}" for child, length in children[node]]
        return f"({','.join(parts)})"

    return _render(roots[0]) + ";"


# ─── PHYLOGENETIC TREE VISUALIZATION ──────────────────────────────────────────

def phylo_to_newick(tree_dict: Dict) -> str:
    """
    Convert tree structure to Newick format (standard phylogenetic format).
    
    Format: (child1:branch1, child2:branch2):parent_branch;
    """
    def _to_newick(node: Dict) -> str:
        if not node.get("children"):
            return f"{node['name']}:{node.get('distance', 0)}"
        
        children_str = ",".join(_to_newick(child) for child in node["children"])
        return f"({children_str}):{node.get('distance', 0)}"
    
    newick = _to_newick(tree_dict) + ";"
    return newick


def newick_to_plotly_tree(newick_str: str) -> Dict:
    """
    Parse Newick format and prepare for Plotly visualization.
    """
    # Simplified Newick parser
    newick_str = newick_str.rstrip(';')
    
    def parse_newick(s: str, idx: int = 0, depth: int = 0) -> Tuple[Dict, int]:
        if idx >= len(s):
            return {}, idx
        
        node = {"name": "", "distance": 0, "children": []}
        
        if s[idx] == '(':
            idx += 1  # skip '('
            while s[idx] != ')':
                if s[idx] == ',':
                    idx += 1
                else:
                    child, idx = parse_newick(s, idx, depth + 1)
                    node["children"].append(child)
            idx += 1  # skip ')'
        
        # Parse name and distance
        name_dist = ""
        while idx < len(s) and s[idx] not in '(),;':
            name_dist += s[idx]
            idx += 1
        
        if ':' in name_dist:
            parts = name_dist.split(':')
            node["name"] = parts[0] or f"Node_{depth}"
            try:
                node["distance"] = float(parts[1])
            except:
                node["distance"] = 0
        else:
            node["name"] = name_dist or f"Node_{depth}"
        
        return node, idx
    
    tree, _ = parse_newick(newick_str)
    return tree
