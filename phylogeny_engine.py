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
    n = len(names)
    
    # Ensure symmetric distance matrix
    dm = distance_matrix.copy()
    dm = (dm + dm.T) / 2  # Make symmetric
    
    # Track which nodes are active
    active_nodes = list(range(n))
    node_names = names.copy()
    node_counter = n
    
    # Store tree edges
    edges = []
    
    while len(active_nodes) > 2:
        # Calculate net divergence (u values)
        u = np.zeros(len(active_nodes))
        for i in range(len(active_nodes)):
            u[i] = dm[active_nodes[i], :].sum() / (len(active_nodes) - 2)
        
        # Find pair with minimum distance (corrected)
        min_dist = float('inf')
        min_i, min_j = 0, 1
        
        for i in range(len(active_nodes)):
            for j in range(i+1, len(active_nodes)):
                corrected_dist = dm[active_nodes[i], active_nodes[j]] - u[i] - u[j]
                if corrected_dist < min_dist:
                    min_dist = corrected_dist
                    min_i, min_j = i, j
        
        # Calculate branch lengths
        branch_i = 0.5 * (dm[active_nodes[min_i], active_nodes[min_j]] + u[min_i] - u[min_j])
        branch_j = dm[active_nodes[min_i], active_nodes[min_j]] - branch_i
        
        # Record edges
        edges.append({
            "parent": f"Node_{node_counter}",
            "child_1": node_names[active_nodes[min_i]],
            "child_2": node_names[active_nodes[min_j]],
            "branch_1": round(branch_i, 4),
            "branch_2": round(branch_j, 4),
        })
        
        # Create new node and update distance matrix
        new_idx = max(active_nodes) + 1
        node_names.append(f"Node_{node_counter}")
        node_counter += 1
        
        # Calculate distances from new node to all remaining
        new_dm_row = []
        for k in range(len(active_nodes)):
            if k != min_i and k != min_j:
                new_dist = 0.5 * (
                    dm[active_nodes[min_i], active_nodes[k]] +
                    dm[active_nodes[min_j], active_nodes[k]] -
                    dm[active_nodes[min_i], active_nodes[min_j]]
                )
                new_dm_row.append(new_dist)
        
        # Update active nodes and distance matrix
        new_active = [active_nodes[k] for k in range(len(active_nodes)) if k != min_i and k != min_j]
        new_active.append(new_idx)
        
        old_size = len(dm)
        new_dm = np.zeros((old_size + 1, old_size + 1))
        new_dm[:old_size, :old_size] = dm
        
        dm_idx = 0
        for k in range(len(active_nodes)):
            if k != min_i and k != min_j:
                new_dm[new_idx, active_nodes[k]] = new_dm_row[dm_idx]
                new_dm[active_nodes[k], new_idx] = new_dm_row[dm_idx]
                dm_idx += 1
        
        dm = new_dm
        active_nodes = new_active
    
    # Final distance between last two nodes
    final_edge = {
        "parent": "Root",
        "child_1": node_names[active_nodes[0]],
        "child_2": node_names[active_nodes[1]],
        "branch_1": round(dm[active_nodes[0], active_nodes[1]] / 2, 4),
        "branch_2": round(dm[active_nodes[0], active_nodes[1]] / 2, 4),
    }
    edges.append(final_edge)
    
    return {
        "algorithm": "Neighbor-Joining",
        "sequence_names": names,
        "edges": edges,
        "tree_type": "Additive (non-clock)",
    }


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
