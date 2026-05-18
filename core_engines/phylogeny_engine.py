"""
Wrapper for phylogeny_engine inside core_engines package.
"""
import phylogeny_engine as _orig

upgma = _orig.upgma
neighbor_joining = _orig.neighbor_joining
phylo_to_newick = _orig.phylo_to_newick
newick_to_plotly_tree = _orig.newick_to_plotly_tree

__all__ = [
    'upgma',
    'neighbor_joining',
    'phylo_to_newick',
    'newick_to_plotly_tree',
]
