"""Standalone tests for phylo_view.py (no dependency on the rest of the app).  Run: pytest test_phylo_view.py"""
import pytest
from phylo_view import (parse_newick, prepare_tree, tree_insights, make_tree_figure, leaf_order,
                        split_label, _leaves, midpoint_root)

A = "'cqi:110688118 | organism=Chenopodium quinoa'"
B = "'cqi:110684354 | organism=Chenopodium quinoa'"
C = "'test copy (A+15) | organism=Chenopodium quinoa'"
D = "'test copy (B+15) | organism=Chenopodium quinoa'"
IQ = f"({A}:0.0000025118,({B}:0.0000020379,{D}:0.0246483335)100:0.0501366272,{C}:0.0246269447):0.0;"
UPGMA = f"(({B}:0.012318,{D}:0.012318):0.024980,({A}:0.012320,{C}:0.012320):0.024977);"
NJ = f"({B}:0.001156,{D}:0.023479,({A}:0.000310,{C}:0.024331):0.049956);"
NJ_NEG = f"(({A}:0.0003,{C}:0.0243):0.025,({B}:-0.0212,{D}:0.0458):0.025);"


def short(names):
    return {split_label(n)[0] for n in names}


def test_quoted_names_and_support():
    root = parse_newick(IQ)
    names = [l.name for l in _leaves(root)]
    assert A.strip("'") in names and C.strip("'") in names
    assert any(c.support == 100 for c in root.children)
    assert parse_newick("((a:1,b:1)95.2/87:1,c:2);").children[0].support == 87


def test_midpoint_root_balances_longest_path():
    root = midpoint_root(parse_newick(NJ))
    assert len(root.children) == 2

    def depth(n):
        return max((c.length + depth(c) for c in n.children), default=0.0)
    left, right = (c.length + depth(c) for c in root.children)
    assert left == pytest.approx(right, abs=1e-9)


def test_insights_groups_and_closest_pair():
    ins = tree_insights(IQ)
    groups, support = ins["split"]
    assert {frozenset(short(g)) for g in groups} == {frozenset({"cqi:110688118", "test copy (A+15)"}),
                                                    frozenset({"cqi:110684354", "test copy (B+15)"})}
    assert support == 100
    # A-C (0.024629) and B-D (0.024650) are tied within 2 %: both must be reported
    pairs = {frozenset(short([a, b])) for _, a, b, *_ in ins["closest"]}
    assert pairs == {frozenset({"cqi:110688118", "test copy (A+15)"}), frozenset({"cqi:110684354", "test copy (B+15)"})}
    assert ins["closest"][0][0] == pytest.approx(0.02463, abs=1e-4)
    assert any("4 sequences" in w for w in ins["warnings"])


def test_upgma_is_ultrametric():
    root, rerooted = prepare_tree(UPGMA)
    assert not rerooted
    xs = [l.x for l in _leaves(root)]
    assert max(xs) - min(xs) < 1e-5 and xs[0] == pytest.approx(0.0373, abs=1e-4)


def test_negative_branch_is_flagged_not_crashing():
    ins = tree_insights(NJ_NEG)
    assert any("Negative branch" in w for w in ins["warnings"])
    make_tree_figure(NJ_NEG)


def test_leaf_order_is_a_permutation_and_figure_builds():
    assert set(leaf_order(IQ)) == {A.strip("'"), B.strip("'"), C.strip("'"), D.strip("'")}
    fig = make_tree_figure(IQ, subtitle="x")
    assert len(fig.data) >= 3


NAMES = [A.strip("'"), B.strip("'"), C.strip("'"), D.strip("'")]
K2P = [[0, .050281, .024641, .074886], [.050281, 0, .076584, .024635],
       [.024641, .076584, 0, .096626], [.074886, .024635, .096626, 0]]


def test_observed_distances_replace_tree_distances():
    ins = tree_insights(UPGMA, dist_names=NAMES, dist_matrix=K2P, method="upgma")
    assert ins["source"] == "observed"
    d, a, b, *_ = ins["farthest"][0]
    assert d == pytest.approx(0.096626) and short([a, b]) == {"test copy (A+15)", "test copy (B+15)"}


def test_upgma_misfit_is_reported():
    ins = tree_insights(UPGMA, dist_names=NAMES, dist_matrix=K2P, method="upgma")
    rel, a, b, dt, do = ins["fit"]
    assert rel == pytest.approx(0.484, abs=0.01) and short([a, b]) == {"cqi:110688118", "cqi:110684354"}
    assert any("does not reproduce the data" in w for w in ins["warnings"])
    assert any("No branch support" in w for w in ins["warnings"])


def test_iqtree_fits_the_data_well():
    ins = tree_insights(IQ, dist_names=NAMES, dist_matrix=K2P, method="iqtree")
    assert ins["fit"][0] < 0.05
    assert not any("does not reproduce" in w for w in ins["warnings"])


def test_heatmap_masks_diagonal():
    from phylo_view import make_distance_heatmap
    fig = make_distance_heatmap(NAMES, K2P, order=leaf_order(IQ))
    z = fig.data[0].z
    assert all(z[i][i] is None for i in range(4))


def test_midpoint_on_binary_rooted_nj_keeps_distances_and_has_no_unary_node():
    """The NJ tree exported by the app has a binary 'anchor' root: midpoint=True must remove it cleanly."""
    from phylo_view import _internal, _patristic
    nj = f"(({A}:0.000310,{C}:0.024331):0.024978,({B}:0.001156,{D}:0.023479):0.024978);"
    root, rerooted = prepare_tree(nj, midpoint=True)
    assert rerooted and len(root.children) == 2
    assert all(len(i.children) != 1 for i in _internal(root))
    d = {frozenset(short([a, b])): v for v, a, b in _patristic(root)}
    assert d[frozenset({"cqi:110688118", "test copy (A+15)"})] == pytest.approx(0.024641, abs=1e-6)
    assert d[frozenset({"test copy (A+15)", "test copy (B+15)"})] == pytest.approx(0.024331 + 0.049956 + 0.023479, abs=1e-6)


def test_negative_branch_is_still_flagged_after_midpoint_rooting():
    ins = tree_insights(NJ_NEG, method="nj", midpoint=True)
    assert any("Negative branch" in w for w in ins["warnings"])


def test_bootstrap_thresholds_differ_from_ufboot():
    from phylo_view import _band
    assert _band(80, "UFBoot") != _band(80, "bootstrap")       # 80 is only 'moderate' for UFBoot, 'reliable' for bootstrap
    ins = tree_insights("((a:1,b:1)75:1,(c:1,d:1)75:1);", support_cut=70.0, method="nj")
    assert ins["weak"] == 0


def test_malformed_newick_raises_clear_errors():
    for bad in ("", "(a:1,b:1", "(a:1,b:1))x;"):
        with pytest.raises(ValueError):
            parse_newick(bad)
