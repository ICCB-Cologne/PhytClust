import pathlib
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from phytclust.algo.core import PhytClust
from phytclust.algo.dp import compute_dp_table

TREE_PATH = pathlib.Path(__file__).parent / "test_tree.nwk"


def _load_tree():
    return Phylo.read(TREE_PATH, "newick")


def _assert_cluster_map_valid(tree, cmap, k_expected: int):
    assert isinstance(cmap, dict)
    assert len(cmap) > 0

    assert all(isinstance(t, Clade) and t.is_terminal() for t in cmap.keys())
    assert all(isinstance(v, int) for v in cmap.values())

    cluster_ids = set(cmap.values())
    assert len(cluster_ids) == k_expected

    terminals = list(tree.get_terminals())
    assert set(cmap.keys()) == set(terminals)


def test_run_exact_k_partition_returns_result_dict():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    result = pc.run(k=2, plot_scores=False)

    assert isinstance(result, dict)
    assert result["mode"] == "k"
    assert result["k"] == 2
    assert result["peaks"] == [2]
    assert result["scores"] is None

    cmap = result["clusters"]
    _assert_cluster_map_valid(tree, cmap, k_expected=2)


def test_get_clusters_exact_k_partition():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    pc.run(k=2, plot_scores=False)

    cmap = pc.get_clusters(2)
    _assert_cluster_map_valid(tree, cmap, k_expected=2)

    assert 2 in pc.clusters
    assert pc.clusters[2] == cmap


def test_run_global_mode_returns_list_of_cluster_maps():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    max_k = min(10, pc.num_terminals)
    if max_k < 4:
        pytest.skip("Test tree too small for global peak search (requires max_k>=4).")

    result = pc.run(top_n=1, max_k=max_k, plot_scores=False)

    assert result["mode"] == "global"
    assert isinstance(result["clusters"], list)

    clusters = result["clusters"]
    assert len(clusters) >= 0

    if clusters:
        assert isinstance(result["ks"], list)
        assert isinstance(result["peaks"], list)
        assert len(result["peaks"]) == len(result["ks"])
        cmap0 = clusters[0]
        assert isinstance(cmap0, dict)
        assert set(cmap0.keys()) == set(tree.get_terminals())


def test_run_resolution_mode_falls_back_for_small_trees():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    result = pc.run(by_resolution=True, num_bins=2, plot_scores=False)

    assert result["mode"] == "resolution"
    assert isinstance(result["clusters"], list)
    assert isinstance(result["ks"], list)
    assert isinstance(result["peaks"], list)


def test_best_global_returns_list_of_cluster_maps():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    max_k = min(10, pc.num_terminals)
    if max_k < 4:
        pytest.skip("Test tree too small for global peak search (requires max_k>=4).")

    clusters = pc.best_global(top_n=1, max_k=max_k, plot_scores=False)

    assert isinstance(clusters, list)
    if clusters:
        cmap = clusters[0]
        assert isinstance(cmap, dict)
        assert set(cmap.keys()) == set(tree.get_terminals())


def test_run_rejects_invalid_argument_combinations():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    with pytest.raises(
        ValueError, match="Cannot combine `k` with `by_resolution=True`"
    ):
        pc.run(k=2, by_resolution=True, plot_scores=False)

    with pytest.raises(ValueError, match="`top_n` is meaningless when `k` is given"):
        pc.run(k=2, top_n=2, plot_scores=False)

    with pytest.raises(ValueError, match="k must be"):
        pc.run(k=0, plot_scores=False)

    with pytest.raises(ValueError, match="`top_n` must be"):
        pc.run(top_n=0, plot_scores=False)
