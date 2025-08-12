import pathlib
from Bio import Phylo
from phytclust.algo.core import PhytClust


def test_exact_k_partition(tmp_path):
    tree_path = pathlib.Path(__file__).parent / "test_tree.nwk"
    tree = Phylo.read(tree_path, "newick")
    pc = PhytClust(tree=tree)
    clusters = pc.get_clusters(2)

    assert isinstance(clusters, dict)
    assert len(set(clusters.values())) == 2


def test_best_global_runs(tmp_path):
    tree_path = pathlib.Path(__file__).parent / "test_tree.nwk"
    tree = Phylo.read(tree_path, "newick")
    pc = PhytClust(tree=tree)
    result = pc.best_global(top_n=1)

    assert result  # should not be empty
    assert isinstance(result[0], dict)
