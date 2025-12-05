import pathlib
from Bio import Phylo
from phytclust.algo.core import PhytClust


def test_exact_k_partition(tmp_path):
    tree_path = pathlib.Path(__file__).parent / "test_tree.nwk"
    tree = Phylo.read(tree_path, "newick")

    pc = PhytClust(tree=tree)
    pc.run(k=2)

    clusters = pc.clusters

    if isinstance(next(iter(clusters.values())), dict):
        cluster_ids = {
            v.get("cluster") or v.get("cluster_id") for v in clusters.values()
        }
    else:
        cluster_ids = set(clusters.values())

    assert len(cluster_ids) == 2


def test_best_global_runs(tmp_path):
    tree_path = pathlib.Path(__file__).parent / "test_tree.nwk"
    tree = Phylo.read(tree_path, "newick")

    pc = PhytClust(tree=tree)
    result = pc.best_global(top_n=1)

    assert result
    assert isinstance(result[0], dict)
