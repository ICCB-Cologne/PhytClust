from ete3 import Tree as EteTree
from Bio import Phylo
from io import StringIO
import re
import random
import numpy as np

def distribute_ntips_into_mparts(
    num_tips, num_groups, min_size, variance, method="random", seed=None
):
    """
    Partition m into n parts using multinomial distribution, avoiding 0s.
    First, give each part 1, then distribute the remaining m-n using multinomial.
    """
    if seed is not None:
        np.random.seed(seed)
    m = num_tips
    n = num_groups
    partition = np.ones(n, dtype=int)
    remaining = m - n

    if remaining > 0:
        partition += np.random.multinomial(remaining, np.ones(n) / n)

    return partition

def convert_ete_to_phylo_newick(ete_tree):
    """
    Convert an ETE3 Tree to a Bio.Phylo Tree using Newick format as intermediary.

    Parameters:
    ete_tree (ete3.Tree): The ETE3 tree to convert.

    Returns:
    Bio.Phylo.BaseTree.Tree: The converted Bio.Phylo tree.
    """
    newick_str = ete_tree.write()
    handle = StringIO(newick_str)
    phylo_tree = Phylo.read(handle, "newick")

    return phylo_tree


def generate_ground_truth(tree):
    clusters = {}
    outlier_counter = 0
    max_cluster_id = 0

    for leaf in tree.get_leaves():
        name = leaf.name
        if not name.startswith("Outlier_"):
            cluster_id = int(re.findall(r"\d+", name.split("_")[1])[0])
            if cluster_id > max_cluster_id:
                max_cluster_id = cluster_id

    for leaf in tree.get_leaves():
        name = leaf.name
        if name.startswith("Outlier_"):
            outlier_cluster_id = max_cluster_id + 1 + outlier_counter
            clusters[outlier_cluster_id] = [name]
            outlier_counter += 1
        else:
            cluster_id = re.findall(r"\d+", name.split("_")[1])[0]
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(name)

    return clusters


def flatten_clusters(clusters):
    labels = []
    for cluster_id, members in clusters.items():
        labels.extend([(member, cluster_id) for member in members])
    return labels
