import numpy as np
import pandas as pd


def pairwise_distances(tree, mode="terminals", as_dataframe=False, mrca=None):
    """
    Calculate pairwise distances between nodes in a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.
        mode (str, optional): Mode to select nodes. Options are "terminals", "nonterminals", "all". Defaults to "terminals".
        as_dataframe (bool, optional): Whether to return the result as a pandas DataFrame. Defaults to False.
        mrca (str, optional): Name of the most recent common ancestor to limit the scope. Defaults to None.

    Returns:
        np.ndarray or pd.DataFrame: Pairwise distance matrix.

    Raises:
        ValueError: If no node is found with the name specified in mrca.
    """
    if mrca:
        clade = next((cl for cl in tree.find_clades() if cl.name == mrca), None)
        if clade is None:
            raise ValueError(f"No node found with name {mrca}")
        terms = list(clade.get_terminals())
    elif mode == "terminals":
        terms = list(tree.get_terminals())
    elif mode == "nonterminals":
        terms = list(tree.get_nonterminals())
    elif mode == "all":
        terms = list(tree.get_terminals()) + list(tree.get_nonterminals())
    else:
        print("Invalid mode")
        return None

    terms = list(terms)
    n = len(terms)
    dist = np.zeros((n, n))

    indices = np.triu_indices(n, k=1)
    for i, j in zip(*indices):
        bl = tree.distance(terms[i], terms[j])
        dist[i, j] = bl
        dist[j, i] = bl

    if as_dataframe:
        dist_df = pd.DataFrame(
            dist,
            index=[term.name for term in terms],
            columns=[term.name for term in terms],
        )
        return dist_df
    else:
        return dist


def get_parent(tree, child):
    # tree = self._no_outgroup_tree if self._no_outgroup_tree else self.tree
    node_path = tree.get_path(child)
    if len(node_path) > 1:
        return node_path[-2]
    else:
        return None


def count_branches_in_clusters(clusters):
    branch_count = 0
    for _, clades in clusters.items():
        N = len(clades)
        if N > 1:
            branch_count += 2 * N - 2
    return branch_count


def find_all_min_indices(arr):
    if len(arr) == 0:
        return []

    min_value = float("inf")
    min_indices = []

    for i, value in enumerate(arr):
        if value < min_value:
            min_value = value
            min_indices = [i]  # Start a new list of indices
        elif value == min_value:
            min_indices.append(i)  # Add index to existing list

    return min_indices, min_value


def rename_internal_nodes(tree):
    """
    Renames all internal nodes of a given Phylo tree object.
    Each internal node will be named as 'internal_X' where X is an incrementing integer.

    :param tree: Bio.Phylo tree object
    """
    internal_node_count = 1
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = f"internal_{internal_node_count}"
            internal_node_count += 1
