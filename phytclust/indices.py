import numpy as np
import statistics
from Bio.Phylo.BaseTree import Tree, Clade
from typing import List, Optional, Tuple, Union


# Colless index
def colless_index_calc(tree: Tree) -> int:
    """
    Calculate the Colless index for a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        int: The Colless index.
    """
    internal_nodes = [node for node in tree.find_clades(terminal=False)]
    colless_sum = 0
    for node in internal_nodes:
        left_size = len(node.clades[0].get_terminals())
        right_size = len(node.clades[1].get_terminals())
        colless_sum += abs(left_size - right_size)
    return colless_sum


def normalized_colless(tree: Tree) -> float:
    """
    Calculate the normalized Colless index for a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The normalized Colless index.
    """
    colless_sum = colless_index_calc(tree)
    n = tree.count_terminals()
    if n <= 2:
        return 0
    normalized_colless = (2 * colless_sum) / ((n - 1) * (n - 2))
    return normalized_colless


# Stemmy/Tippy
def calculate_internal_terminal_ratio(tree: Tree) -> float:
    """
    Calculate the ratio of internal to terminal branch lengths.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The ratio of internal to terminal branch lengths.
    """
    internal_length_sum = 0
    terminal_length_sum = 0
    for node in tree.find_clades():
        if node.is_terminal():
            terminal_length_sum += node.branch_length or 0
        else:
            internal_length_sum += node.branch_length or 0

    return (
        internal_length_sum / terminal_length_sum
        if terminal_length_sum != 0
        else float("inf")
    )


def calculate_int_term_ratio(tree: Tree) -> float:
    """
    Calculate the normalized ratio of internal to terminal branch lengths.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The normalized ratio of internal to terminal branch lengths.
    """
    ratio = calculate_internal_terminal_ratio(tree)
    num_terminals = tree.count_terminals()
    normalized_ratio = ratio * (num_terminals / ((2 * num_terminals) - 2))
    return normalized_ratio


# Branch length variance
def collect_branch_lengths(node: Clade) -> List[float]:
    """
    Collect all branch lengths from a given node.

    Args:
        node (Clade): The starting node.

    Returns:
        List[float]: A list of branch lengths.
    """
    lengths = []
    if node.clades:
        for child in node.clades:
            if child.branch_length is not None:
                lengths.append(child.branch_length)
            lengths.extend(collect_branch_lengths(child))
    return lengths


def calculate_variance_branch_length(tree: Tree) -> float:
    """
    Calculate the variance of branch lengths in a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The variance of branch lengths.
    """
    branch_lengths = collect_branch_lengths(tree.root)
    return np.std(branch_lengths) if branch_lengths else 0


# Variance in total branch length (node to root) - punctuated vs gradual evolution
def total_branch_lengths(tree: Tree) -> List[float]:
    """
    Calculate the total branch lengths from root to each terminal node.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        List[float]: A list of total branch lengths.
    """
    return [tree.distance(terminal) for terminal in tree.get_terminals()]


def calculate_total_length_variation(tree: Tree) -> float:
    """
    Calculate the variance of total branch lengths in a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The variance of total branch lengths.
    """
    branch_lengths = total_branch_lengths(tree)
    return np.std(branch_lengths) if branch_lengths else 0


# Variance in internal nodes
def calculate_internal_variance(tree: Tree) -> float:
    """
    Calculate the variance of internal branch lengths in a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The variance of internal branch lengths.
    """
    internal_branch_lengths = [
        node.branch_length
        for node in tree.get_nonterminals()
        if node.branch_length is not None
    ]
    return np.std(internal_branch_lengths) if internal_branch_lengths else 0


# Variance in terminal nodes
def calculate_terminal_variance(tree: Tree) -> float:
    """
    Calculate the variance of terminal branch lengths in a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The variance of terminal branch lengths.
    """
    terminal_branch_lengths = [
        node.branch_length
        for node in tree.get_terminals()
        if node.branch_length is not None
    ]
    return np.std(terminal_branch_lengths) if terminal_branch_lengths else 0


# Ratio of internal to terminal branch length variance
def variation_ratio(tree: Tree) -> float:
    """
    Calculate the ratio of internal to terminal branch length variance.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The ratio of internal to terminal branch length variance.
    """
    internal_variance = calculate_internal_variance(tree)
    terminal_variance = calculate_terminal_variance(tree)
    return internal_variance / terminal_variance if terminal_variance else float("inf")


# Contribution of terminal nodes to total branch length
def calculate_terminal_contributions(tree: Tree) -> float:
    """
    Calculate the percentage contribution of terminal branch lengths to the total branch length.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The percentage contribution of terminal branch lengths.
    """
    total_branch_length = sum(
        node.branch_length for node in tree.find_clades() if node.branch_length
    )
    terminal_branch_length = sum(
        leaf.branch_length for leaf in tree.get_terminals() if leaf.branch_length
    )
    return (
        (terminal_branch_length / total_branch_length) * 100
        if total_branch_length
        else 0
    )


# Coefficient of variation for sibling distances
def find_siblings(tree: Tree) -> List[float]:
    """
    Find the distances between sibling terminal nodes.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        List[float]: A list of distances between sibling terminal nodes.
    """
    sibling_distances = []
    for clade in tree.find_clades():
        if clade.is_preterminal():
            terminals = clade.get_terminals()
            if len(terminals) == 2:
                sibling_distances.append(tree.distance(terminals[0], terminals[1]))
    return sibling_distances


def calculate_variance_of_distances(distances: List[float]) -> Optional[float]:
    """
    Calculate the variance of a list of distances.

    Args:
        distances (List[float]): A list of distances.

    Returns:
        Optional[float]: The variance of the distances, or None if not enough data.
    """
    if distances and len(distances) > 1:
        return statistics.variance(distances)
    return None


def calculate_coefficient_of_variation(distances: List[float]) -> Optional[float]:
    """
    Calculate the coefficient of variation for a list of distances.

    Args:
        distances (List[float]): A list of distances.

    Returns:
        Optional[float]: The coefficient of variation, or None if not enough data.
    """
    if distances and len(distances) > 1:
        mean = statistics.mean(distances)
        if mean != 0:
            return statistics.stdev(distances) / mean
    return None


# Gini coefficient
def calculate_proportions(tree: Tree, sibling_distances: List[float]) -> List[float]:
    """
    Calculate the proportions of sibling distances relative to the total branch length.

    Args:
        tree (Tree): The phylogenetic tree.
        sibling_distances (List[float]): A list of sibling distances.

    Returns:
        List[float]: A list of proportions.
    """
    total_branch_length = sum(
        node.branch_length for node in tree.find_clades() if node.branch_length
    )
    return [distance / total_branch_length for distance in sibling_distances]


def gini_coefficient(proportions: List[float]) -> float:
    """
    Calculate the Gini coefficient for a list of proportions.

    Args:
        proportions (List[float]): A list of proportions.

    Returns:
        float: The Gini coefficient.
    """
    n = len(proportions)
    sorted_proportions = sorted(proportions)
    numerator = sum(
        (2 * i - n - 1) * x for i, x in enumerate(sorted_proportions, start=1)
    )
    total = sum(sorted_proportions)
    return numerator / (n * total) if total != 0 else 0
