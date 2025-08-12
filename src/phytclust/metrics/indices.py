from __future__ import annotations

import statistics
from typing import List, Optional
import numpy as np
from Bio.Phylo.BaseTree import Tree, Clade


# Colless index (binary only)
def colless_index_calc(tree: Tree) -> int:
    """
    Colless index for *binary* trees: sum over internal nodes of
    |#leaves(left) - #leaves(right)|.
    """
    colless_sum = 0
    for node in tree.find_clades(terminal=False):
        # guard: Colless is defined for bifurcating nodes
        if len(node.clades) != 2:
            # You could raise or skip; here we skip non-binary nodes
            # to keep behavior tolerant.
            continue
        left_size = len(node.clades[0].get_terminals())
        right_size = len(node.clades[1].get_terminals())
        colless_sum += abs(left_size - right_size)
    return colless_sum


def normalized_colless(tree: Tree) -> float:
    """
    Normalized Colless (simple normalization). For non-binary trees,
    we skip non-binary nodes as above.
    """
    colless_sum = colless_index_calc(tree)
    n = tree.count_terminals()
    if n <= 2:
        return 0.0
    return (2.0 * colless_sum) / ((n - 1) * (n - 2))


# Stemmy / Tippy ratio
def calculate_internal_terminal_ratio(tree: Tree) -> float:
    """
    Ratio of internal to terminal branch lengths.
    Treats missing branch lengths as 0.
    """
    internal_len = 0.0
    terminal_len = 0.0
    for node in tree.find_clades():
        bl = node.branch_length or 0.0
        if node.is_terminal():
            terminal_len += bl
        else:
            internal_len += bl
    return (internal_len / terminal_len) if terminal_len != 0 else float("inf")


def calculate_int_term_ratio(tree: Tree) -> float:
    """
    Normalized internal:terminal branch length ratio (simple scaling).
    """
    ratio = calculate_internal_terminal_ratio(tree)
    n = tree.count_terminals()
    return ratio * (n / ((2 * n) - 2)) if n > 1 else float("inf")


# Branch-length variance utilities
def collect_branch_lengths(node: Clade) -> List[float]:
    lengths: List[float] = []
    for child in getattr(node, "clades", []):
        if child.branch_length is not None:
            lengths.append(float(child.branch_length))
        lengths.extend(collect_branch_lengths(child))
    return lengths


def calculate_variance_branch_length(tree: Tree) -> float:
    """
    Standard deviation of all branch lengths in the tree.
    """
    bl = collect_branch_lengths(tree.root)
    return float(np.std(bl)) if bl else 0.0


# Root-to-tip (total) distances
def total_branch_lengths(tree: Tree) -> List[float]:
    """Root-to-tip path length for each terminal."""
    return [float(tree.distance(t)) for t in tree.get_terminals()]


def calculate_total_length_variation(tree: Tree) -> float:
    """Std. dev. of root-to-tip distances."""
    vals = total_branch_lengths(tree)
    return float(np.std(vals)) if vals else 0.0


# Variance on internal vs terminal branches
def calculate_internal_variance(tree: Tree) -> float:
    vals = [
        float(node.branch_length)
        for node in tree.get_nonterminals()
        if node.branch_length is not None
    ]
    return float(np.std(vals)) if vals else 0.0


def calculate_terminal_variance(tree: Tree) -> float:
    vals = [
        float(node.branch_length)
        for node in tree.get_terminals()
        if node.branch_length is not None
    ]
    return float(np.std(vals)) if vals else 0.0


def variation_ratio(tree: Tree) -> float:
    """Internal variance / terminal variance."""
    v_int = calculate_internal_variance(tree)
    v_term = calculate_terminal_variance(tree)
    return (v_int / v_term) if v_term else float("inf")


# Terminal contribution to total length
def calculate_terminal_contributions(tree: Tree) -> float:
    """
    Percentage of total branch length contributed by terminal edges.
    """
    total = sum((node.branch_length or 0.0) for node in tree.find_clades())
    term = sum((leaf.branch_length or 0.0) for leaf in tree.get_terminals())
    return (term / total) * 100.0 if total else 0.0


# Sibling distances (for CV / variance)
def find_siblings(tree: Tree) -> List[float]:
    """
    Distances between siblings in every preterminal clade with exactly 2 terminals.
    """
    dists: List[float] = []
    for clade in tree.find_clades():
        if clade.is_preterminal():
            terminals = clade.get_terminals()
            if len(terminals) == 2:
                dists.append(float(tree.distance(terminals[0], terminals[1])))
    return dists


def calculate_variance_of_distances(distances: List[float]) -> Optional[float]:
    if distances and len(distances) > 1:
        return float(statistics.variance(distances))
    return None


def calculate_coefficient_of_variation(distances: List[float]) -> Optional[float]:
    if distances and len(distances) > 1:
        mean = statistics.mean(distances)
        if mean != 0:
            return float(statistics.stdev(distances) / mean)
    return None


# Gini coefficient on proportions
def calculate_proportions(tree: Tree, sibling_distances: List[float]) -> List[float]:
    total = sum((node.branch_length or 0.0) for node in tree.find_clades())
    return [(d / total) if total else 0.0 for d in sibling_distances]


def gini_coefficient(proportions: List[float]) -> float:
    """
    Gini coefficient of a vector of non-negative proportions.
    Returns 0 if empty or sum is zero.
    """
    n = len(proportions)
    if n == 0:
        return 0.0
    sorted_p = sorted(proportions)
    total = sum(sorted_p)
    if total == 0:
        return 0.0
    numerator = sum((2 * i - n - 1) * x for i, x in enumerate(sorted_p, start=1))
    return float(numerator / (n * total))
