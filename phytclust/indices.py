import numpy as np
import statistics


# colless index
def colless_index_calc(tree):
    internal_nodes = [node for node in tree.find_clades(terminal=False)]
    colless_sum = 0
    for node in internal_nodes:
        left_size = len(node.clades[0].get_terminals())
        right_size = len(node.clades[1].get_terminals())
        colless_sum += abs(left_size - right_size)
    return colless_sum


def normalized_colless(tree):
    colless_sum = colless_index_calc(tree)
    n = tree.count_terminals()
    normalized_colless = colless_sum * (2 / (n * (n - 1)))
    return normalized_colless


# stemmy/tippy
def calculate_internal_terminal_ratio(tree):
    internal_length_sum = 0
    terminal_length_sum = 0
    for node in tree.find_clades():
        if node.is_terminal():
            terminal_length_sum += (
                node.branch_length if node.branch_length is not None else 0
            )
        else:
            internal_length_sum += (
                node.branch_length if node.branch_length is not None else 0
            )

    if terminal_length_sum != 0:
        ratio = internal_length_sum / terminal_length_sum
    else:
        ratio = float("inf")
    return ratio


def calculate_int_term_ratio(tree):
    ratio = calculate_internal_terminal_ratio(tree)
    num_terminals = tree.count_terminals()
    normalized_ratio = ratio * (num_terminals / ((2 * num_terminals) - 2))
    return normalized_ratio


# branch length variance
def collect_branch_lengths(node):
    lengths = []
    if node.clades:
        for child in node.clades:
            if child.branch_length is not None:
                lengths.append(child.branch_length)
            lengths.extend(collect_branch_lengths(child))
    return lengths


def calculate_variance_branch_length(tree):
    branch_lengths = collect_branch_lengths(tree.root)
    if branch_lengths:
        return np.std(branch_lengths)
    else:
        return 0


# variance in total branch length (node to root) - punctuated vs gradual evolution
def total_branch_lengths(tree):
    total_bl = []
    for terminal in tree.get_terminals():
        branch_length = tree.distance(terminal)
        total_bl.append(branch_length)
    return total_bl


def calculate_total_length_variation(tree):
    branch_lengths = total_branch_lengths(tree)
    if branch_lengths:
        return np.std(branch_lengths)
    else:
        return 0


# variance = internal nodes
def calculate_internal_variance(tree):
    internal_branch_lengths = [
        node.branch_length
        for node in tree.get_nonterminals()
        if node.branch_length is not None
    ]
    internal_bl_variation = (
        np.std(internal_branch_lengths) if internal_branch_lengths else 0
    )
    return internal_bl_variation


# terminal
def calculate_terminal_variance(tree):
    terminal_branch_lengths = [
        node.branch_length
        for node in tree.get_terminals()
        if node.branch_length is not None
    ]
    terminal_bl_variation = (
        np.std(terminal_branch_lengths) if terminal_branch_lengths else 0
    )

    return terminal_bl_variation


# ratio int_variation/terminal variation
def variation_ratio(tree):
    internal_branch_lengths = [
        node.branch_length
        for node in tree.get_nonterminals()
        if node.branch_length is not None
    ]
    terminal_branch_lengths = [
        node.branch_length
        for node in tree.get_terminals()
        if node.branch_length is not None
    ]
    internal_bl_variation = (
        np.std(internal_branch_lengths) if internal_branch_lengths else 0
    )
    terminal_bl_variation = (
        np.std(terminal_branch_lengths) if terminal_branch_lengths else 0
    )
    variation_ratio = (
        (internal_bl_variation / terminal_bl_variation)
        if terminal_bl_variation
        else float("inf")
    )

    return variation_ratio


# how much temrinal nodes contribute to total length
def calculate_terminal_contributions(tree):
    total_branch_length = sum(
        node.branch_length for node in tree.find_clades() if node.branch_length
    )
    terminal_branch_length = sum(
        leaf.branch_length for leaf in tree.get_terminals() if leaf.branch_length
    )
    percentage_contribution = (terminal_branch_length / total_branch_length) * 100
    return percentage_contribution


# coeff of var for
def find_siblings(tree):
    sibling_distances = []
    for clade in tree.find_clades():
        if clade.is_preterminal():
            terminals = clade.get_terminals()
            if len(terminals) == 2:
                dist = tree.distance(terminals[0], terminals[1])
                sibling_distances.append(dist)
    return sibling_distances


def calculate_variance_of_distances(distances):
    if distances and len(distances) > 1:
        std_dev = statistics.stdev(distances)
        mean = statistics.mean(distances)
        variance = std_dev**2
        return variance
    else:
        return None


def calculate_coefficient_of_variation(distances):
    if distances and len(distances) > 1:
        std_dev = statistics.stdev(distances)
        mean = statistics.mean(distances)
        if mean != 0:  # Avoid division by zero
            cv = std_dev / mean
            return cv
        else:
            return None  # Mean is zero, CV is not defined
    else:
        return None  # Not enough data to calculate CV


# gini coeff
def calculate_proportions(tree, sibling_distances):
    total_branch_length = sum(
        node.branch_length for node in tree.find_clades() if node.branch_length
    )
    proportions = [distance / total_branch_length for distance in sibling_distances]
    return proportions


def gini_coefficient(proportions):
    n = len(proportions)
    sorted_proportions = sorted(proportions)
    numerator = sum(
        (2 * i - n - 1) * x for i, x in enumerate(sorted_proportions, start=1)
    )
    total = sum(sorted_proportions)
    gini = numerator / (n * total)
    return gini

