import io
import os
from Bio import Phylo
import random
from ete3 import Tree, TreeStyle
import string
import random


def total_branch_length(tree):
    """
    Calculate the total of all path lengths from each terminal node to the root.

    Parameters:
    tree (Tree): A phylogenetic tree object from the ETE toolkit.

    Returns:
    float: The total of all path lengths from each terminal node to the root.
    """
    total_path_length = 0.0
    for leaf in tree.iter_leaves():
        current_node = leaf
        while not current_node.is_root():
            total_path_length += current_node.dist
            current_node = current_node.up

    return total_path_length


def distribute_subtree_sums(
    total_length, num_subtrees, num_tips_list, mode="equal", seed=None
):
    if seed is not None:
        random.seed(seed)

    if mode == "equal":
        return [total_length / num_subtrees]

    elif mode == "random":
        min_length = total_length * 0.1
        remaining_length = total_length - (min_length * num_subtrees)
        lengths = [min_length] * num_subtrees

        for _ in range(int(remaining_length)):
            lengths[random.randint(0, num_subtrees - 1)] += 1

        return lengths

    elif mode == "proportional":
        total_tips = sum(num_tips_list)
        lengths = [(total_length / total_tips) * tips for tips in num_tips_list]
        return lengths

    else:
        raise ValueError(
            "Invalid mode. Choose from 'equal', 'random', or 'proportional'."
        )


def assign_branch_lengths_based_on_leaves(tree, total_branch_length):
    """
    Assign branch lengths based on the number of terminal nodes (leaves) under each branch.

    Parameters:
    tree (Tree): A phylogenetic tree object.
    total_branch_length (float): Total desired sum of branch lengths.

    Returns:
    Tree: The same tree but with branch lengths assigned.
    """
    num_leaves = {}
    for node in tree.traverse():
        num_leaves[node] = len(node.get_leaves())
    total_leaves = sum(
        num_leaves[node] for node in tree.traverse() if not node.is_root()
    )

    for node in tree.traverse():
        if not node.is_root():
            # Proportional length assignment
            node.dist = total_branch_length / total_leaves

    return tree

# 1. Outliers: add them by % of total terminals, branch lengths drawn
# from the same min max of og tree, and attachment is random, keeping distances between OG clusters constant


def split_number(n, seed=None):
    if seed is not None:
        random.seed(seed)
    epsilon = 0.1
    first_part = round(random.uniform(epsilon, n - epsilon), 1)
    return (first_part, round(n - first_part, 1))


def random_supertree_internal_node(tree, seed=None):
    if seed is not None:
        random.seed(seed)
    internal_nodes = [
        node
        for node in tree.traverse()
        if not node.is_leaf()
        and hasattr(node, "is_cluster")
        and node.is_cluster is False
        and not hasattr(node, "has_outlier")
    ]
    if not internal_nodes:
        raise ValueError("No suitable internal nodes found in the tree.")
    return random.choice(internal_nodes)


counter = 0


def add_outlier(tree, parent_node, min_branch_length=1, max_branch_length=10, seed=None):
    global counter
    if seed is not None:
        random.seed(seed)
    counter += 1
    child_node = random.choice(parent_node.children)
    outlier_branch_length = random.randint(min_branch_length, max_branch_length)
    original_parent_child_distance = child_node.dist

    split_distance = split_number(original_parent_child_distance)
    length_parent_new_child = split_distance[0]
    length_new_child_old_child = split_distance[1]

    child_node.detach()

    intermediate_name = f"Intermediate_{counter}"
    new_node = parent_node.add_child(
        name=intermediate_name, dist=length_parent_new_child
    )
    new_node.is_cluster = False

    new_node.add_child(child_node, dist=length_new_child_old_child)

    outlier_name = f"Outlier_{counter}"
    outlier_node = new_node.add_child(name=outlier_name, dist=outlier_branch_length)
    outlier_node.is_outlier = True

    parent_node.has_outlier = True

# 2. Sampling. we take the OG tree and sample a % of the terminal nodes,
# ensuring a binary structure, keeping length constant


def sample_tree(original_tree, percentage, seed=None):
    if seed is not None:
        random.seed(seed)
    leaves = original_tree.get_leaves()
    total_leaves = len(leaves)
    num_to_sample = int(total_leaves * percentage)

    sampled_leaves = random.sample(leaves, num_to_sample)
    sampled_leaf_names = set(leaf.name for leaf in sampled_leaves)
    pruned_tree = original_tree.copy()
    pruned_tree.prune(sampled_leaf_names, preserve_branch_length=True)

    return pruned_tree

# 3. Branch Length noise, apply small distribution noise to all branch
# lengths, extent determined by some factor alpha


def add_noise_to_branch_lengths(tree, noise_level=0.1,seed=None):
    if seed is not None:
        random.seed(seed)
    for node in tree.traverse():
        if not node.is_root():
            noise = random.gauss(
                0, noise_level
            )
            scaled_noise = noise * node.dist * noise_level
            node.dist += scaled_noise
            if node.dist < 0:
                node.dist = 0 
    return tree
