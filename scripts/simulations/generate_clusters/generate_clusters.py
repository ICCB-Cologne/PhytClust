import io
import os
from Bio import Phylo
import random
from ete3 import Tree, TreeStyle
import string


def generate_random_tree(ntips, min_branch_length=1, max_branch_length=10,seed=None):
    """
    Generate a random binary phylogenetic tree with specified number of tips and random integer branch lengths.

    Parameters:
    ntips (int): Number of terminal nodes (tips) in the tree.
    min_branch_length (int): Minimum branch length (inclusive).
    max_branch_length (int): Maximum branch length (inclusive).

    Returns:
    Tree: A randomly generated binary phylogenetic tree.
    """
    if seed is not None:
        random.seed(seed)
    if ntips <= 0:
        raise ValueError("Number of tips (ntips) must be a positive integer.")

    tree = Tree()
    tree.populate(ntips)

    leaf_counter = 1
    internal_counter = 1
    for node in tree.traverse():
        node.is_cluster = False
        if not node.is_root():
            node.dist = random.randint(min_branch_length, max_branch_length)
        if node.is_leaf():
            node.name = f"terminal_{leaf_counter}"
            leaf_counter += 1
        if not node.is_leaf():
            node.name = f"internal_{internal_counter}"
            internal_counter += 1

    return tree


def replace_terminals_with_subtrees(super_tree, subtrees):
    """
    Replace each terminal node in the super tree with a corresponding subtree,
    preserving the original branch lengths to the new subtree roots, and
    marking all nodes in the attached subtrees with an 'in_cluster' attribute set to True.

    Parameters:
    super_tree (Tree): The super tree where terminal nodes will be replaced.
    subtrees (list of Tree): A list of subtrees that will replace each terminal node in the super tree.
    """
    terminal_nodes = list(
        super_tree.iter_leaves()
    )

    if len(terminal_nodes) != len(subtrees):
        raise ValueError(
            "The number of subtrees must match the number of terminal nodes in the super tree."
        )

    for terminal, subtree in zip(terminal_nodes, subtrees):
        original_branch_length = terminal.dist

        parent = terminal.up
        terminal.detach()

        subtree_root = subtree.detach()
        subtree_root.dist = original_branch_length

        parent.add_child(subtree_root)

        for node in subtree_root.traverse():
            node.in_cluster = True

    return super_tree
