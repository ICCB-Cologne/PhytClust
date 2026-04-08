from __future__ import annotations

import logging
import string
from collections import deque
from typing import Any, Optional, Tuple

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from .exceptions import InvalidTreeError

logger = logging.getLogger(__name__)


def root_tree_at_taxon(tree: Tree, root_taxon: Optional[str]) -> Tree:
    """
    Root tree at a specified taxon or use midpoint rooting if root_taxon is "midpoint".

    Parameters
    ----------
    tree : Tree
        The phylogenetic tree (modified in-place).
    root_taxon : str, optional
        Name of taxon to root on. If "midpoint", use midpoint rooting.
        If None, tree is used as-is (assumed already rooted).

    Returns
    -------
    tree : Tree
        The rooted tree.
    """
    if root_taxon is None:
        return tree

    if root_taxon == "midpoint":
        tree.root_at_midpoint()
        return tree

    # Root at specified taxon
    if not is_outgroup_valid(tree, root_taxon):
        raise InvalidTreeError(f"Root taxon '{root_taxon}' not found in the tree.")
    tree.root_with_outgroup(tree.find_clades(name=root_taxon).__next__())
    return tree


def validate_and_set_outgroup(
    tree: Tree,
    outgroup: Optional[str],
    root_taxon: Optional[str] = None,
) -> Tuple[Tree, Optional[str]]:
    """
    Optionally root the tree, validate presence of outgroup (if provided),
    and normalize node names/branch lengths.
    Returns (tree, outgroup).
    """
    tree = root_tree_at_taxon(tree, root_taxon)

    if outgroup and not is_outgroup_valid(tree, outgroup):
        raise InvalidTreeError(f"Outgroup '{outgroup}' not found in the tree.")
    validate_tree(tree, outgroup)
    rename_nodes(tree, outgroup)
    ensure_branch_lengths(tree)
    return tree, outgroup


def prune_outgroup(
    tree: Tree, outgroup: Optional[str]
) -> Tuple[dict[Any, list[Any]], dict[Any, int]]:
    """
    Return mappings after pruning the outgroup from a COPY of the tree.
    If the outgroup is at root with two children, we keep the sibling as the new root.
    Otherwise we simply prune the outgroup clade.
    """
    if outgroup is None:
        node_terminals = {node: node.get_terminals() for node in tree.find_clades()}
        terminal_count = {node: len(terms) for node, terms in node_terminals.items()}
        return node_terminals, terminal_count

    outgroup_clade = next((cl for cl in tree.find_clades(name=outgroup)), None)
    if outgroup_clade is None:
        raise InvalidTreeError(
            f"Outgroup '{outgroup}' not found during prune_outgroup()."
        )

    # If outgroup is a direct child of root and root is bifurcating,
    # set sibling as new root - else, prune normally.
    if tree.root and len(tree.root.clades) == 2 and outgroup_clade in tree.root.clades:
        sibling = (
            tree.root.clades[0]
            if tree.root.clades[1] is outgroup_clade
            else tree.root.clades[1]
        )
        tree.root = sibling
    else:
        tree.prune(outgroup_clade)

    node_terminals = {node: node.get_terminals() for node in tree.find_clades()}
    terminal_count = {node: len(terms) for node, terms in node_terminals.items()}
    return node_terminals, terminal_count


def is_outgroup_valid(tree: Tree, outgroup: str) -> bool:
    """True if any clade in the tree has name == outgroup."""
    return next(tree.find_clades(name=outgroup), None) is not None


def validate_tree(tree: Tree, outgroup: Optional[str] = None) -> None:
    """
    Collapse single-child chains and log polytomous nodes (handled natively by DP).
    """
    # collapse single-child chains first
    merge_single_child_clades(tree)

    polytomy_nodes = []
    for node in tree.get_nonterminals():
        children = node.clades
        if len(children) > 2 and all(
            getattr(c, "name", None) != outgroup for c in children
        ):
            polytomy_nodes.append(node)

    if polytomy_nodes:
        logger.info(
            "Polytomous nodes (handled natively by DP convolution): "
            + ", ".join(
                f"{node.name or '?'} (children={len(node.clades)})"
                for node in polytomy_nodes
            )
        )


def rename_nodes(tree: Tree, outgroup: Optional[str] = None) -> None:
    """
    Give unique names to all nodes; keep the outgroup name as-is if unique.
    """
    if outgroup:
        outgroup_clades = list(tree.find_clades(name=outgroup))
        outgroup_count = len(outgroup_clades)
    else:
        outgroup_count = 0

    node_names = set([outgroup]) if outgroup else set()
    internal_node_counter = 0

    for node in tree.get_nonterminals() + tree.get_terminals():
        name = getattr(node, "name", None)

        if not name or (outgroup and name == outgroup and outgroup_count > 1):
            while True:
                new_name = "internal_node_" + (
                    string.ascii_uppercase[internal_node_counter % 26]
                    + str(internal_node_counter // 26)
                )
                internal_node_counter += 1
                if new_name not in node_names:
                    node.name = new_name
                    break

        elif (not outgroup or name != outgroup) and name in node_names:
            suffix = 1
            base = name
            new_name = f"{base}_{suffix}"
            while new_name in node_names:
                suffix += 1
                new_name = f"{base}_{suffix}"
            logger.info(f"Node name '{base}' is duplicated. Renaming to '{new_name}'")
            node.name = new_name

        node_names.add(node.name)


def merge_single_child_clades(tree: Tree) -> None:
    """
    Collapse chains of single-child clades by summing branch lengths.
    """
    queue: deque[Clade] = deque([tree.root])
    while queue:
        clade = queue.popleft()
        while len(clade.clades) == 1:
            child = clade.clades[0]
            clade.name = getattr(child, "name", clade.name)
            clade.branch_length = (clade.branch_length or 0.0) + (
                child.branch_length or 0.0
            )
            clade.clades = child.clades
        queue.extend(clade.clades)


def resolve_polytomies(tree: Tree) -> Tree:
    """
    Break nodes with >2 children by inserting 0-length dummy internal nodes until binary.
    """
    to_visit: deque[Clade] = deque([tree.root])
    while to_visit:
        node = to_visit.popleft()
        while len(node.clades) > 2:
            new_clade = Clade(
                branch_length=0.0, clades=[node.clades.pop(0), node.clades.pop(0)]
            )
            new_clade.comment = "DUMMY_NODE"
            node.clades.append(new_clade)
            to_visit.append(new_clade)
        to_visit.extend(node.clades)
    return tree


def ensure_branch_lengths(tree: Tree) -> None:
    """
    If the tree has no positive branch lengths (i.e. all are None or 0),
    set all non-root branches to length 1.0
    """
    clades = [cl for cl in tree.find_clades() if cl is not tree.root]
    has_positive_len = any((cl.branch_length or 0.0) > 0.0 for cl in clades)

    if has_positive_len:
        return

    logger.warning(
        "Tree has no branch lengths. "
        "PhytClust will assume all branches have length 1.0."
    )

    for cl in clades:
        if cl.branch_length is None or cl.branch_length == 0.0:
            cl.branch_length = 1.0
