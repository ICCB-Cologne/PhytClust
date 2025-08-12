from __future__ import annotations

import logging
import string
from collections import deque
from typing import Any, Optional, Tuple, Dict, List

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

logger = logging.getLogger(__name__)


def validate_and_set_outgroup(
    tree: Tree, outgroup: Optional[str]
) -> Tuple[Tree, Optional[str]]:
    """
    Validate presence of outgroup (if provided), normalize names, resolve polytomies if needed.
    Returns (tree, outgroup).
    """
    if outgroup and not is_outgroup_valid(tree, outgroup):
        raise ValueError(f"Outgroup '{outgroup}' not found in the tree.")
    validate_tree(tree, outgroup)
    rename_nodes(tree, outgroup)
    return tree, outgroup


def prune_outgroup(
    tree: Tree, outgroup: Optional[str]
) -> Tuple[Dict[Any, List[Any]], Dict[Any, int]]:
    """
    Return mappings after pruning the outgroup from a COPY of the tree.
    If the outgroup is at root with two children, we keep the sibling as the new root.
    Otherwise we simply prune the outgroup clade.
    """
    if outgroup is None:
        # nothing to prune
        node_terminals = {node: node.get_terminals() for node in tree.find_clades()}
        terminal_count = {node: len(terms) for node, terms in node_terminals.items()}
        return node_terminals, terminal_count

    # locate the clade by name
    outgroup_clade = next((cl for cl in tree.find_clades(name=outgroup)), None)
    if outgroup_clade is None:
        raise ValueError(f"Outgroup '{outgroup}' not found during prune_outgroup().")

    # If outgroup is a direct child of root and root is bifurcating,
    # set sibling as new root; else, prune normally.
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
    Ensure binary branching except possibly at the outgroup, and collapse single-child chains.
    """
    invalid_nodes = [
        node
        for node in tree.get_nonterminals()
        if len(node.clades) != 2
        and all(getattr(c, "name", None) != outgroup for c in node.clades)
    ]
    # extra guard at root children
    for clade in getattr(tree.root, "clades", []):
        if len(clade.clades) != 2 and all(
            getattr(child, "name", None) != outgroup for child in clade.clades
        ):
            invalid_nodes.append(clade)

    if invalid_nodes:
        logger.info(
            "Nodes with > 2 children (excluding specified outgroup): "
            + ", ".join(
                [
                    f"{getattr(node, 'name', '?')} (children={len(node.clades)})"
                    for node in invalid_nodes
                ]
            )
        )
        logger.info("Resolving polytomies and merging single-child clades…")
        merge_single_child_clades(tree)
        resolve_polytomies(tree)


def rename_nodes(tree: Tree, outgroup: Optional[str] = None) -> None:
    """
    Give unique names to all nodes; keep the outgroup name as-is if unique.
    """
    node_names = set([outgroup]) if outgroup else set()
    internal_node_counter = 0

    for node in tree.get_nonterminals() + tree.get_terminals():
        if not getattr(node, "name", None) or (
            outgroup
            and node.name == outgroup
            and len(list(tree.find_clades(name=outgroup))) > 1
        ):
            while True:
                new_name = "internal_node_" + (
                    string.ascii_uppercase[internal_node_counter % 26]
                    + str(internal_node_counter // 26)
                )
                internal_node_counter += 1
                if new_name not in node_names:
                    node.name = new_name
                    break
        elif (not outgroup or node.name != outgroup) and node.name in node_names:
            suffix = 1
            base = node.name
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
