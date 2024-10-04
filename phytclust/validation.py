import string
import logging
from Bio import Phylo
from collections import deque
from Bio.Phylo.BaseTree import Clade
from typing import Any, Optional, List

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_and_set_outgroup(tree: Any, outgroup: Optional[Any]) -> Optional[Any]:
    """
    Validates and sets the outgroup for the tree.

    Parameters:
    tree (Any): The phylogenetic tree.
    outgroup (Optional[Any]): The outgroup to be set.

    Raises:
    ValueError: If the outgroup is not valid.

    Returns:
    Optional[Any]: The validated outgroup.
    """
    if outgroup and not is_outgroup_valid(tree, outgroup):
        raise ValueError("Outgroup not found, please check input.")
    validate_tree(tree, outgroup)
    rename_nodes(tree, outgroup)
    return tree, outgroup


def prune_outgroup(tree: Any, outgroup: Optional[Any]) -> None:
    """
    Prunes the outgroup from the tree.

    Parameters:
    tree (Any): The phylogenetic tree.
    outgroup (Optional[Any]): The outgroup to be pruned.

    Returns:
    None
    """
    outgroup_clade: Optional[Any] = next(tree.find_clades(outgroup), None)
    if len(tree.root.clades) > 2:
        tree.prune(outgroup_clade)
    elif len(tree.root.clades) == 2:
        outgroup_clade = next((clade for clade in tree.root.clades if clade.name == outgroup), None)
        sibling_clade: Optional[Any] = next((clade for clade in tree.root.clades if clade.name != outgroup), None)
        if sibling_clade:
            sibling_clade.name = outgroup
        tree.root = sibling_clade

    node_terminals = {
        node: node.get_terminals() for node in tree.find_clades()
    }
    terminal_count = {
        node: len(terminals)
        for node, terminals in node_terminals.items()
        if node != outgroup_clade
    }
    return node_terminals, terminal_count


def is_outgroup_valid(tree: Any, outgroup: Any) -> bool:
    """
    Checks if outgroup exists in the tree.

    Parameters:
    tree (Any): The phylogenetic tree.
    outgroup (Any): The outgroup to be checked.

    Returns:
    bool: True if the outgroup exists, False otherwise.
    """
    return any([clade.name == outgroup for clade in tree.root.clades])


def validate_tree(tree: Any, outgroup: Optional[Any] = None) -> None:
    """
    Validate the tree structure and resolve polytomies.

    Checks for nodes with more than two children (excluding the specified outgroup)
    and resolves polytomies by merging single child clades and creating new internal nodes with branch length 0.

    Parameters:
    tree (Any): The phylogenetic tree.
    outgroup (Optional[Any]): The name of the outgroup. Defaults to None.

    Returns:
    None
    """
    invalid_nodes = [
        node
        for node in tree.get_nonterminals()
        if len(node.clades) != 2
        and all(clade.name != outgroup for clade in node.clades)
    ]
    for clade in tree.root.clades:
        if len(clade.clades) != 2 and all(
            child.name != outgroup for child in clade.clades
        ):
            invalid_nodes.append(clade)

    if invalid_nodes:
        logger.info(
            "Nodes with > 2 children (excluding specified outgroup): "
            + ", ".join(
                [
                    f"Node: {node.name}, Children: {len(node.clades)}"
                    for node in invalid_nodes
                ]
            )
        )

        logger.info("Resolving polytomies and merging single child clades...")
        merge_single_child_clades(tree)


def rename_nodes(tree: Any, outgroup: Optional[Any] = None) -> None:
    """
    Rename nodes in the tree to ensure unique names.

    This function renames internal nodes and resolves name conflicts by adding suffixes.

    Parameters:
    tree (Any): The phylogenetic tree.
    outgroup (Optional[Any]): The name of the outgroup. Defaults to None.

    Returns:
    None
    """
    node_names = set([outgroup]) if outgroup else set()
    internal_node_counter = 0

    for node in tree.get_nonterminals() + tree.get_terminals():
        if not node.name or (
            node.name == outgroup and len(list(tree.find_clades(outgroup))) > 1
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
        elif node.name != outgroup and node.name in node_names:
            suffix = 1
            new_name = f"{node.name}_{suffix}"
            while new_name in node_names:
                suffix += 1
                new_name = f"{node.name}_{suffix}"
            logger.info(
                f"Node name '{node.name}' is duplicated. Adding suffix to make it '{new_name}'"
            )
            node.name = new_name

        node_names.add(node.name)


def merge_single_child_clades(tree: Any) -> None:
    """
    Merge clades that have only one child recursively, combining their branch lengths.

    Parameters:
    tree (Any): The phylogenetic tree.

    Returns:
    None
    """
    queue = deque([tree.root])
    while queue:
        clade = queue.popleft()
        while len(clade.clades) == 1:
            single_child = clade.clades[0]
            clade.name = single_child.name
            clade.branch_length = (clade.branch_length or 0) + (
                single_child.branch_length or 0
            )
            clade.clades = single_child.clades
        queue.extend(clade.clades)


def resolve_polytomies(tree: Any) -> Any:
    """
    Resolve polytomies in the tree by creating new internal nodes until each node has at most two children.

    Parameters:
    tree (Any): The phylogenetic tree.

    Returns:
    Any: The modified tree with resolved polytomies.
    """
    to_visit = deque([tree.root])
    while to_visit:
        node = to_visit.popleft()
        while len(node.clades) > 2:
            new_clade = Clade(
                branch_length=0.0, clades=[node.clades.pop(0), node.clades.pop(0)]
            )
            node.clades.append(new_clade)
            to_visit.append(new_clade)
        to_visit.extend(node.clades)
    return tree
