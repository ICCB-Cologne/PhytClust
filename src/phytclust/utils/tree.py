from __future__ import annotations

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from Bio.Phylo.BaseTree import Tree, Clade


def get_pairwise_distances(
    tree: Tree,
    mode: str = "terminals",
    as_dataframe: bool = False,
    mrca: Optional[str] = None,
) -> Union[np.ndarray, pd.DataFrame]:
    """
    Pairwise distances for a subset of nodes.
    mode ∈ {"terminals","nonterminals","all"}.
    If mrca is given, restrict to descendants of that named clade.
    """
    if mrca:
        clade = next((cl for cl in tree.find_clades(name=mrca)), None)
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
        raise ValueError("mode must be one of {'terminals','nonterminals','all'}")

    n = len(terms)
    dist = np.zeros((n, n), dtype=float)

    iu, ju = np.triu_indices(n, k=1)
    for i, j in zip(iu, ju):
        bl = float(tree.distance(terms[i], terms[j]))
        dist[i, j] = bl
        dist[j, i] = bl

    if as_dataframe:
        names = [getattr(t, "name", f"node_{i}") for i, t in enumerate(terms)]
        return pd.DataFrame(dist, index=names, columns=names)
    return dist


def get_parent(tree: Tree, child: Clade) -> Optional[Clade]:
    path = tree.get_path(child)
    return path[-2] if len(path) > 1 else None


def count_branches_in_clusters(clusters: dict) -> int:
    """
    clusters: {cluster_id: [clades]}
    Total branches ~ 2*N-2 per cluster with N>1.
    """
    branch_count = 0
    for clades in clusters.values():
        N = len(clades)
        if N > 1:
            branch_count += 2 * N - 2
    return branch_count


def find_all_min_indices(arr: List[float]) -> Tuple[List[int], float]:
    if not arr:
        return [], float("inf")
    min_value = float("inf")
    min_indices: List[int] = []
    for i, value in enumerate(arr):
        if value < min_value:
            min_value = value
            min_indices = [i]
        elif value == min_value:
            min_indices.append(i)
    return min_indices, min_value


def rename_internal_nodes(tree: Tree) -> None:
    """
    Rename every internal node to 'internal_X' with a running counter.

    Notes
    -----
    This is a generic helper for ad hoc renaming and is NOT part of the
    main PhytClust preprocessing / DP pipeline. The official renaming
    used during clustering is in ``phytclust.validation.rename_nodes``.
    """
    internal_node_count = 1
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = f"internal_{internal_node_count}"
            internal_node_count += 1
