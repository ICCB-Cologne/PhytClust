import numpy as np
from math import ceil
from typing import Any, Dict, List, Optional

from ..validation import (
    validate_and_set_outgroup,
    prune_outgroup,
    resolve_polytomies,
)


def validate_args(pc) -> None:
    if pc.k is not None and pc.k < 1:
        raise ValueError("k must be ≥ 1 if provided.")
    if not 0 < pc.max_k_limit <= 1:
        raise ValueError("max_k_limit must be between 0 and 1")
    if pc.num_peaks < 1:
        raise ValueError("num_peaks must be ≥ 1")


def prepare_tree(pc) -> None:
    pc.tree, pc.outgroup = validate_and_set_outgroup(pc.tree, pc.outgroup)
    resolve_polytomies(pc.tree)
    pc.name_leaves_per_node = {n: n.get_terminals() for n in pc.tree.find_clades()}
    pc.num_leaves_per_node = {n: len(t) for n, t in pc.name_leaves_per_node.items()}
    root_cnt = pc.num_leaves_per_node[pc.tree.root]
    pc.num_terminals = root_cnt - 1 if pc.outgroup else root_cnt

    pc.max_k = (
        ceil(pc.num_terminals * pc.max_k_limit)
        if pc.k is None and pc.max_k is None
        else pc.max_k if pc.k is None else None
    )

    pc._tree_wo_outgroup = None
    if pc.outgroup:
        import copy

        pc._tree_wo_outgroup = copy.deepcopy(pc.tree)
        pc.name_leaves_per_node, pc.num_leaves_per_node = prune_outgroup(
            pc._tree_wo_outgroup, pc.outgroup
        )


def compute_dp_table(pc) -> None:
    tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    nodes = list(tree.find_clades(order="postorder"))
    pc.postorder_nodes = nodes
    num_nodes = len(nodes)
    pc.node_to_id = {node: i for i, node in enumerate(nodes)}

    pc.dp_table: List[Optional[np.ndarray]] = [None] * num_nodes
    pc.backptr: List[Optional[np.ndarray]] = [None] * num_nodes
    dtype = np.float64 if pc.num_terminals > 80_000 else np.float32

    for node in nodes:
        node_id = pc.node_to_id[node]
        nleaf = pc.num_leaves_per_node[node]
        dp_array = np.full(nleaf + 1, np.inf, dtype=dtype)
        backptr_array = np.full((2, nleaf + 1), -1, dtype=np.int32)

        if node.is_terminal():
            dp_array[0] = 0.0
            pc.dp_table[node_id] = dp_array
            pc.backptr[node_id] = backptr_array
            continue

        left, right = node.clades
        left_dp = pc.dp_table[pc.node_to_id[left]]
        right_dp = pc.dp_table[pc.node_to_id[right]]

        cost_one_cluster = (
            left_dp[0]
            + right_dp[0]
            + pc.num_leaves_per_node[left] * (left.branch_length or 0)
            + pc.num_leaves_per_node[right] * (right.branch_length or 0)
        )
        if pc.use_branch_support:
            support = (node.confidence or 100) / 100
            cost_one_cluster /= support

        dp_array[0] = cost_one_cluster
        backptr_array[0, 0] = 0
        backptr_array[1, 0] = 0

        for k in range(1, nleaf + 1):
            i_vals = np.arange(k)
            valid_mask = (i_vals < len(left_dp)) & ((k - 1 - i_vals) < len(right_dp))
            i_vals = i_vals[valid_mask]
            j_vals = k - 1 - i_vals

            if i_vals.size == 0:
                continue

            scores = left_dp[i_vals] + right_dp[j_vals]
            best = np.argmin(scores)
            dp_array[k] = scores[best]
            backptr_array[0, k] = i_vals[best]
            backptr_array[1, k] = j_vals[best]

        pc.dp_table[node_id] = dp_array
        pc.backptr[node_id] = backptr_array
        pc.dp_table[pc.node_to_id[left]] = None
        pc.dp_table[pc.node_to_id[right]] = None

    pc.postorder_nodes = nodes


def backtrack(pc, k: int, *, verbose: bool = False) -> Dict[Any, int]:
    if k is None:
        raise ValueError("value of k is missing.")
    if k <= 0:
        raise ValueError("k must be a positive integer.")

    cluster_index = k - 1
    active_tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    root = active_tree.root
    root_id = pc.node_to_id[root]

    root_dp = pc.dp_table[root_id]
    if k - 1 >= len(root_dp) or np.isinf(root_dp[k - 1]):
        raise ValueError(f"Cannot form {k} clusters on this tree.")

    clusters = {}
    current_cluster_id = 0
    stack = [(root_id, cluster_index)]

    while stack:
        node_id, c_index = stack.pop()
        node = pc.postorder_nodes[node_id]

        if verbose:
            print(f"Visiting node {getattr(node, 'name', '')} with c_index={c_index}")

        if c_index == 0:
            for t in pc.name_leaves_per_node[node]:
                clusters[t] = current_cluster_id
            current_cluster_id += 1
        else:
            left_k = pc.backptr[node_id][0, c_index]
            right_k = pc.backptr[node_id][1, c_index]

            if left_k < 0 or right_k < 0:
                raise RuntimeError("Back-pointer missing - fatal error. Check DP")

            left, right = node.clades
            stack.append((pc.node_to_id[right], right_k))
            stack.append((pc.node_to_id[left], left_k))

    if current_cluster_id != k:
        raise ValueError(
            f"Number of clusters found: {current_cluster_id}, expected: {k}"
        )

    return clusters


def cluster_map(pc, k: int) -> Optional[Dict[Any, int]]:
    if pc.clusters is None:
        pc.clusters = {}

    cmap = pc.clusters.get(k)
    if cmap is not None:
        return cmap

    if pc._dp_ready:
        return pc.get_clusters(k)

    return None
