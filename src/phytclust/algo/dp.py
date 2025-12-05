import numpy as np
from math import ceil
from typing import Any, Dict, List, Optional
import logging

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


def _eff_length(node, pc):
    """
    If `pc.use_branch_support` is True, we augment the branch length
    with a penalty derived from the support value:

        eff_len = branch_length + support_weight * (-log(p))

    where p = max(confidence/100, pc.min_support).
    """
    base = node.branch_length or 0.0
    if not pc.use_branch_support:
        return base

    raw = node.confidence if node.confidence is not None else 100.0
    p = max(raw / 100.0, pc.min_support)
    penalty = -np.log(p)
    return base + pc.support_weight * penalty


def prepare_tree(pc) -> None:
    pc.tree, pc.outgroup = validate_and_set_outgroup(pc.tree, pc.outgroup)
    resolve_polytomies(pc.tree)
    pc.name_leaves_per_node = {n: n.get_terminals() for n in pc.tree.find_clades()}
    pc.num_leaves_per_node = {}
    for n in pc.tree.find_clades(order="postorder"):
        if n.is_terminal():
            pc.num_leaves_per_node[n] = 1
        else:
            left, right = n.clades
            pc.num_leaves_per_node[n] = (
                pc.num_leaves_per_node[left] + pc.num_leaves_per_node[right]
            )

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
    """
    Build DP tables with a minimum cluster size constraint.

    - pc.min_cluster_size (int >= 1): minimal number of leaves allowed in any final cluster.
    - pc.max_k (optional): cap on maximum number of clusters.
    """
    tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    nodes = list(tree.find_clades(order="postorder"))
    pc.postorder_nodes = nodes
    num_nodes = len(nodes)
    pc.node_to_id = {node: i for i, node in enumerate(nodes)}

    pc.dp_table = [None] * num_nodes
    pc.backptr = [None] * num_nodes
    pc.cluster_cost = {}

    dtype = np.float64 if pc.num_terminals > 80_000 else np.float32

    max_states_global = pc.max_k if pc.max_k is not None else pc.num_terminals
    if max_states_global < 1:
        raise ValueError("max_k (or implied max_states_global) must be ≥ 1.")

    bp_dtype = np.int16 if max_states_global <= 32767 else np.int32

    min_cluster_size = getattr(pc, "min_cluster_size", 1)
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be ≥ 1.")

    outlier_thresh = getattr(pc, "outlier_size_threshold", None)
    outlier_penalty = getattr(pc, "outlier_penalty", 0.0)

    for node in nodes:
        node_id = pc.node_to_id[node]
        n_leaves = pc.num_leaves_per_node[node]

        n_states = min(n_leaves, max_states_global)

        dp_array = np.full(n_states + 1, np.inf, dtype=dtype)
        backptr_array = np.full((2, n_states + 1), -1, dtype=bp_dtype)

        if node.is_terminal():
            pc.cluster_cost[node] = 0.0

            if n_leaves >= min_cluster_size:
                dp_array[0] = 0.0
            else:
                dp_array[0] = np.inf

            if (
                outlier_thresh is not None
                and outlier_penalty > 0
                and n_leaves < outlier_thresh
                and np.isfinite(dp_array[0])
            ):
                dp_array[0] += outlier_penalty

            pc.dp_table[node_id] = dp_array
            pc.backptr[node_id] = backptr_array
            continue

        left, right = node.clades
        left_id = pc.node_to_id[left]
        right_id = pc.node_to_id[right]

        left_dp = pc.dp_table[left_id]
        right_dp = pc.dp_table[right_id]

        if left_dp is None or right_dp is None:
            raise RuntimeError("Child DP table missing – compute_dp_table order bug.")

        n_left = pc.num_leaves_per_node[left]
        n_right = pc.num_leaves_per_node[right]

        len_left = _eff_length(left, pc)
        len_right = _eff_length(right, pc)

        cost_one_cluster = (
            left_dp[0] + right_dp[0] + n_left * len_left + n_right * len_right
        )

        if getattr(pc, "use_branch_support", False):
            raw_support = node.confidence if node.confidence is not None else 100.0
            support = max(raw_support / 100.0, pc.min_support)
            cost_one_cluster /= support

        pc.cluster_cost[node] = float(cost_one_cluster)

        if n_leaves >= min_cluster_size:
            dp_array[0] = cost_one_cluster
        else:
            dp_array[0] = np.inf

        if (
            outlier_thresh is not None
            and outlier_penalty > 0
            and n_leaves < outlier_thresh
            and np.isfinite(dp_array[0])
        ):
            dp_array[0] += outlier_penalty

        backptr_array[0, 0] = 0
        backptr_array[1, 0] = 0

        for k in range(1, n_states + 1):
            max_i = min(k, len(left_dp) - 1)
            min_i = max(0, k - 1 - (len(right_dp) - 1))

            if min_i > max_i:
                continue

            i_vals = np.arange(min_i, max_i + 1)
            j_vals = k - 1 - i_vals

            scores = left_dp[i_vals] + right_dp[j_vals]
            best = np.argmin(scores)
            dp_array[k] = scores[best]
            backptr_array[0, k] = i_vals[best]
            backptr_array[1, k] = j_vals[best]

        if getattr(node, "comment", None) == "DUMMY_NODE":
            dp_array[1:] = np.inf

        pc.dp_table[node_id] = dp_array
        pc.backptr[node_id] = backptr_array

        pc.dp_table[left_id] = None
        pc.dp_table[right_id] = None

    pc.postorder_nodes = nodes


def backtrack(pc, k: int, *, verbose: bool = False) -> Dict[Any, int]:
    if k is None:
        raise ValueError("value of k is missing.")
    if k <= 0:
        raise ValueError("k must be a positive integer.")
    if not getattr(pc, "_dp_ready", False):
        raise RuntimeError("DP table not computed. Call compute_dp_table(pc) first.")

    active_tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    root = active_tree.root
    root_id = pc.node_to_id[root]

    root_dp = pc.dp_table[root_id]
    if root_dp is None:
        raise RuntimeError(
            "DP table at root is missing. Did compute_dp_table(pc) fail?"
        )

    cluster_index = k - 1

    if cluster_index >= len(root_dp) or np.isinf(root_dp[cluster_index]):
        raise ValueError(
            f"No feasible partition into {k} clusters with "
            f"min_cluster_size={getattr(pc, 'min_cluster_size', 1)}."
        )

    clusters: Dict[Any, int] = {}
    current_cluster_id = 0
    stack = [(root_id, cluster_index)]

    while stack:
        node_id, c_index = stack.pop()
        node = pc.postorder_nodes[node_id]

        if verbose:
            print(f"Visiting node {getattr(node, 'name', '')} with c_index={c_index}")

        if c_index == 0:
            for t in node.get_terminals():
                clusters[t] = current_cluster_id
            current_cluster_id += 1
        else:
            left_k = pc.backptr[node_id][0, c_index]
            right_k = pc.backptr[node_id][1, c_index]

            if left_k < 0 or right_k < 0:
                raise RuntimeError("Back-pointer missing - fatal error. Check DP.")

            left, right = node.clades
            stack.append((pc.node_to_id[right], int(right_k)))
            stack.append((pc.node_to_id[left], int(left_k)))

    if current_cluster_id != k:
        raise ValueError(
            f"Number of clusters found: {current_cluster_id}, expected: {k}"
        )

    return clusters


def cluster_map(pc, k: int) -> Dict[Any, int]:
    if not getattr(pc, "_dp_ready", False):
        validate_args(pc)
        prepare_tree(pc)
        compute_dp_table(pc)

    pc._dp_ready = True
    if k <= 0:
        raise ValueError("k must be a positive integer.")

    if pc.clusters is None:
        pc.clusters = {}

    if k in pc.clusters:
        return pc.clusters[k]

    if not getattr(pc, "_dp_ready", False):
        compute_dp_table(pc)
        pc._dp_ready = True

    cmap = pc.get_clusters(k)
    pc.clusters[k] = cmap
    return cmap
