import numpy as np
from math import ceil
from typing import Any, Optional
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


def _compute_polytomy_dp(node, pc, max_states_global, min_cluster_size,
                          outlier_thresh, outlier_penalty, dtype, bp_dtype):
    """
    Compute DP table for a polytomous (k>2 children) node by sequentially
    convolving child DP tables. Equivalent to trying all (2k-3)!! binary
    resolutions and keeping the best at each k — but O(n_children * max_k^2).

    Parameters
    ----------
    node : Clade
        The polytomous node (has > 2 children).
    pc : PhytClust
        The clustering object.
    max_states_global : int
        Maximum k+1 clusters to compute for.
    min_cluster_size : int
        Minimum leaves per cluster.
    outlier_thresh : int or None
        Threshold for outlier penalty.
    outlier_penalty : float
        Cost to add when isolating a small cluster.
    dtype : numpy dtype
        Data type for DP table.
    bp_dtype : numpy dtype
        Data type for backpointer array.

    Returns
    -------
    dp_array : numpy array
        DP table for this node.
    step_backptrs : list of numpy arrays
        Backpointers from each convolution step.
    """
    children = node.clades
    n_leaves = pc.num_leaves_per_node[node]
    n_states = min(n_leaves, max_states_global)

    # Start: combined = first child's DP table
    first_child = children[0]
    first_id = pc.node_to_id[first_child]
    combined = pc.dp_table[first_id].copy()
    n_combined = pc.num_leaves_per_node[first_child]

    step_backptrs = []  # step_backptrs[i][k] = clusters from first i+1 children

    for child in children[1:]:
        child_id = pc.node_to_id[child]
        child_dp = pc.dp_table[child_id]
        n_child = pc.num_leaves_per_node[child]
        n_total = n_combined + n_child
        n_new_states = min(n_total, max_states_global)

        new_combined = np.full(n_new_states + 1, np.inf, dtype=dtype)
        step_bp = np.full(n_new_states + 1, -1, dtype=bp_dtype)

        # k=0: all in one cluster (no edge cost added for polytomous branch)
        if n_total >= min_cluster_size:
            new_combined[0] = combined[0] + child_dp[0]
        step_bp[0] = 0

        # k>=1: k+1 total clusters via convolution
        for k in range(1, n_new_states + 1):
            max_i = min(k, len(combined) - 1)
            min_i = max(0, k - 1 - (len(child_dp) - 1))
            if min_i > max_i:
                continue
            i_vals = np.arange(min_i, max_i + 1)
            j_vals = k - 1 - i_vals
            scores = combined[i_vals] + child_dp[j_vals]

            # Outlier penalty: when combined (or child) is cut as single cluster
            if outlier_thresh is not None and outlier_penalty > 0:
                if i_vals[0] == 0 and n_combined < outlier_thresh:
                    scores[0] = scores[0] + outlier_penalty
                if j_vals[-1] == 0 and n_child < outlier_thresh:
                    scores[-1] = scores[-1] + outlier_penalty

            best = int(np.argmin(scores))
            new_combined[k] = scores[best]
            step_bp[k] = i_vals[best]

        step_backptrs.append(step_bp)
        combined = new_combined
        n_combined = n_total

    # Apply min_cluster_size at node level
    if n_leaves < min_cluster_size:
        combined[0] = np.inf

    return combined, step_backptrs


def prepare_tree(pc) -> None:
    pc.tree, pc.outgroup = validate_and_set_outgroup(
        pc.tree, pc.outgroup,
        optimize_polytomies=getattr(pc, 'optimize_polytomies', True),
        root_taxon=getattr(pc, 'root_taxon', None)
    )
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
    pc.polytomy_backptr = [None] * num_nodes  # for polytomous nodes only
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

            pc.dp_table[node_id] = dp_array
            pc.backptr[node_id] = backptr_array
            continue

        # Handle polytomous vs binary nodes
        if len(node.clades) > 2 and getattr(pc, 'optimize_polytomies', True):
            # Polytomy: use sequential convolution
            dp_array, step_backptrs = _compute_polytomy_dp(
                node, pc, max_states_global, min_cluster_size,
                outlier_thresh, outlier_penalty, dtype, bp_dtype
            )
            pc.dp_table[node_id] = dp_array
            pc.polytomy_backptr[node_id] = step_backptrs
            # Free children DP tables
            for child in node.clades:
                pc.dp_table[pc.node_to_id[child]] = None
        else:
            # Binary node (or optimize_polytomies=False): existing logic
            left, right = node.clades[0], node.clades[1]
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

            # Check if both branches are zero-length and no_split_zero_length is enabled
            no_split_zero = getattr(pc, "no_split_zero_length", False)
            both_zero_length = (len_left == 0 and len_right == 0) and no_split_zero

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

            backptr_array[0, 0] = 0
            backptr_array[1, 0] = 0

            for k in range(1, n_states + 1):
                # If no_split_zero_length is enabled and both branches are zero-length,
                # only k=0 (no split) is allowed; all k >= 1 remain infeasible
                if both_zero_length:
                    dp_array[k] = np.inf
                    backptr_array[0, k] = -1
                    backptr_array[1, k] = -1
                    continue

                max_i = min(k, len(left_dp) - 1)
                min_i = max(0, k - 1 - (len(right_dp) - 1))

                if min_i > max_i:
                    continue

                i_vals = np.arange(min_i, max_i + 1)
                j_vals = k - 1 - i_vals

                scores = left_dp[i_vals] + right_dp[j_vals]

                # Apply outlier penalty at split point: when left or right is cut as a single cluster
                if outlier_thresh is not None and outlier_penalty > 0:
                    # Left being cut off as a single cluster (i=0) and it's small
                    if i_vals[0] == 0 and n_left < outlier_thresh:
                        scores[0] = scores[0] + outlier_penalty
                    # Right being cut off as a single cluster (j=0) and it's small
                    if j_vals[-1] == 0 and n_right < outlier_thresh:
                        scores[-1] = scores[-1] + outlier_penalty

                best = np.argmin(scores)
                dp_array[k] = scores[best]
                backptr_array[0, k] = i_vals[best]
                backptr_array[1, k] = j_vals[best]

            pc.dp_table[node_id] = dp_array
            pc.backptr[node_id] = backptr_array

            pc.dp_table[left_id] = None
            pc.dp_table[right_id] = None

    pc.postorder_nodes = nodes


def backtrack(pc, k: int, *, verbose: bool = False) -> dict[Any, int]:
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

    if cluster_index >= len(root_dp):
        raise ValueError(
            f"Cannot partition into {k} clusters — exceeds tree capacity."
        )

    # Warn if solution violates min_cluster_size constraint
    min_cluster_size = getattr(pc, "min_cluster_size", 1)
    if np.isinf(root_dp[cluster_index]):
        import warnings
        warnings.warn(
            f"Partition into {k} clusters violates min_cluster_size={min_cluster_size}. "
            f"Using best available solution (some clusters < {min_cluster_size}).",
            UserWarning
        )

    clusters: dict[Any, int] = {}
    current_cluster_id = 0
    stack = [(root_id, cluster_index)]

    while stack:
        node_id, c_index = stack.pop()
        node = pc.postorder_nodes[node_id]

        if verbose:
            print(f"Visiting node {getattr(node, 'name', '')} with c_index={c_index}")

        if c_index == 0:
            # All leaves in this clade form one cluster
            for t in node.get_terminals():
                clusters[t] = current_cluster_id
            current_cluster_id += 1
        elif pc.polytomy_backptr[node_id] is not None:
            # Polytomous node: unwind the convolution steps in reverse
            step_backptrs = pc.polytomy_backptr[node_id]
            n = c_index
            child_k_indices = []
            # step_backptrs[s][n] = clusters from first s+1 children
            for step_bp in reversed(step_backptrs):
                j = int(step_bp[n])        # clusters from combined part
                child_k_indices.append(n - j)  # this child gets n-j clusters
                n = j
            child_k_indices.append(n)      # first child gets n clusters
            child_k_indices.reverse()
            for child, ci in zip(node.clades, child_k_indices):
                stack.append((pc.node_to_id[child], ci))
        else:
            # Binary node: unchanged
            left_k = pc.backptr[node_id][0, c_index]
            right_k = pc.backptr[node_id][1, c_index]

            if left_k < 0 or right_k < 0:
                raise RuntimeError("Back-pointer missing - fatal error. Check DP.")

            left, right = node.clades[0], node.clades[1]
            stack.append((pc.node_to_id[right], int(right_k)))
            stack.append((pc.node_to_id[left], int(left_k)))

    if current_cluster_id != k:
        raise ValueError(
            f"Number of clusters found: {current_cluster_id}, expected: {k}"
        )

    return clusters


def cluster_map(pc, k: int) -> dict[Any, int]:
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
