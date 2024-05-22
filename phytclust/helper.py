def get_parent(tree, child):
    # tree = self._no_outgroup_tree if self._no_outgroup_tree else self.tree
    node_path = tree.get_path(child)
    if len(node_path) > 1:
        return node_path[-2]
    else:
        return None


def count_branches_in_clusters(clusters):
    branch_count = 0
    for _, clades in clusters.items():
        N = len(clades)
        if N > 1:
            branch_count += 2 * N - 2
    return branch_count


def find_all_min_indices(arr):
    if len(arr) == 0:
        return []

    min_value = float("inf")
    min_indices = []

    for i, value in enumerate(arr):
        if value < min_value:
            min_value = value
            min_indices = [i]  # Start a new list of indices
        elif value == min_value:
            min_indices.append(i)  # Add index to existing list

    return min_indices, min_value


def rename_internal_nodes(tree):
    """
    Renames all internal nodes of a given Phylo tree object.
    Each internal node will be named as 'internal_X' where X is an incrementing integer.

    :param tree: Bio.Phylo tree object
    """
    internal_node_count = 1
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = f"internal_{internal_node_count}"
            internal_node_count += 1
