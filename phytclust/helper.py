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
