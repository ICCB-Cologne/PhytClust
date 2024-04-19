import string


def is_outgroup_valid(tree, outgroup):
    return any([clade.name == outgroup for clade in tree.root.clades])


def validate_tree(tree, outgroup=None):
    invalid_nodes = [
        node
        for node in tree.get_nonterminals()
        if len(node.clades) != 2 and node.name != outgroup
    ]
    assert not invalid_nodes, f"Nodes must have 2 children. Violating nodes: {invalid_nodes}"


def rename_nodes(tree, outgroup=None):
    node_names = set([outgroup]) if outgroup else set()
    internal_node_counter = 0
    for node in tree.get_nonterminals() + tree.get_terminals():
        if not node.name or (
            node.name == outgroup
            and len(list(tree.find_clades(outgroup))) > 1
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
            print(
                f"Node name '{node.name}' is duplicated. Adding suffix to make it '{new_name}'"
            )
            node.name = new_name

        node_names.add(node.name)
