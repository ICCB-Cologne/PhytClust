import string


def is_outgroup_valid(tree, outgroup):
    for clade in tree.root.clades:
        if clade.name == outgroup:
            return True
    return False


def validate_tree(tree, outgroup=None):
    invalid_nodes = [
        node
        for node in tree.get_nonterminals()
        if len(node.clades) != 2 and node.name != outgroup
    ]

    if invalid_nodes:
        print("Nodes with > 2 children:(excluding specified outgroup)")
        for node in invalid_nodes:
            print(f"Node: {node.name}, Children: {len(node.clades)}")

        raise AssertionError("All internal nodes should have 2 children")

    node_names = set([outgroup]) if outgroup else set()
    internal_node_counter = 0
    for node in tree.get_nonterminals() + tree.get_terminals():
        if not node.name or (
            node.name == outgroup
            and len([n for n in tree.find_clades() if n.name == outgroup]) > 1
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
