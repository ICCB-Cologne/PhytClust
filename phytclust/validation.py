import string
from Bio import Phylo

from Bio.Phylo.BaseTree import Clade


def is_outgroup_valid(tree, outgroup):
    return any([clade.name == outgroup for clade in tree.root.clades])


def validate_tree(tree, outgroup=None):
    invalid_nodes = [
        node
        for node in tree.get_nonterminals()
        if len(node.clades) != 2
        and all(clade.name != outgroup for clade in node.clades)
    ]
    assert not invalid_nodes, f"Nodes must have 2 children. Violating nodes: {invalid_nodes}"

        print("Resolving polytomies and merging single child clades...")
        merge_single_child_clades(tree)
        resolve_polytomies(tree)



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


def merge_single_child_clades(tree):
    # Recursively merge clades with a single child
    def merge_clades(clade):
        if len(clade.clades) == 1:
            # Merge with the single child
            single_child = clade.clades[0]
            clade.name = single_child.name
            clade.branch_length = (clade.branch_length or 0) + (
                single_child.branch_length or 0
            )
            clade.clades = single_child.clades
            # Continue merging for the new clade
            merge_clades(clade)
        else:
            # Recursively merge children
            for child in clade.clades:
                merge_clades(child)

    # Start merging from the root
    merge_clades(tree.root)


def resolve_polytomies(tree):
    stack = [tree.root]

    while stack:
        clade = stack.pop()

        # While there are more than 2 children, we need to resolve it
        while len(clade.clades) > 2:
            # Select two children to create a new internal node
            selected_clades = [clade.clades.pop(0), clade.clades.pop(0)]

            # Create a new internal node with zero branch length
            new_internal_node = Clade(branch_length=0.0)
            new_internal_node.clades.extend(selected_clades)

            # Add the new internal node back to the current clade
            clade.clades.append(new_internal_node)

        # Add children to stack to resolve their polytomies
        stack.extend(clade.clades)

    return tree
