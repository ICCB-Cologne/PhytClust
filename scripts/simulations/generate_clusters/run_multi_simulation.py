import os
import random
import math
import string
import json
from ete3 import Tree
import numpy as np
from copy import deepcopy
from utils import (
    distribute_ntips_into_mparts,
    convert_ete_to_phylo_newick,
    generate_ground_truth,
    flatten_clusters,
)

from generate_clusters import generate_random_tree, replace_terminals_with_subtrees
from tree_operations import (
    add_noise_to_branch_lengths,
    sample_tree,
    add_outlier,
    random_supertree_internal_node,
    split_number,
    assign_branch_lengths_based_on_leaves,
    total_branch_length,
    distribute_subtree_sums,
)


def generate_label(index):
    """Generate a label in the sequence A, B, ..., Z, AA, AB, ..., ZZ, AAA, ..."""
    label = ""
    while index >= 0:
        label = string.ascii_uppercase[index % 26] + label
        index = index // 26 - 1
    return label


def get_clusters(nodes, lookup):
    clusters = {}
    cluster_number = 1

    def assign_cluster(node, cluster_number):
        clusters[node.node.name] = cluster_number

        def assign_to_descendants(sub_node):
            sub_node.node.add_feature("in_cluster", True)
            if sub_node.left:
                assign_to_descendants(sub_node.left)
            if sub_node.right:
                assign_to_descendants(sub_node.right)

        if node.left:
            assign_to_descendants(node.left)
        if node.right:
            assign_to_descendants(node.right)

    for node_name in nodes:
        node = lookup[node_name]
        assign_cluster(node, cluster_number)
        cluster_number += 1

    return clusters


def add_outliers(supertree, num_outliers, seed=42):
    for _ in range(num_outliers):
        parent = random_supertree_internal_node(supertree)
        add_outlier(supertree, parent, seed=seed)


def convert_clade_to_phylonode(node):
    """Recursively convert ete3 TreeNode objects to custom PhyloNode objects."""
    return PhyloNode(node)


def set_branch_lengths_to_one(tree,n=1):
    if n is None:
        n = 1
    for node in tree.traverse():
        node.dist = n


def get_random_int(min_val, max_val):
    return random.randint(min_val, max_val)


class PhyloNode:
    def __init__(self, node):
        self.node = node
        self.left = PhyloNode(node.children[0]) if len(node.children) > 0 else None
        self.right = PhyloNode(node.children[1]) if len(node.children) > 1 else None
        self.name = node.name
    @property
    def terminal_nodes(self):
        if not self.left and not self.right:
            return 1  # This is a terminal node
        else:
            return (self.left.terminal_nodes if self.left else 0) + (
                self.right.terminal_nodes if self.right else 0
            )

    def distribute(self, value, list_nodes=None):
        if list_nodes is None:
            list_nodes = []
        self.node.add_feature("confidence", value)

        if value == 1:
            list_nodes.append(self.node.name)
            return

        if not self.left and not self.right:
            return list_nodes

        left_terminals = self.left.terminal_nodes if self.left else 0
        right_terminals = self.right.terminal_nodes if self.right else 0

        if random.choice([True, False]):
            max_left_value = min(left_terminals, value - 1)
            min_left_value = max(1, value - right_terminals)
            left_value = get_random_int(min_left_value, max(1, max_left_value))
            right_value = value - left_value
        else:
            max_right_value = min(right_terminals, value - 1)
            min_right_value = max(1, value - left_terminals)
            right_value = get_random_int(min_right_value, max(1, max_right_value))
            left_value = value - right_value

        if left_value > 1:
            self.left.distribute(left_value, list_nodes)
        else:
            list_nodes.append(self.left.node.name)

        if right_value > 1:
            self.right.distribute(right_value, list_nodes)
        else:
            list_nodes.append(self.right.node.name)

        return list_nodes


def build_lookup_dict(node, lookup):
    """Build a lookup dictionary for node names to PhyloNode objects."""
    lookup[node.node.name] = node
    if node.left:
        build_lookup_dict(node.left, lookup)
    if node.right:
        build_lookup_dict(node.right, lookup)


def wrapper_simulation(config):
    # Load configuration
    ntips = config["ntips"]
    total_terminals = config["total_terminals"]
    wss_bss_ratios = config.get("wss_bss_ratio", None)
    increase_factor = config.get("increase_factor", None)
    num_increases = config.get("num_increases", 1)
    output_dir_template = config["output_dir"]
    base_seed = config.get("seed", 42)
    outlier_steps = config["outlier_steps"]
    sampling_steps = config["sampling_steps"]
    noise_steps = config["noise_steps"]
    num_trees = config.get("num_trees", 1)

    for tree_index in range(1, num_trees + 1):
        seed = base_seed + tree_index
        random.seed(seed)
        np.random.seed(seed)

        random_tree = generate_random_tree(total_terminals)
        root_node = convert_clade_to_phylonode(random_tree.get_tree_root())
        node_counter = 0
        for node in random_tree.iter_descendants():
            if not node.is_leaf() and (not node.name or node.name == 1):
                node.name = f"node_{node_counter}"
                node_counter += 1
        nodes = root_node.distribute(ntips)
        lookup = {}
        build_lookup_dict(root_node, lookup)
        clusters = get_clusters(nodes, lookup)

        leaf_counter = 0
        for node_name, cluster_number in clusters.items():
            node = lookup[node_name]
            for leaf in node.node.get_leaves():
                label = generate_label(leaf_counter)
                leaf.name = f"set_{cluster_number}{label}"
                leaf_counter += 1

        set_branch_lengths_to_one(random_tree, n=0)
        for node in random_tree.traverse():
            if hasattr(node, "in_cluster") and node.in_cluster:
                node.dist = 1

        output_dir = output_dir_template.format(
            ntips=ntips,
            total_terminals=total_terminals,
            wss_bss_ratio="original",
        )
        os.makedirs(output_dir, exist_ok=True)
        tree_output_dir = os.path.join(output_dir, f"tree_{tree_index}")
        os.makedirs(tree_output_dir, exist_ok=True)

        ground_truth_clusters = generate_ground_truth(random_tree)
        ground_truth_labels = flatten_clusters(ground_truth_clusters)
        with open(
            os.path.join(
                tree_output_dir,
                "ground_truth_labels.txt",
            ),
            "w",
        ) as f:
            for label in ground_truth_labels:
                if isinstance(label, (list, tuple)) and len(label) >= 2:
                    f.write(f"{label[0]}\t{label[1]}\n")
                else:
                    print(f"Unexpected label format: {label}")
        random_tree.write(
            outfile=os.path.join(tree_output_dir, "tree_increase_0.nw")
        )
        if increase_factor is not None and wss_bss_ratios is None:
            tree = deepcopy(random_tree)

            for i in range(1, num_increases + 1):
                new_branch_length = increase_factor * i
                for node in tree.traverse():
                    if not hasattr(node, "in_cluster") or not node.in_cluster:
                        node.dist = new_branch_length

                filename_suffix = f"increase_{i}"
                tree.write(
                    outfile=os.path.join(tree_output_dir, f"tree_{filename_suffix}.nw")
                )
        else:
            for wss_bss_ratio in wss_bss_ratios:
                bss = total_branch_length(random_tree)
                wss = wss_bss_ratio * bss

                new_branch_length = wss / subtrees_sum
                for subtree in subtrees:
                    for node in subtree.traverse():
                        node.dist = new_branch_length

                combined_tree = deepcopy(random_tree)
                combined_tree = replace_terminals_with_subtrees(combined_tree, subtrees)
                combined_tree.write(
                    outfile=os.path.join(tree_output_dir, f"tree_wss_{wss_bss_ratio}.nw")
                )

                for outlier_fraction in outlier_steps:
                    for sampling_fraction in sampling_steps:
                        for noise_level in noise_steps:
                            if (
                                outlier_fraction > 0
                                or sampling_fraction < 1
                                or noise_level > 0
                            ):
                                current_tree = deepcopy(combined_tree)
                                if outlier_fraction > 0:
                                    num_outliers = math.floor(
                                        total_terminals * outlier_fraction
                                    )
                                    add_outliers(current_tree, num_outliers)
                                if sampling_fraction < 1:
                                    current_tree = sample_tree(
                                        current_tree, sampling_fraction
                                    )
                                if noise_level > 0:
                                    current_tree = add_noise_to_branch_lengths(
                                        current_tree, noise_level
                                    )

                                filename_suffix = f"outliers_{int(outlier_fraction*100)}_sampling_{int(sampling_fraction*100)}_noise_{int(noise_level*100)}_wss_{wss_bss_ratio}"
                                current_tree.write(
                                    outfile=os.path.join(
                                        tree_output_dir, f"tree_{filename_suffix}.nw")
                                )
                                ground_truth_clusters = generate_ground_truth(current_tree)
                                ground_truth_labels = flatten_clusters(
                                    ground_truth_clusters
                                )
                                with open(
                                    os.path.join(
                                        tree_output_dir,
                                        f"ground_truth_labels_{filename_suffix}.txt",
                                    ),
                                    "w",
                                ) as f:
                                    for label in ground_truth_labels:
                                        f.write(f"{label[0]}\t{label[1]}\n")
        print(f"Processed tree {tree_index} with seed {seed}")


def main():
    with open("config_multi.json", "r") as config_file:
        config = json.load(config_file)
    wrapper_simulation(config)


if __name__ == "__main__":
    main()
