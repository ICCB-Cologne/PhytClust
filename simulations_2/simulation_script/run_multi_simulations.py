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


def add_outliers(supertree, num_outliers, seed=42):
    for _ in range(num_outliers):
        parent = random_supertree_internal_node(supertree)
        add_outlier(supertree, parent, seed=seed)


def set_branch_lengths_to_one(tree,n=1):
    if n is None:
        n = 1
    for node in tree.traverse():
        node.dist = n


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
        # new_ntips =
        supertree = generate_random_tree(ntips)
        set_branch_lengths_to_one(supertree, n=0)
        # supertree_sum = total_branch_length(supertree)

        output_dir = output_dir_template.format(
            ntips=ntips,
            total_terminals=total_terminals,
            wss_bss_ratio="original",
        )
        os.makedirs(output_dir, exist_ok=True)
        tree_output_dir = os.path.join(output_dir, f"tree_{tree_index}")
        os.makedirs(tree_output_dir, exist_ok=True)
        variance = random.uniform(0, 1)
        distributed_terminals = distribute_ntips_into_mparts(total_terminals, ntips, 1, variance, method="random",seed=seed)
        print(distributed_terminals)
        subtrees = []
        subtrees_sum = 0
        for index, terminals in enumerate(distributed_terminals, start=1):
            subtree = Tree()
            subtree.populate(terminals, random_branches=False)
            for node_index, leaf in enumerate(subtree.iter_leaves(), start=1):
                leaf.name = f"set_{index}{generate_label(node_index - 1)}"
                leaf.is_cluster = True
            letter_counter = 0
            for node in subtree.traverse():
                if not node.is_leaf():
                    node.name = f"internal_{index}{generate_label(letter_counter)}"
                    letter_counter += 1
                    node.is_cluster = True
                node.dist = 1  # Set branch length to 1
            subtrees.append(subtree)
            subtrees_sum += total_branch_length(subtree)
        random.shuffle(subtrees)

        combined_tree = deepcopy(supertree)
        combined_tree = replace_terminals_with_subtrees(combined_tree, subtrees)

        ground_truth_clusters = generate_ground_truth(combined_tree)
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

        if increase_factor is not None and wss_bss_ratios is None:
            # Set all branch lengths to one
            for subtree in subtrees:
                set_branch_lengths_to_one(subtree, n=1)
            set_branch_lengths_to_one(supertree, n=0)

            # Save the initial tree (increase_0)
            combined_tree = deepcopy(supertree)
            combined_tree = replace_terminals_with_subtrees(combined_tree, subtrees)
            combined_tree.write(
                outfile=os.path.join(tree_output_dir, "tree_increase_0.nw")
            )

            for i in range(1, num_increases + 1):
                # Increase branch lengths of the supertree only
                for node in supertree.traverse():
                    node.dist += increase_factor
                for subtree in subtrees:
                    set_branch_lengths_to_one(subtree, n=1)

                # Combine the supertree with subtrees
                combined_tree = deepcopy(supertree)
                combined_tree = replace_terminals_with_subtrees(combined_tree, subtrees)

                # Save the tree
                filename_suffix = f"increase_{i}"
                combined_tree.write(
                    outfile=os.path.join(tree_output_dir, f"tree_{filename_suffix}.nw")
                )
        else:
            for wss_bss_ratio in wss_bss_ratios:
                bss = supertree_sum
                wss = wss_bss_ratio * bss

                # Calculate the new branch length for each branch in the subtrees
                new_branch_length = wss / subtrees_sum
                for subtree in subtrees:
                    for node in subtree.traverse():
                        node.dist = new_branch_length

                combined_tree = deepcopy(supertree)
                combined_tree = replace_terminals_with_subtrees(combined_tree, subtrees)
                # Save or process the combined tree as needed
                combined_tree.write(
                    outfile=os.path.join(tree_output_dir, f"tree_wss_{wss_bss_ratio}.nw")
                )

                # Apply combinations of transformations
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
                                        tree_output_dir, f"tree_{filename_suffix}.nw"
                                    )
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
    # Load configuration from file
    with open("config_multi.json", "r") as config_file:
        config = json.load(config_file)

    # Run the wrapper simulation with the loaded configuration
    wrapper_simulation(config)


if __name__ == "__main__":
    main()
