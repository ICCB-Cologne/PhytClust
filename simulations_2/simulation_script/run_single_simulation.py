from utils import (
    distribute_ntips_into_mparts,
    convert_ete_to_phylo_newick,
    generate_ground_truth,
    flatten_clusters,
)
import numpy as np
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
import os
import random
import math
import string
import json
from io import StringIO
from ete3 import Tree


def run_simulation(
    ntips,
    total_terminals,
    wss_bss_ratio,
    noise_level=0.1,
    percentage_outliers=0.05,
    sample_fraction=0.1,
    seed=42,
):
    # Generate random supertree
    supertree = generate_random_tree(ntips, seed=seed)

    # Distribute terminals among ntips
    distributed_terminals = distribute_ntips_into_mparts(
        total_terminals, ntips, 2, 0.2, method="random", seed=seed
    )

    trees = []

    for index, terminals in enumerate(distributed_terminals, start=1):
        subtree = Tree()
        subtree.populate(terminals, random_branches=False)
        # Set unique names to each leaf node
        for node_index, leaf in enumerate(subtree.iter_leaves(), start=1):
            leaf.name = f"set_{index}{string.ascii_uppercase[node_index - 1]}"
            leaf.is_cluster = True
        letter_counter = 0
        for node in subtree.traverse():
            if not node.is_leaf():
                node.name = f"internal_{index}{string.ascii_uppercase[letter_counter]}"
                letter_counter += 1
                node.is_cluster = True
            # Set branch lengths to 0 initially
            node.dist = 0
        trees.append(subtree)

    # Attach subtrees to the supertree without assigning branch lengths
    supertree = replace_terminals_with_subtrees(supertree, trees)

    # Calculate the total branch length of the combined tree
    combined_tree_sum = total_branch_length(supertree)
    total_subtree_sum = combined_tree_sum * wss_bss_ratio

    # Assign branch lengths to subtrees based on the total branch length
    num_tips_list = [terminals for terminals in distributed_terminals]
    subtree_sums = distribute_subtree_sums(
        total_subtree_sum, ntips, num_tips_list, mode="proportional", seed=seed
    )

    for subtree, subtree_sum in zip(trees, subtree_sums):
        subtree = assign_branch_lengths_based_on_leaves(subtree, subtree_sum)

    # Replace the terminals with subtrees again to update branch lengths
    supertree = replace_terminals_with_subtrees(supertree, trees)

    # Add outliers
    num_outliers = math.floor(total_terminals * percentage_outliers)
    add_outliers(supertree, num_outliers, seed=seed)

    # Sample the tree
    sampled_tree = sample_tree(supertree, sample_fraction, seed=seed)

    # Add noise to the tree
    noisy_tree = add_noise_to_branch_lengths(sampled_tree, noise_level, seed=seed)

    # Generate ground truth clusters
    ground_truth_clusters = generate_ground_truth(supertree)
    ground_truth_labels = flatten_clusters(ground_truth_clusters)

    return noisy_tree, ground_truth_labels


# def run_simulation(
#     ntips,
#     total_terminals,
#     wss_bss_ratio,
#     noise_level=0.1,
#     percentage_outliers=0.05,
#     sample_fraction=0.1,
#     seed=42,
# ):
#     # Generate random supertree
#     supertree = generate_random_tree(ntips, seed=seed)

#     # Distribute terminals among ntips
#     distributed_terminals = distribute_ntips_into_mparts(
#         total_terminals, ntips, 2, 0.2, method="random", seed=seed
#     )
#     supertree_sum = total_branch_length(supertree)
#     total_subtree_sum = supertree_sum * wss_bss_ratio
#     num_tips_list = [terminals for terminals in distributed_terminals]
#     subtree_sums = distribute_subtree_sums(
#         total_subtree_sum, ntips, num_tips_list, mode="proportional", seed=seed
#     )

#     trees = []

#     for index, (terminals, subtree_sum) in enumerate(
#         zip(distributed_terminals, subtree_sums), start=1
#     ):
#         subtree = Tree()
#         subtree.populate(terminals, random_branches=False)
#         # Set unique names to each leaf node
#         for node_index, leaf in enumerate(subtree.iter_leaves(), start=1):
#             leaf.name = f"set_{index}{string.ascii_uppercase[node_index - 1]}"
#             leaf.is_cluster = True
#         letter_counter = 0
#         for node in subtree.traverse():
#             if not node.is_leaf():
#                 node.name = f"internal_{index}{string.ascii_uppercase[letter_counter]}"
#                 letter_counter += 1
#                 node.is_cluster = True
#         subtree = assign_branch_lengths_based_on_leaves(subtree, subtree_sum)
#         trees.append(subtree)

#     supertree = replace_terminals_with_subtrees(supertree, trees)

#     # Add outliers
#     num_outliers = math.floor(total_terminals * percentage_outliers)
#     add_outliers(supertree, num_outliers, seed=seed)

#     # Sample the tree
#     sampled_tree = sample_tree(supertree, sample_fraction, seed=seed)

#     # Add noise to the tree
#     noisy_tree = add_noise_to_branch_lengths(sampled_tree, noise_level, seed=seed)

#     # Generate ground truth clusters
#     ground_truth_clusters = generate_ground_truth(supertree)
#     ground_truth_labels = flatten_clusters(ground_truth_clusters)

#     return noisy_tree, ground_truth_labels


def add_outliers(supertree, num_outliers, seed=42):
    for _ in range(num_outliers):
        parent = random_supertree_internal_node(supertree)
        add_outlier(supertree, parent, seed=seed)


if __name__ == "__main__":
    # Load configuration
    with open("config_single.json", "r") as config_file:
        config = json.load(config_file)

    # Set the seed for reproducibility
    seed = config.get("seed", 42)
    random.seed(seed)
    np.random.seed(seed)

    # Ensure the output directory exists
    output_dir = config["output_dir"]
    os.makedirs(output_dir, exist_ok=True)

    # Get wss_bss_ratios from config and ensure it's a list
    wss_bss_ratios = config["wss_bss_ratio"]
    if not isinstance(wss_bss_ratios, list):
        wss_bss_ratios = [wss_bss_ratios]

    # Iterate over the list of wss_bss_ratios
    for wss_bss_ratio in wss_bss_ratios:
        # Run the simulation with parameters from the config file
        noisy_tree, ground_truth_labels = run_simulation(
            config["ntips"],
            config["total_terminals"],
            wss_bss_ratio,
            config["noise_level"],
            config["percentage_outliers"],
            config["sample_fraction"],
            seed=seed,
        )

        # Save the noisy tree and ground truth labels
        ratio_str = str(wss_bss_ratio).replace(".", "_")
        noisy_tree.write(outfile=os.path.join(output_dir, f"noisy_tree_{ratio_str}.nw"))
        with open(
            os.path.join(output_dir, f"noisy_tree_labels_{ratio_str}.txt"), "w"
        ) as f:
            for label in ground_truth_labels:
                f.write(f"{label[0]}\t{label[1]}\n")
