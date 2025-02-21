from utils import (
    distribute_ntips_into_mparts,
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
    assign_branch_lengths_based_on_leaves,
    total_branch_length,
    distribute_subtree_sums,
)
import os
import random
import math
import string
import json
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
    supertree = generate_random_tree(ntips, seed=seed)

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
            node.dist = 0
        trees.append(subtree)

    supertree = replace_terminals_with_subtrees(supertree, trees)

    combined_tree_sum = total_branch_length(supertree)
    total_subtree_sum = combined_tree_sum * wss_bss_ratio

    num_tips_list = [terminals for terminals in distributed_terminals]
    subtree_sums = distribute_subtree_sums(
        total_subtree_sum, ntips, num_tips_list, mode="proportional", seed=seed
    )

    for subtree, subtree_sum in zip(trees, subtree_sums):
        subtree = assign_branch_lengths_based_on_leaves(subtree, subtree_sum)

    supertree = replace_terminals_with_subtrees(supertree, trees)

    num_outliers = math.floor(total_terminals * percentage_outliers)
    add_outliers(supertree, num_outliers, seed=seed)

    sampled_tree = sample_tree(supertree, sample_fraction, seed=seed)

    noisy_tree = add_noise_to_branch_lengths(sampled_tree, noise_level, seed=seed)

    ground_truth_clusters = generate_ground_truth(supertree)
    ground_truth_labels = flatten_clusters(ground_truth_clusters)

    return noisy_tree, ground_truth_labels

def add_outliers(supertree, num_outliers, seed=42):
    for _ in range(num_outliers):
        parent = random_supertree_internal_node(supertree)
        add_outlier(supertree, parent, seed=seed)


if __name__ == "__main__":
    with open("config_single.json", "r") as config_file:
        config = json.load(config_file)

    seed = config.get("seed", 42)
    random.seed(seed)
    np.random.seed(seed)

    output_dir = config["output_dir"]
    os.makedirs(output_dir, exist_ok=True)

    wss_bss_ratios = config["wss_bss_ratio"]
    if not isinstance(wss_bss_ratios, list):
        wss_bss_ratios = [wss_bss_ratios]

    for wss_bss_ratio in wss_bss_ratios:
        noisy_tree, ground_truth_labels = run_simulation(
            config["ntips"],
            config["total_terminals"],
            wss_bss_ratio,
            config["noise_level"],
            config["percentage_outliers"],
            config["sample_fraction"],
            seed=seed,
        )

        ratio_str = str(wss_bss_ratio).replace(".", "_")
        noisy_tree.write(outfile=os.path.join(output_dir, f"noisy_tree_{ratio_str}.nw"))
        with open(
            os.path.join(output_dir, f"noisy_tree_labels_{ratio_str}.txt"), "w"
        ) as f:
            for label in ground_truth_labels:
                f.write(f"{label[0]}\t{label[1]}\n")
