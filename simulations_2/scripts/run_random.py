import random
from sklearn.metrics import adjusted_rand_score, v_measure_score
import numpy as np
import os
import csv
from Bio import Phylo
from phytclust import PhytClust
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_prominences
from scipy.interpolate import UnivariateSpline
import json
from collections import defaultdict


def load_ground_truth_labels(tree_folder):
    ground_truth = {}
    for file_name in os.listdir(tree_folder):
        if file_name.startswith("ground_truth_labels") and file_name.endswith(".txt"):
            # Extract parameters from the filename
            params = file_name.replace("ground_truth_labels", "").replace(".txt", "")
            ground_truth[params] = {}
            with open(os.path.join(tree_folder, file_name), "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    ground_truth[params][row[0]] = row[1]

    # If params is an empty string, return the dictionary inside it
    if "" in ground_truth:
        return ground_truth[""]

    return ground_truth


def assign_names_to_internal_nodes(tree):
    counter = 1
    for clade in tree.find_clades():
        if not clade.is_terminal() and not clade.name:
            clade.name = f"Internal_{counter}"
            counter += 1


def get_random_int(min_val, max_val):
    """Generate a random integer between min_val and max_val (inclusive)."""
    return random.randint(min_val, max_val)


def convert_ground_truth_to_names(ground_truth_dict):
    """Convert ground truth dictionary to use clade names as keys."""
    return {
        clade.name: cluster_value for clade, cluster_value in ground_truth_dict.items()
    }


def calculate_ari(clustering_dict, ground_truth_dict):
    # Ensure both dictionaries have the same keys
    clustering_keys = set(clustering_dict.keys())
    ground_truth_keys = set(ground_truth_dict.keys())

    if clustering_keys != ground_truth_keys:
        assert (
            clustering_keys == ground_truth_keys
        ), "Keys in both dictionaries must match."

    # Extract the cluster labels
    clustering_labels = [clustering_dict[key] for key in clustering_dict.keys()]
    ground_truth_labels = [ground_truth_dict[key] for key in ground_truth_dict.keys()]

    # Calculate the Adjusted Rand Index
    ari = v_measure_score(ground_truth_labels, clustering_labels)
    return ari


def convert_clade_to_phylonode(clade):
    """Recursively convert Biopython Clade objects to custom PhyloNode objects."""
    return PhyloNode(clade)


def build_lookup_dict(node, lookup):
    """Build a lookup dictionary for clade names to PhyloNode objects."""
    lookup[node.clade.name] = node
    if node.left:
        build_lookup_dict(node.left, lookup)
    if node.right:
        build_lookup_dict(node.right, lookup)


def get_clusters(nodes, lookup):
    clusters = {}
    cluster_number = 1

    def assign_cluster(node, cluster_number):
        if not node.left and not node.right:
            # This is a terminal node
            clusters[node.clade.name] = cluster_number
        else:
            # Recursively assign cluster number to terminal nodes
            if node.left:
                assign_cluster(node.left, cluster_number)
            if node.right:
                assign_cluster(node.right, cluster_number)

    for node_name in nodes:
        node = lookup[node_name]
        if not node.left and not node.right:
            # Node is a terminal node, forms its own unique cluster
            clusters[node.clade.name] = cluster_number
        else:
            # Assign cluster number to all terminal nodes under this node
            assign_cluster(node, cluster_number)
        cluster_number += 1

    return clusters


class PhyloNode:
    def __init__(self, clade):
        self.clade = clade
        self.left = PhyloNode(clade.clades[0]) if len(clade.clades) > 0 else None
        self.right = PhyloNode(clade.clades[1]) if len(clade.clades) > 1 else None

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
        self.clade.confidence = value

        if value == 1:
            list_nodes.append(self.clade.name)
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
            list_nodes.append(self.left.clade.name)

        if right_value > 1:
            self.right.distribute(right_value, list_nodes)
        else:
            list_nodes.append(self.right.clade.name)

        return list_nodes


def run_random_clustering(tree, ground_truth_dict, num_runs=100):
    root_node = convert_clade_to_phylonode(tree.root)
    terminal_nodes_count = root_node.terminal_nodes
    print(f"Total terminal nodes: {terminal_nodes_count}")

    ari_values = []

    for run in range(num_runs):
        # Random number of clusters between 2 and the number of terminal nodes
        num_clusters = get_random_int(2, terminal_nodes_count)
        # print(f"Run {run + 1}: Number of clusters = {num_clusters}")

        # Distribute the value across the phylogenetic tree
        nodes = root_node.distribute(num_clusters)
        lookup = {}
        build_lookup_dict(root_node, lookup)

        # Get random clusters
        random_clusters = get_clusters(nodes, lookup)

        # Calculate ARI
        ari = calculate_ari(random_clusters, ground_truth_dict)
        ari_values.append(ari)
        # print(f"ARI for run {run + 1}: {ari}")

    # Calculate average and standard deviation of ARI
    average_ari = np.mean(ari_values)
    std_ari = np.std(ari_values)

    print(f"Average ARI over {num_runs} runs: {average_ari}")
    print(f"Standard deviation of ARI: {std_ari}")

    return average_ari, std_ari


def process_trees(case_folder, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    ari_results = []

    for tree_dir in os.listdir(case_folder):
        tree_subdir_path = os.path.join(case_folder, tree_dir)
        if os.path.isdir(tree_subdir_path):
            tree_output_dir = os.path.join(output_dir, tree_dir)
            os.makedirs(tree_output_dir, exist_ok=True)

            ground_truth = load_ground_truth_labels(tree_subdir_path)

            for tree_file in os.listdir(tree_subdir_path):
                if tree_file.endswith(".nw"):
                    tree_path = os.path.join(tree_subdir_path, tree_file)
                    params = tree_file.replace("tree_", "").replace(".nw", "")
                    tree = Phylo.read(tree_path, "newick")
                    assign_names_to_internal_nodes(tree)
                    average_ari, std_ari = run_random_clustering(
                        tree, ground_truth, num_runs=500
                    )
                    ari_results.append(
                        {"Tree": params, "ARI": average_ari, "Std_ARI": std_ari}
                    )

            # Save ARI results
            ari_results_path = os.path.join(tree_output_dir, "comparison_results.json")
            with open(ari_results_path, "w") as jsonfile:
                json.dump(ari_results, jsonfile, indent=4)


if __name__ == "__main__":
    case_folder = (
        "/home/ganesank/project/phytclust/simulations_2/data_new/N_100_K_5"
    )
    output_dir = "/home/ganesank/project/phytclust/simulations_2/results_fmi/random_all_k/N_100_K_5"
    process_trees(case_folder, output_dir)
    print("Processing completed.")


# def run_random_clustering_until_unique_saturates(
#     tree, ground_truth_dict, max_runs=500, max_no_new_solutions=5
# ):
#     root_node = convert_clade_to_phylonode(tree.root)
#     terminal_nodes_count = root_node.terminal_nodes
#     print(f"Total terminal nodes: {terminal_nodes_count}")

#     ari_results = []

#     for k in range(2, terminal_nodes_count + 1):
#         print(f"\nRunning for k = {k} clusters")

#         unique_clusters = set()
#         ari_values = []
#         no_new_solution_runs = 0

#         for run in range(max_runs):
#             # Distribute the k clusters across the phylogenetic tree
#             nodes = root_node.distribute(k)
#             lookup = {}
#             build_lookup_dict(root_node, lookup)

#             # Get random clusters
#             random_clusters = get_clusters(nodes, lookup)
#             random_clusters_tuple = tuple(
#                 sorted(random_clusters.items())
#             )  # Make hashable for uniqueness check

#             # Check if this clustering is unique
#             if random_clusters_tuple in unique_clusters:
#                 no_new_solution_runs += 1
#             else:
#                 unique_clusters.add(random_clusters_tuple)
#                 no_new_solution_runs = 0  # Reset since we found a new unique solution

#             # Calculate ARI
#             ari = calculate_ari(random_clusters, ground_truth_dict)
#             ari_values.append(ari)

#             # If no new unique solutions found in consecutive runs, stop for this k
#             if no_new_solution_runs >= max_no_new_solutions:
#                 print(
#                     f"No new unique solutions found after {max_no_new_solutions} runs. Stopping for k = {k}"
#                 )
#                 break

#         # Calculate average and standard deviation of ARI for this k
#         average_ari = np.mean(ari_values)
#         std_ari = np.std(ari_values)
#         ari_results.append(
#             {
#                 "k": k,
#                 "unique_solutions": len(unique_clusters),
#                 "Average_ARI": average_ari,
#                 "Std_ARI": std_ari,
#                 "Runs": len(ari_values),
#             }
#         )

#         print(
#             f"Completed k = {k}: Average ARI = {average_ari}, Std ARI = {std_ari}, Unique Solutions = {len(unique_clusters)}"
#         )

#     return ari_results


# def process_trees_varying_k(case_folder, output_dir):
#     os.makedirs(output_dir, exist_ok=True)

#     for tree_dir in os.listdir(case_folder):
#         tree_subdir_path = os.path.join(case_folder, tree_dir)
#         if os.path.isdir(tree_subdir_path):
#             tree_output_dir = os.path.join(output_dir, tree_dir)
#             os.makedirs(tree_output_dir, exist_ok=True)

#             ground_truth = load_ground_truth_labels(tree_subdir_path)

#             for tree_file in os.listdir(tree_subdir_path):
#                 if tree_file.endswith(".nw"):
#                     tree_path = os.path.join(tree_subdir_path, tree_file)
#                     params = tree_file.replace("tree_", "").replace(".nw", "")
#                     tree = Phylo.read(tree_path, "newick")
#                     assign_names_to_internal_nodes(tree)

#                     # Run random clustering with varying k
#                     ari_results = run_random_clustering_until_unique_saturates(
#                         tree, ground_truth, max_runs=100, max_no_new_solutions=5
#                     )

#                     # Save ARI results
#                     ari_results_path = os.path.join(
#                         tree_output_dir, "comparison_results_varying_k.json"
#                     )
#                     with open(ari_results_path, "w") as jsonfile:
#                         json.dump(ari_results, jsonfile, indent=4)


# if __name__ == "__main__":
#     case_folder = (
#         "/home/ganesank/project/phytclust/simulations_2/unbalanced_trees/N_100_K_10"
#     )
#     output_dir = "/home/ganesank/project/phytclust/simulations_2/unbalanced_trees_results/random_all_k/N_100_K_10"
#     process_trees_varying_k(case_folder, output_dir)
#     print("Processing completed.")
