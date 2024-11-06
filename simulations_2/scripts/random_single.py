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

num_clusters = 20


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


# def get_random_int(min_val, max_val):
#     """Generate a random integer between min_val and max_val (inclusive)."""
#     return random.randint(min_val, max_val)


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


def get_random_int(min_val, max_val):
    """Generate a random integer between min_val and max_val (inclusive)."""
    result = random.randint(min_val, max_val)
    # print(f"Random integer chosen between {min_val} and {max_val}: {result}")
    return result


# class PhyloNode:
#     def __init__(self, clade):
#         self.clade = clade
#         self.left = PhyloNode(clade.clades[0]) if len(clade.clades) > 0 else None
#         self.right = PhyloNode(clade.clades[1]) if len(clade.clades) > 1 else None

#     @property
#     def cherries(self):
#         if not self.left and not self.right:
#             return 0  # This is a terminal node, not a cherry
#         elif (
#             self.left
#             and self.left.is_terminal()
#             and self.right
#             and self.right.is_terminal()
#         ):
#             return 1  # This is a cherry
#         else:
#             return (self.left.cherries if self.left else 0) + (
#                 self.right.cherries if self.right else 0
#             )

#     def distribute(self, value, list_nodes=None):
#         if list_nodes is None:
#             list_nodes = []
#         self.clade.confidence = value

#         # Stop recursion if the value is 1
#         if value == 1:
#             list_nodes.append(self.clade.name)
#             return

#         if not self.left and not self.right:
#             return list_nodes

#         if not self.left:
#             return self.right.distribute(value, list_nodes)
#         if not self.right:
#             return self.left.distribute(value, list_nodes)

#         # Get the number of cherries in both subtrees
#         left_cherries = self.left.cherries
#         right_cherries = self.right.cherries

#         # Ensure each subtree gets at least 2 cherries
#         if random.choice([True, False]):
#             # Left subtree first
#             max_left_value = min(left_cherries, value - 1)
#             min_left_value = max(1, value - right_cherries)
#             left_value = get_random_int(min_left_value, max(1, max_left_value))
#             right_value = value - left_value
#         else:
#             # Right subtree first
#             max_right_value = min(right_cherries, value - 1)
#             min_right_value = max(1, value - left_cherries)
#             right_value = get_random_int(min_right_value, max(1, max_right_value))
#             left_value = value - right_value

#         # Recursively distribute values
#         if left_value > 1:
#             self.left.distribute(left_value, list_nodes)
#         else:
#             list_nodes.append(self.left.clade.name)

#         if right_value > 1:
#             self.right.distribute(right_value, list_nodes)
#         else:
#             list_nodes.append(self.right.clade.name)

#         return list_nodes

#     def is_terminal(self):
#         return not self.left and not self.right


def collect_nodes_with_value_1(node, nodes_with_value_1):
    """Collect all nodes with confidence value of 1."""
    if node.clade.confidence == 1:
        print(f"Node {node.clade.name} has value 1")
        nodes_with_value_1.append(node.clade.name)

    # Recursively check left and right children
    if node.left:
        collect_nodes_with_value_1(node.left, nodes_with_value_1)
    if node.right:
        collect_nodes_with_value_1(node.right, nodes_with_value_1)


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


def count_cherries(tree):
    cherries_count = 0

    # Traverse the tree to find cherries
    for clade in tree.find_clades():
        if len(clade.clades) == 2:
            if clade.clades[0].is_terminal() and clade.clades[1].is_terminal():
                cherries_count += 1

    return cherries_count


def run_random_clustering(tree, ground_truth_dict, num_runs=100):
    root_node = convert_clade_to_phylonode(tree.root)
    # terminal_nodes_count = root_node.terminal_nodes
    # print(f"Total terminal nodes: {terminal_nodes_count}")

    ari_values = []

    for run in range(num_runs):
        # Random number of clusters between 2 and the number of terminal nodes
        # num_clusters = 5  # get_random_int(2, terminal_nodes_count)
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


# def run_random_clustering(tree, ground_truth_dict, num_runs=100):
#     root_node = convert_clade_to_phylonode(tree.root)
#     terminal_nodes_count = root_node.terminal_nodes
#     print(f"Total terminal nodes: {terminal_nodes_count}")

#     ari_values = []

#     for run in range(num_runs):
#         # Random number of clusters between 2 and the number of terminal nodes
#         num_clusters = 50  #get_random_int(2, terminal_nodes_count)
#         # print(f"Run {run + 1}: Number of clusters = {num_clusters}")

#         # Distribute the value across the phylogenetic tree
#         nodes = root_node.distribute(num_clusters)
#         lookup = {}
#         build_lookup_dict(root_node, lookup)

#         # Get random clusters
#         random_clusters = get_clusters(nodes, lookup)

#         # Calculate ARI
#         ari = calculate_ari(random_clusters, ground_truth_dict)
#         ari_values.append(ari)
#         # print(f"ARI for run {run + 1}: {ari}")

#     # Calculate average and standard deviation of ARI
#     average_ari = np.mean(ari_values)
#     std_ari = np.std(ari_values)

#     print(f"Average ARI over {num_runs} runs: {average_ari}")
#     print(f"Standard deviation of ARI: {std_ari}")

#     return average_ari, std_ari


def process_trees(case_folder, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for tree_dir in os.listdir(case_folder):
        tree_subdir_path = os.path.join(case_folder, tree_dir)
        if os.path.isdir(tree_subdir_path):
            tree_output_dir = os.path.join(output_dir, tree_dir)
            os.makedirs(tree_output_dir, exist_ok=True)

            ground_truth = load_ground_truth_labels(tree_subdir_path)

            increase_ari_dict = {}
            ari_results = []

            for tree_file in os.listdir(tree_subdir_path):
                if tree_file.endswith(".nw"):
                    tree_path = os.path.join(tree_subdir_path, tree_file)
                    params = tree_file.replace("tree_", "").replace(".nw", "")
                    increase_value = int(
                        params.split("_")[1]
                    )  # Assuming the increase value is the second part
                    tree = Phylo.read(tree_path, "newick")
                    assign_names_to_internal_nodes(tree)
                    average_ari, std_ari = run_random_clustering(
                        tree, ground_truth, num_runs=100
                    )
                    if increase_value not in increase_ari_dict:
                        increase_ari_dict[increase_value] = []
                    increase_ari_dict[increase_value].append((average_ari, std_ari))

            for increase, ari_list in increase_ari_dict.items():
                avg_ari = np.mean([ari for ari, _ in ari_list])
                std_ari = np.mean([std for _, std in ari_list])
                ari_results.append(
                    {"Tree": f"increase_{increase}", "ARI": avg_ari, "Std_ARI": std_ari}
                )

            # Save ARI results
            ari_results_path = os.path.join(tree_output_dir, "ari_results.json")
            with open(ari_results_path, "w") as jsonfile:
                json.dump(ari_results, jsonfile, indent=4)


if __name__ == "__main__":
    case_folder = f"/home/ganesank/project/phytclust/simulations_2/data_new/N_100_K_{num_clusters}"
    output_dir = f"/home/ganesank/project/phytclust/simulations_2/results_fmi/random_single_k/N_100_K_{num_clusters}"
    process_trees(case_folder, output_dir)
    print("Processing completed.")
