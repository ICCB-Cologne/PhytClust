import os
import csv
from Bio import Phylo
from phytclust import PhytClust
from scipy.ndimage import gaussian_filter1d
import numpy as np
from scipy.signal import find_peaks, peak_prominences
from scipy.interpolate import UnivariateSpline
import json
from sklearn.metrics import (
    v_measure_score,
    fowlkes_mallows_score,
    adjusted_mutual_info_score,
)
from collections import defaultdict
from math import log


# Colless index
def colless_index_calc(tree):
    """
    Calculate the Colless index for a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        int: The Colless index.
    """
    internal_nodes = [node for node in tree.find_clades(terminal=False)]
    colless_sum = 0
    for node in internal_nodes:
        left_size = len(node.clades[0].get_terminals())
        right_size = len(node.clades[1].get_terminals())
        colless_sum += abs(left_size - right_size)
    return colless_sum


def normalized_colless(tree):
    """
    Calculate the normalized Colless index for a phylogenetic tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: The normalized Colless index.
    """
    colless_sum = colless_index_calc(tree)
    n = tree.count_terminals()
    if n <= 2:
        return 0
    normalized_colless = (2 * colless_sum) / ((n - 1) * (n - 2))
    return normalized_colless


def variation_of_information(true_labels, predicted_labels):
    # Convert label lists to clusters
    true_clusters = [
        set(np.where(np.array(true_labels) == label)[0]) for label in set(true_labels)
    ]
    predicted_clusters = [
        set(np.where(np.array(predicted_labels) == label)[0])
        for label in set(predicted_labels)
    ]

    # Now calculate VI based on clusters
    n = float(sum([len(x) for x in true_clusters]))  # Total number of elements
    sigma = 0.0
    for x in true_clusters:
        p = len(x) / n  # Probability of cluster x
        for y in predicted_clusters:
            q = len(y) / n  # Probability of cluster y
            r = len(set(x) & set(y)) / n  # Joint probability of x and y overlap
            if r > 0.0:
                sigma += r * (log(r / p, 2) + log(r / q, 2))
    return abs(sigma)


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
    return ground_truth


def assign_names_to_internal_nodes(tree):
    counter = 1
    for clade in tree.find_clades():
        if not clade.is_terminal() and not clade.name:
            clade.name = f"Internal_{counter}"
            counter += 1


def calculate_total_branch_length_to_root(tree):
    total_length = 0
    for terminal in tree.get_terminals():
        path = tree.get_path(terminal)
        total_length += sum(
            clade.branch_length for clade in path if clade.branch_length
        )
    return total_length


def calculate_within_cluster_branch_length(tree, cluster):
    mrca = tree.common_ancestor(cluster)
    total_length = 0
    for terminal in cluster:
        distance = tree.distance(mrca, terminal)
        total_length += distance
    return total_length


def process_trees(case_folder, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    vmeasure_results = []

    for tree_dir in os.listdir(case_folder):
        tree_subdir_path = os.path.join(case_folder, tree_dir)
        if os.path.isdir(tree_subdir_path):
            tree_output_dir = os.path.join(output_dir, tree_dir)
            os.makedirs(tree_output_dir, exist_ok=True)

            ground_truth = load_ground_truth_labels(tree_subdir_path)

            vmeasure_results = []

            for tree_file in os.listdir(tree_subdir_path):
                if tree_file.endswith(".nw"):
                    tree_path = os.path.join(tree_subdir_path, tree_file)
                    params = tree_file.replace("tree_", "").replace(".nw", "")
                    tree = Phylo.read(tree_path, "newick")
                    clust_obj = PhytClust(
                        tree, should_plot_scores=False, num_peaks=1000
                    )
                    actual_scores = clust_obj.scores

                    # peaks, _ = find_peaks(actual_scores)
                    if len(clust_obj.best_peaks) == 0:
                        clusters = clust_obj.clusters[
                            clust_obj.tree.count_terminals() - 1
                        ]
                    else:
                        chosen_peak = clust_obj.best_peaks[0] - 1
                        clusters = clust_obj.clusters[chosen_peak]

                    results = {}
                    for clade, label in clusters.items():
                        clade_name = clade.name if clade.name else "Unnamed_Clade"
                        if clade_name not in results:
                            results[clade_name] = {}
                        results[clade_name]["ALG_Label"] = label

                    output_csv_path = os.path.join(
                        tree_output_dir, f"{params}_clustering_results.csv"
                    )
                    with open(output_csv_path, "w", newline="") as csvfile:
                        fieldnames = ["ID", "ALG_Label"]
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        for clade_name, labels in results.items():
                            row = {"ID": clade_name}
                            row.update(labels)
                            writer.writerow(row)

                    if params in ground_truth:
                        ground_truth_labels = []
                        predicted_labels = []
                        for clade_name, labels in results.items():
                            ground_truth_label = ground_truth[params].get(
                                clade_name, "Unknown"
                            )
                            ground_truth_labels.append(ground_truth_label)
                            predicted_labels.append(labels["ALG_Label"])
                            vmeasure_results.append(
                                {
                                    "Tree": params,
                                    "Clade": clade_name,
                                    "Ground_Truth": ground_truth_label,
                                    "Predicted_Labels": labels["ALG_Label"],
                                }
                            )

                        vmeasure = v_measure_score(
                            ground_truth_labels, predicted_labels
                        )

                        # Calculate the ratio
                        ground_truth_clusters = defaultdict(list)
                        for clade_name, cluster_id in ground_truth[params].items():
                            ground_truth_clusters[int(cluster_id)].append(clade_name)

                        total_branch_length = calculate_total_branch_length_to_root(
                            tree
                        )
                        assign_names_to_internal_nodes(tree)
                        # Calculate the sum of within-cluster branch lengths
                        within_cluster_branch_lengths = sum(
                            calculate_within_cluster_branch_length(tree, cluster)
                            for cluster in ground_truth_clusters.values()
                        )
                        ratio = (
                            total_branch_length - within_cluster_branch_lengths
                        ) / within_cluster_branch_lengths

                        vmeasure_results.append(
                            {
                                "Tree": params,
                                "vmeasure": vmeasure,
                                "Ratio": ratio,
                            }
                        )
                    else:
                        nested_ground_truth = ground_truth.get("", {})
                        ground_truth_labels = []
                        predicted_labels = [
                            labels["ALG_Label"] for labels in results.values()
                        ]
                        for clade_name, predicted_label in results.items():
                            ground_truth_label = nested_ground_truth.get(
                                clade_name, "Unknown"
                            )
                            ground_truth_labels.append(ground_truth_label)

                        vmeasure = v_measure_score(
                            ground_truth_labels, predicted_labels
                        )
                        ground_truth_clusters = defaultdict(list)
                        for clade_name, cluster_id in nested_ground_truth.items():
                            ground_truth_clusters[int(cluster_id)].append(clade_name)

                        total_branch_length = calculate_total_branch_length_to_root(
                            tree
                        )
                        assign_names_to_internal_nodes(tree)
                        # Calculate the sum of within-cluster branch lengths
                        within_cluster_branch_lengths = sum(
                            calculate_within_cluster_branch_length(tree, cluster)
                            for cluster in ground_truth_clusters.values()
                        )
                        ratio = (
                            total_branch_length - within_cluster_branch_lengths
                        ) / within_cluster_branch_lengths
                        colless = normalized_colless(tree)

                        # vmeasure_results.append(
                        #     {
                        #         "Tree": params,
                        #         "VMeasure": vmeasure,
                        #         "Ratio":ratio
                        #     }
                        # )
                    vmeasure_results.append(
                        {
                            "Tree": params,
                            "VMeasure": vmeasure,
                            "Ratio": ratio,
                            "Colless": colless,
                        }
                    )

            vmeasure_output_path = os.path.join(
                tree_output_dir, "vmeasure_results.json"
            )
            with open(vmeasure_output_path, "w") as jsonfile:
                json.dump(vmeasure_results, jsonfile, indent=4)


if __name__ == "__main__":
    case_folder = "/home/ganesank/project/phytclust/simulations_2/data/N_1000_K_5"
    output_dir = "/home/ganesank/project/phytclust/simulations_2/results/phytclust_results_all_k/N_1000_K_5"
    process_trees(case_folder, output_dir)
    print("Processing completed.")
