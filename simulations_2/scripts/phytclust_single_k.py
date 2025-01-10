import os
import csv
from Bio import Phylo
from phytclust import PhytClust
from scipy.ndimage import gaussian_filter1d
import numpy as np
from scipy.signal import find_peaks, peak_prominences
from scipy.interpolate import UnivariateSpline
import json
from sklearn.metrics import adjusted_rand_score, v_measure_score

num_clusters = 5


def load_ground_truth_labels(tree_folder):
    ground_truth = {}
    ground_truth_file = os.path.join(tree_folder, "ground_truth_labels.txt")
    if os.path.isfile(ground_truth_file):
        with open(ground_truth_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                ground_truth[row[0]] = row[1]
    return ground_truth


def process_trees(case_folder, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each tree subdirectory
    for tree_dir in os.listdir(case_folder):
        tree_subdir_path = os.path.join(case_folder, tree_dir)
        if os.path.isdir(tree_subdir_path):
            tree_output_dir = os.path.join(output_dir, tree_dir)
            os.makedirs(tree_output_dir, exist_ok=True)

            ground_truth = load_ground_truth_labels(tree_subdir_path)

            comparison_results = []
            # Iterate over each tree file in the subdirectory
            for tree_file in os.listdir(tree_subdir_path):
                if tree_file.endswith(".nw"):
                    tree_path = os.path.join(tree_subdir_path, tree_file)
                    params = tree_file.replace("tree_", "").replace(".nw", "")
                    tree = Phylo.read(tree_path, "newick")

                    # Determine the number of clusters in the ground truth
                    ground_truth_labels = list(ground_truth.values())
                    k = num_clusters
                    # k = len(set(ground_truth_labels))  # Number of unique clusters
                    # print(k)

                    # Create the PhytClust object with the determined k
                    clust_obj = PhytClust(
                        tree, should_plot_scores=False, num_peaks=1000, k=k
                    )
                    clusters = clust_obj.clusters
                    results = {}  # Store results for each tree here
                    for clade, label in clusters.items():
                        clade_name = clade.name if clade.name else "Unnamed_Clade"
                        if clade_name not in results:
                            results[clade_name] = {}
                        results[clade_name]["ALG_Label"] = label

                    # Save results to CSV
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

                    # Compare with ground truth and save comparison results
                    ground_truth_labels = []
                    predicted_labels = []
                    for clade_name, labels in results.items():
                        ground_truth_label = ground_truth.get(clade_name, "Unknown")
                        ground_truth_labels.append(ground_truth_label)
                        predicted_labels.append(labels["ALG_Label"])
                        comparison_results.append(
                            {
                                "Tree": params,
                                "Clade": clade_name,
                                "Ground_Truth": ground_truth_label,
                                "Predicted_Labels": labels["ALG_Label"],
                            }
                        )

                    # Calculate ARI
                    ari = v_measure_score(ground_truth_labels, predicted_labels)
                    comparison_results.append(
                        {
                            "Tree": params,
                            "ARI": ari,
                        }
                    )

            # Save comparison results to a JSON file for each tree
            comparison_output_path = os.path.join(
                tree_output_dir, "comparison_results.json"
            )
            with open(comparison_output_path, "w") as jsonfile:
                json.dump(comparison_results, jsonfile, indent=4)


if __name__ == "__main__":
    case_folder = f"/home/ganesank/project/phytclust/simulations_2/data/N_100_K_{num_clusters}_short"
    output_dir = f"/home/ganesank/project/phytclust/simulations_2/results/phytclust_single_k/N_100_K_{num_clusters}_short"
    process_trees(case_folder, output_dir)
    print("Processing completed.")
