import os
import json
import glob
import pandas as pd
from sklearn.metrics import adjusted_rand_score, silhouette_score, v_measure_score
from ete3 import Tree
from pathlib import Path
import numpy as np
import subprocess
import argparse
from Bio import Phylo


# Read ground truth labels
def read_ground_truth(file_path):
    ground_truth = {}
    with open(file_path, "r") as f:
        for line in f:
            sample, cluster = line.strip().split("\t")
            ground_truth[sample] = int(cluster)
    return ground_truth


# Read cluster results and assign unique integer labels to decimal clusters
def read_cluster_results(file_path):
    cluster_results = pd.read_csv(file_path)
    return cluster_results


def assign_clusters_2(cluster_results):
    clusters = {}
    for _, row in cluster_results.iterrows():
        cluster = row["clustername"]
        leaf = row["leafname"]
        if cluster == 0:
            clusters[leaf] = leaf  # Assign each leaf to its own cluster
        else:
            clusters[leaf] = cluster
    return clusters


def calculate_vmeasure(ground_truth, clusters):
    common_samples = set(ground_truth.keys()).intersection(clusters.keys())
    ground_truth_labels = [ground_truth[sample] for sample in common_samples]
    cluster_labels = [clusters[sample] for sample in common_samples]
    return v_measure_score(ground_truth_labels, cluster_labels)


# Add bootstrap values to the tree
def add_bootstrap_values(tree_file, output_file):
    tree = Tree(tree_file)
    for node in tree.traverse():
        if not node.is_leaf():
            node.support = 100  # Assign a default bootstrap value
    tree.write(outfile=output_file)
    print(f"Bootstrap values added to {output_file}")


# Validate Newick file
def validate_newick_file(tree_file):
    try:
        tree = Tree(tree_file)
        print(f"{tree_file} is a valid Newick file.")
        return True
    except Exception as e:
        print(f"Error in {tree_file}: {e}")
        return False


# Clean Newick file
def clean_newick_file(tree_file):
    with open(tree_file, "r") as file:
        content = file.read()
    # Remove any multiple decimal points in branch lengths
    cleaned_content = content.replace("..", ".")
    with open(tree_file, "w") as file:
        file.write(cleaned_content)
    print(f"Cleaned {tree_file}")


# Process trees and add bootstrap values
def process_trees(tree_dir, phylopart_input_dir):
    for root, _, files in os.walk(tree_dir):
        for file in files:
            if file.endswith(".nw"):
                tree_file = os.path.join(root, file)
                # Create corresponding path in phylopart_input_dir
                relative_path = os.path.relpath(tree_file, tree_dir)
                output_tree_file = os.path.join(phylopart_input_dir, relative_path)
                os.makedirs(os.path.dirname(output_tree_file), exist_ok=True)
                # Add bootstrap values to the new file
                add_bootstrap_values(tree_file, output_tree_file)
                # Validate and clean the Newick file
                if not validate_newick_file(output_tree_file):
                    clean_newick_file(output_tree_file)
                    validate_newick_file(output_tree_file)


# Function to run PhyloPart with the calculated threshold
def run_phylopart(tree_file, threshold, output_file, phylopart_jar):
    cmd = [
        "java",
        "-Xmx4G",
        "-jar",
        phylopart_jar,
        tree_file,
        str(threshold),
        f"-o{output_file}",
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running PhyloPart for {tree_file} with threshold {threshold}")
        print(e.stdout)
        print(e.stderr)


# Process trees and calculate global distribution
def process_trees_and_calculate_distribution(tree_dir, output_base_dir, phylopart_jar):
    thresholds = np.arange(
        0.05, 0.45, 0.05
    )  # Thresholds from 0.05 to 0.4 with step size 0.05
    for root, _, files in os.walk(tree_dir):
        for file in files:
            if file.endswith(".nw"):
                tree_file = os.path.join(root, file)
                relative_path = os.path.relpath(tree_file, tree_dir)
                tree_output_dir = os.path.join(
                    output_base_dir,
                    os.path.dirname(relative_path),
                    os.path.splitext(file)[0],
                )
                # Ensure the output directory structure is created
                os.makedirs(tree_output_dir, exist_ok=True)
                # Run PhyloPart on the file for each threshold
                for threshold in thresholds:
                    output_file = os.path.join(
                        tree_output_dir,
                        f"{os.path.splitext(file)[0]}_partition_output_{threshold}.txt",
                    )
                    run_phylopart(tree_file, threshold, output_file, phylopart_jar)


# Load tree
def load_tree(tree_file):
    try:
        print(f"Loading tree from file: {tree_file}")
        tree = Phylo.read(tree_file, "newick")
        return tree
    except Exception as e:
        print(f"Error loading tree from file {tree_file}: {e}")
        return None


# Calculate distance matrix
def calculate_distance_matrix(tree):
    leaves = tree.get_terminals()
    num_leaves = len(leaves)
    distance_matrix = np.zeros((num_leaves, num_leaves))
    for i, leaf1 in enumerate(leaves):
        for j, leaf2 in enumerate(leaves):
            distance_matrix[i, j] = tree.distance(leaf1, leaf2)
    return distance_matrix, [leaf.name for leaf in leaves]


# Load clusters
def load_clusters(cluster_file):
    df = pd.read_csv(cluster_file)
    clusters = df.groupby("clustername")["leafname"].apply(list).to_dict()
    return clusters


# Assign clusters
def assign_clusters(leaves, clusters):
    leaf_to_cluster = {
        leaf_name: cluster
        for cluster, leaf_names in clusters.items()
        for leaf_name in leaf_names
    }
    labels = [leaf_to_cluster.get(leaf, -1) for leaf in leaves]
    return labels


# Calculate silhouette score
def calculate_silhouette(distance_matrix, labels):
    return silhouette_score(distance_matrix, labels, metric="precomputed")


# Process files to calculate silhouette scores
def process_files(tree_file, cluster_files, output_file):
    tree = load_tree(tree_file)
    if tree is None:
        print(f"Skipping processing for {tree_file} due to loading error.")
        return
    distance_matrix, leaves = calculate_distance_matrix(tree)
    with open(output_file, "w") as f:
        f.write("Threshold,Silhouette Score\n")
        for cluster_file in cluster_files:
            clusters = load_clusters(cluster_file)
            labels = assign_clusters(leaves, clusters)
            # Check if there are at least two unique labels
            unique_labels = set(labels)
            if len(unique_labels) < 2:
                silhouette_avg = 0
                print(
                    f"Only one unique label in {cluster_file}. Setting silhouette score to 0."
                )
            else:
                silhouette_avg = calculate_silhouette(distance_matrix, labels)
            threshold = os.path.splitext(os.path.basename(cluster_file))[0].split("_")[
                -1
            ]
            f.write(f"{threshold},{silhouette_avg}\n")
            print(f"Silhouette score for {cluster_file}: {silhouette_avg}")


# Process directories to calculate silhouette scores
def process_directories(output_range_dir, data_dir):
    for output_dir in Path(output_range_dir).iterdir():
        if output_dir.is_dir():
            for sub_dir in output_dir.iterdir():
                if sub_dir.is_dir():
                    # Construct the tree file name based on the sub-directory name
                    tree_file = Path(data_dir) / output_dir.name / f"{sub_dir.name}.nw"
                    tree_file = tree_file.resolve()
                    print(tree_file)
                    if tree_file.exists():
                        cluster_files = list(sub_dir.glob("*.txt"))
                        output_file = sub_dir / "silhouette_scores.csv"
                        process_files(tree_file, cluster_files, output_file)
                    else:
                        print(f"No tree file found for {tree_file}")


def calculate_vmeasure_results(data_dir, cluster_results_base_dir):
    for tree_dir in os.listdir(data_dir):
        tree_path = os.path.join(data_dir, tree_dir)
        if not os.path.isdir(tree_path):
            continue
        print(f"Processing tree directory: {tree_dir}")
        ground_truth_files = glob.glob(
            os.path.join(tree_path, "ground_truth_labels*.txt")
        )
        if not ground_truth_files:
            print(f"No ground truth files found for {tree_dir}")
            continue
        cluster_results_dir = os.path.join(cluster_results_base_dir, tree_dir)
        if not os.path.exists(cluster_results_dir):
            print(f"No cluster results found for {tree_dir}")
            continue
        aggregated_results = []
        for ground_truth_path in ground_truth_files:
            print(f"Processing ground truth file: {ground_truth_path}")
            ground_truth = read_ground_truth(ground_truth_path)
            ground_truth_filename = os.path.basename(ground_truth_path)
            # Find the corresponding folder for the ground truth file
            ground_truth_suffix = ground_truth_filename.split("ground_truth_labels")[
                1
            ].replace(".txt", "")
            cluster_result_folder_pattern = os.path.join(
                cluster_results_dir, f"*{ground_truth_suffix}*"
            )
            matching_folders = [
                folder
                for folder in glob.glob(cluster_result_folder_pattern)
                if os.path.isdir(folder)
            ]
            print(f"Matching folders for {ground_truth_suffix}: {matching_folders}")
            if not matching_folders:
                print(f"No matching folder found for {ground_truth_filename}")
                continue
            for folder in matching_folders:
                # Read silhouette scores from the CSV file in the results directory
                silhouette_scores_file = os.path.join(folder, "silhouette_scores.csv")
                if not os.path.exists(silhouette_scores_file):
                    print(f"No silhouette scores file found in {folder}")
                    continue
                silhouette_scores = pd.read_csv(silhouette_scores_file)
                best_solution = silhouette_scores.loc[
                    silhouette_scores["Silhouette Score"].idxmax()
                ]
                best_threshold = round(best_solution["Threshold"], 2)
                best_solution_pattern = os.path.join(
                    folder,
                    f"tree{ground_truth_suffix}*partition_output_{best_threshold}*.txt",
                )
                best_solution_files = glob.glob(best_solution_pattern)
                if not best_solution_files:
                    print(
                        f"No matching solution file found for pattern: {best_solution_pattern}"
                    )
                    continue
                best_solution_file = best_solution_files[0]
                print(f"Processing best solution file: {best_solution_file}")
                cluster_results = read_cluster_results(best_solution_file)
                clusters = assign_clusters_2(cluster_results)
                vmeasure = calculate_vmeasure(ground_truth, clusters)
                tree_name = os.path.basename(best_solution_file).split(
                    "_partition_output"
                )[0]
                aggregated_results.append({"Tree": tree_name, "VMeasure": vmeasure})
        # Save the aggregated results for the current tree directory
        output_file_path = os.path.join(cluster_results_dir, "vmeasure_results.json")
        try:
            with open(output_file_path, "w") as f:
                json.dump(aggregated_results, f, indent=4)
            print(f"Saved aggregated ARI results to {output_file_path}")
        except Exception as e:
            print(f"Failed to save aggregated ARI results to {output_file_path}: {e}")


# Main function
def main():
    parser = argparse.ArgumentParser(
        description="Process Newick files, run PhyloPart, calculate silhouette scores, and calculate ARI."
    )
    parser.add_argument("--input", required=True, help="Input base directory")
    parser.add_argument("--output", required=True, help="Output base directory")
    parser.add_argument(
        "--phylopart_jar", required=True, help="Path to the PhyloPart JAR file"
    )
    args = parser.parse_args()

    input_base_dir = Path(args.input)
    result_base_dir = Path(args.output)
    phylopart_jar = args.phylopart_jar

    # Process trees and add bootstrap values
    process_trees(input_base_dir, result_base_dir)

    # Run PhyloPart and calculate global distribution
    process_trees_and_calculate_distribution(
        result_base_dir, result_base_dir, phylopart_jar
    )

    # Calculate silhouette scores
    process_directories(result_base_dir, input_base_dir)

    # Calculate ARI results
    calculate_vmeasure_results(input_base_dir, result_base_dir)


if __name__ == "__main__":
    main()
