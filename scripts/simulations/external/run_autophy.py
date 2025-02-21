import os
import json
import glob
import subprocess
from pathlib import Path
import pandas as pd
from sklearn.metrics import v_measure_score
from Bio import Phylo
import argparse


# Read ground truth labels
def read_ground_truth(file_path):
    ground_truth = {}
    with open(file_path, "r") as f:
        for line in f:
            taxa, cluster = line.strip().split()
            ground_truth[taxa] = int(cluster)
    return ground_truth


# Read cluster results and assign unique integer labels to decimal clusters
def read_cluster(file_path):
    df = pd.read_csv(file_path)
    clusters = {}
    cluster_mapping = {}
    next_cluster_id = 0

    for _, row in df.iterrows():
        taxa = row["Label"]
        cluster = float(row["Cluster"])
        if cluster not in cluster_mapping:
            cluster_mapping[cluster] = next_cluster_id
            next_cluster_id += 1
        clusters[taxa] = cluster_mapping[cluster]

    return clusters


# Calculate ARI
def calculate_vmeasure(ground_truth, clusters):
    all_taxa = set(ground_truth.keys()).union(set(clusters.keys()))
    ground_truth_labels = []
    cluster_labels = []
    for taxa in all_taxa:
        ground_truth_labels.append(
            ground_truth.get(taxa, max(ground_truth.values()) + 1)
        )
        cluster_labels.append(clusters.get(taxa, max(clusters.values()) + 1))
    return v_measure_score(ground_truth_labels, cluster_labels)


# Process Newick files, convert to Nexus, and run AutoPhy
def process_newick_files(input_base_dir, result_base_dir):
    print(f"Starting traversal from base directory: {input_base_dir}")

    # Loop through each tree directory
    for tree_dir in input_base_dir.iterdir():
        if tree_dir.is_dir() and tree_dir.name.startswith("tree_"):
            print(f"Checking tree directory: {tree_dir}")
            # Loop through each .nwk file within the tree directory
            for newick_file in tree_dir.iterdir():
                if newick_file.is_file() and newick_file.suffix == ".nw":
                    print(f"Found file: {newick_file}")
                    newick_filename = (
                        newick_file.stem
                    )  # Get the base filename without extension
                    tree_name = tree_dir.name

                    # Create output directory structure with tree and iteration subfolders
                    supertree_output_subdir = (
                        result_base_dir / tree_name / f"{newick_filename}_output"
                    )
                    os.makedirs(supertree_output_subdir, exist_ok=True)

                    nexus_file_path = (
                        supertree_output_subdir / f"{newick_filename}.nexus"
                    )

                    # Check if the Nexus file already exists to avoid reprocessing
                    if nexus_file_path.exists():
                        print(f"Nexus file already exists, skipping: {nexus_file_path}")
                        continue

                    # Convert Newick to Nexus
                    print(f"Processing Newick file: {newick_file}")
                    try:
                        tree = Phylo.read(newick_file, "newick")
                        Phylo.write(tree, nexus_file_path, "nexus")
                        print(f"Converted {newick_file} to {nexus_file_path}")
                    except Exception as e:
                        print(f"Error converting {newick_file} to Nexus format: {e}")
                        continue

                    # Construct the command
                    command = [
                        "/home/ganesank/ENTER/envs/autophy/bin/autophy",
                        "-t",
                        str(nexus_file_path),
                    ]

                    # Change the current working directory to the supertree output subdirectory
                    current_dir = os.getcwd()
                    try:
                        os.chdir(supertree_output_subdir)
                        # Execute the command
                        print(
                            f"Running command in directory: {supertree_output_subdir}"
                        )
                        print(f"Full command: {' '.join(command)}")
                        result = subprocess.run(command, capture_output=True, text=True)
                        result.check_returncode()  # Check for errors
                        print(f"Command output: {result.stdout}")
                    except subprocess.CalledProcessError as e:
                        print(
                            f"Error processing {newick_file} with AutoPhy: {e.stderr}"
                        )
                    except Exception as e:
                        print(f"Unexpected error: {e}")
                    finally:
                        os.chdir(current_dir)  # Change back to the original directory

                    print(
                        f"Processed {newick_file} and saved results in {supertree_output_subdir}"
                    )

    print("All .nw files processed.")


def calculate_vmeasure_results(input_base_dir, result_base_dir):
    for tree_dir in os.listdir(input_base_dir):
        tree_path = os.path.join(input_base_dir, tree_dir)
        if not os.path.isdir(tree_path):
            continue

        # Assuming there's only one ground truth file per tree folder
        ground_truth_path = os.path.join(tree_path, "ground_truth_labels.txt")
        if not os.path.isfile(ground_truth_path):
            print(f"No ground truth file found for {tree_dir}")
            continue

        print(f"Processing ground truth file: {ground_truth_path}")
        ground_truth = read_ground_truth(ground_truth_path)
        print(f"Ground truth data: {ground_truth}")

        results = []

        # Iterate through each subfolder inside the corresponding tree folder in result_base_dir
        result_tree_path = os.path.join(result_base_dir, tree_dir)
        for subfolder in os.listdir(result_tree_path):
            subfolder_path = os.path.join(result_tree_path, subfolder)
            if not os.path.isdir(subfolder_path):
                continue

            result_dir = os.path.join(subfolder_path, "output")
            if not os.path.exists(result_dir):
                print(f"No result directory found for {subfolder}")
                continue

            # Check if the output file already exists
            output_files = glob.glob(os.path.join(result_dir, "*_RF.csv"))
            # if output_files:
            #     print(f"Output file already exists, skipping: {output_files[0]}")
            #     continue

            for file_name in os.listdir(result_dir):
                if file_name.endswith("RF.csv"):
                    file_path = os.path.join(result_dir, file_name)
                    clusters = read_cluster(file_path)
                    vmeasure = calculate_vmeasure(ground_truth, clusters)
                    parent_folder_name = os.path.basename(subfolder_path)
                    results.append({"Tree": parent_folder_name, "VMeasure": vmeasure})

        # Save the results for the current tree directory
        output_json_file = os.path.join(result_base_dir, tree_dir, "validation_results.json")
        with open(output_json_file, "w") as f:
            json.dump(results, f, indent=4)
        print(f"Saved ARI results to {output_json_file}")


# Main function
def main():
    parser = argparse.ArgumentParser(
        description="Process Newick files and calculate ARI."
    )
    parser.add_argument("--input", required=True, help="Input base directory")
    parser.add_argument("--output", required=True, help="Output base directory")
    args = parser.parse_args()

    input_base_dir = Path(args.input)
    result_base_dir = Path(args.output)

    process_newick_files(input_base_dir, result_base_dir)
    calculate_vmeasure_results(input_base_dir, result_base_dir)


if __name__ == "__main__":
    main()
