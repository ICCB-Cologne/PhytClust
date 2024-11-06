import os
import shutil

base_dir = "/home/ganesank/project/phytclust/simulations_2/data/N_100_K_80"
target_dir = "/home/ganesank/project/phytclust/simulations_2/data/N_100_K_80_short"

# Define the increase values you want to copy
increase_values = list(range(0, 101, 5)) + [1]  # 0, 20, 40, ..., 200, and 1

# Loop through each tree folder in the specified directory
for root, dirs, files in os.walk(base_dir):
    if os.path.basename(root).startswith("tree_"):
        for file in files:
            filename = "ground_truth_labels.txt"
            if file == filename:
                # Construct the relative path inside the base directory
                relative_path = os.path.relpath(root, base_dir)

                # Construct the source and target paths
                source_file = os.path.join(root, file)
                target_folder = os.path.join(target_dir, relative_path)
                target_file = os.path.join(target_folder, file)

                # Create the target directory if it doesn't exist
                os.makedirs(target_folder, exist_ok=True)

                # Copy the file to the target directory
                shutil.copy2(source_file, target_file)
                print(f"Copied {source_file} to {target_file}")
            for increase in increase_values:

                filename = f"tree_increase_{increase}.nw"
                if file == filename:
                    # Construct the relative path inside the base directory
                    relative_path = os.path.relpath(root, base_dir)

                    # Construct the source and target paths
                    source_file = os.path.join(root, file)
                    target_folder = os.path.join(target_dir, relative_path)
                    target_file = os.path.join(target_folder, file)

                    # Create the target directory if it doesn't exist
                    os.makedirs(target_folder, exist_ok=True)

                    # Copy the file to the target directory
                    shutil.copy2(source_file, target_file)
                    print(f"Copied {source_file} to {target_file}")
