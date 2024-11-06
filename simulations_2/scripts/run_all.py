import subprocess

# input and output paths
input_base_dir = "/home/ganesank/project/phytclust/simulations_2/data/N_100_K_80_short"
result_base_dir = "/home/ganesank/project/phytclust/simulations_2/results"
phylopart_jar = "/home/ganesank/project/phytclust/simulations_2/software/PhyloPart_v2.1/PhyloPart_v2.1.jar"
sub_dir_autophy = "autophy_results/N_100_K_80_short"
sub_dir_phylopart = "phylopart_results/N_100_K_80_short"

# AutoPhy

autophy_command = [
    "python",
    "/home/ganesank/project/phytclust/simulations_2/scripts/run_autophy.py",
    "--input",
    input_base_dir,
    "--output",
    f"{result_base_dir}/{sub_dir_autophy}",
]

subprocess.run(autophy_command)

# Phylopart

phylopart_command = [
    "python",
    "/home/ganesank/project/phytclust/simulations_2/scripts/run_phylopart.py",
    "--input",
    input_base_dir,
    "--output",
    f"{result_base_dir}/{sub_dir_phylopart}",
    "--phylopart_jar",
    phylopart_jar,
]
subprocess.run(phylopart_command)
