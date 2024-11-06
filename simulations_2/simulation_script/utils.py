from ete3 import Tree as EteTree
from Bio import Phylo
from io import StringIO
import re
import random
import numpy as np


# def distribute_ntips_into_mparts(
#     num_tips, num_groups, min_size, variance, method="random", seed=None
# ):
#     """
#     Distribute n (int) into m partitions, adjusted by minimum set size and how varied the distribution is.
#     Supports "equal" and "random" distribution methods, with the ability to control variance.
#     """
#     if seed is not None:
#         np.random.seed(seed)  # Set the seed for reproducibility

#     if num_groups * min_size > num_tips:
#         raise ValueError("Minimum size too large for this distribution.")

#     while True:
#         tips_per_group = np.full(num_groups, min_size)
#         remaining_tips = num_tips - (num_groups * min_size)
#         mean_num_tips = remaining_tips / num_groups

#         # Adjust variance based on the number of tips
#         adjusted_variance = variance * (num_tips / num_groups)

#         random_distribution = np.random.normal(
#             loc=mean_num_tips, scale=adjusted_variance, size=num_groups
#         )
#         random_distribution = np.clip(random_distribution, 0, remaining_tips)

#         if np.sum(random_distribution) == 0:
#             continue  # Retry if the sum of random_distribution is zero

#         # Directly add the integer part of the random distribution to tips_per_group
#         tips_per_group += random_distribution.astype(int)

#         # Calculate the remaining tips after adding the integer parts
#         remaining_tips = num_tips - np.sum(tips_per_group)

#         if remaining_tips == 0:
#             break
#         if remaining_tips > 0 and remaining_tips <= num_groups:
#             tips_per_group[
#                 np.random.choice(num_groups, remaining_tips, replace=False)
#             ] += 1
#             break

#     return tips_per_group


def distribute_ntips_into_mparts(
    num_tips, num_groups, min_size, variance, method="random", seed=None
):
    """
    Partition m into n parts using multinomial distribution, avoiding 0s.
    First, give each part 1, then distribute the remaining m-n using multinomial.
    """
    if seed is not None:
        np.random.seed(seed)
    m = num_tips
    n = num_groups
    # Start by giving 1 to each part to avoid 0s
    partition = np.ones(n, dtype=int)
    remaining = m - n  # Remaining value to partition

    if remaining > 0:
        # Use multinomial to distribute the remaining value
        partition += np.random.multinomial(remaining, np.ones(n) / n)

    return partition


# def distribute_ntips_into_mparts(
#     num_tips, num_groups, min_size, variance, method="random", seed=None
# ):
#     """
#     Distribute n (int) into m partitions, adjusted by minimum set size and how varied the distribution is.
#     Supports "equal", "random", and "extreme" distribution methods, with the ability to control variance.
#     """
#     if seed is not None:
#         np.random.seed(seed)  # Set the seed for reproducibility

#     if num_groups * min_size > num_tips:
#         raise ValueError("Minimum size too large for this distribution.")

#     while True:
#         tips_per_group = np.full(num_groups, min_size)
#         remaining_tips = num_tips - (num_groups * min_size)
#         mean_num_tips = remaining_tips / num_groups

#         # Adjust variance based on the number of tips
#         adjusted_variance = variance * (num_tips / num_groups)

#         if method == "random":
#             scale = adjusted_variance * mean_num_tips
#             random_distribution = np.random.normal(
#                 loc=mean_num_tips, scale=scale, size=num_groups
#             )
#         elif method == "extreme":
#             scale_factor = (1 - variance) * mean_num_tips  # Smaller scale for more extremeness
#             random_distribution = np.random.exponential(
#                 scale=scale_factor, size=num_groups
#             )
#             # random_distribution = np.random.normal(
#             # loc=mean_num_tips, scale=adjusted_variance*100, size=num_groups
#             # )
#         elif method == "equal":
#             random_distribution = np.full(num_groups, mean_num_tips)
#         else:
#             raise ValueError("Unsupported method. Use 'random', 'equal', or 'extreme'.")

#         # Ensure no negative values
#         random_distribution = np.maximum(random_distribution, 0)

#         total_random = np.sum(random_distribution)

#         # If the sum of the distribution is 0, we'll retry
#         if total_random == 0:
#             continue  # Retry if the sum of random_distribution is zero

#         # Scale the distribution to sum up to remaining_tips without removing the randomness
#         scaled_distribution = random_distribution * (remaining_tips / total_random)

#         # Convert to integers and adjust remaining tips
#         integer_distribution = np.floor(scaled_distribution).astype(int)
#         tips_per_group += integer_distribution
#         remaining_tips -= np.sum(integer_distribution)

#         # Distribute any remaining tips randomly to account for rounding down
#         if remaining_tips == 0:
#             break
#         if remaining_tips > 0 and remaining_tips <= num_groups:
#             tips_per_group[
#                 np.random.choice(num_groups, remaining_tips, replace=False)
#             ] += 1
#             break

#     return tips_per_group


# def distribute_ntips_into_mparts(n, m, min_size, variance_scale, method="random", seed=None):
#     """
#     Distribute n (int) into m partitions, adjusted by minimum set size and how varied the distribution is.
#     Supports "equal" and "random" distribution methods.
#     """
#     if seed is not None:
#         random.seed(seed)  # Set the seed for reproducibility

#     if m * min_size > n:
#         raise ValueError("Minimum size too large for this distribution.")

#     remaining = n - m * min_size
#     groups = [min_size] * m  # Start with minimum size in each group

#     if method == "random":
#         # Randomly distribute the remaining elements
#         for _ in range(remaining):
#             groups[random.randint(0, m - 1)] += 1
#     else:
#         extras = [0] * m
#         increment_factor = variance_scale * (
#             remaining / m
#         )  # Scale how much extra each gets

#         # Determine distribution based on variance_scale
#         for i in range(m):
#             if variance_scale < 1:
#                 # Distribute incrementally more to earlier groups as an example
#                 extra = increment_factor * (1 - variance_scale) * (m - i) / m
#                 extras[i] = int(extra)
#             else:
#                 # Allocate all remaining to the first group for maximum variance
#                 extras[0] = remaining
#                 break

#         remaining -= sum(extras)

#         # Evenly distribute any small remainder
#         i = 0
#         while remaining > 0:
#             extras[i % m] += 1
#             i += 1
#             remaining -= 1

#         # Combine minimums with extras
#         groups = [groups[i] + extras[i] for i in range(m)]

#     return groups


def convert_ete_to_phylo_newick(ete_tree):
    """
    Convert an ETE3 Tree to a Bio.Phylo Tree using Newick format as intermediary.

    Parameters:
    ete_tree (ete3.Tree): The ETE3 tree to convert.

    Returns:
    Bio.Phylo.BaseTree.Tree: The converted Bio.Phylo tree.
    """
    # Export the ETE tree to a Newick string
    newick_str = ete_tree.write()

    # Use StringIO to simulate a file-like object containing the Newick string
    handle = StringIO(newick_str)

    # Read the Newick string into a Bio.Phylo tree
    phylo_tree = Phylo.read(handle, "newick")

    return phylo_tree


def generate_ground_truth(tree):
    clusters = {}
    outlier_counter = 0
    max_cluster_id = 0

    # First pass to determine the maximum cluster ID
    for leaf in tree.get_leaves():
        name = leaf.name
        if not name.startswith("Outlier_"):
            cluster_id = int(re.findall(r"\d+", name.split("_")[1])[0])
            if cluster_id > max_cluster_id:
                max_cluster_id = cluster_id

    # Second pass to assign clusters
    for leaf in tree.get_leaves():
        name = leaf.name
        if name.startswith("Outlier_"):
            outlier_cluster_id = max_cluster_id + 1 + outlier_counter
            clusters[outlier_cluster_id] = [name]
            outlier_counter += 1
        else:
            cluster_id = re.findall(r"\d+", name.split("_")[1])[0]
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(name)

    return clusters


def flatten_clusters(clusters):
    labels = []
    for cluster_id, members in clusters.items():
        labels.extend([(member, cluster_id) for member in members])
    return labels
