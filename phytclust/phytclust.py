import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from Bio import Phylo
import io
import os
import string
import random
import matplotlib as mpl
from phytclust.plotting import plot_tree
import matplotlib.colors as mcolors

plt.rcParams["axes.prop_cycle"] = plt.cycler(
    "color", plt.cm.tab20.colors
)  # matplotlib style


class PhytClust:
    """
    This class uses a dynamic programming algorithm to cluster nodes on a phylogenetic tree to annotate monophyletic clades

        -   tree: a newick file

        -   k: number of clusters you would like

        -   max_k: if you'd like to find optimal clusters for a range of 1 to max_k clusters, if k is not given, Default = number of total terminal nodes in a tree
    """

    def __init__(self, tree, k=None, max_k=None):
        self.tree = tree  # Phylo.read(io.StringIO(tree), "newick")
        self.validate_tree()  # to check if everything is in its place, no planted root, internal node names assigned
        self.k = k
        if k is None:
            if max_k is None or max_k == "max":
                self.max_k = len(self.tree.get_terminals())
            else:
                self.max_k = max_k
        else:
            self.max_k = None
        self.clusters = {}
        self.score_type = "CH_score"
        self.dpTable = {}
        self.backtrack = {}
        self.scores = []
        self.tree_clust_dpTable()

    def validate_tree(self):
        assert all(
            [len(node.clades) == 2 for node in self.tree.get_nonterminals()]
        ), "All internal nodes should have 2 children"
        node_names = set()
        internal_node_counter = 0
        for node in self.tree.get_nonterminals() + self.tree.get_terminals():
            if not node.name:
                while True:
                    new_name = "internal_node_" + (
                        string.ascii_uppercase[internal_node_counter % 26]
                        + str(internal_node_counter // 26)
                    )
                    internal_node_counter += 1
                    if new_name not in node_names:
                        node.name = new_name
                        break
            elif node.name in node_names:
                suffix = 1
                new_name = f"{node.name}_{suffix}"
                while new_name in node_names:
                    suffix += 1
                    new_name = f"{node.name}_{suffix}"
                print(
                    f"Node name '{node.name} is duplicated. Adding suffix to make it '{new_name}'"
                )
                node.name = new_name

            node_names.add(node.name)

    def tree_clust_dpTable(self):  # dpTable algorithm
        n = self.k if self.k else self.max_k

        for clade in self.tree.find_clades(order="postorder"):
            scores = np.full(n, np.inf)
            back = np.full((2, n), np.inf)

            if clade.is_terminal():
                scores[0] = 0
            else:
                left, right = clade.clades
                left_size = left.count_terminals()
                right_size = right.count_terminals()

                scores[0] = (
                    self.dpTable[left][0]
                    + self.dpTable[right][0]
                    + left.branch_length * left_size
                    + right.branch_length * right_size
                )
                back[:, 0] = (0, 0)  ## (left, right)

                for i in range(1, n):
                    sub_scores = [
                        self.dpTable[left][j] + self.dpTable[right][i - j - 1]
                        for j in range(i)
                    ]
                    min_score = np.min(sub_scores)
                    scores[i] = min_score

                    if not np.isinf(min_score):
                        min_index = np.argmin(sub_scores)
                        back[:, i] = (min_index, i - min_index - 1)

            self.dpTable[clade], self.backtrack[clade] = scores, back

        # return self.dpTable, self.backtrack

    def tree_clust_backtrack(self, k=None, verbose=False):
        if k is None:
            raise ValueError("value of k is missing.")

        ptr = k - 1
        stack = [(self.tree.clade, 0, ptr)]
        clusters = {}
        current_cluster = 0

        while stack:
            node, counter, ptr = stack.pop()  # Get the last node from the stack
            ptr = int(ptr)
            if verbose:
                print("Visiting node %s (count %d):" % (node.name, counter))

            if ptr == 0:
                clusters.update({x: current_cluster for x in node.get_terminals()})
                current_cluster += 1
            else:
                if len(node.clades) > counter:
                    stack.append((node, counter + 1, ptr))

                    next_ptr = self.backtrack[node][counter, ptr]
                    if not 0 <= next_ptr < k:
                        raise ValueError(f"Pointer out of bounds: {next_ptr}")

                    stack.append((node.clades[counter], 0, int(next_ptr)))
        if current_cluster != k:
            raise ValueError(f"Expected {k} clusters, but found {current_cluster}.")

        return clusters

    def run(self, verbose=False):
        # print(f"self.k = {self.k}, self.max_k = {self.max_k}")  # Debug print

        if self.max_k is None and self.k:
            self.clusters = self.tree_clust_backtrack(self.k, verbose=verbose)
        elif self.k is None:
            if self.max_k is None or self.max_k == "max":
                self.max_k = self.tree.get_terminals()
                self.clusters = [
                    self.tree_clust_backtrack(i, verbose=verbose)
                    for i in range(1, self.max_k + 1)
                ]
            elif isinstance(self.max_k, int) and self.max_k > 0:
                self.clusters = [
                    self.tree_clust_backtrack(i, verbose=verbose)
                    for i in range(1, self.max_k + 1)
                ]
            else:
                raise ValueError("invalid cluster number")
        else:
            raise ValueError("invalid cluster number. Please choose either k or max_k")

        # return self.clusters

    def cluster_score(self, clusters=None, score_type=None, output_all=False):
        clusters = clusters if clusters is not None else self.clusters
        score_type = score_type if score_type is not None else self.score_type

        if not clusters:
            raise ValueError("Clusters not found")
        if score_type not in ["CH_score", "Alt"]:
            raise ValueError(f"Invalid score type: {score_type}")

        # Inverting the clusters dictionary
        clust_inv = {}
        for k, v in clusters.items():
            clust_inv.setdefault(v, []).append(k)

        num_clusters = len(clust_inv)
        num_terminals = self.tree.count_terminals()

        # Check if the number of clusters is 1 or equal to the number of terminals
        if num_clusters == 1 or num_clusters == num_terminals:
            return float("inf")

        # Calculating the left hand side (lhs) and right hand side (rhs) of the score equation
        lhs_list = [
            len(taxa)
            * self.tree.distance(self.tree.common_ancestor(taxa), self.tree.root)
            for taxa in clust_inv.values()
        ]
        rhs_list = [
            np.sum(
                [self.tree.distance(self.tree.common_ancestor(taxa), t) for t in taxa]
            )
            for taxa in clust_inv.values()
        ]

        # Normalizing factors
        lhs_normal = max(1, len(clust_inv) - 1)
        rhs_normal = max(1, num_terminals - len(clust_inv))

        # Computing the score
        lhs = np.sum(lhs_list) / lhs_normal
        rhs = np.sum(rhs_list) / rhs_normal
        score = 0  # Default score value

        if score_type == "CH_score":
            score = lhs / rhs if rhs != 0 else float("nan")

        elif score_type == "Alt":
            tree_depth = (
                max(self.tree.depths().values()) if self.tree.depths().values() else 1
            )
            lambda_clusters = 1 / tree_depth
            lambda_size = 1 / num_terminals

            variance_penalty = np.mean([len(clust) for clust in clust_inv.values()])
            regularization = (
                lambda_clusters * len(clust_inv) + lambda_size * variance_penalty
            )

            # Calculating the final score
            score = lhs / rhs - regularization if rhs != 0 else 0

        # Assuming self.scores is a list attribute where scores are stored
        return (score, lhs, rhs) if output_all else score

    def score_calc(self, plot=True):
        self.scores = []
        if self.k is None:
            for cluster in self.clusters:
                score = self.cluster_score(clusters=cluster)
                self.scores.append(score)
        else:
            score = self.cluster_score(clusters=self.clusters)
            print(f"Single Cluster: {self.clusters}, Score: {score}")  # For debugging
            self.scores.append(score)

        self.scores = np.array(self.scores)

        if plot and self.k is None:
            plt.plot(range(1, len(self.scores) + 1), self.scores, "o-")
            plt.xlabel("No. of Clusters")
            plt.ylabel("Scores")
            plt.title("Scores for Different No. of Clusters")
            plt.show()

        # return self.scores

    def find_peaks(self, n=3, plot=False):
        if self.scores is None or len(self.scores) == 0:
            print("Please calculate scores first")
            return []
        if self.max_k is None:
            print("max_k is not given")
        peaks, _ = find_peaks(self.scores, distance=1)
        if peaks.size < 1:
            print("No peaks found")
            return []

        sorted_peaks = peaks[np.argsort(-self.scores[peaks])]
        top_peaks = sorted_peaks[: min(n, len(sorted_peaks))]

        if plot:
            plt.plot(self.scores)
            plt.plot(
                top_peaks, self.scores[top_peaks], "x", markersize=10, label="Top Peaks"
            )
            plt.legend()
            plt.show()

        if len(top_peaks) < n:
            print(f"Found only {len(top_peaks)} peak(s)")

        self.top_peaks = list(top_peaks)

        return [peak + 1 for peak in top_peaks]

    def plot(
        self,
        results_dir=None,
        top_n=1,
        n=None,
        outlier=True,
        save=False,
        filename=None,
        hide_internal_nodes=True,
        width_scale=2,
        height_scale=0.5,
        label_func=lambda x: None,
        show_branch_lengths=False,
        marker_size=50,
        **kwargs,
    ):
        if self.scores is None:
            print("Please calculate scores first")
            return
        clusters_to_plot = []

        if n is not None:
            clusters_to_plot.append((n - 1, self.clusters[n - 1]))
        elif hasattr(self, "top_peaks") and self.top_peaks is not None:
            top_peaks = self.top_peaks[:top_n]
            clusters_to_plot.extend([(peak, self.clusters[peak]) for peak in top_peaks])
        else:
            peaks, _ = find_peaks(self.scores, distance=1)
            if peaks.size < 1:
                print("No peaks found")
                return
            top_peaks = peaks[np.argsort(-self.scores[peaks])][:top_n]
            clusters_to_plot.extend([(peak, self.clusters[peak]) for peak in top_peaks])

        cmap = plt.get_cmap("tab20")

        # Plot clusters corresponding to the top peaks
        for i, (peak, cluster) in enumerate(clusters_to_plot):
            cluster = self.clusters[peak]
            cluster_sizes = {}
            cluster_number = peak + 1

            for clade, cluster_id in cluster.items():
                cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1

            clumap = {}
            for clade in cluster:
                cluster_id = cluster[clade]
                if outlier and cluster_sizes[cluster_id] == 1:
                    clumap[clade.name] = "black"
                else:
                    color_index = cluster_id % len(cmap.colors)
                    clumap[clade.name] = cmap.colors[color_index]

            plot_tree(
                self.tree,
                title=f"No. of clusters: {cluster_number}, Score: {self.scores[peak]:.4f}",
                label_colors=clumap,
                hide_internal_nodes=hide_internal_nodes,
                width_scale=width_scale,
                height_scale=height_scale,
                label_func=label_func,
                show_branch_lengths=show_branch_lengths,
                marker_size=marker_size,
                **kwargs,
            )

            if save:
                file_to_save = f"{filename}_{i}.png" if filename else f"tree_{peak}.png"
                full_path = os.path.join(results_dir, file_to_save)
                plt.savefig(full_path)
            plt.show()

    def save(
        self, results_dir, top_n=1, filename="phyclust_results.csv", outliers=True
    ):
        if self.scores is None or len(self.scores) == 0:
            print("Please calculate scores first.")
            return

        # Find peaks in the scores
        peaks, _ = find_peaks(self.scores, distance=1)
        if peaks.size < 1:
            print("No peaks found.")
            return

        # Sort peaks from highest score to lowest and select top_n peaks
        top_peaks = peaks[np.argsort(-self.scores[peaks])][:top_n]

        # Get the cluster configurations corresponding to the top peaks
        selected_clusters = [self.clusters[i] for i in top_peaks]

        # Prepare data for CSV
        data = []
        all_cluster_ids = set()  # To keep track of all non-outlier cluster IDs
        for cluster in selected_clusters:
            cluster_ids = list(cluster.values())
            outlier_clusters = {
                cluster_id
                for cluster_id in cluster_ids
                if cluster_ids.count(cluster_id) == 1
            }
            non_outlier_clusters = set(cluster_ids) - outlier_clusters

            # Create a mapping from old cluster IDs to new consecutive cluster IDs starting from 0
            new_cluster_id_map = {
                old_id: new_id
                for new_id, old_id in enumerate(sorted(non_outlier_clusters))
            }

            for node, cluster_id in cluster.items():
                # Nodes that are alone in a cluster are marked as outliers
                if outliers and cluster_id in outlier_clusters:
                    cluster_id = -1
                else:
                    # Update cluster ID to new consecutive ID
                    cluster_id = new_cluster_id_map[cluster_id]

                all_cluster_ids.add(cluster_id)
                data.append({"Node Name": node, "Cluster Number": cluster_id})

        # Create a DataFrame
        df = pd.DataFrame(data)
        full_path = os.path.join(results_dir, filename)
        # Save the DataFrame as a CSV file
        df.to_csv(full_path, index=False)
        print(f"Cluster data saved to {filename}")
