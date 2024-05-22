import time
from collections import defaultdict
from copy import deepcopy
import string
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib.ticker import MaxNLocator
from phytclust.find_peaks import find_plateau_edges, select_representative_edges
from phytclust.plotting import plot_cluster, plot_peaks
from phytclust.save import save_clusters
from phytclust.validation import is_outgroup_valid, validate_tree, rename_nodes
from phytclust.helper import find_all_min_indices

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

    def __init__(self, tree, num_peaks=3, plot_scores=False, k=None, max_k=None, outgroup=None):
        self.logger = logging.getLogger("phytclust") #logging
        self.tree = tree
        self.validate_and_set_outgroup(outgroup)
        validate_tree(self.tree, self.outgroup)
        rename_nodes(self.tree, self.outgroup)
        self.node_terminals = {
            node: node.get_terminals() for node in self.tree.find_clades()
        }
        self.terminal_count = {
            node: len(terminals) for node, terminals in self.node_terminals.items()
        }
        self.num_terminals = (
            (self.terminal_count[self.tree.root] - 1)
            if self.outgroup
            else self.terminal_count[self.tree.root]
        )

        self.k = k
        if k is None:
            self.max_k = self.num_terminals if max_k is None else max_k
        else:
            self.max_k = None

        self.clusters = {}
        self.dp_table = {}
        self.backtrack = {}
        self.scores = []
        self.top_peaks = []
        self.total_runtime = 0.0
        self.runtimes = {}
        self._no_outgroup_tree = deepcopy(tree) if self.outgroup else None
        if self.outgroup:
            self.prune_outgroup()

        self.node_terminals = {node: node.get_terminals() for node in self.tree.find_clades()}
        self.terminal_count = {node: len(terminals) for node, terminals in self.node_terminals.items()}

        self.tree_clust_dp_table()
        self.num = []
        self.den = []

        self.tree_clust_dp_table()
        self.run()
        self.score_calc(plot=False)
        self.find_peaks(n=num_peaks, plot=plot_scores)

    def validate_and_set_outgroup(self, outgroup):
        if outgroup and not is_outgroup_valid(self.tree, outgroup):
            raise ValueError("Outgroup not found, please check input.")
        self.outgroup = outgroup

    def prune_outgroup(self):
        outgroup_clade = next(
            self._no_outgroup_tree.find_clades(self.outgroup), None)
        if len(self._no_outgroup_tree.root.clades) > 2:

            self._no_outgroup_tree.prune(outgroup_clade)
        elif len(self._no_outgroup_tree.root.clades) == 2:
            outgroup_clade = None
            sibling_clade = None
            for clade in self.tree.root.clades:
                if clade.name == self.outgroup:
                    outgroup_clade = clade
                else:
                    sibling_clade = clade
            sibling_clade.name = self.outgroup

            self._no_outgroup_tree.root = sibling_clade

    # def tree_clust_dp_table(self):
    #     start_time = time.time()

    #     n = self.k or self.max_k
    #     active_tree = self._no_outgroup_tree if self.outgroup else self.tree

    #     # Prepare dictionaries to store the terminal counts and branch scores
    #     terminal_count = {}
    #     branch_score = {}
    #     dp_table = {}
    #     backtrack = {}

    #     # Postorder traversal to compute terminal counts and initial DP values
    #     for clade in active_tree.find_clades(order="postorder"):
    #         scores = np.full(n, np.inf)
    #         back = np.full((2, n), np.inf)
    #         if clade.is_terminal():
    #             terminal_count[clade] = 1
    #             scores[0] = 0
    #         else:
    #             left, right = clade.clades
    #             left_size = terminal_count[left]
    #             right_size = terminal_count[right]
    #             total_terminals = left_size + right_size
    #             terminal_count[clade] = total_terminals

    #             left_branch_score = left.branch_length * left_size
    #             right_branch_score = right.branch_length * right_size
    #             branch_score[clade] = left_branch_score + right_branch_score

    #             limit = min(n, total_terminals)

    #             scores[0] = (
    #                 dp_table[left][0]
    #                 + dp_table[right][0]
    #                 + left_branch_score
    #                 + right_branch_score
    #             )
    #             back[:, 0] = (0, 0)

    #             left_scores = dp_table[left][:limit]
    #             right_scores_reversed = dp_table[right][:limit][::-1]

    #             for i in range(1, limit):
    #                 left_sub_scores = left_scores[:i]
    #                 right_sub_scores = right_scores_reversed[:i]

    #                 # Broadcasting to compute sub_scores
    #                 sub_scores = (
    #                     left_sub_scores[:, np.newaxis] + right_sub_scores[np.newaxis, :]
    #                 )

    #                 # Get the index of the minimum value in the flattened sub_scores array
    #                 flat_min_index = np.argmin(sub_scores)
    #                 # Convert the flat index to 2D indices
    #                 left_index, right_index = np.unravel_index(
    #                     flat_min_index, sub_scores.shape
    #                 )
    #                 min_score = sub_scores[left_index, right_index]

    #                 scores[i] = min_score

    #                 if not np.isinf(min_score):
    #                     back[:, i] = (left_index, right_index)

    #         dp_table[clade], backtrack[clade] = scores, back

    #     self.dp_table = dp_table
    #     self.backtrack = backtrack

    #     end_time = time.time()
    #     self._track_runtime("dp_table", start_time, end_time)

    def tree_clust_dp_table(self):
        start_time = time.time()

        n = self.k or self.max_k

        active_tree = self._no_outgroup_tree if self.outgroup else self.tree

        for clade in active_tree.find_clades(order="postorder"):
            scores = np.full(n, np.inf)
            back = np.full((2, n), np.inf)
            if clade.is_terminal():
                scores[0] = 0
            else:
                left, right = clade.clades
                left_size = self.terminal_count[left]
                right_size = self.terminal_count[right]
                total_terminals = left_size + right_size

                limit = min(n, total_terminals)

                scores[0] = (
                    self.dp_table[left][0]
                    + self.dp_table[right][0]
                    + (left.branch_length * (left_size))
                    + (right.branch_length * (right_size))
                )
                back[:, 0] = (0, 0)

                for i in range(1, limit):
                    sub_scores = self.dp_table[left][:i] + self.dp_table[right][:i][::-1]

                    min_index = np.argmin(sub_scores)
                    min_score = sub_scores[min_index]
                    scores[i] = min_score

                    if not np.isinf(min_score):
                        back[:, i] = (min_index, i - min_index - 1)
                    # min_indices, min_score = find_all_min_indices(sub_scores)
                    # scores[i] = min_score

                    # if not np.isinf(min_score):
                    #     if len(min_indices) > 1:
                    #         print(f"Multiple equal solutions detected at {clade.name}")
                    #     min_index = min_indices[0]
                    #     back[:, i] = (min_index, i - min_index - 1)

            self.dp_table[clade], self.backtrack[clade] = scores, back
        end_time = time.time()
        self._track_runtime("dp_table", start_time, end_time)

    def tree_clust_backtrack(self, k=None, verbose=False):
        if k is None:
            raise ValueError("value of k is missing.")

        active_tree = self._no_outgroup_tree if self.outgroup else self.tree
        ptr = k - 1
        stack = [(active_tree.clade, 0, ptr)]
        clusters = {}
        current_cluster = 0

        while stack:
            node, counter, ptr = stack.pop()

            if verbose:
                print("Visiting node %s (count %d):" % (node.name, counter))

            if ptr == 0:
                clusters.update({x: current_cluster for x in self.node_terminals[node]})
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

    def def_clusters(self, verbose):
        self.clusters = [
            self.tree_clust_backtrack(i, verbose=verbose)
            for i in range(1, self.max_k + 1)
        ]

    def run(self, verbose=False):
        # print(f"self.k = {self.k}, self.max_k = {self.max_k}")  # Debug
        start_time = time.time()
        if self.max_k is None and self.k:
            self.clusters = self.tree_clust_backtrack(k=self.k, verbose=verbose)

        elif self.k is None:
            if self.max_k is None:
                self.max_k = self.num_terminals
                self.def_clusters(verbose)

            elif isinstance(self.max_k, int) and self.max_k > 0:
                self.def_clusters(verbose)

            else:
                raise ValueError("invalid cluster number")
        else:
            raise ValueError("invalid cluster number. Please choose either k or max_k")

        end_time = time.time()
        self._track_runtime("run", start_time, end_time)

    # def calculate_coefficient_of_variation(self):
    #     # Extract branch lengths
    #     branch_lengths = [
    #         branch.branch_length
    #         for branch in self.tree.get_terminals()
    #         if branch.branch_length is not None
    #     ]

    #     # Calculate mean and standard deviation
    #     mean_length = np.mean(branch_lengths)
    #     std_deviation = np.std(branch_lengths)

    #     # Calculate coefficient of variation
    #     cv = std_deviation / mean_length

    #     return cv
    def calculate_coefficient_of_variation(self, branch_lengths):
        if len(branch_lengths) < 2:
            return None  # Not enough data to calculate CV
        mean_length = np.mean(branch_lengths)
        std_deviation = np.std(branch_lengths)
        cv = (std_deviation / mean_length) if mean_length > 0 else float('inf')
        return cv

    def cluster_score(self, clusters=None, output_all=False):
        clusters = clusters or self.clusters

        if not clusters:
            raise ValueError("Clusters not found")

        active_tree = self._no_outgroup_tree if self.outgroup else self.tree
        clust_inv = {}
        for k, v in clusters.items():
            clust_inv.setdefault(v, []).append(k)

        num_clusters = len(clust_inv)
        num_terminals = self.num_terminals

        root_clade = active_tree.root
        dp_row = self.dp_table.get(root_clade, None)
        if dp_row is None:
            raise ValueError("Root not found in dp_table")

        beta_1 = dp_row[0]
        beta = dp_row[num_clusters - 1]

        num = (beta_1 - beta) / (beta) if beta else float("inf")
        den = (
            # ((num_terminals - num_clusters)) / (num_clusters - 1)
            ((num_terminals - num_clusters)) / ((num_clusters - 1))
            if (num_clusters - 1)
            else float("nan")
        )
        score = (num * den) if den not in [float("inf")] else float("inf")

        return num, (beta / beta_1), score

    def score_calc(self, plot=False, output_all=False):
        start_time = time.time()
        if self.k is None:
            results = [self.cluster_score(clusters=cluster, output_all=output_all) for cluster in self.clusters]
            num_list, den_list, scores = zip(*results)
        else:
            num, den, score = self.cluster_score(
                clusters=self.clusters, output_all=output_all
            )
            print(f"Single Cluster: {self.clusters}, Score: {score}")
            scores, num_list, den_list = [score], [num], [den]

        self.scores = np.array(scores)
        self.num = np.array(num_list)
        self.den = np.array(den_list)

        if plot and self.k is None:
            self.plot_scores()

        end_time = time.time()
        self._track_runtime("score_calc", start_time, end_time)

    def plot_scores(self, scores_subset, peaks, k_start):
        fig, ax = plt.subplots()
        # Adjust x_indices to reflect the original full list indexing
        x_indices = np.arange(
            k_start + 1, k_start + len(scores_subset) + 1
        )  # k_start + 1 for 1-based label
        ax.plot(x_indices, scores_subset, "o-", label="Scores", markersize=2)

        # Highlight peaks if provided
        if peaks is not None:
            # Adjust peak x-indices to match the 1-based plot index
            peak_x_indices = [k_start + 1 + p for p in peaks]
            peak_scores = [scores_subset[p] for p in peaks]
            ax.plot(peak_x_indices, peak_scores, "rx", label="Peaks")

        plt.xlabel("No. of Clusters")
        plt.ylabel("Scores")
        plt.yscale("log")
        plt.title("Scores for Different No. of Clusters")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.show()
        return fig

    def find_peaks(
        self,
        n=3,
        plot=True,
        k_start=0,
        k_end=None,
    ):
        """
        Finds and plots the top n maxima in the scores within a specified range of k.

        Parameters:
        - n: int, the number of top peaks to find.
        - plot: bool, whether to plot the score graph with the peaks highlighted.
        - method: str, 'default' or 'alt' for different peak finding methods.
        - prominence: tuple or None, the prominence parameter for peak finding, applicable in 'alt' method.
        - k_start: int or None, the starting index for the range of k values to consider.
        - k_end: int or None, the ending index for the range of k values to consider.
        """
        start_time = time.time()
        if not hasattr(self, "scores") or self.scores is None or len(self.scores) == 0:
            print("Please calculate scores first")
            return []

        # Adjust k_start for 0-based indexing within Python
        k_start = max(k_start - 1, 0)
        k_end = k_end or len(self.scores)

        if k_end > len(self.scores) or k_end < k_start:
            raise ValueError(f"Invalid k_end value. Max allowed value is {len(self.scores)}")

        scores_subset = self.scores[k_start:k_end]
        scores_subset = np.where(np.isinf(scores_subset), np.nan, scores_subset)

        #####
        smoothing = max(1, int(self.num_terminals * 0.1))
        edges, plateau_sizes, dF, d2F = find_plateau_edges(
            scores_subset, smoothing)

        # Add edge cases to the edges list
        if scores_subset[0] > scores_subset[1]:
            edges.append(k_start)

        # Find the last non-NaN index
        last_valid_index = len(scores_subset) - np.isnan(scores_subset[::-1]).argmax() - 1

        # Check if the last non-NaN value is greater than the one before it
        if scores_subset[last_valid_index] > scores_subset[last_valid_index - 1]:
            edges.append(k_start + last_valid_index)

        # If there's only one edge, it's the only peak
        if len(edges) == 1:
            self.top_peaks = list(peak + k_start + 1 for peak in edges)
            if plot:
                self.plot_scores(scores_subset, list(edges), k_start)
    
            return self.top_peaks

        edges = list(set(edges))

        representatives, representative_plateau_sizes = select_representative_edges(
            edges, plateau_sizes
        )

        # Sort representatives and their plateau sizes by plateau size
        sorted_indices = np.argsort(representative_plateau_sizes)[::-1]
        representatives = np.array(representatives)[sorted_indices]
        representative_plateau_sizes = np.array(representative_plateau_sizes)[sorted_indices]

        top_peaks = [representatives[: min(n, len(representatives))]]

        if plot:
            self.plot_scores(scores_subset, top_peaks, k_start)

        if len(top_peaks) < n:
            print(f"Found only {len(top_peaks)} peak(s)")

        # Adjust top_peaks to the original indexing context
        self.top_peaks = list(peak + k_start + 1 for peak in top_peaks)
        end_time = time.time()
        self._track_runtime("find_peaks", start_time, end_time)

        # Return the peaks adjusted to 1-based indexing
        return [peak + 1 for peak in top_peaks]

    def plot(
        self,
        results_dir=None,
        top_n=1,
        n=None,
        cmap=plt.get_cmap("tab20"),
        show_terminal_labels=False,
        outlier=False,
        save=False,
        filename=None,
        hide_internal_nodes=True,
        width_scale=2,
        height_scale=0.1,
        label_func=lambda x: None,
        show_branch_lengths=False,
        marker_size=40,
        **kwargs,
    ):
        if not hasattr(self, "scores") or self.scores is None:
            print("Please calculate scores first")
            return

        clusters_to_plot = []
        # Check if k is specified and directly use self.clusters for plotting
        if self.k is not None:
            clusters_to_plot.append((self.k, self.clusters))
        elif n is not None:
            clusters_to_plot.append((n, self.clusters[n - 1]))
        else:
            # Existing logic for when k is not specified
            if hasattr(self, "top_peaks") and self.top_peaks:
                top_peaks = self.top_peaks[:top_n]
                clusters_to_plot.extend(
                    [(peak, self.clusters[peak]) for peak in top_peaks]
                )
            else:
                peaks, _ = find_peaks(self.scores, distance=1)
                if peaks.size < 1:
                    print("No peaks found")
                    return
                top_peaks = peaks[np.argsort(-self.scores[peaks])][:top_n]
                clusters_to_plot.extend(
                    [(peak, self.clusters[peak]) for peak in top_peaks]
                )
            # Plotting for specified cluster or top peaks/scores-based clusters
        for cluster_number, cluster in clusters_to_plot:
            plot_cluster(
                cluster=cluster,
                cluster_number=cluster_number,
                tree=self.tree,
                cmap=cmap,
                outlier=outlier,
                save=save,
                filename=filename,
                hide_internal_nodes=hide_internal_nodes,
                show_terminal_labels=show_terminal_labels,
                width_scale=width_scale,
                height_scale=height_scale,
                label_func=label_func,
                show_branch_lengths=show_branch_lengths,
                marker_size=marker_size,
                outgroup=self.outgroup,
                results_dir=results_dir,
                **kwargs,
            )

    def save(
        self,
        results_dir,
        top_n=1,
        filename="phyclust_results.csv",
        outlier=True,
        n=None,
    ):
        save_clusters(
            scores=self.scores,
            clusters=self.clusters,
            results_dir=results_dir,
            top_n=top_n,
            filename=filename,
            outlier=outlier,
            selected_peaks=self.top_peaks,
            n=n,
        )

    def _track_runtime(self, function_name, start, end):
        runtime = end - start
        self.total_runtime += runtime
        if function_name in self.runtimes:
            self.runtimes[function_name].append(runtime)
        else:
            self.runtimes[function_name] = [runtime]
