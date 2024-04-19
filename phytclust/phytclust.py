import time
from collections import defaultdict
from copy import deepcopy
import string
import logging

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib.ticker import MaxNLocator

from phytclust.plotting import plot_cluster, plot_peaks
from phytclust.save import save_clusters
from phytclust.validation import is_outgroup_valid, validate_tree, rename_nodes


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

    def __init__(self, tree, k=None, max_k=None, score_type=None, outgroup=None):
        self.logger = logging.getLogger('phytclust')
        self.tree = tree
        self.validate_and_set_outgroup(outgroup)
        self.num_terminals = self.calculate_num_terminals()
        validate_tree(self.tree, self.outgroup)
        rename_nodes(self.tree, self.outgroup)
        self.k = k
        self.max_k = (
            self.num_terminals if max_k is None else max_k if k is None else None
        )
        self.clusters = {}
        self.score_type = score_type
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

    def calculate_num_terminals(self):
        if self.outgroup:
            return len(self.tree.get_terminals()) - 1
        return len(self.tree.get_terminals())

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
                # scores = np.full(limit, np.inf)
                scores[0] = (
                    self.dp_table[left][0]
                    + self.dp_table[right][0]
                    + (left.branch_length * (left_size))
                    + (right.branch_length * (right_size))
                )
                back[:, 0] = (0, 0)  # (left, right)

                for i in range(1, limit):
                    sub_scores = self.dp_table[left][:i] + self.dp_table[right][:i][::-1]

                    min_index = np.argmin(sub_scores)
                    min_score = sub_scores[min_index]
                    scores[i] = min_score

                    if not np.isinf(min_score):
                        back[:, i] = (min_index, i - min_index - 1)

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
            self.clusters = self.tree_clust_backtrack(self.k, verbose=verbose)

        elif self.k is None:
            if self.max_k is None or self.max_k == "max":
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

    def cluster_score(self, clusters=None, score_type=None, output_all=False):
        clusters = clusters or self.clusters
        score_type = self.score_type if self.score_type is not None else "CH_score"

        if not clusters:
            raise ValueError("Clusters not found")
        if score_type not in ["CH_score", "PDI", "New_Score"]:
            raise ValueError(f"Invalid score type: {score_type}")

        active_tree = self._no_outgroup_tree if self.outgroup else self.tree
        clust_inv = {}
        for k, v in clusters.items():
            clust_inv.setdefault(v, []).append(k)

        num_clusters = len(clust_inv)
        num_terminals = self.num_terminals

        score = 0
        if score_type == "PDI":
            score = self.pdi_index(clusters)
        elif score_type == "CH_score":
            root_clade = active_tree.root
            dp_row = self.dp_table.get(root_clade, None)
            if dp_row is None:
                raise ValueError("Root not found in dp_table")

            beta_1 = dp_row[0]
            beta = dp_row[num_clusters - 1]

            num = (beta_1 - beta) / (beta + (beta_1 * 0.01))
            den = (
                (num_terminals - num_clusters)
                / (num_clusters - 1)
                if (num_clusters - 1)
                else float("nan")
            )
            score = num * den if den not in [float("inf")] else float("inf")

        return num, den, score

    def score_calc(self, plot=True, output_all=False):
        start_time = time.time()
        scores = []
        num_list = []
        den_list = []
        if self.k is None:
            for cluster in self.clusters:
                num, den, score = self.cluster_score(
                    clusters=cluster, output_all=output_all
                )
                scores.append(score)
                num_list.append(num)
                den_list.append(den)
        else:
            num, den, score = self.cluster_score(
                clusters=self.clusters, output_all=output_all
            )
            print(f"Single Cluster: {self.clusters}, Score: {score}")
            scores.append(score)
            num_list.append(num)
            den_list.append(den)

        self.scores = np.array(scores)
        self.num = np.array(num_list)
        self.den = np.array(den_list)

        if plot and self.k is None:
            self.plot_scores()

        end_time = time.time()
        self._track_runtime("score_calc", start_time, end_time)

    def plot_scores(self):
        fig = plt.figure()
        ax = plt.gca()
        plt.plot(range(1, len(self.scores) + 1), self.scores, "o-")
        plt.xlabel("No. of Clusters")
        plt.ylabel("Scores")
        plt.title("Scores for Different No. of Clusters")

        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        return fig

    def find_peaks(
        self,
        n=3,
        plot=True,
        method="default",
        prominence=None,
        k_start=None,
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

        k_start = k_start or 0
        k_end = k_end or len(self.scores)

        if k_end > len(self.scores) or k_end < k_start:
            raise ValueError(f"Invalid k_end value. Max allowed value is {len(self.scores)}")

        scores_subset = self.scores[k_start:k_end]
        scores_subset = np.where(np.isinf(scores_subset), np.nan, scores_subset)

        if method == "default":
            peaks, _ = find_peaks(scores_subset, prominence=prominence)
        elif method == "alt":
            peaks, _ = find_peaks(scores_subset, prominence=prominence)

        global_maxima_index = np.nanargmax(scores_subset)
        global_maxima_index += k_start
        if global_maxima_index not in peaks:
            peaks = np.append(peaks, global_maxima_index)

        if peaks.size < 1:
            print("No peaks found")
            return []

        sorted_peaks = peaks[np.argsort(-scores_subset[peaks])]
        top_peaks = sorted_peaks[: min(n, len(sorted_peaks))]

        if plot:
            plot_peaks(scores_subset, top_peaks, k_start, k_end=k_end)

        if len(top_peaks) < n:
            print(f"Found only {len(top_peaks)} peak(s)")

        self.top_peaks = list(top_peaks + k_start)
        end_time = time.time()
        self._track_runtime("find_peaks", start_time, end_time)

        return [peak + 1 + k_start for peak in top_peaks]

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
                **kwargs
            )

    def save(
        self, results_dir, top_n=1, filename="phyclust_results.csv", outlier=True, n=None
    ):
        save_clusters(
            scores=self.scores,
            clusters=self.clusters,
            results_dir=results_dir,
            top_n=top_n,
            filename=filename,
            outlier=outlier,
            selected_peaks=self.top_peaks,
            n=n
        )

    def _track_runtime(self, function_name, start, end):
        runtime = end - start
        self.total_runtime += runtime
        if function_name in self.runtimes:
            self.runtimes[function_name].append(runtime)
        else:
            self.runtimes[function_name] = [runtime]
