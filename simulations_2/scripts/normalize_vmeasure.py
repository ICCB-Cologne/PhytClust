import logging
import copy
from math import log
from dataclasses import dataclass, field
from typing import Any, Optional, Dict, List, Tuple
from functools import lru_cache


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib.ticker import MaxNLocator
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
import matplotlib.colors as mcolors
import random
from phytclust.find_peaks import (
    find_plateau_edges,
    select_representative_edges,
    normalize,
    elbow_point,
)
from phytclust.plotting import plot_cluster
from phytclust.save import save_clusters
import phytclust.greedy_alg as greedy_alg
from phytclust.validation import (
    validate_and_set_outgroup,
    prune_outgroup,
    resolve_polytomies,
)

plt.rcParams["axes.prop_cycle"] = plt.cycler(
    "color", plt.cm.tab20.colors
)  # matplotlib style

num_clusters = 80


@dataclass
class PhytClust_worst:
    """
    This class uses a dynamic programming algorithm to cluster nodes on a phylogenetic tree to annotate monophyletic clades

    Parameters:
    -   tree: a newick file
    -   k: number of clusters you would like
    -   max_k: if you'd like to find optimal clusters for a range of 1 to max_k clusters, if k is not given, Default = number of total terminal nodes in a tree
    """

    tree: Any
    num_peaks: int = 3
    should_plot_scores: bool = False
    resolution_on: bool = True
    num_bins: int = 3
    k: Optional[int] = None
    max_k: Optional[int] = None
    outgroup: Optional[Any] = None
    method: Optional[str] = None
    logger: logging.Logger = field(init=False, default=logging.getLogger("phytclust"))
    node_terminals: Dict[Any, List[Any]] = field(init=False, default_factory=dict)
    terminal_count: Dict[Any, int] = field(init=False, default_factory=dict)
    num_terminals: int = field(init=False)
    beta_1: Optional[Any] = field(init=False, default=None)
    terminal_child_count_cache: Dict[Any, int] = field(init=False, default_factory=dict)
    clusters: Dict[Any, Any] = field(init=False, default_factory=dict)
    dp_table: Dict[Any, Any] = field(init=False, default_factory=dict)
    backtrack: Dict[Any, Any] = field(init=False, default_factory=dict)
    scores: List[float] = field(init=False, default_factory=list)
    best_peaks: List[float] = field(init=False, default_factory=list)
    highest_peaks: List[float] = field(init=False, default_factory=list)
    states: List[Any] = field(init=False, default_factory=list)
    _no_outgroup_tree: Optional[Any] = field(init=False, default=None)
    num: List[float] = field(init=False, default_factory=list)
    den: List[float] = field(init=False, default_factory=list)

    def __post_init__(self):
        self.tree = copy.deepcopy(self.tree)
        self.initialize()

    def initialize(self) -> None:
        self.tree, self.outgroup = validate_and_set_outgroup(self.tree, self.outgroup)
        resolve_polytomies(self.tree)
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

        self.max_k = (
            self.num_terminals
            if self.k is None and self.max_k is None
            else self.max_k if self.k is None else None
        )
        self._no_outgroup_tree = copy.deepcopy(self.tree) if self.outgroup else None
        if self.outgroup:
            self.node_terminals, self.terminal_count = prune_outgroup(
                self._no_outgroup_tree, self.outgroup
            )
        self.method = self.method if self.method is not None else "default"
        if self.method == "default":
            self.run_dp_method(
                num_peaks=self.num_peaks,
                should_plot_scores=self.should_plot_scores,
                resolution_on=self.resolution_on,
                num_bins=self.num_bins,
            )
        elif self.method == "greedy":
            self.run_greedy_alg()  # tbc

    def run_dp_method(
        self,
        num_peaks: int = 3,
        should_plot_scores: bool = False,
        resolution_on: bool = True,
        num_bins: int = 3,
    ) -> None:
        self.tree_clust_dp_table()
        self.run(
            num_peaks=num_peaks,
            should_plot_scores=should_plot_scores,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )

    def tree_clust_dp_table(self) -> None:
        n = self.k or self.max_k
        active_tree = self._no_outgroup_tree if self.outgroup else self.tree

        for clade in active_tree.find_clades(order="postorder"):
            scores = np.full(n, -np.inf)
            back = np.full((2, n), -np.inf)
            if clade.is_terminal():
                scores[0] = 0
            else:
                left, right = clade.clades
                left_size = self.terminal_count[left]
                right_size = self.terminal_count[right]
                total_terminals = left_size + right_size

                limit = min(n, total_terminals)

                left_dp = self.dp_table[left]
                right_dp = self.dp_table[right]
                left_branch_length = left.branch_length
                right_branch_length = right.branch_length

                scores[0] = (
                    left_dp[0]
                    + right_dp[0]
                    + (left_branch_length * left_size)
                    + (right_branch_length * right_size)
                )
                back[:, 0] = (0, 0)

                for i in range(1, limit):
                    left_indices = np.arange(i)
                    right_indices = i - left_indices - 1

                    left_dp_values = left_dp[left_indices]
                    right_dp_values = right_dp[right_indices]

                    sub_scores = left_dp_values + right_dp_values

                    valid_indices = ~np.isinf(sub_scores)

                    if np.any(valid_indices):
                        valid_sub_scores = sub_scores[valid_indices]
                        max_indices = np.where(valid_sub_scores == np.max(valid_sub_scores))[0]
                        chosen_index = np.random.choice(max_indices)
                        max_score = valid_sub_scores[chosen_index]

                        scores[i] = max_score

                        if not np.isinf(max_score):
                            back[:, i] = (
                                left_indices[valid_indices][chosen_index],
                                right_indices[valid_indices][chosen_index],
                            )

            self.dp_table[clade], self.backtrack[clade] = scores, back

    def tree_clust_backtrack(
        self, k: Optional[int] = None, verbose: bool = False
    ) -> Dict[Any, int]:
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
                node_terminals = self.node_terminals[node]
                clusters.update(
                    {terminal: current_cluster for terminal in node_terminals}
                )
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

    def extract_clusters(self, verbose: bool) -> None:
        """Define clusters using tree_clust_backtrack."""
        self.clusters = [
            self.tree_clust_backtrack(i, verbose=verbose)
            for i in range(1, self.max_k + 1)
        ]

    def run(
        self,
        verbose: bool = False,
        num_peaks: int = 3,
        should_plot_scores: bool = False,
        resolution_on: bool = False,
        num_bins: int = 3,
    ) -> None:
        """Run the clustering algorithm."""
        if self.max_k is None and self.k:
            self.clusters = self.tree_clust_backtrack(k=self.k, verbose=verbose)
        elif self.k is None:
            self._run_with_max_k(
                verbose, num_peaks, should_plot_scores, resolution_on, num_bins
            )
        else:
            raise ValueError("Invalid cluster number. Please choose either k or max_k")

    def _run_with_max_k(
        self,
        verbose: bool,
        num_peaks: int,
        should_plot_scores: bool,
        resolution_on: bool = False,
        num_bins: int = 3,
    ) -> None:
        """Helper method to handle the case when k is None."""
        if self.max_k is None:
            self.max_k = self.num_terminals
        if isinstance(self.max_k, int) and self.max_k > 0:
            self.extract_clusters(verbose)
            self.score_calc(plot=should_plot_scores)
            self.find_peaks(
                n=num_peaks,
                plot=should_plot_scores,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )
        else:
            raise ValueError("Invalid max_k value. It should be a positive integer.")

    def cluster_score(
        self, clusters: Optional[Dict[Any, Any]] = None, output_all: bool = False
    ) -> Tuple[float, float, float]:
        """Calculate the cluster score."""
        clusters = clusters or self.clusters

        if not clusters:
            raise ValueError("Clusters not found")

        active_tree = self._no_outgroup_tree if self.outgroup else self.tree
        clust_inv = {}
        for k, v in clusters.items():
            clust_inv.setdefault(v, set()).add(k)

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
            (num_terminals - num_clusters) / (num_clusters)
            if (num_clusters)
            else float("nan")
        )

        score = (num * den) if den not in [float("inf")] else float("inf")
        # score = np.log((beta_1/(num_terminals**2))**1/2) + np.log(num_clusters)
        return beta, score, score

    def score_calc(self, plot: bool = False, output_all: bool = False) -> None:
        """Calculate scores for all clusters."""
        self.scores = np.zeros(len(self.clusters))
        self.num = np.zeros(len(self.clusters))
        self.den = np.zeros(len(self.clusters))
        if self.k is None:
            results = [
                self.cluster_score(clusters=cluster, output_all=output_all)
                for cluster in self.clusters
            ]
            beta_values, den_list, scores = zip(*results)
        else:
            beta_values, den, scores = self.cluster_score(
                clusters=self.clusters, output_all=output_all
            )
            scores, beta_values, den_list = [scores], [beta_values], [den]

        # adding elbow index
        scores = np.array(scores)
        scores[scores < 0] = 0

        beta_values = np.array(beta_values)
        beta_values[beta_values < 0] = 0
        beta_values = np.nan_to_num(beta_values, nan=0.0, posinf=0.0, neginf=0.0)

        elbow_scores = []
        n = len(beta_values)
        for k in range(n - 1):
            if beta_values[k] - beta_values[k + 1] != 0:
                elbow_score = (
                    (beta_values[k - 1] - beta_values[k])
                    / (beta_values[k] - beta_values[k + 1])
                    if k > 0
                    else 0
                )
            else:
                elbow_score = 0
            elbow_scores.append(elbow_score)
        elbow_scores.append(0)

        min_invalid_index = None

        invalid_new_scores = np.where(np.isnan(elbow_scores) | np.isinf(elbow_scores))[
            0
        ]
        invalid_clust_scores = np.where(np.isnan(scores) | np.isinf(scores))[0]

        if len(invalid_new_scores) > 0:
            min_invalid_index = invalid_new_scores[0]
        if len(invalid_clust_scores) > 0:
            if min_invalid_index is not None:
                min_invalid_index = min(min_invalid_index, invalid_clust_scores[0])
            else:
                min_invalid_index = invalid_clust_scores[0]
        if min_invalid_index is None:
            min_invalid_index = len(elbow_scores)

        elbow_scores = elbow_scores[:min_invalid_index]
        scores = scores[:min_invalid_index]

        def min_max_normalize(arr):
            min_val = np.min(arr)
            max_val = np.max(arr)
            if max_val - min_val != 0:
                return (arr - min_val) / (max_val - min_val)
            else:
                return arr

        normalized_elbow_scores = min_max_normalize(elbow_scores)
        normalized_initial_scores = min_max_normalize(scores)
        combined_scores = normalized_elbow_scores * normalized_initial_scores

        combined_scores = np.nan_to_num(
            combined_scores, nan=0.0, posinf=0.0, neginf=0.0
        )
        self.scores = np.array(combined_scores)
        self.beta_values = np.array(beta_values)
        self.den = np.array(den_list)

    def plot_scores(
        self,
        scores_subset: np.ndarray,
        k_start: int,
        k_end: Optional[int] = None,
        peaks: Optional[List[int]] = None,
        resolution_on=True,
        num_bins=3,
        fig_width=14,
        fig_height=6,
        log_base=None,
    ) -> plt.Figure:
        """Plot scores and highlight peaks if provided."""
        if scores_subset is None:
            scores_subset = self.scores

        # Normalize scores to a scale of 0 to 1
        scores_subset = (scores_subset - np.min(scores_subset)) / (
            np.max(scores_subset) - np.min(scores_subset)
        )

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        if k_end is None:
            k_end = k_start + len(scores_subset)

        scores_subset_slice = scores_subset[k_start:k_end]
        x_indices = np.arange(k_start + 1, k_end + 1)
        ax.plot(x_indices, scores_subset_slice, "o-", label="Scores", markersize=2)

        if peaks is not None:
            peak_x_indices = [peak for peak in peaks]
            peak_scores = [scores_subset_slice[peak - k_start - 1] for peak in peaks]
            ax.plot(peak_x_indices, peak_scores, "rx", label="Peaks", markersize=8)
            for peak in peaks:
                ax.text(
                    peak,
                    scores_subset_slice[peak - k_start - 1]
                    + 0.02,  # Offset to avoid overlap
                    str(peak),
                    fontsize=12,
                    color="red",
                    ha="center",
                )

        if resolution_on:
            xmin = 1
            xmax = len(scores_subset_slice)

            # Calculate the base for the logarithm of bins
            if log_base is None:
                log_base = (xmax / xmin) ** (1 / num_bins)

            # Define bins using the specified logarithmic base
            bins = xmin * (log_base ** np.arange(num_bins + 1))
            bins = np.ceil(bins).astype(int) + 1

            if num_bins == 3:
                colors = ["#5DADE2", "#58D68D", "#EC7063"]
                darker_colors = ["#2874A6", "#239B56", "#C0392B"]
                labels = ["Low", "Intermediate", "High"]
            elif num_bins == 2:
                colors = ["#5DADE2", "#EC7063"]
                darker_colors = ["#2874A6", "#C0392B"]
                labels = ["Low", "High"]
            else:
                # Generate random light colors excluding white
                light_colors = [
                    color
                    for color in mcolors.CSS4_COLORS.values()
                    if not color.lower() in ["white", "#ffffff"]
                ]
                colors = random.sample(light_colors, num_bins)
                labels = [f"Bin {i + 1}" for i in range(num_bins)]

            # Add the shaded areas for each bin
            for i in range(len(bins) - 1):
                plt.axvspan(
                    bins[i],
                    bins[i + 1],
                    color=colors[i],
                    alpha=0.2,
                    label=(
                        f"Bin {i + 1}: {bins[i]} - {bins[i + 1]}"
                        if num_bins in [2, 3]
                        else None
                    ),
                )

            # Add vertical lines for each bin edge
            for bin_edge in bins:
                plt.axvline(x=bin_edge, color="grey", linestyle="--", linewidth=1)

            # Add vertical labels within the plot if num_bins is 2 or 3
            if num_bins in [2, 3]:
                for i, label in enumerate(labels):
                    plt.text(
                        bins[i] + (bins[i + 1] - bins[i]) / 2,
                        0.5,  # Position within the plot
                        label,
                        verticalalignment="center",
                        horizontalalignment="center",
                        fontsize=12,
                        color=darker_colors[i],
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                        rotation="horizontal",
                        clip_on=True,  #
                    )

        plt.xlabel("No. of Clusters", fontsize=14, labelpad=10)
        plt.ylabel("Scores (log)", fontsize=14, labelpad=10)
        plt.yscale("log")
        plt.title("Scores", fontsize=16, pad=20)
        plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tick_params(axis="both", which="major", labelsize=12)
        plt.tick_params(axis="both", which="minor", labelsize=10)
        plt.tight_layout()
        plt.show()
        return fig

    def find_peaks(
        self,
        scores: Optional[np.ndarray] = None,
        n: int = 3,
        plot: bool = True,
        k_start: Optional[int] = None,
        k_end: Optional[int] = None,
        resolution_on: bool = True,
        num_bins: int = 3,
    ) -> List[int]:
        """
        Finds and plots the top n prominent peaks in the scores within a specified range of k.

        Parameters:
        - n: int, the number of top peaks to find.
        - plot: bool, whether to plot the score graph with the peaks highlighted.
        - k_start: int or None, the starting index for the range of k values to consider.
        - k_end: int or None, the ending index for the range of k values to consider.
        - resolution_on: bool, whether to use resolution levels.
        - num_bins: int, the number of bins for resolution levels.
        """
        if scores is None:
            scores = self.scores  # Use the instance attribute if no argument was passed

        if len(scores) == 0:
            print("Please calculate scores first")
            return []

        # Adjust k_start for 0-based indexing within Python
        k_start = k_start if k_start is not None else 0
        k_end = k_end or len(scores)

        if k_end > len(scores) or k_end < k_start:
            raise ValueError(
                f"Invalid k_end value. Max allowed value is {len(self.scores)}"
            )

        scores_subset = scores[k_start:k_end]
        scores_subset = np.where(np.isinf(scores_subset), np.nan, scores_subset)

        # Transform scores to log scale to make the method independent of scale
        log_scores_subset = np.log(
            scores_subset + 1e-10
        )  # Adding a small value to avoid log(0)

        peaks, properties = find_peaks(log_scores_subset, prominence=1e-10)

        # Extract prominences and sort by them in descending order
        prominences = properties["prominences"]
        sorted_indices = np.argsort(prominences)[::-1]
        # Retrieve the top n peaks (sorted by prominence)
        top_peaks_indices = sorted_indices[:n]
        top_peaks = peaks[top_peaks_indices]
        self.best_peaks = [peak + k_start + 1 for peak in top_peaks]

        if plot:
            self.plot_scores(
                scores_subset=scores_subset,
                peaks=self.best_peaks,
                k_start=k_start,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )

        return self.best_peaks

    def plot(
        self,
        results_dir: Optional[str] = None,
        top_n: int = 1,
        n: Optional[int] = None,
        cmap: plt.cm = plt.get_cmap("tab20"),
        show_terminal_labels: bool = False,
        outlier: bool = False,
        save: bool = False,
        filename: Optional[str] = None,
        hide_internal_nodes: bool = True,
        width_scale: float = 2,
        height_scale: float = 0.1,
        label_func: Optional[callable] = None,
        show_branch_lengths: bool = False,
        marker_size: int = 40,
        **kwargs,
    ) -> None:
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
            if hasattr(self, "ranked_peaks") and self.ranked_peaks:
                num_peaks = len(self.ranked_peaks)
                top_peaks = [
                    self.ranked_peaks[i][1] for i in range(min(top_n, num_peaks))
                ]
                clusters_to_plot.extend(
                    [(peak, self.clusters[peak - 1]) for peak in top_peaks]
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
        for _, cluster in clusters_to_plot:
            plot_cluster(
                cluster=cluster,
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
        results_dir: str,
        top_n: int = 1,
        filename: str = "phyclust_results.csv",
        outlier: bool = True,
        n: Optional[int] = None,
    ) -> None:
        save_clusters(
            scores=self.scores,
            clusters=self.clusters,
            results_dir=results_dir,
            top_n=top_n,
            filename=filename,
            outlier=outlier,
            selected_peaks=self.best_peaks,
            n=n,
        )

    def adjust_for_elbow(self):
        normalized_scores = normalize(self.scores)
        x = np.arange(len(self.den)).reshape(-1, 1)
        y = self.den

        # Assuming the elbow point is calculated directly from y
        max_product_index = elbow_point(y)
        distances = np.abs(np.arange(len(normalized_scores)) - max_product_index)

        normalized_distances = distances / np.max(distances)
        reciprocal_distances = 1 / (
            1 + normalized_distances
        )  # Smoothing to avoid division by zero

        self.adjusted_scores = normalized_scores * reciprocal_distance

    def rank_peaks(self) -> None:
        normalized_prominences = normalize(self.peak_prominence)
        normalized_prominences = np.array(normalized_prominences)
        x = np.arange(len(self.den)).reshape(-1, 1)
        y = self.den
        model = LinearRegression()
        model.fit(x, y)
        predictions = model.predict(x)
        r_squared = r2_score(y, predictions)

        indices = [peak - 1 for peak in self.best_peaks]
        max_product_index = elbow_point(y)
        distances = [abs(i - max_product_index) for i in indices]

        normalized_distances = distances / np.max(distances)
        reciprocal_distances = 1 / (1 + normalized_distances)
        prominence_ranks = np.argsort(np.argsort(-normalized_prominences[indices])) + 1
        distance_ranks = np.argsort(np.argsort(-reciprocal_distances)) + 1

        # Combine rankings by summing: Lower scores are better
        combined_ranks = prominence_ranks + distance_ranks

        # Sort by combined ranks in ascending order (lower is better)
        sorted_indices = np.argsort(combined_ranks)
        indices = np.array(indices)
        sorted_indices = (
            sorted_indices  # This sorts in ascending order, lowest rank to highest
        )

        # Use the sorted indices to sort the indices and ranks
        final_indices = indices[sorted_indices]
        final_ranks = combined_ranks[sorted_indices]

        # Calculate the rankings based on sorted combined ranks
        rankings = np.arange(1, len(final_ranks) + 1)

        # Create a list of tuples (ranking, index, combined rank)
        ranked_data = [
            (rank, idx + 1, score)
            for rank, idx, score in zip(rankings, final_indices, final_ranks)
        ]

        # Assign ranked data to class attribute
        self.ranked_peaks = ranked_data
        self.best_peaks = [tup[1] for tup in ranked_data]

    def run_greedy_alg(self, num_peaks=3, should_plot_scores=False):
        for clade in self.tree.find_clades(order="postorder"):
            clade.children = list(
                clade.clades
            )  # Ensuring children are a list, not a generator
            clade.terminal_child_count = self.calculate_terminal_child_count(clade)

            # Calculate self_cost efficiently
            clade.self_cost = sum(
                child.terminal_child_count * (child.branch_length or 0)
                for child in clade.children
            )

            # Calculate total_recursive_cost efficiently
            if (
                not clade.children
            ):  # If leaf node, total_recursive_cost is just self_cost
                clade.total_recursive_cost = clade.self_cost
            else:
                clade.total_recursive_cost = clade.self_cost + sum(
                    child.total_recursive_cost for child in clade.children
                )

        if self.k is not None:
            # Use the provided k value to get the terminal clusters for the last split
            _, states, _ = greedy_alg.split_tree(self.tree.root, max_splits=self.k - 1)
            self.states = states
            clades = self.node_terminals
            # Create a dictionary for clades by name
            clades_by_name = {clade.name: clade for clade in clades}

            # all_terminal_clusters = []

            for node_list in states:
                terminal_clusters = {}
                cluster_number = 0  # Reset cluster number for each node_list
                for internal_node in node_list:
                    # Find the clade with the matching name using the dictionary
                    clade = clades_by_name.get(internal_node)
                    if clade:
                        terminal_nodes = clades[clade]
                        # Assign the current cluster number to these terminal nodes
                        terminal_clusters.update(
                            {
                                terminal_node: cluster_number
                                for terminal_node in terminal_nodes
                            }
                        )
                    # Increment the cluster number counter for the next internal node
                    cluster_number += 1
                # Append the dictionary to the list of mappings
                self.clusters = terminal_clusters

        else:
            # Use the default num_terminals value
            scores, states, total = greedy_alg.split_tree(
                self.tree.root, max_splits=self.max_k
            )

            beta_scores, self.num, self.den = greedy_alg.calculate_beta_scores(
                scores, self.num_terminals
            )
            self.states = states
            # self.beta_list = self.num
            clades = self.node_terminals
            # Create a dictionary for clades by name
            clades_by_name = {clade.name: clade for clade in clades}

            all_terminal_clusters = []

            for node_list in states:
                terminal_clusters = {}
                cluster_number = 0  # Reset cluster number for each node_list
                for internal_node in node_list:
                    # Find the clade with the matching name using the dictionary
                    clade = clades_by_name.get(internal_node)
                    if clade:
                        terminal_nodes = clades[clade]
                        # Assign the current cluster number to these terminal nodes
                        terminal_clusters.update(
                            {
                                terminal_node: cluster_number
                                for terminal_node in terminal_nodes
                            }
                        )
                    # Increment the cluster number counter for the next internal node
                    cluster_number += 1
                # Append the dictionary to the list of mappings
                all_terminal_clusters.append(terminal_clusters)

            self.clusters = all_terminal_clusters
            self.scores = np.array(beta_scores)
            self.find_peaks(n=1000, plot=should_plot_scores)
            self.rank_peaks()

    # @staticmethod
    def calculate_terminal_child_count(self, clade):
        if clade in self.terminal_child_count_cache:
            return self.terminal_child_count_cache[clade]

        if not clade.clades:
            self.terminal_child_count_cache[clade] = 1
        else:
            self.terminal_child_count_cache[clade] = sum(
                self.calculate_terminal_child_count(child) for child in clade.clades
            )

        return self.terminal_child_count_cache[clade]

    def count_leaves(self, tree, node):
        # Recursively count the number of leaves in the subtree rooted at `node`
        if node.is_terminal():
            return 1
        return sum(self.count_leaves(tree, child) for child in node.clades)

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
                    clust_obj = PhytClust_worst(
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
    output_dir = f"/home/ganesank/project/phytclust/simulations_2/results/phytclust_single_k_worst/N_100_K_{num_clusters}_short"
    process_trees(case_folder, output_dir)
    print("Processing completed.")
