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

from phytclust.find_peaks import find_plateau_edges, select_representative_edges, normalize, elbow_point
from phytclust.plotting import plot_cluster
from phytclust.save import save_clusters
import phytclust.greedy_alg as greedy_alg
from phytclust.validation import (
    validate_and_set_outgroup, prune_outgroup, resolve_polytomies
)
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    "color", plt.cm.tab20.colors
)  # matplotlib style


@dataclass
class PhytClust:
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

        self.max_k = self.num_terminals if self.k is None and self.max_k is None else self.max_k if self.k is None else None
        self._no_outgroup_tree = copy.deepcopy(self.tree) if self.outgroup else None
        if self.outgroup:
            self.node_terminals, self.terminal_count = prune_outgroup(
                self.tree, self.outgroup
            )
        self.method = self.method if self.method is not None else "default"
        if self.method == "default":
            self.run_dp_method(num_peaks=self.num_peaks, should_plot_scores=self.should_plot_scores)
        elif self.method == "greedy":
            self.run_greedy_alg()  # tbc
        elif self.method == "shannon":
            self.run_shannon_method()
        elif self.method == "stairs2":
            self.run_stairs2_method()

    def run_dp_method(self, num_peaks: int = 3, should_plot_scores: bool = False) -> None:
        self.tree_clust_dp_table()
        self.run(num_peaks=num_peaks, should_plot_scores=should_plot_scores)

    def tree_clust_dp_table(self) -> None:
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
                    sub_scores = left_dp[:i] + right_dp[:i][::-1]

                    min_index = np.argmin(sub_scores)
                    min_score = sub_scores[min_index]
                    scores[i] = min_score

                    if not np.isinf(min_score):
                        back[:, i] = (min_index, i - min_index - 1)

            self.dp_table[clade], self.backtrack[clade] = scores, back

    def tree_clust_backtrack(self, k: Optional[int] = None, verbose: bool = False) -> Dict[Any, int]:
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

    def run(self, verbose: bool = False, num_peaks: int = 3, should_plot_scores: bool = False) -> None:
        """Run the clustering algorithm."""
        if self.max_k is None and self.k:
            self.clusters = self.tree_clust_backtrack(k=self.k, verbose=verbose)
        elif self.k is None:
            self._run_with_max_k(verbose, num_peaks, should_plot_scores)
        else:
            raise ValueError("Invalid cluster number. Please choose either k or max_k")

    def _run_with_max_k(self, verbose: bool, num_peaks: int, should_plot_scores: bool) -> None:
        """Helper method to handle the case when k is None."""
        if self.max_k is None:
            self.max_k = self.num_terminals
        if isinstance(self.max_k, int) and self.max_k > 0:
            self.extract_clusters(verbose)
            self.score_calc(plot=False)
            self.find_peaks(n=num_peaks, plot=should_plot_scores)
            self.rank_peaks()
        else:
            raise ValueError("Invalid cluster number")

    def cluster_score(self, clusters: Optional[Dict[Any, Any]] = None, output_all: bool = False) -> Tuple[float, float, float]:
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

        score = (num*den) if den not in [float("inf")] else float("inf")
        # score = np.log((beta_1/(num_terminals**2))**1/2) + np.log(num_clusters)
        return beta/num_clusters, (beta / beta_1), score/beta_1

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
            num_list, den_list, scores = zip(*results)
        else:
            num, den, score = self.cluster_score(clusters=self.clusters, output_all=output_all)
            print(f"Single Cluster: {self.clusters}, Score: {score}")
            scores, num_list, den_list = [score], [num], [den]

        self.scores = np.array(scores)
        self.num = np.array(num_list)
        self.den = np.array(den_list)

        if plot and self.k is None:
            self.plot_scores(scores_subset=self.scores, peaks=self.peaks, k_start=0)

    def plot_scores(self, scores_subset: np.ndarray, k_start: int, k_end: Optional[int] = None, peaks: Optional[List[int]] = None) -> plt.Figure:
        """Plot scores and highlight peaks if provided."""
        fig, ax = plt.subplots()
        if k_end is None:
            k_end = k_start + len(scores_subset)

        scores_subset_slice = scores_subset[k_start:k_end]
        x_indices = np.arange(k_start + 1, k_end + 1)
        ax.plot(x_indices, scores_subset_slice, "o-", label="Scores", markersize=2)
        if peaks is None:
            peaks, _ = find_peaks(scores_subset_slice)
        if peaks is not None:
            peak_x_indices = [k_start + 1 + p for p in peaks]
            peak_scores = [scores_subset_slice[p] for p in peaks]
            ax.plot(peak_x_indices, peak_scores, "rx", label="Peaks")

        plt.xlabel("No. of Clusters")
        plt.ylabel("Scores")
        plt.title("Scores for Different No. of Clusters")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.show()
        return fig

    def find_peaks(
        self,
        scores: Optional[np.ndarray] = None,
        n: int = 3,
        plot: bool = True,
        k_start: Optional[int] = None,
        k_end: Optional[int] = None,
    ) -> List[int]:
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

        #####
        smoothing = max(1, int(self.num_terminals * 0.1))
        edges, prominence_edges, prominences = find_plateau_edges(
            scores_subset, k_start, smoothing
        )

        self.peak_prominence = prominences

        if len(edges) == 1:
            self.best_peaks = list(peak + k_start + 1 for peak in edges)
            print(f"Found only {len(self.best_peaks)} peak(s)")
            if plot:
                self.plot_scores(scores_subset=scores_subset, peaks=list(edges), k_start=k_start)

            return self.best_peaks

        # edges = list(set(edges))

        representatives, representative_plateau_sizes = select_representative_edges(
            edges, prominence_edges
        )

        # Sort representatives and their plateau sizes by plateau size
        # sorted_indices = np.argsort(representative_plateau_sizes)[::-1]
        representatives = np.array(representatives)  # [sorted_indices]
        representative_plateau_sizes = np.array(
            representative_plateau_sizes
        )  # [sorted_indices]

        top_peaks = list(representatives[: min(n, len(representatives))])
        top_peak_plateau_sizes = list(
            representative_plateau_sizes[: min(n, len(representative_plateau_sizes))]
        )
        self.top_peak_plateau_sizes = top_peak_plateau_sizes

        if plot:
            self.plot_scores(scores_subset=scores_subset, peaks=top_peaks, k_start=k_start)

        # Adjust top_peaks to the original indexing context
        self.best_peaks = list(peak + k_start + 1 for peak in top_peaks)

        if len(self.best_peaks) < n:
            print(f"Found only {len(self.best_peaks)} peak(s)")
        # Return the peaks adjusted to 1-based indexing
        return [peak + 1 for peak in top_peaks]

    def rank_peaks_alt(peak_prominence, best_peaks, den):
        normalized_prominences = normalize(peak_prominence)
        normalized_prominences = np.array(normalized_prominences)
        x = np.arange(len(den)).reshape(-1, 1)
        y = den
        model = LinearRegression()
        model.fit(x, y)
        predictions = model.predict(x)
        r_squared = r2_score(y, predictions)

        indices = [peak - 1 for peak in best_peaks]
        max_product_index = elbow_point(y)
        distances = [abs(i - max_product_index) for i in indices]

        normalized_distances = distances / np.max(distances)
        reciprocal_distances = 1 / (1 + normalized_distances)

        prominence_ranks = np.argsort(np.argsort(-normalized_prominences[indices])) + 1
        distance_ranks = np.argsort(np.argsort(-reciprocal_distances)) + 1

        combined_ranks = prominence_ranks + distance_ranks

        sorted_indices = np.argsort(combined_ranks)
        indices = np.array(indices)
        sorted_indices = sorted_indices

        final_indices = indices[sorted_indices]
        final_ranks = combined_ranks[sorted_indices]

        rankings = np.arange(1, len(final_ranks) + 1)

        ranked_data = [
            (rank, idx + 1, score)
            for rank, idx, score in zip(rankings, final_indices, final_ranks)
        ]

        return ranked_data, [tup[1] for tup in ranked_data]

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
        reciprocal_distances = 1 / (1 + normalized_distances)  # Smoothing to avoid division by zero

        self.adjusted_scores = normalized_scores * reciprocal_distances

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

        # scores = ((1-r_squared) * (reciprocal_distances)) + ((r_squared) * normalized_prominences[indices])

        # sorted_indices = np.argsort(scores)[::-1]
        # indices = np.array(indices)
        # # Use the sorted indices to sort the scores and indices
        # scores = scores[sorted_indices]
        # indices = indices[sorted_indices]

        # # Calculate the rankings
        # rankings = np.arange(1, len(scores) + 1)

        # # Create a list of tuples (ranking, index, score)
        # ranked_data = [
        #     (rank, (idx + 1), score) for rank, idx, score in zip(rankings, indices, scores)
        # ]

        # self.ranked_peaks = ranked_data
        # Rank normalized prominences and reciprocal distances
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

    def run_shannon_method(self, num_peaks=3, should_plot_scores=False):
        self.tree_clust_dp_table()
        self.run(num_peaks=num_peaks, should_plot_scores=should_plot_scores)
        self.extract_clusters_shannon(verbose=False)
        entropy_sums = []

        for clade_dict in self.mrcas:
            total_entropy = 0
            for clade in clade_dict.keys():  # Iterate over the keys, which are Clade objects
                # Calculate the Shannon entropy for the clade
                total_entropy += self.shannon_entropy(clade)
            entropy_sums.append(total_entropy)
        self.shannon_score = entropy_sums
        scores_array = np.array(self.scores)
        factors_array = np.array(entropy_sums)
        new_scores = scores_array * (1+factors_array)
        positions = np.arange(1, len(new_scores) + 1)
        self.adjusted_scores = new_scores / (positions)
        # Find peaks and their properties
        peaks, properties = find_peaks(self.adjusted_scores, prominence=0.0001)

        # Filter out peaks with negative values in adjusted_scores
        positive_peaks = [
            peak for peak in peaks
            if (peak > 0 and self.adjusted_scores[peak - 1] > 0) and 
            (peak < len(self.adjusted_scores) - 1 and self.adjusted_scores[peak + 1] > 0)
        ]
        positive_peaks = np.array(positive_peaks)
        positive_prominences = properties['prominences'][np.isin(peaks, positive_peaks)]

        # Sort peaks by prominence
        sorted_indices = np.argsort(positive_prominences)[::-1]
        sorted_peaks = positive_peaks[sorted_indices]
        # Ensure sorted_peaks is an integer array
        if len(sorted_peaks) == 1:
            sorted_peaks = np.array([sorted_peaks[0]], dtype=int)
        else:
            sorted_peaks = sorted_peaks.astype(int)

        # Save the sorted peaks and their values
        self.best_peaks = sorted_peaks + 1
        self.best_peak_values = self.adjusted_scores[sorted_peaks] 
        self.scores = self.adjusted_scores
        self.rank_peaks()
        if should_plot_scores:
            plt.figure(figsize=(10, 6))
            x_values = np.arange(1, len(self.adjusted_scores) + 1)
            plt.plot(x_values, self.adjusted_scores, label="Scores", color="blue")

            # Highlight the peaks
            plt.scatter(self.best_peaks, self.best_peak_values, color='red', marker='o', label='Peaks')

            # Add labels and title
            plt.xlabel('Position')
            plt.ylabel('Score')
            plt.title('Scores with Identified Peaks')
            plt.legend()

            # Show the plot
            plt.show()

    def run_stairs2_method(self, num_peaks=3, should_plot_scores=False):
        self.tree_clust_dp_table()
        self.run(num_peaks=num_peaks, should_plot_scores=should_plot_scores)
        self.extract_clusters_shannon(verbose=False)
        self.stairs2()
        scores_array = np.array(self.scores)
        factors_array = np.array(self.ratio_sums)
        new_scores = scores_array * (1 + factors_array)
        positions = np.arange(1, len(new_scores) + 1)
        self.adjusted_scores = new_scores / (positions)
        # Find peaks and their properties
        peaks, properties = find_peaks(self.adjusted_scores, prominence=0.0001)

        # Filter out peaks with negative values in adjusted_scores
        positive_peaks = [
            peak for peak in peaks
            if (peak > 0 and self.adjusted_scores[peak - 1] > 0) and 
            (peak < len(self.adjusted_scores) - 1 and self.adjusted_scores[peak + 1] > 0)
        ]
        positive_peaks = np.array(positive_peaks)
        positive_prominences = properties['prominences'][np.isin(peaks, positive_peaks)]

        # Sort peaks by prominence
        sorted_indices = np.argsort(positive_prominences)[::-1]
        sorted_peaks = positive_peaks[sorted_indices]
        self.scores = self.adjusted_scores
        # Ensure sorted_peaks is an integer array
        if len(sorted_peaks) == 1:
            sorted_peaks = np.array([sorted_peaks[0]], dtype=int)
        else:
            sorted_peaks = sorted_peaks.astype(int)
        # Save the sorted peaks and their values
        self.best_peaks = sorted_peaks + 1
        self.best_peak_values = self.adjusted_scores[sorted_peaks] 
        self.rank_peaks()

        if should_plot_scores:
            plt.figure(figsize=(10, 6))
            x_values = np.arange(1, len(self.adjusted_scores) + 1)
            plt.plot(x_values, self.adjusted_scores, label="Scores", color="blue")

            # Highlight the peaks
            plt.scatter(self.best_peaks, self.best_peak_values, color='red', marker='o', label='Peaks')

            # Add labels and title
            plt.xlabel('Position')
            plt.ylabel('Score')
            plt.title('Scores with Identified Peaks')
            plt.legend()

            # Show the plot
            plt.show()

    def tree_clust_backtrack_shannon(self, k=None, verbose=False):
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
                # node_terminals = self.node_terminals[node]
                clusters[node] = current_cluster
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

    def extract_clusters_shannon(self, verbose):
        self.mrcas = [
            self.tree_clust_backtrack_shannon(i, verbose=verbose)
            for i in range(1, self.max_k + 1)
        ]

    def run_greedy_alg(self, num_peaks=3, should_plot_scores=False):
        for clade in self.tree.find_clades(order="postorder"):
            clade.children = list(clade.clades)  # Ensuring children are a list, not a generator
            clade.terminal_child_count = self.calculate_terminal_child_count(clade)

            # Calculate self_cost efficiently
            clade.self_cost = sum(
                child.terminal_child_count * (child.branch_length or 0) for child in clade.children
            )

            # Calculate total_recursive_cost efficiently
            if not clade.children:  # If leaf node, total_recursive_cost is just self_cost
                clade.total_recursive_cost = clade.self_cost
            else:
                clade.total_recursive_cost = clade.self_cost + sum(
                    child.total_recursive_cost for child in clade.children
                )

        scores, states, total = greedy_alg.split_tree(
            self.tree.root, max_splits=self.num_terminals
        )
        beta_scores, self.num, self.den = greedy_alg.calculate_beta_scores(
            scores, self.num_terminals
        )
        self.states = states
        clades = self.node_terminals   
        # Create a dictionary for clades by name
        clades_by_name = {clade.name: clade for clade in clades}

        all_terminal_clusters = []

        for node_list in states:
            terminal_clusters = {}
            cluster_number = 0  # Reset cluster number for each node_list
            for internal_node in node_list:
                # Find the clade with the matching name using the dictionary
                if internal_node in clades_by_name:
                    clade = clades_by_name[internal_node]
                    terminal_nodes = clades[clade]
                    # Assign the current cluster number to these terminal nodes
                    terminal_clusters.update({terminal_node: cluster_number for terminal_node in terminal_nodes})
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

    def shannon_entropy(self, node):
        # Find all child nodes of the given node
        children = [child for child in node.clades]

        # Calculate the size of each subtree rooted at each child node by counting leaves
        sizes = [self.count_leaves(self.tree, child) for child in children]
        total_size = sum(sizes)
        if total_size == 0:
            return 0  # To handle cases where total size is 0 to avoid division by zero

        # Calculate probabilities
        probabilities = [size / total_size for size in sizes]

        # Calculate Shannon's entropy
        entropy = -sum(
            p * log(p, 2) for p in probabilities if p > 0
        )  # Ensure p > 0 to avoid log(0)

        # Multiplying entropy by total size of the node (optional, depends on the interpretation)
        return entropy * 1/total_size

    def stairs2(self):
        ratio_sums = []
        terminal_count = self.terminal_count
        for mrcas_dict in self.mrcas:
            sum_ratios = 0

            for clade in mrcas_dict.keys():
                clade_name = clade.name
            corresponding_clade = next(
                (c for c in terminal_count.keys() if c.name == clade_name), None
            )
            # print(f"Corresponding clade: {corresponding_clade}")

            if corresponding_clade and terminal_count[corresponding_clade] > 1:
                children = corresponding_clade.clades
                if len(children) == 2:
                    child1, child2 = children
                    count1 = terminal_count.get(child1, 0)
                    count2 = terminal_count.get(child2, 0)

                    # print(f"Child1: {child1}, Count1: {count1}")
                    # print(f"Child2: {child2}, Count2: {count2}")

                    if count1 > 0 and count2 > 0:
                        lower_count = min(count1, count2)
                        higher_count = max(count1, count2)
                        ratio = lower_count / higher_count
                        sum_ratios += ratio
            ratio_sums.append(sum_ratios)
        self.ratio_sums = ratio_sums
