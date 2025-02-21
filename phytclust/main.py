import os
import re
import copy
import random
import logging
from math import log, ceil
from functools import lru_cache
from dataclasses import dataclass, field
from typing import Any, Optional, Dict, List, Tuple
import multiprocessing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LogLocator, LogFormatter
from scipy.signal import find_peaks
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.colors as mcolors


from phytclust.find_peaks import (
    find_plateau_edges,
    select_representative_edges,
    normalize,
    elbow_point,
)
from phytclust.plotting import plot_cluster
from phytclust.save import save_clusters
from phytclust.validation import (
    validate_and_set_outgroup,
    prune_outgroup,
    resolve_polytomies,
)

plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors)


@dataclass
class PhytClust:
    """
    Implements a dynamic programming-based phylogenetic clustering algorithm.

    Args:
    tree: a newick file
    k: number of clusters you would like
    max_k: if you'd like to find optimal clusters for a range of 1 to max_k clusters, if k is not given, Default = 90% of no. of total terminal nodes in a tree
    """

    tree: Any
    use_branch_support_values: bool = False
    num_peaks: int = 3
    should_plot_scores: bool = False
    compute_all_clusters: bool = (False,)
    resolution_on: bool = False
    max_k_limit: int = 0.9
    num_bins: int = 3
    k: Optional[int] = None
    max_k: Optional[int] = None
    outgroup: Optional[Any] = None
    method: Optional[str] = None
    logger: logging.Logger = field(init=False, default=logging.getLogger("phytclust"))

    name_leaves_per_node: Dict[Any, List[Any]] = field(init=False, default_factory=dict)
    num_leaves_per_node: Dict[Any, int] = field(init=False, default_factory=dict)
    num_terminals: int = field(init=False)  # N

    dp_table: Dict[Any, Any] = field(init=False, default_factory=dict)
    backtrack: Dict[Any, Any] = field(init=False, default_factory=dict)

    postorder_nodes: List[Any] = field(init=False, default_factory=list)
    node_to_id: Dict[Any, int] = field(init=False, default_factory=dict)

    beta_1: Optional[Any] = field(init=False, default=None)
    terminal_child_count_cache: Dict[Any, int] = field(init=False, default_factory=dict)
    clusters: Dict[Any, Any] = field(init=False, default_factory=dict)
    scores: List[float] = field(init=False, default_factory=list)
    peaks_by_rank: List[float] = field(init=False, default_factory=list)
    states: List[Any] = field(init=False, default_factory=list)
    _no_outgroup_tree: Optional[Any] = field(init=False, default=None)
    beta_ratios: List[float] = field(init=False, default_factory=list)
    norm_ratios: List[float] = field(init=False, default_factory=list)

    def __post_init__(self):
        # self.tree = copy.deepcopy(self.tree)
        self.initialize()

    def initialize(self) -> None:
        self.tree, self.outgroup = validate_and_set_outgroup(self.tree, self.outgroup)
        resolve_polytomies(self.tree)
        self.name_leaves_per_node = {
            node: node.get_terminals() for node in self.tree.find_clades()
        }
        self.num_leaves_per_node = {
            node: len(terminals)
            for node, terminals in self.name_leaves_per_node.items()
        }
        self.num_terminals = (
            (self.num_leaves_per_node[self.tree.root] - 1)
            if self.outgroup
            else self.num_leaves_per_node[self.tree.root]
        )

        self.max_k = (
            ceil(self.num_terminals * self.max_k_limit)
            if self.k is None and self.max_k is None
            else self.max_k if self.k is None else None
        )
        self._no_outgroup_tree = copy.deepcopy(self.tree) if self.outgroup else None
        if self.outgroup:
            self.name_leaves_per_node, self.num_leaves_per_node = prune_outgroup(
                self._no_outgroup_tree, self.outgroup
            )
        self.method = self.method if self.method is not None else "default"
        if self.method == "default":
            self.run_dp_clustering(
                num_peaks=self.num_peaks,
                should_plot_scores=self.should_plot_scores,
                resolution_on=self.resolution_on,
                num_bins=self.num_bins,
            )
        elif self.method == "greedy":
            self.run_greedy_alg()  # tbc

    def run_dp_clustering(
        self,
        num_peaks: int = 3,
        should_plot_scores: bool = False,
        resolution_on: bool = True,
        num_bins: int = 3,
    ) -> None:
        self.compute_dp_table()
        self.execute_clustering(
            num_peaks=num_peaks,
            should_plot_scores=should_plot_scores,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )

    def compute_dp_table(self) -> None:
        n = self.k or self.max_k
        active_tree = self._no_outgroup_tree if self.outgroup else self.tree
        self.postorder_nodes = list(active_tree.find_clades(order="postorder"))
        num_nodes = len(self.postorder_nodes)

        self.node_to_id = {node: i for i, node in enumerate(self.postorder_nodes)}

        self.dp_table = np.full((num_nodes, n), np.inf, dtype=np.float32)
        self.backtrack = np.full((num_nodes, 2, n), -1, dtype=np.int32)

        for node in self.postorder_nodes:
            node_id = self.node_to_id[node]
            if node.is_terminal():
                self.dp_table[node_id, 0] = 0.0
                continue
            left, right = node.clades
            left_id = self.node_to_id[left]
            right_id = self.node_to_id[right]
            left_size = self.num_leaves_per_node[left]
            right_size = self.num_leaves_per_node[right]

            left_dp = self.dp_table[left_id]
            right_dp = self.dp_table[right_id]
            left_branch_length = left.branch_length or 0
            right_branch_length = right.branch_length or 0

            cost_one_cluster = (
                left_dp[0]
                + right_dp[0]
                + left_size * left_branch_length
                + right_size * right_branch_length
            )

            if self.use_branch_support_values:
                branch_support = (node.confidence / 100) if node.confidence else 1
                cost_one_cluster /= branch_support

            self.dp_table[node_id, 0] = cost_one_cluster
            self.backtrack[node_id, 0, 0] = 0
            self.backtrack[node_id, 1, 0] = 0

            total_terminals = left_size + right_size
            limit = min(n, total_terminals)

            for cluster_size in range(1, limit):
                left_indices = np.arange(cluster_size)
                right_indices = cluster_size - left_indices - 1
                left_dp_values = left_dp[left_indices]
                right_dp_values = right_dp[right_indices]
                sub_scores = left_dp_values + right_dp_values
                min_index = np.argmin(sub_scores)
                min_score = sub_scores[min_index]
                # if not np.isinf(min_score):
                #     back[:, cluster_size] = (left_indices[min_index], right_indices[min_index])

                self.dp_table[node_id, cluster_size] = min_score
                self.backtrack[node_id, 0, cluster_size] = left_indices[min_index]
                self.backtrack[node_id, 1, cluster_size] = right_indices[min_index]

    def backtrack_dp_assignments(
        self, k: Optional[int] = None, verbose: bool = False
    ) -> Dict[Any, int]:
        if k is None:
            raise ValueError("value of k is missing.")
        cluster_index = k - 1
        active_tree = self._no_outgroup_tree if self.outgroup else self.tree
        root = active_tree.root
        root_id = self.node_to_id[root]

        clusters = {}
        current_cluster_id = 0

        stack = [(root_id, cluster_index)]

        while stack:
            node_id, c_index = stack.pop()
            node = self.postorder_nodes[node_id]

            if verbose:
                print(f"Visiting node {node.name} with c_index={c_index}")
            if c_index == 0:
                for t in self.name_leaves_per_node[node]:
                    clusters[t] = current_cluster_id
                current_cluster_id += 1
            else:
                left_k = self.backtrack[node_id, 0, c_index]
                right_k = self.backtrack[node_id, 1, c_index]

                right = node.clades[1]
                right_id = self.node_to_id[right]
                stack.append((right_id, right_k))

                left = node.clades[0]
                left_id = self.node_to_id[left]
                stack.append((left_id, left_k))

        if current_cluster_id != k:
            raise ValueError(
                f"Number of clusters found: {current_cluster_id}, expected: {k}"
            )

        return clusters

    def extract_clusters(self, verbose: bool) -> None:
        """Define clusters using backtrack_dp_assignments."""
        self.clusters = [
            self.backtrack_dp_assignments(i, verbose=verbose)
            for i in range(1, self.max_k + 1)
        ]

    def execute_clustering(
        self,
        verbose: bool = False,
        num_peaks: int = 3,
        should_plot_scores: bool = False,
        resolution_on: bool = False,
        num_bins: int = 3,
        compute_all_clusters: bool = False,
    ) -> None:
        """Run the clustering algorithm."""
        if self.max_k is None and self.k:
            self.clusters = self.backtrack_dp_assignments(k=self.k, verbose=verbose)
        elif self.k is None:
            self._execute_clustering_for_all_k(
                verbose,
                num_peaks,
                should_plot_scores,
                resolution_on,
                num_bins,
                compute_all_clusters,
            )
        else:
            raise ValueError("Invalid cluster number. Please choose either k or max_k")

    def _execute_clustering_for_all_k(
        self,
        verbose: bool,
        num_peaks: int,
        should_plot_scores: bool,
        resolution_on: bool = False,
        num_bins: int = 3,
        compute_all_clusters: bool = False,
    ) -> None:
        """Helper method to handle the case when k is None."""
        if self.max_k is None:
            self.max_k = ceil(self.num_terminals * self.max_k_limit)
        if not isinstance(self.max_k, int) or self.max_k <= 0:
            raise ValueError("Invalid max_k value. It should be a positive integer.")

        self.calculate_scores(plot=should_plot_scores)
        self.find_score_peaks(
            global_peaks=num_peaks,
            plot=should_plot_scores,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )
        self.clusters = [None] * self.max_k

        if compute_all_clusters:
            for k_val in range(1, self.max_k + 1):
                c = self.backtrack_dp_assignments(k=k_val, verbose=verbose)
                self.clusters[k_val - 1] = c
        else:
            for peak_k in self.peaks_by_rank:
                if 1 <= peak_k <= self.max_k:
                    c = self.backtrack_dp_assignments(k=peak_k, verbose=verbose)
                    self.clusters[peak_k - 1] = c
                else:
                    pass

    def compute_cluster_score(
        self,
        clusters: Optional[Dict[Any, Any]] = None,
        k: Optional[int] = None,
        output_all: bool = False,
    ) -> Tuple[float, float, float]:
        """
        Calculate the cluster score in one of two ways:
        1) If 'clusters' is provided (dict of {terminal->cluster_id}):
            - We determine num_clusters by counting unique cluster IDs.
            - Then read dp_row[num_clusters-1].
        2) If 'k' is provided (int):
            - We skip the dictionary. We just read dp_row[k-1].
        3) If both are None, raise an error.

        Returns: (beta, beta_ratios, score)

        Where 'beta' = dp_row[num_clusters-1],
            'beta_ratios' = (beta_1 - beta) / beta,
            'score' = beta_ratios * ((num_terminals - num_clusters)/num_clusters)
        """

        if not self.max_k or self.max_k <= 0:
            raise ValueError(
                "max_k must be set and positive to compute cluster scores."
            )

        active_tree = self._no_outgroup_tree if self.outgroup else self.tree
        root = active_tree.root
        root_id = self.node_to_id[root]
        dp_row = self.dp_table[root_id, :]

        if dp_row is None:
            raise ValueError("Root not found in dp_table")

        self.beta_1 = dp_row[0]
        num_terminals = self.num_terminals

        if clusters is not None:
            clust_inv = {}
            for leaf, cluster_id in clusters.items():
                clust_inv.setdefault(cluster_id, set()).add(leaf)

            num_clusters = len(clust_inv)
            if num_clusters < 1 or num_clusters > self.max_k:
                return (float("inf"), float("inf"), float("inf"))

            beta = dp_row[num_clusters - 1]

        elif k is not None:
            if k < 1 or k > self.max_k:
                return (float("inf"), float("inf"), float("inf"))
            num_clusters = k
            beta = dp_row[k - 1]

        else:
            raise ValueError(
                "Either 'clusters' or 'k' must be provided to compute the score."
            )

        if np.isinf(beta):
            return (beta, float("inf"), 0.0)

        beta_ratios = (self.beta_1 - beta) / beta if beta != 0 else float("inf")
        if num_clusters == 0:
            norm_ratios = float("inf")
        else:
            norm_ratios = (num_terminals - num_clusters) / float(num_clusters)

        if norm_ratios in [float("inf"), float("nan")] or np.isinf(beta_ratios):
            score = float("inf")
        else:
            score = beta_ratios * norm_ratios

        return (beta, beta_ratios, score)

    def calculate_scores(self, plot: bool = False, output_all: bool = False) -> None:
        """Calculate scores for all clusters."""
        if self.k is not None:
            if isinstance(self.clusters, dict):
                beta_val, norm_val, sc = self.compute_cluster_score(
                    clusters=self.clusters, output_all=output_all
                )
            else:
                beta_val, norm_val, sc = self.compute_cluster_score(
                    k=self.k, output_all=output_all
                )
            beta_values = np.array([beta_val])
            den_list = np.array([norm_val])
            scores = np.array([sc])

        elif isinstance(self.clusters, list) and len(self.clusters) > 0:
            results = []
            for cluster_dict in self.clusters:
                if cluster_dict is None:
                    results.append((np.inf, np.inf, 0.0))
                else:
                    results.append(
                        self.compute_cluster_score(
                            clusters=cluster_dict, output_all=output_all
                        )
                    )
            beta_values, den_list, scores = zip(*results)
            beta_values = np.array(beta_values)
            den_list = np.array(den_list)
            scores = np.array(scores)

        else:
            if not self.max_k or self.max_k <= 0:
                raise ValueError(
                    "max_k must be set and positive to compute DP-based scores."
                )

            results = []
            for k_val in range(1, self.max_k + 1):
                b, r, s = self.compute_cluster_score(k=k_val, output_all=output_all)
                results.append((b, r, s))

            beta_values, den_list, scores = zip(*results)
            beta_values = np.array(beta_values)
            den_list = np.array(den_list)
            scores = np.array(scores)

        scores[scores < 0] = 0

        beta_values[beta_values < 0] = 0
        beta_values = np.nan_to_num(beta_values, nan=0.0, posinf=0.0, neginf=0.0)

        elbow_scores = []
        n = len(beta_values)
        for i in range(n - 1):
            if beta_values[i] - beta_values[i + 1] != 0:
                elbow_score = (
                    (beta_values[i - 1] - beta_values[i])
                    / (beta_values[i] - beta_values[i + 1])
                    if i > 0
                    else 0
                )
            else:
                elbow_score = 0
            elbow_scores.append(elbow_score)
        elbow_scores.append(0)

        min_invalid_index = None
        invalid_elbow = np.where(np.isnan(elbow_scores) | np.isinf(elbow_scores))[0]
        invalid_scores = np.where(np.isnan(scores) | np.isinf(scores))[0]

        if len(invalid_elbow) > 0:
            min_invalid_index = invalid_elbow[0]
        if len(invalid_scores) > 0:
            if min_invalid_index is not None:
                min_invalid_index = min(min_invalid_index, invalid_scores[0])
            else:
                min_invalid_index = invalid_scores[0]

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

        normalized_elbow = min_max_normalize(elbow_scores)
        normalized_scores = min_max_normalize(scores)
        combined_scores = normalized_elbow * normalized_scores
        combined_scores = np.nan_to_num(
            combined_scores, nan=0.0, posinf=0.0, neginf=0.0
        )

        self.scores = combined_scores
        self.beta_values = beta_values
        self.norm_ratios = den_list

        if plot:
            pass

    def define_bins(
        self, num_bins: int = 3, log_base: int = None
    ) -> List[Tuple[int, int]]:

        if not self.num_terminals:
            raise ValueError("Number of terminals must be set to define bins.")

        k_min = 1
        k_max = self.num_terminals

        if log_base is None:
            log_base = (k_max / k_min) ** (1 / num_bins)

        bin_edges = k_min * (log_base ** np.arange(num_bins + 1))
        bin_edges = np.ceil(bin_edges).astype(int)  # + 1

        if bin_edges[-1] < k_max:
            bin_edges[-1] = k_max

        bin_ranges = []
        for i in range(len(bin_edges) - 1):
            start = bin_edges[i]
            end = bin_edges[i + 1]
            bin_ranges.append((start, end))

        return bin_ranges

    def plot_scores(
        self,
        scores_subset: Optional[np.ndarray] = None,
        k_start: int = 0,
        k_end: Optional[int] = None,
        peaks: Optional[List[int]] = None,
        resolution_on: bool = False,
        num_bins: int = 3,
        fig_width: int = 18,
        fig_height: int = 6,
        log_scale_x: bool = True,
        log_scale_y: bool = True,
        log_base: Optional[float] = None,
    ) -> plt.Figure:

        if scores_subset is None:
            scores_subset = self.scores
        if len(scores_subset) == 0:
            raise ValueError("Scores are empty, please calculate scores first.")

        min_val, max_val = np.min(scores_subset), np.max(scores_subset)
        if max_val - min_val == 0:
            scores_subset = np.zeros_like(scores_subset)
        else:
            scores_subset = (scores_subset - min_val) / (max_val - min_val)

        if k_end is None:
            k_end = k_start + len(scores_subset)

        scores_slice = scores_subset[k_start:k_end]
        x_indices = np.arange(k_start + 1, k_end + 1)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.plot(x_indices, scores_slice, "o-", label="Scores", markersize=2)

        if peaks is not None and len(peaks) > 0:
            peak_x = [p for p in peaks]
            peak_scores = []
            for p in peaks:
                idx_in_slice = p - (k_start + 1)
                if 0 <= idx_in_slice < len(scores_slice):
                    peak_scores.append(scores_slice[idx_in_slice])
                else:
                    peak_scores.append(np.nan)  # out of slice range

            ax.plot(peak_x, peak_scores, "rx", label="Peaks", markersize=8)
            for p, sc in zip(peaks, peak_scores):
                ax.text(p, sc + 0.02, str(p), fontsize=10, color="red", ha="center")

        if resolution_on:
            xmin = 1
            xmax = len(scores_slice)

            if log_base is None:
                log_base = (xmax / xmin) ** (1 / num_bins)

            bin_edges = xmin * (log_base ** np.arange(num_bins + 1))
            bin_edges = np.ceil(bin_edges).astype(int)

            bin_edges[-1] = max(bin_edges[-1], xmax)

            if num_bins == 3:
                colors = ["#5DADE2", "#58D68D", "#EC7063"]
                darker_colors = ["#2874A6", "#239B56", "#C0392B"]
                labels = ["Low", "Intermediate", "High"]
            elif num_bins == 2:
                colors = ["#5DADE2", "#EC7063"]
                darker_colors = ["#2874A6", "#C0392B"]
                labels = ["Low", "High"]
            else:
                light_colors = [
                    c
                    for c in mcolors.CSS4_COLORS.values()
                    if c.lower() not in ["white", "#ffffff"]
                ]
                colors = random.sample(light_colors, num_bins)
                labels = [f"Bin {i+1}" for i in range(num_bins)]
                darker_colors = colors

            for i in range(len(bin_edges) - 1):
                left, right = bin_edges[i], bin_edges[i + 1]
                plt.axvspan(
                    left,
                    right,
                    color=colors[i],
                    alpha=0.2,
                    label=(
                        f"Bin {i+1}: {left}-{right}" if num_bins in [2, 3] else None
                    ),
                )

            for edge in bin_edges:
                plt.axvline(x=edge, color="grey", linestyle="--", linewidth=1)

            if num_bins in [2, 3]:
                for i, lab in enumerate(labels):
                    left, right = bin_edges[i], bin_edges[i + 1]
                    mid = left + (right - left) / 2.0
                    plt.text(
                        mid,
                        0.5,
                        lab,
                        ha="center",
                        va="center",
                        fontsize=12,
                        color=darker_colors[i],
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                    )

        plt.xlabel("No. of Clusters (log)", fontsize=14)
        plt.ylabel("Scores (log)", fontsize=14)

        if log_scale_x:
            plt.xscale("log")
            ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
            ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs="auto", numticks=10))
            ax.xaxis.set_major_formatter(LogFormatter(base=10.0, labelOnlyBase=False))

        if log_scale_y:
            plt.yscale("log")

        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs="auto", numticks=10))
        ax.xaxis.set_major_formatter(LogFormatter(base=10.0, labelOnlyBase=False))
        plt.title("Scores", fontsize=16)
        plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tick_params(axis="both", which="major", labelsize=12, rotation=45)
        plt.tick_params(axis="both", which="minor", labelsize=12)
        plt.tight_layout()
        plt.show()
        return fig

    def find_score_peaks(
        self,
        scores: Optional[np.ndarray] = None,
        global_peaks: int = 3,
        peaks_per_bin: int = 1,
        resolution_on: bool = False,
        num_bins: int = 3,
        min_k: int = 3,
        k_start: Optional[int] = None,
        k_end: Optional[int] = None,
        plot: bool = True,
    ) -> List[int]:
        if scores is None:
            scores = self.scores
        if len(scores) == 0:
            raise ValueError("Please calculate scores first")
            return []
        k_start = k_start if k_start is not None else 0
        k_end = k_end if k_end is not None else len(scores)
        if k_end > len(scores) or k_end < k_start:
            raise ValueError(f"Invalid k_end value. Must be <= {len(scores)}")

        scores_subset = scores[k_start:k_end]
        scores_subset = np.where(np.isinf(scores_subset), np.nan, scores_subset)

        if len(scores_subset) > 1:
            if scores_subset[-1] <= 0:
                log_scores = np.log(scores_subset[:-1] + 1e-10)
            else:
                log_scores = np.log(scores_subset + 1e-10)
        else:
            log_scores = np.log(scores_subset + 1e-10)

        peaks, props = find_peaks(log_scores, prominence=1e-10)
        prominences = props["prominences"]
        adjusted_prom = prominences / (peaks + k_start + 1)

        sorted_indices = np.argsort(adjusted_prom)[::-1]
        peaks_sorted = peaks[sorted_indices]
        prom_sorted = adjusted_prom[sorted_indices]

        peaks_converted = [p + k_start + 1 for p in peaks_sorted]

        filtered_peaks = []
        filtered_prom = []
        for pk, pr in zip(peaks_converted, prom_sorted):
            if pk >= min_k:
                filtered_peaks.append(pk)
                filtered_prom.append(pr)

        final_peaks = []
        if not resolution_on:
            self.resolution_info = None
            final_peaks = filtered_peaks[:global_peaks]
        else:
            bin_ranges = self.define_bins(num_bins=num_bins)
            self.resolution_info = {}
            peak_to_prom = {pk: pr for pk, pr in zip(filtered_peaks, filtered_prom)}

            for i, (start_kv, end_kv) in enumerate(bin_ranges, start=1):
                bin_label = f"Bin {i}: {start_kv}-{end_kv}"
                cands = [pk for pk in filtered_peaks if start_kv <= pk <= end_kv]
                cands.sort(key=lambda x: peak_to_prom[x], reverse=True)
                chosen = cands[:peaks_per_bin]
            self.resolution_info[bin_label] = chosen
            final_peaks.extend(chosen)

        self.peaks_by_rank = sorted(final_peaks)

        if plot:
            self.plot_of_scores = self.plot_scores(
                scores_subset=scores_subset,
                peaks=self.peaks_by_rank,
                k_start=k_start,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )

        return self.peaks_by_rank

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
        if self.k is not None:
            clusters_to_plot.append((self.k, self.clusters))
        elif n is not None:
            clusters_to_plot.append((n, self.clusters[n - 1]))
        else:
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
        output_all: bool = False,
    ) -> None:

        os.makedirs(results_dir, exist_ok=True)

        if hasattr(self, "plot_of_scores") and self.scores is not None:
            fig = self.plot_of_scores
            fig.savefig(os.path.join(results_dir, "scores.png"))

        # save_clusters(
        #     scores=self.scores,
        #     clusters=self.clusters,
        #     results_dir=results_dir,
        #     top_n=top_n,
        #     filename=filename,
        #     outlier=outlier,
        #     selected_peaks=self.peaks_by_rank,
        #     n=n,
        #     output_all=output_all,
        # )

        peaks_by_rank_file = os.path.join(results_dir, "peaks_by_rank.txt")
        with open(peaks_by_rank_file, "w") as f:
            for rank,kval in enumerate(self.peaks_by_rank, start=1):
                # i = 1
                f.write(f"Rank {rank}: {kval} clusters\n")
                # i+=1

        if getattr(self, "resolution_info", None) is not None:
            resolution_file = os.path.join(results_dir, "resolution_bins.txt")
            with open(resolution_file, "w") as f:
                for bin_label, pk_list in self.resolution_info.items():
                    f.write(f"{bin_label}: {pk_list}\n")
        print(f"Saved results in {results_dir}")
