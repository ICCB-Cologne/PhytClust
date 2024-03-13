import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import string
from phytclust.plotting import plot_cluster, plot_peaks
from phytclust.save import save_clusters
import time

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

    def __init__(self, tree, k=None, max_k=None, score_type=None):
        self.tree = tree
        self.validate_tree()  # to check if everything is in its place, no planted root, internal node names assigned
        self.k = k

        if k is None:
            self.max_k = len(self.tree.get_terminals()) if max_k is None else max_k
        else:
            self.max_k = None
        self.clusters = {}
        self.score_type = score_type
        self.dpTable = {}
        self.backtrack = {}
        self.scores = []
        self.top_peaks = []
        self.total_runtime = 0.0
        self.runtimes = {}
        self.tree_clust_dpTable()

    def validate_tree(self):
        invalid_nodes = [
            node for node in self.tree.get_nonterminals() if len(node.clades) != 2
        ]

        if invalid_nodes:
            print("Nodes with > 2 children:")
            for node in invalid_nodes:
                print(f"Node: {node.name}, Children: {len(node.clades)}")

            raise AssertionError("All internal nodes should have 2 children")

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

    def tree_clust_dpTable(self):
        start_time = time.time()
        n = self.k or self.max_k

        for clade in self.tree.find_clades(order="postorder"):

            scores = np.full(n, np.inf)
            back = np.full((2, n), np.inf)

            if clade.is_terminal():
                scores[0] = 0
            else:
                left, right = clade.clades
                left_size = left.count_terminals()
                right_size = right.count_terminals()
                total_terminals = left_size + right_size
                limit = min(n, total_terminals)

                scores[0] = (
                    self.dpTable[left][0]
                    + self.dpTable[right][0]
                    + left.branch_length * left_size
                    + right.branch_length * right_size
                )
                back[:, 0] = (0, 0)  # (left, right)

                for i in range(1, limit):
                    sub_scores = [
                        self.dpTable[left][j] + self.dpTable[right][i - j - 1]
                        for j in range(i)
                        if 0 <= j < len(self.dpTable[left])
                        and 0 <= i - j - 1 < len(self.dpTable[right])
                    ]
                    min_score = np.min(sub_scores)
                    scores[i] = min_score

                    if not np.isinf(min_score):
                        min_index = np.argmin(sub_scores)
                        back[:, i] = (min_index, i - min_index - 1)

            self.dpTable[clade], self.backtrack[clade] = scores, back
        end_time = time.time()
        self._track_runtime("dpTable", start_time, end_time)

    def tree_clust_backtrack(self, k=None, verbose=False):
        if k is None:
            raise ValueError("value of k is missing.")

        ptr = k - 1
        stack = [(self.tree.clade, 0, ptr)]
        clusters = {}
        current_cluster = 0

        while stack:
            node, counter, ptr = stack.pop()

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
                self.max_k = self.tree.get_terminals()
                self.def_clusters(verbose)

            elif isinstance(self.max_k, int) and self.max_k > 0:
                self.def_clusters(verbose)

            else:
                raise ValueError("invalid cluster number")
        else:
            raise ValueError("invalid cluster number. Please choose either k or max_k")

        end_time = time.time()
        self._track_runtime("run", start_time, end_time)

    def calinski_harabasz_index(self, clusters):
        clust_inv = {}
        for k, v in clusters.items():
            clust_inv[v] = clust_inv.get(v, []) + [k]
        num_clusters = len(clust_inv)
        num_terminals = self.tree.count_terminals()

        if num_clusters <= 1:
            return 0  # Not enough clusters to calculate the index

        lhs = np.sum(
            [
                self.tree.distance(self.tree.common_ancestor(taxa), self.tree.root)
                for k, taxa in clust_inv.items()
            ]
        )
        lhs /= num_clusters - 1

        rhs = np.sum(
            [
                np.sum(
                    [
                        self.tree.distance(self.tree.common_ancestor(taxa), t)
                        for t in taxa
                    ]
                )
                for k, taxa in clust_inv.items()
            ]
        )
        if num_terminals - num_clusters > 0:
            rhs /= num_terminals - num_clusters
        else:
            return 0  # Avoid division by zero

        chi = lhs / rhs if rhs != 0 else 0  # Avoid division by zero
        return chi

    def cluster_score(self, clusters=None, score_type=None, output_all=False):
        clusters = clusters or self.clusters
        score_type = self.score_type if self.score_type is not None else "CH_score"

        if not clusters:
            raise ValueError("Clusters not found")
        if score_type not in ["CH_score", "Alt", "New_Score"]:
            raise ValueError(f"Invalid score type: {score_type}")

        clust_inv = {}
        for k, v in clusters.items():
            clust_inv.setdefault(v, []).append(k)

        num_clusters = len(clust_inv)
        num_terminals = self.tree.count_terminals()

        if num_clusters == num_terminals or num_clusters == 1:
            return float("inf")

        score = 0
        if score_type == "Alt":
            score = self.calinski_harabasz_index(clusters)
        elif score_type == "CH_score":
            root_clade = self.tree.root
            dp_row = self.dpTable.get(root_clade, None)
            if dp_row is None:
                raise ValueError("Root not found in dpTable")

            beta_1 = dp_row[0]
            beta = dp_row[num_clusters - 1]

            num = (beta_1 - beta) / (num_clusters - 1)
            den = beta / (num_terminals - num_clusters)
            score = num / den

        return (score, num, den) if output_all else score

    def score_calc(self, plot=True):
        start_time = time.time()
        self.scores = []
        if self.k is None:
            for cluster in self.clusters:
                self.scores.append(self.cluster_score(clusters=cluster))
        else:
            score = self.cluster_score(clusters=self.clusters)
            print(f"Single Cluster: {self.clusters}, Score: {score}")
            self.scores.append(score)

        self.scores = np.array(self.scores)

        if plot and self.k is None:
            plt.plot(range(1, len(self.scores) + 1), self.scores, "o-")
            plt.xlabel("No. of Clusters")
            plt.ylabel("Scores")
            plt.title("Scores for Different No. of Clusters")
            plt.show()

        end_time = time.time()
        self._track_runtime("score_calc", start_time, end_time)

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
        Finds and plots the top n peaks in the scores within a specified range of k.

        Parameters:
        - n: int, the number of top peaks to find.
        - plot: bool, whether to plot the score graph with the peaks highlighted.
        - method: str, 'default' or 'alt' for different peak finding methods.
        - prominence: tuple or None, the prominence parameter for peak finding, applicable in 'alt' method.
        - k_start: int or None, the starting index for the range of k values to consider.
        - k_end: int or None, the ending index for the range of k values to consider.
        """
        start_time = time.time()
        if self.scores is None or len(self.scores) == 0:
            print("Please calculate scores first")
            return []

        if k_end is not None and (k_end > len(self.scores) or k_end < k_start):
            raise ValueError("Invalid k_end value.")

        k_start = k_start or 0
        k_end = k_end or len(self.scores)

        scores_subset = self.scores[k_start:k_end]

        if method == "default":
            peaks, _ = find_peaks(scores_subset, prominence=(None, None))
        elif method == "alt":
            peaks, _ = find_peaks(scores_subset, prominence=prominence)

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
        outlier=True,
        save=False,
        filename=None,
        hide_internal_nodes=True,
        width_scale=2,
        height_scale=0.1,
        label_func=lambda x: None,
        show_branch_lengths=False,
        marker_size=50,
        **kwargs,
    ):
        if self.scores is None:
            print("Please calculate scores first")
            return

        # Check if k is specified and directly use self.clusters for plotting
        if self.k is not None:
            cluster = self.clusters  # Directly use the mapping
            cluster_number = self.k
            plot_cluster(
                cluster,
                cluster_number,
                self.tree,
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
                **kwargs,
            )
        elif n is not None:
            cluster = self.clusters[n - 1]
            cluster_number = n
            plot_cluster(
                cluster,
                cluster_number,
                self.tree,
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
                **kwargs,
            )
        else:
            # Existing logic for when k is not specified
            clusters_to_plot = []
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

            # Plotting for top peaks or scores-based clusters
            for i, (peak, cluster) in enumerate(clusters_to_plot):
                plot_cluster(
                    cluster,
                    peak + 1,
                    self.tree,
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
                    **kwargs,
                )

    def save(
        self, results_dir, top_n=1, filename="phyclust_results.csv", outliers=True
    ):
        save_clusters(
            self.scores,
            self.clusters,
            results_dir=results_dir,
            top_n=top_n,
            filename=filename,
            outliers=outliers,
        )

    def _track_runtime(self, function_name, start, end):
        runtime = end - start
        self.total_runtime += runtime
        if function_name in self.runtimes:
            self.runtimes[function_name].append(runtime)
        else:
            self.runtimes[function_name] = [runtime]
