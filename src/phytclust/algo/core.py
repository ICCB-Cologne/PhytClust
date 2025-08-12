import logging
from dataclasses import dataclass, field
from math import ceil
from typing import Any, Optional, Dict, List

import numpy as np

from ..algo.dp import (
    validate_args,
    prepare_tree,
    compute_dp_table,
    backtrack,
    cluster_map,
)
from ..algo.scoring import calculate_scores, find_score_peaks
from ..viz.cluster_plot import plot_clusters
from ..io.save_results import save_clusters

logger = logging.getLogger("phytclust")

IntMap = Dict[Any, int]


@dataclass
class PhytClust:
    """
    Dynamic-programming phylogenetic clustering.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The input tree (pass a Newick string / file if you prefer; it is
        parsed upstream).
    outgroup : str | None, default=None
        Taxon to exclude from all clusters.
    min_cluster_size : int, default=1
        Ignore clusters smaller than this.

    After construction
    ------------------
    • `pc.get_clusters(k)`         – exact k-way partition
    • `pc.best_global(top_n=3)`    – CalBow peak search
    • `pc.best_by_resolution()`    – one peak per log-bin

    Results in **`pc.clusters`**
    dict keyed by *k*, each value a `{terminal → cluster_id}` mapping.
    """

    tree: Any
    # _: KW_ONLY

    outgroup: Optional[str] = None
    min_cluster_size: int = 1
    k: Optional[int] = None
    max_k: Optional[int] = None
    max_k_limit: float = 0.9
    should_plot_scores: bool = True
    resolution_on: bool = False
    num_peaks: int = 3
    num_bins: int = 3

    # tunables
    use_branch_support: bool = False
    compute_all_clusters: bool = False
    drop_outliers: bool = False

    # internal caches/state
    _dp_ready: bool = field(init=False, default=False, repr=False)
    name_leaves_per_node: Dict[Any, List[Any]] = field(
        init=False, repr=False, default_factory=dict
    )
    num_leaves_per_node: Dict[Any, int] = field(
        init=False, repr=False, default_factory=dict
    )
    backptr: Dict[int, np.ndarray] = field(init=False, repr=False, default_factory=dict)

    clusters: Optional[Dict[int, IntMap]] = field(init=False, default=None)
    scores: Optional[np.ndarray] = field(init=False, default=None)
    peaks_by_rank: Optional[List[int]] = field(init=False, default=None)

    def __post_init__(self) -> None:
        if (
            getattr(self, "_dp_ready", False)
            and getattr(self, "_prepared_tree", None) is self.tree
        ):
            return
        validate_args(self)
        prepare_tree(self)
        compute_dp_table(self)
        setattr(self.tree, "_phytclust_done", True)
        self._dp_ready = True

    # explicit k
    def get_clusters(self, k: int, *, verbose: bool = False) -> Dict[Any, int]:
        if not self._dp_ready:
            raise RuntimeError("DP table not built. Check for error")
        if k is None:
            raise ValueError("Please provide k")
        self.k = k
        cmap = backtrack(self, k, verbose=verbose)
        if self.clusters is None:
            self.clusters = {}
        self.clusters[k] = cmap
        return cmap

    # optimise globally
    def best_global(
        self,
        *,
        top_n: int = 1,
        max_k: Optional[int] = None,
        max_k_limit: Optional[float] = None,
        plot_scores: bool = True,
        compute_all_clusters: bool = False,
    ) -> List[Dict[Any, int]]:
        self.max_k_limit = max_k_limit if max_k_limit else 0.9
        self.max_k = max_k or max(2, ceil(self.num_terminals * self.max_k_limit))
        calculate_scores(self, plot=plot_scores)
        if self.scores is None or len(self.scores) == 0:
            try:
                cmap = self.get_clusters(2)
                self.clusters = {2: cmap}
                self.peaks_by_rank = [2]
                self.k = None
                return [cmap]
            except Exception:
                self.clusters = {}
                self.peaks_by_rank = []
                self.k = None
                return []

        score_len = min(self.max_k, len(self.scores))
        find_score_peaks(
            self,
            global_peaks=top_n,
            resolution_on=False,
            k_start=2,
            k_end=score_len,
            plot=plot_scores,
        )
        self.clusters = {}
        if compute_all_clusters:
            for k_val in range(1, self.max_k + 1):
                self.get_clusters(k_val)
        else:
            for k_val in self.peaks_by_rank or []:
                self.get_clusters(k_val)
        self.k = None
        return [
            self.clusters[kv]
            for kv in (self.peaks_by_rank or [])
            if kv in self.clusters
        ]

    def best_by_resolution(
        self,
        *,
        num_bins: int = 3,
        max_k: Optional[int] = None,
        plot_scores: bool = True,
    ) -> List[Dict[Any, int]]:
        self.max_k = max_k or ceil(self.num_terminals * self.max_k_limit)
        calculate_scores(self, plot=plot_scores)
        score_len = min(self.max_k, len(self.scores)) - 1
        find_score_peaks(
            self,
            resolution_on=True,
            num_bins=num_bins,
            peaks_per_bin=1,
            k_start=2,
            k_end=score_len,
            plot=plot_scores,
        )
        self.clusters = {}
        for k_val in self.peaks_by_rank or []:
            self.get_clusters(k_val)
        self.k = None
        return [
            self.clusters[kv]
            for kv in (self.peaks_by_rank or [])
            if kv in self.clusters
        ]

    # thin wrappers to viz / io
    def plot(self, results_dir: Optional[str] = None, **kwargs) -> None:
        plot_clusters(self, results_dir=results_dir, **kwargs)

    def save(
        self,
        results_dir: str,
        top_n: int = 1,
        filename: str = "phyclust_results.csv",
        outlier: bool = True,
        n: Optional[int] = None,
        output_all: bool = False,
    ) -> None:
        save_clusters(
            self,
            results_dir=results_dir,
            top_n=top_n,
            filename=filename,
            outlier=outlier,
            n=n,
            output_all=output_all,
        )
