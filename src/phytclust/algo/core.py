import logging
from dataclasses import dataclass, field
from math import ceil
from typing import Any, Optional, Dict, List

import numpy as np
from pathlib import Path
from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from ..algo.dp import (
    validate_args,
    prepare_tree,
    compute_dp_table,
    backtrack,
)
from ..algo.scoring import calculate_scores, find_score_peaks
from ..viz.cluster_plot import plot_clusters
from ..io.save_results import save_clusters

logger = logging.getLogger("phytclust")

IntMap = Dict[Any, int]


def _coerce_to_tree(obj: Any) -> Tree:
    """
    Accepts:
      - Bio.Phylo.BaseTree.Tree  → returned as-is
      - pathlib.Path             → read as Newick
      - str:
          * if it looks like a file path and exists → read as Newick file
          * otherwise → treat as a Newick string

    Raises TypeError if the object cannot be interpreted as a tree.
    """
    if isinstance(obj, Tree):
        return obj

    if isinstance(obj, Path):
        return Phylo.read(str(obj), "newick")

    if isinstance(obj, str):
        candidate = Path(obj)
        if candidate.exists() and not any(ch in obj for ch in "() ;"):
            return Phylo.read(str(candidate), "newick")

        handle = StringIO(obj)
        return Phylo.read(handle, "newick")

    raise TypeError(
        f"Unsupported tree input type: {type(obj)!r}. "
        "Expected a Bio.Phylo Tree, a Newick string, or a path to a Newick file."
    )


@dataclass
class PhytClust:
    """
    Dynamic-programming phylogenetic clustering.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The input tree (you can pass a Newick string/file upstream).
    outgroup : str | None, default=None
        Taxon to exclude from all clusters (treated as outgroup).
    min_cluster_size : int, default=1
        Hard constraint: final clusters smaller than this are disallowed.
    k : int | None, default=None
        Fixed number of clusters (only used when you call `run(k=...)`
        or `get_clusters(k)` explicitly).
    max_k : int | None, default=None
        Upper bound on k for scoring / peak search. If None, derived
        from `max_k_limit * num_terminals`.
    max_k_limit : float, default=0.9
        When `max_k` is not set, `max_k` = ceil(max_k_limit * num_terminals).
    resolution_on : bool, default=False
        Internal flag for multi-resolution peak search.
    num_peaks : int, default=3
        Number of global peaks to report in best_global.
    num_bins : int, default=3
        Number of log-resolution bins for best_by_resolution.

    Support / branch-length tuning
    ------------------------------
    use_branch_support : bool, default=False
        If True, branch supports are incorporated into effective branch
        lengths and internal split penalties.
    min_support : float, default=0.05
        Minimal support used when normalizing (avoid division by 0).
    support_weight : float, default=1.0
        Weight of the support-derived penalty in `_eff_length`.

    Outlier handling (soft constraint)
    ----------------------------------
    outlier_size_threshold : int | None, default=None
        If set, clusters with size < threshold incur `outlier_penalty`
        added to their cost (soft penalty, not a hard constraint).
    outlier_penalty : float, default=0.0
        Additive penalty for clusters smaller than outlier_size_threshold.

    Other flags
    -----------
    compute_all_clusters : bool, default=False
        If True in best_global, compute and cache all k clusterings up to max_k.
    drop_outliers : bool, default=False
        Reserved flag (used downstream in plotting/saving).
    """

    tree: Any
    outgroup: Optional[str] = None
    min_cluster_size: int = 1
    k: Optional[int] = None
    max_k: Optional[int] = None
    max_k_limit: float = 0.9
    resolution_on: bool = False
    num_peaks: int = 3
    num_bins: int = 3

    # tunables
    use_branch_support: bool = False
    min_support: float = 0.05
    support_weight: float = 1.0

    # outlier penalties
    outlier_size_threshold: Optional[int] = None
    outlier_penalty: float = 0.0

    compute_all_clusters: bool = False
    drop_outliers: bool = False

    def __post_init__(self) -> None:
        self.tree = _coerce_to_tree(self.tree)

        self.name_leaves_per_node = {}
        self.num_leaves_per_node = {}
        self.backptr = {}
        self.dp_table = None
        self.postorder_nodes = None
        self.node_to_id = None
        self.num_terminals = 0
        self._tree_wo_outgroup = None

        self.scores = None
        self.peaks_by_rank = None

        self._dp_ready = False
        self._tree_hash: Optional[int] = None
        self.clusters: Dict[int, IntMap] = {}
        self._last_result: Optional[Dict[str, Any]] = None

    def _hash_tree(self) -> int:
        """tree fingerprint, used to detect modifications."""
        try:
            return hash(self.tree.format("newick"))
        except Exception:
            return hash(repr(self.tree))

    def _ensure_dp(self) -> None:
        current = self._hash_tree()

        if self._dp_ready and self._tree_hash == current:
            return

        validate_args(self)
        prepare_tree(self)
        compute_dp_table(self)

        self._dp_ready = True
        self._tree_hash = current

        self.clusters = {}
        self.scores = None
        self.peaks_by_rank = None

        if self.max_k is None or self.max_k < 1:
            self.max_k = max(2, ceil(self.num_terminals * self.max_k_limit))

    # explicit k
    def get_clusters(self, k: int, *, verbose: bool = False) -> Dict[Any, int]:
        """
        Exact k-cluster partition
        """
        if k is None:
            raise ValueError("Please provide k")
        if k < 1:
            raise ValueError("k must be ≥ 1")

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
        alpha: float = 0.7,
    ) -> List[Dict[Any, int]]:
        """
        cluster-validity index-based global peak search.

        Returns a list of cluster maps in peak-rank order.
        """
        self._ensure_dp()

        self.max_k_limit = max_k_limit if max_k_limit is not None else 0.9
        self.max_k = max_k or max(2, ceil(self.num_terminals * self.max_k_limit))

        calculate_scores(self, plot=plot_scores)
        logger.debug(
            f"score vector length = {0 if self.scores is None else len(self.scores)}"
        )
        logger.debug(f"scores = {self.scores}")

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
            alpha=alpha,
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

    # choose best solutions for different clade levels
    def best_by_resolution(
        self,
        *,
        num_bins: int = 3,
        max_k: Optional[int] = None,
        plot_scores: bool = True,
        alpha: float = 0.7,
    ) -> List[Dict[Any, int]]:
        """
        Multi-resolution mode: one peak per log-bin.

        Returns a list of cluster maps, one per selected k.
        """
        self._ensure_dp()

        self.max_k = max_k or ceil(self.num_terminals * self.max_k_limit)
        calculate_scores(self, plot=plot_scores)

        score_k_count = min(self.max_k, len(self.scores))

        if score_k_count < 50:
            logger.info(
                "Tree too small for resolution mode, finding globally optimal solution instead."
            )
            return self.best_global(
                top_n=min(self.num_peaks, score_k_count - 1),
                max_k=self.max_k,
                max_k_limit=self.max_k_limit,
                plot_scores=plot_scores,
                compute_all_clusters=False,
                alpha=alpha,
            )

        score_len = score_k_count - 1

        find_score_peaks(
            self,
            resolution_on=True,
            num_bins=num_bins,
            peaks_per_bin=1,
            k_start=2,
            k_end=score_len,
            plot=plot_scores,
            alpha=alpha,
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

    def run(
        self,
        *,
        k: Optional[int] = None,
        top_n: int = 1,
        by_resolution: bool = False,
        num_bins: Optional[int] = None,
        max_k: Optional[int] = None,
        max_k_limit: Optional[float] = None,
        plot_scores: bool = True,
        alpha: float = 0.7,
    ) -> Dict[str, Any]:
        """
        Unified high-level entry point.

        Modes
        -----
        1. Exact k:
            pc.run(k=5)

        2. Global peaks (CalBow):
            pc.run(top_n=3)

        3. Multi-resolution peaks (one per log-bin):
            pc.run(by_resolution=True, num_bins=3)
        """
        # validation
        k_val = k if k is not None else self.k
        if k_val is not None:
            if k_val < 1:
                raise ValueError("k must be ≥ 1.")
            if by_resolution:
                raise ValueError("Cannot combine `k` with `by_resolution=True`.")
            if top_n != 1:
                raise ValueError("`top_n` is meaningless when `k` is given.")
            if num_bins is not None:
                logger.warning("`num_bins` ignored when `k` is specified.")

            # exact-k mode
            self._ensure_dp()
            cmap = self.get_clusters(k_val)
            result = {
                "mode": "k",
                "k": k_val,
                "clusters": cmap,
                "scores": None,
                "peaks": [k_val],
            }
            self._last_result = result
            return result

        if top_n < 1:
            raise ValueError("`top_n` must be ≥ 1.")

        if by_resolution:
            if num_bins is None:
                num_bins = self.num_bins
            if max_k_limit is not None:
                self.max_k_limit = max_k_limit

            clusters = self.best_by_resolution(
                num_bins=num_bins,
                max_k=max_k,
                plot_scores=plot_scores,
                alpha=alpha,
            )
            result = {
                "mode": "resolution",
                "ks": list(self.peaks_by_rank or []),
                "clusters": clusters,
                "scores": None if self.scores is None else self.scores.copy(),
                "peaks": list(self.peaks_by_rank or []),
            }
            self._last_result = result
            return result

        # global peak mode
        if max_k_limit is None:
            max_k_limit = self.max_k_limit

        clusters = self.best_global(
            top_n=top_n,
            max_k=max_k,
            max_k_limit=max_k_limit,
            plot_scores=plot_scores,
            compute_all_clusters=False,
            alpha=alpha,
        )
        result = {
            "mode": "global",
            "ks": list(self.peaks_by_rank or []),
            "clusters": clusters,
            "scores": None if self.scores is None else self.scores.copy(),
            "peaks": list(self.peaks_by_rank or []),
        }
        self._last_result = result
        return result

    def plot(self, results_dir: Optional[str] = None, **kwargs) -> None:
        plot_clusters(self, results_dir=results_dir, **kwargs)

    def save(
        self,
        results_dir: str,
        top_n: int = 1,
        filename: str = "phytclust_results.tsv",
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
