import os
from typing import Optional
import matplotlib.pyplot as plt

from ..algo.dp import cluster_map
from .plotting import plot_cluster


def plot_clusters(
    pc,
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
    if pc.clusters is None:
        pc.clusters = {}

    if pc.k is None and (pc.scores is None):
        print("Scores not available – continuing without a score plot.")

    def _get_map(k: int):
        try:
            get_clusters = getattr(pc, "get_clusters", None)
            if callable(get_clusters):
                return get_clusters(int(k))
        except Exception:
            pass
        return cluster_map(pc, int(k))

    clusters_to_plot: list[tuple[int, dict]] = []

    if pc.k is not None:
        k_val = int(pc.k)
        clmap = _get_map(k_val)
        if clmap is not None:
            clusters_to_plot.append((k_val, clmap))
    elif n is not None:
        k_val = int(n)
        clmap = _get_map(k_val)
        if clmap is not None:
            clusters_to_plot.append((k_val, clmap))
    else:
        for k_val in (pc.peaks_by_rank or [])[:top_n]:
            kv = int(k_val)
            clmap = _get_map(kv)
            if clmap is not None:
                clusters_to_plot.append((kv, clmap))

    if not clusters_to_plot:
        print("No clusters to plot – check your arguments.")
        return

    if (save or results_dir) and results_dir is not None:
        os.makedirs(results_dir, exist_ok=True)

    for k_val, clmap in clusters_to_plot:
        fig = plot_cluster(
            cluster=clmap,
            tree=pc.tree,
            cmap=cmap,
            outlier=outlier,
            hide_internal_nodes=hide_internal_nodes,
            show_terminal_labels=show_terminal_labels,
            width_scale=width_scale,
            height_scale=height_scale,
            label_func=label_func,
            show_branch_lengths=show_branch_lengths,
            marker_size=marker_size,
            outgroup=pc.outgroup,
            results_dir=None,
            **kwargs,
        )

        if save or results_dir:
            out_dir = results_dir or "."
            out_name = filename or f"tree_k{k_val}.png"
            out_path = os.path.join(out_dir, out_name)
            fig.savefig(out_path, bbox_inches="tight")

        plt.close(fig)
