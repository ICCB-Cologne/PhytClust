import os
from typing import Optional, Dict, Any, List, Tuple
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
        print("No clusters to plot – run get_clusters / best_* first.")
        return

    if pc.k is None and (pc.scores is None):
        print("Scores not available – continuing without a score plot.")

    clusters_to_plot: list[tuple[int, dict]] = []

    if pc.k is not None:
        clmap = cluster_map(pc, pc.k)
        if clmap is not None:
            clusters_to_plot.append((pc.k, clmap))
    elif n is not None:
        clmap = cluster_map(pc, n)
        if clmap is not None:
            clusters_to_plot.append((n, clmap))
    else:
        for k_val in (pc.peaks_by_rank or [])[:top_n]:
            clmap = cluster_map(pc, k_val)
            if clmap is not None:
                clusters_to_plot.append((k_val, clmap))

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
        print(f"Plotted clusters for k={k_val}.")

        if save or results_dir:
            out_dir = results_dir or "."
            out_name = filename or f"tree_k{k_val}.png"
            fig.savefig(os.path.join(out_dir, out_name), bbox_inches="tight")
            print(f"Saved figure → {os.path.join(out_dir, out_name)}")
