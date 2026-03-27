import os
import numpy as np
from typing import Optional, Callable, Tuple, List
import matplotlib.pyplot as plt

from ..algo.dp import cluster_map
from .plotting import plot_cluster


def plot_clusters(
    pc,
    results_dir: Optional[str] = None,
    top_n: int = 1,
    n: Optional[int] = None,
    k: Optional[int] = None,
    cmap: plt.cm = plt.get_cmap("tab20"),
    show_terminal_labels: bool = False,
    outlier: bool = False,
    save: bool = False,
    filename: Optional[str] = None,
    hide_internal_nodes: bool = True,
    width_scale: float = 2,
    height_scale: float = 0.1,
    label_func: Optional[Callable[[any], Tuple[float, str]]] = None,
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

    # Priority: explicit k param > n param > pc.k > peaks
    if k is not None:
        k_val = int(k)
        clmap = _get_map(k_val)
        if clmap is not None:
            clusters_to_plot.append((k_val, clmap))
    elif n is not None:
        k_val = int(n)
        clmap = _get_map(k_val)
        if clmap is not None:
            clusters_to_plot.append((k_val, clmap))
    elif pc.k is not None:
        k_val = int(pc.k)
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
        else:
            plt.show()


def plot_multiple_k(
    pc,
    k_values: Optional[List[int]] = None,
    results_dir: Optional[str] = None,
    top_n: int = 3,
    cmap: plt.cm = plt.get_cmap("tab20"),
    show_terminal_labels: bool = False,
    hide_internal_nodes: bool = True,
    width_scale: float = 1.5,
    height_scale: float = 0.08,
    label_func: Optional[Callable[[any], Tuple[float, str]]] = None,
    show_branch_lengths: bool = False,
    marker_size: int = 30,
    save: bool = False,
    **kwargs,
) -> None:
    """
    Plot multiple cluster solutions (different k values) in a grid for comparison.

    Parameters
    ----------
    pc : PhytClust
        The PhytClust object with computed clusters.
    k_values : list[int], optional
        Specific k values to plot. If None, uses top_n peaks.
    results_dir : str, optional
        Directory for saving output.
    top_n : int, optional
        Number of peaks to use if k_values not specified (default: 3).
    cmap : colormap, optional
        Matplotlib colormap for clusters (default: "tab20").
    show_terminal_labels : bool, optional
        Show leaf names (default: False).
    hide_internal_nodes : bool, optional
        Hide internal node markers (default: True).
    width_scale : float, optional
        Scale factor for subplot width (default: 1.5).
    height_scale : float, optional
        Scale factor for subplot height (default: 0.08).
    label_func : callable, optional
        Function to format labels.
    show_branch_lengths : bool, optional
        Show branch length labels (default: False).
    marker_size : int, optional
        Size of leaf markers (default: 30).
    save : bool, optional
        Save figures to files named tree_k{k}.png (default: False).
    **kwargs : dict
        Additional arguments passed to plot_cluster.
    """
    if pc.clusters is None:
        pc.clusters = {}

    def _get_map(k: int):
        try:
            get_clusters = getattr(pc, "get_clusters", None)
            if callable(get_clusters):
                return get_clusters(int(k))
        except Exception:
            pass
        return cluster_map(pc, int(k))

    # Determine which k values to plot
    clusters_to_plot: List[Tuple[int, dict]] = []

    if k_values is not None:
        for k_val in k_values:
            try:
                clmap = _get_map(k_val)
                if clmap is not None:
                    clusters_to_plot.append((int(k_val), clmap))
            except Exception:
                print(f"Warning: Could not compute clusters for k={k_val}")
    else:
        for k_val in (pc.peaks_by_rank or [])[:top_n]:
            try:
                clmap = _get_map(int(k_val))
                if clmap is not None:
                    clusters_to_plot.append((int(k_val), clmap))
            except Exception:
                print(f"Warning: Could not compute clusters for k={k_val}")

    if not clusters_to_plot:
        print("No clusters to plot – check your arguments.")
        return

    # Create grid of subplots
    n_plots = len(clusters_to_plot)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    n_leaves = len(list(pc.tree.get_terminals()))
    plot_height = max(1.5, height_scale * n_leaves * n_rows * 0.25)
    plot_width = width_scale * n_cols + 2
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(min(250, plot_width), min(250, plot_height)),
        squeeze=False
    )
    axes_flat = axes.flatten()

    # Hide unused subplots
    for i in range(n_plots, len(axes_flat)):
        axes_flat[i].set_visible(False)

    if (save or results_dir) and results_dir is not None:
        os.makedirs(results_dir, exist_ok=True)

    # Plot each cluster solution (saves individually or displays)
    for k_val, clmap in clusters_to_plot:
        try:
            fig_single = plot_cluster(
                cluster=clmap,
                tree=pc.tree,
                cmap=cmap,
                outlier=False,
                hide_internal_nodes=hide_internal_nodes,
                show_terminal_labels=show_terminal_labels,
                width_scale=width_scale,
                height_scale=height_scale,
                label_func=label_func,
                show_branch_lengths=show_branch_lengths,
                marker_size=marker_size,
                outgroup=pc.outgroup,
                **kwargs,
            )
            if save or results_dir:
                out_dir = results_dir or "."
                out_name = f"tree_k{k_val}.png"
                out_path = os.path.join(out_dir, out_name)
                fig_single.savefig(out_path, bbox_inches="tight", dpi=150)
                plt.close(fig_single)
            else:
                plt.show()
        except Exception as e:
            print(f"Error plotting k={k_val}: {e}")

    # Close the empty grid figure
    plt.close(fig)
