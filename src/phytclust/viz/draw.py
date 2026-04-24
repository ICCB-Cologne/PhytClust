from __future__ import annotations

import os
import logging
from typing import Any, Callable, Optional, Union

import matplotlib as mpl
import matplotlib.collections as mpcollections

from ..exceptions import ValidationError
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

COLORS = {
    "allele_a": mpl.colors.to_rgba("orange"),
    "allele_b": mpl.colors.to_rgba("teal"),
    "clonal": mpl.colors.to_rgba("lightgrey"),
    "normal": mpl.colors.to_rgba("dimgray"),
    "gain": mpl.colors.to_rgba("red"),
    "wgd": mpl.colors.to_rgba("green"),
    "loss": mpl.colors.to_rgba("blue"),
    "chr_label": mpl.colors.to_rgba("grey"),
    "vlines": "#1f77b4",
    "marker_internal": "#1f77b4",
    "marker_terminal": "black",
    "marker_normal": "green",
    "summary_label": "grey",
    "background": "white",
    "background_hatch": "lightgray",
    "patch_background": "white",
}

LINEWIDTHS = {"copy_numbers": 2, "chr_boundary": 1, "segment_boundary": 0.5}
ALPHAS = {"patches": 0.15, "patches_wgd": 0.3, "clonal": 0.3}
SIZES = {
    "tree_marker": 40,
    "ylabel_font": 8,
    "ylabel_tick": 6,
    "xlabel_font": 10,
    "xlabel_tick": 8,
    "chr_label": 8,
}


class PlotError(Exception):
    pass


def _value_to_str(value: Optional[float]) -> Optional[str]:
    """Format a float branch-length value for display; returns None for zero/None."""
    if value is None or value == 0:
        return None
    return str(int(value)) if int(value) == value else str(value)


def _make_branch_label_func(
    branch_labels: Optional[Any],
    show_branch_lengths: bool,
) -> Callable[[Any], Optional[str]]:
    """Return the branch-label callable appropriate for the given configuration."""
    if not branch_labels:
        if show_branch_lengths:

            def _format_with_lengths(clade: Any) -> Optional[str]:
                if getattr(clade, "name", None) in (None, "root"):
                    return None
                bl = getattr(clade, "branch_length", None)
                if bl is None:
                    return None
                return _value_to_str(np.round(bl, 1))

            return _format_with_lengths
        return lambda _: None

    if isinstance(branch_labels, dict):
        return lambda clade: branch_labels.get(clade)

    if not callable(branch_labels):
        raise ValidationError("branch_labels must be either a dict or a callable")

    return lambda clade: _value_to_str(branch_labels(clade))


def plot_tree(
    input_tree: Any,
    label_func: Optional[Callable[[Any], str]] = None,
    title: str = "",
    ax: Optional[plt.Axes] = None,
    output_name: Optional[str] = None,
    outgroup: Optional[str] = None,
    width_scale: float = 1.0,
    height_scale: float = 1.0,
    show_terminal_labels: bool = False,  # currently unused, but kept for API compatibility
    show_branch_lengths: bool = True,
    show_branch_support: bool = False,
    show_events: bool = False,  # unused; kept for compatibility
    branch_labels: Optional[Union[dict[Any, str], Callable[[Any], str]]] = None,
    label_colors: Optional[Union[dict[str, str], Callable[[str], str]]] = None,
    hide_internal_nodes: bool = True,
    marker_size: Optional[int] = None,
    line_width: Optional[float] = None,
    **kwargs: Any,
) -> plt.Figure:
    """
    Minimal, fast Matplotlib phylogram for Bio.Phylo trees.
    Returns a matplotlib Figure.
    """
    marker_size = marker_size or SIZES["tree_marker"]
    line_width = line_width or LINEWIDTHS.get("segment_boundary", 1.0)
    label_func = label_func or (
        lambda x: x if isinstance(x, str) else getattr(x, "name", str(x))
    )

    horizontal_lines: list[list[tuple[float, float]]] = []
    vertical_lines: list[list[tuple[float, float]]] = []
    horizontal_colors: list[Any] = []
    vertical_colors: list[Any] = []
    horizontal_lws: list[float] = []
    vertical_lws: list[float] = []

    marker_x: list[float] = []
    marker_y: list[float] = []
    marker_sizes: list[float] = []
    marker_colors: list[Any] = []

    text_x: list[float] = []
    text_y: list[float] = []
    texts: list[str] = []
    text_colors: list[Any] = []

    if ax is None:
        nsamp = len(list(input_tree.find_clades()))
        plot_height = max(1.5, height_scale * nsamp * 0.25)
        plot_width = 5 + width_scale
        fig, ax = plt.subplots(figsize=(min(250, plot_width), min(250, plot_height)))
    else:
        fig = ax.figure

    if label_colors is None:
        clade_colors: dict[str, Any] = {}
        for clade in input_tree.find_clades():
            name = getattr(clade, "name", None)
            if not name:
                continue
            is_term = clade.is_terminal()
            clade_colors[name] = (
                COLORS["marker_terminal"] if is_term else COLORS["marker_internal"]
            )
            if outgroup is not None and name == outgroup:
                clade_colors[name] = COLORS["marker_normal"]

        get_label_color = lambda label: clade_colors.get(label, "black")
    else:
        get_label_color = (
            label_colors
            if callable(label_colors)
            else (lambda label: label_colors.get(label, "black"))
        )

    marker_func = lambda node: (
        (marker_size, get_label_color(getattr(node, "name", "")))
        if getattr(node, "name", None)
        else None
    )

    # setup axes
    ax.axes.get_yaxis().set_visible(False)
    for spine in ("right", "left", "top"):
        ax.spines[spine].set_visible(False)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune=None))
    ax.xaxis.set_tick_params(labelsize=SIZES["xlabel_tick"])
    ax.xaxis.label.set_size(SIZES["xlabel_font"])
    ax.set_title(
        title,
        x=0.01,
        y=1.0,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=16,
        zorder=10,
    )

    x_posns = _get_x_positions(input_tree)
    y_posns = _get_y_positions(
        input_tree, adjust=not hide_internal_nodes, outgroup=outgroup
    )

    xmax = max(x_posns.values()) if x_posns else 1.0
    ax.set_xlim(-0.05 * xmax, 1.05 * xmax)
    top_margin = 0.5
    ymax = (max(y_posns.values()) if y_posns else 0.0) + top_margin
    ax.set_ylim(ymax, -0.5)
    ax_scale = ax.get_xlim()[1] - ax.get_xlim()[0]

    format_branch_label = _make_branch_label_func(branch_labels, show_branch_lengths)

    def draw_clade_lines(
        *,
        use_linecollection: bool,
        orientation: str,
        y_here: float = 0.0,
        x_start: float = 0.0,
        x_here: float = 0.0,
        y_bot: float = 0.0,
        y_top: float = 0.0,
        color: Any = "black",
        lw: float = 0.1,
    ) -> None:
        if use_linecollection and orientation == "horizontal":
            horizontal_lines.append([(x_start, y_here), (x_here, y_here)])
            horizontal_colors.append(color)
            horizontal_lws.append(lw)
        elif use_linecollection and orientation == "vertical":
            vertical_lines.append([(x_here, y_bot), (x_here, y_top)])
            vertical_colors.append(color)
            vertical_lws.append(lw)

    def draw_clade(clade: Any, x_start: float, color: Any, lw: float) -> None:
        x_here = x_posns.get(clade, 0.0)
        y_here = y_posns.get(clade, 0.0)

        if hasattr(clade, "color") and clade.color is not None:
            try:
                color = clade.color.to_hex()
            except Exception:
                color = clade.color

        if hasattr(clade, "width") and clade.width is not None:
            lw = float(clade.width) * float(plt.rcParams["lines.linewidth"])

        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )

        if marker_func is not None and not (
            hide_internal_nodes and not clade.is_terminal()
        ):
            marker = marker_func(clade)
            if marker is not None:
                m_size, m_color = marker
                marker_x.append(x_here)
                marker_y.append(y_here)
                marker_sizes.append(m_size)
                marker_colors.append(m_color)

        lab = (
            label_func(clade)
            if callable(label_func)
            else str(getattr(clade, "name", ""))
        )
        if lab not in (None, clade.__class__.__name__) and not (
            hide_internal_nodes and not clade.is_terminal()
        ):
            text_x.append(x_here + min(0.02 * ax_scale, 1.0))
            text_y.append(y_here)
            texts.append(f" {lab}")
            text_colors.append(get_label_color(lab))

        if clade.clades:
            y_top = y_posns.get(clade.clades[0], y_here)
            y_bot = y_posns.get(clade.clades[-1], y_here)
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            for child in clade:
                draw_clade(child, x_here, color, lw)

    line_width = float(
        line_width if line_width is not None else plt.rcParams["lines.linewidth"]
    )
    draw_clade(input_tree.root, 0.0, "k", line_width)

    if horizontal_lines:
        h = mpcollections.LineCollection(
            horizontal_lines, colors=horizontal_colors, linewidths=horizontal_lws
        )
        ax.add_collection(h)
    if vertical_lines:
        v = mpcollections.LineCollection(
            vertical_lines, colors=vertical_colors, linewidths=vertical_lws
        )
        ax.add_collection(v)

    if marker_x:
        ax.scatter(marker_x, marker_y, s=marker_sizes, c=marker_colors, zorder=3)

    for x, y, text, color in zip(text_x, text_y, texts, text_colors):
        ax.text(x, y, text, va="center", color=color)

    ax.set_xlabel("branch length")
    ax.set_ylabel("taxa")

    # pass-through pyplot ops, e.g. axvline={'x':...}
    for key, value in kwargs.items():
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif isinstance(value, tuple) and value and isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))
        else:
            getattr(plt, str(key))(*value)

    if output_name is not None:
        plt.savefig(output_name + ".png", bbox_inches="tight")

    return fig


def _get_x_positions(tree: Any) -> dict[Any, float]:
    """
    Get x-coordinates for tree nodes based on branch lengths.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The phylogenetic tree.

    Returns
    -------
    dict[Any, float]
        Mapping of clade objects to x-coordinates (branch lengths from root).
    """
    depths = tree.depths()
    if not depths or not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths


def _get_y_positions(
    tree: Any,
    adjust: bool = False,
    outgroup: Optional[str] = None,
) -> dict[Any, float]:
    """
    Get y-coordinates for tree nodes.

    Maps each clade to a y-coordinate for drawing. Terminal leaves are assigned integer
    coordinates in reverse order. Internal nodes are positioned at the midpoint of their
    children. If an outgroup is specified, it is placed at the top.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The phylogenetic tree.
    adjust : bool, optional
        If True, adjust positions to account for variable numbers of children.
    outgroup : str, optional
        Name of outgroup clade to place at top. If provided, must match a clade name.

    Returns
    -------
    dict[Any, float]
        Mapping of clade objects to y-coordinates (vertical position).
    """
    maxheight = tree.count_terminals()
    terms = [
        x for x in tree.get_terminals() if (outgroup is None or x.name != outgroup)
    ]
    heights: dict[Any, float] = {
        tip: maxheight - 1 - i for i, tip in enumerate(reversed(terms))
    }

    if outgroup is not None:
        # Bio.Phylo find_clades(name=...)
        normal_clades = list(tree.find_clades(name=outgroup))
        if not normal_clades:
            raise PlotError(f"Outgroup '{outgroup}' not found in tree")
        heights[normal_clades[0]] = maxheight

    def calc_row(clade: Any) -> None:
        for sub in clade:
            if sub not in heights:
                calc_row(sub)
        if clade.clades:
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0

    if tree.root.clades:
        calc_row(tree.root)

    if adjust:
        sorted_clades = sorted(heights, key=heights.get)
        count = 0
        new_heights: dict[Any, float] = {}
        for cl in sorted_clades:
            if cl != tree.root:
                count += 1
            new_heights[cl] = count
        heights = new_heights

    return heights


def plot_peaks(
    scores_subset: list[float],
    peaks: list[int],
    k_start: int,
    k_end: Optional[int] = None,
    fig_width: int = 10,
    fig_height: int = 10,
    log_scale_x: bool = False,
    log_scale_y: bool = False,
    show_plot: bool = True,
) -> plt.Figure:
    """
    Plot score peaks over a range of k values.

    Creates a line plot of scores with highlighted peak points, useful for visualizing
    the optimal clustering solutions across different k values.

    Parameters
    ----------
    scores_subset : list[float]
        Array of scores for each k value.
    peaks : list[int]
        Indices of peak positions within scores_subset.
    k_start : int
        Starting k value (displayed on x-axis).
    k_end : int, optional
        Ending k value. If None, uses the full range.
    fig_width : int, optional
        Figure width in inches (default: 10).
    fig_height : int, optional
        Figure height in inches (default: 10).
    log_scale_x : bool, optional
        Apply log scale to x-axis (default: False).
    log_scale_y : bool, optional
        Apply log scale to y-axis (default: False).
    show_plot : bool, optional
        Display plot immediately (default: True).

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """
    title_fontsize, label_fontsize = 16, 14
    tick_labelsize, legend_fontsize, peak_labelsize = 12, 12, 12

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    x_vals = np.arange(len(scores_subset)) + (k_start + 1)
    ax.plot(x_vals, scores_subset, label="Scores")

    peak_x_vals = np.array(peaks) + (k_start + 1)
    peak_y_vals = np.array([scores_subset[p] for p in peaks])
    ax.plot(
        peak_x_vals, peak_y_vals, "x", markersize=10, color="red", label="Top Peaks"
    )

    ax.set_title("Top Peaks for the given Score Range", fontsize=title_fontsize)
    ax.set_xlabel("k", fontsize=label_fontsize)
    ax.set_ylabel("Score", fontsize=label_fontsize)
    ax.tick_params(axis="both", labelsize=tick_labelsize)

    if log_scale_x:
        ax.set_xscale("log")
    if log_scale_y:
        ax.set_yscale("log")

    if k_end is not None:
        ax.set_xlim(k_start, k_end + 1)
        ax.set_title(
            f"Top Peaks for k between {k_start} and {k_end}", fontsize=title_fontsize
        )

    data_min, data_max = min(scores_subset), max(scores_subset)
    offset_y = 0.02 if data_max == data_min else 0.02 * (data_max - data_min)

    for i, (px, py) in enumerate(zip(peak_x_vals, peak_y_vals)):
        ax.text(
            px,
            py + offset_y,
            str(peaks[i]),
            fontsize=peak_labelsize,
            ha="center",
            va="bottom",
        )

    ax.legend(fontsize=legend_fontsize)
    if show_plot:
        plt.show()
    return fig


def plot_cluster(
    cluster: dict[Any, int],
    tree: Any,
    *,
    cmap: str | mcolors.Colormap = "phytclust",
    save: bool = False,
    filename: str | None = None,
    results_dir: str | None = None,
    outlier: bool = False,
    hide_internal_nodes: bool = True,
    show_terminal_labels: bool = False,  # passed to plot_tree (currently not used inside)
    width_scale: float = 2.0,
    height_scale: float = 0.4,
    label_func: Callable[[Any], str] | None = None,
    show_branch_lengths: bool = False,
    marker_size: int = 50,
    outgroup: str | None = None,
    scores: list[float] | np.ndarray | None = None,
    **kwargs: Any,
) -> plt.Figure:
    """
    Plot a single clustering solution on the phylogenetic tree.

    Displays clusters as colored leaf nodes, with each cluster assigned a distinct color
    from the colormap. Can optionally save to file and annotate with clustering scores.

    Parameters
    ----------
    cluster : dict[Any, int]
        Mapping of leaf clades to cluster IDs.
    tree : Bio.Phylo.BaseTree.Tree
        The phylogenetic tree.
    cmap : str | matplotlib.colors.Colormap, optional
        Colormap name or object (default: "tab20").
    save : bool, optional
        Save figure to file (default: False).
    filename : str, optional
        Output filename (default: auto-generated as "tree_k{n}.png").
    results_dir : str, optional
        Directory for saving results (default: current directory).
    outlier : bool, optional
        Highlight singleton clusters differently (currently unused).
    hide_internal_nodes : bool, optional
        Hide internal node markers (default: True).
    show_terminal_labels : bool, optional
        Display terminal labels (default: False).
    width_scale : float, optional
        Scale factor for figure width (default: 2.0).
    height_scale : float, optional
        Scale factor for figure height (default: 0.4).
    label_func : Callable, optional
        Function to format node labels.
    show_branch_lengths : bool, optional
        Display branch length labels (default: False).
    marker_size : int, optional
        Size of leaf node markers (default: 50).
    outgroup : str, optional
        Name of outgroup (displayed in grey).
    scores : list[float] | np.ndarray, optional
        Clustering scores for title annotation.
    **kwargs : Any
        Additional arguments passed to plot_tree.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """
    n_unique = len(set(cluster.values()))
    if isinstance(cmap, str) and cmap == "phytclust":
        from .palette import expand_palette

        palette = expand_palette(max(n_unique, 8))
    elif isinstance(cmap, str):
        cmap_obj = plt.get_cmap(cmap)
        palette = list(
            cmap_obj.colors
            if hasattr(cmap_obj, "colors")
            else cmap_obj(np.linspace(0, 1, getattr(cmap_obj, "N", 20)))
        )
    else:
        palette = list(
            cmap.colors
            if hasattr(cmap, "colors")
            else cmap(np.linspace(0, 1, getattr(cmap, "N", 20)))
        )

    # Ensure we never repeat colors: extend palette if needed
    if n_unique > len(palette):
        from .palette import expand_palette

        palette = expand_palette(n_unique)

    # Keep RGBA channels (including alpha) and map actual cluster IDs to unique
    # palette indices. This avoids modulo collisions for sparse/non-contiguous IDs.
    palette = [tuple(col[:4]) for col in palette]

    ids = np.fromiter(cluster.values(), dtype=int)
    unique_ids = sorted(set(int(v) for v in ids.tolist()))
    id_to_idx = {cid: i for i, cid in enumerate(unique_ids)}
    colours: list[Any] = [palette[id_to_idx[int(cid)] % len(palette)] for cid in ids]

    # outgroup → grey
    if outgroup is not None:
        for i, leaf in enumerate(cluster.keys()):
            if getattr(leaf, "name", None) == outgroup:
                colours[i] = "grey"

    unique, counts = np.unique(ids, return_counts=True)

    leaf_names = [getattr(leaf, "name", str(leaf)) for leaf in cluster.keys()]
    clumap = dict(zip(leaf_names, colours))
    n_clusters = len(unique)

    title = ""
    if scores is not None and (n_clusters - 1) < len(scores):
        try:
            title = f"clusters={n_clusters}, score={float(scores[n_clusters-1]):.4f}"
        except Exception:
            title = f"clusters={n_clusters}"

    fig = plot_tree(
        tree,
        title=title,
        label_colors=clumap,
        hide_internal_nodes=hide_internal_nodes,
        show_terminal_labels=show_terminal_labels,
        width_scale=width_scale,
        height_scale=height_scale,
        label_func=label_func,
        show_branch_lengths=show_branch_lengths,
        marker_size=marker_size,
        outgroup=outgroup,
        **kwargs,
    )

    if save:
        results_dir = results_dir or "."
        os.makedirs(results_dir, exist_ok=True)
        filename = filename or f"tree_k{n_clusters}.png"
        fig.savefig(os.path.join(results_dir, filename), bbox_inches="tight")

    return fig


def plot_multiple_clusters(
    input_df: pd.DataFrame,
    final_tree: Optional[Any] = None,
    y_posns: Optional[dict[str, int]] = None,
    cmax: Optional[int] = None,
    tree_width_ratio: float = 1.0,
    cbar_width_ratio: float = 0.05,
    figsize: tuple[int, int] = (20, 10),
    tree_marker_size: int = 0,
    show_internal_nodes: bool = False,
    title: str = "",
    tree_label_func: Optional[Callable[[Any], str]] = None,
    cmap: str = "tab20b",
    outgroup: str = "diploid",
    fixed_x_range: tuple[int, int] = (10000, 50000),
) -> plt.Figure:
    """
    Plot multiple clustering solutions simultaneously.

    Creates a heatmap-style visualization showing how cluster assignments change
    across multiple clustering solutions, displayed alongside the phylogenetic tree.

    Parameters
    ----------
    input_df : pd.DataFrame
        DataFrame with species names as index and cluster assignments as columns,
        where each column represents a different clustering solution.
    final_tree : Bio.Phylo.BaseTree.Tree, optional
        The phylogenetic tree to display alongside clusters.
    y_posns : dict[str, int], optional
        Pre-computed y-positions for species (default: auto-computed from tree).
    cmax : int, optional
        Maximum cluster ID for colormap scaling (default: auto).
    tree_width_ratio : float, optional
        Relative width of tree to heatmap (default: 1.0).
    cbar_width_ratio : float, optional
        Width of colorbar as fraction of heatmap (default: 0.05).
    figsize : tuple[int, int], optional
        Figure dimensions in inches (default: (20, 10)).
    tree_marker_size : int, optional
        Size of tree node markers (default: 0).
    show_internal_nodes : bool, optional
        Display internal node labels (default: False).
    title : str, optional
        Figure title (default: empty).
    tree_label_func : Callable, optional
        Function to format tree node labels.
    cmap : str, optional
        Colormap name for clusters (default: "tab20b").
    outgroup : str, optional
        Name of outgroup species (default: "diploid").
    fixed_x_range : tuple[int, int], optional
        Fixed x-axis range for tree (default: (10000, 50000)).

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """
    cmax = cmax or int(np.max(input_df.values.astype(int)))
    sample_labels = input_df.index.get_level_values("leaf_name").unique()
    input_df = input_df.sort_index(level="comparison_IDs")

    if final_tree is None:
        fig, ax = plt.subplots(
            figsize=figsize, ncols=1, sharey=False, gridspec_kw={"width_ratios": [1]}
        )
        if not show_internal_nodes:
            logger.warning(
                'No tree provided, so "show_internal_nodes=False" is ignored'
            )
        y_posns = y_posns or {s: i for i, s in enumerate(sample_labels)}
        ax.set_title(
            title,
            x=0,
            y=1,
            ha="left",
            va="bottom",
            pad=20,
            fontweight="bold",
            fontsize=16,
            zorder=10,
        )
    else:
        fig, axs = plt.subplots(
            figsize=(figsize[0] * 0.5, figsize[1]),
            ncols=2,
            sharey=False,
            gridspec_kw={"width_ratios": [tree_width_ratio * 0.2, 0.5], "wspace": 0.05},
        )
        tree_ax, ax = axs
        y_posns = {
            k.name: v
            for k, v in _get_y_positions(
                final_tree, adjust=show_internal_nodes, outgroup=outgroup
            ).items()
        }
        plot_tree(
            final_tree,
            ax=tree_ax,
            outgroup=outgroup,
            label_func=tree_label_func or (lambda _: ""),
            hide_internal_nodes=True,
            show_branch_lengths=False,
            show_events=False,
            line_width=0.5,
            marker_size=tree_marker_size,
            title=title,
            label_colors=None,
        )
        tree_ax.set_axis_off()
        fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0.0, wspace=150)

    # The heatmap axis: when a tree is present it's the second subplot,
    # otherwise it's the single axis returned by subplots.
    heat_ax = ax if final_tree is None else axs[1]

    # reorder labels to match y positions
    ind = [y_posns.get(x, -1) for x in sample_labels]
    sample_labels = sample_labels[np.argsort(ind)]
    color_norm = mcolors.Normalize(vmin=0, vmax=cmax)

    solution_ends = input_df.loc[sample_labels[0]].copy()
    solution_ends["end_pos"] = np.cumsum([1] * len(solution_ends))
    solution_ends = (
        solution_ends.reset_index()
        .groupby("comparison_IDs")
        .max()["end_pos"]
        .dropna()
        .astype(int)
    )
    x_pos = np.linspace(
        fixed_x_range[0],
        fixed_x_range[1],
        len(input_df.loc[sample_labels].astype(int).unstack("leaf_name")) + 1,
    )
    y_pos = np.arange(len(sample_labels) + 1) + 0.5

    # heatmap data assumes columns MultiIndex (..., 'cluster_ID')
    data = (
        input_df.loc[sample_labels]
        .astype(int)
        .unstack("leaf_name")
        .loc[:, "cluster_ID"]
        .loc[:, sample_labels]
        .values.T
    )
    heat_ax.pcolormesh(x_pos, y_pos, data, cmap=cmap, norm=color_norm)

    # vertical separators between solutions
    for line in solution_ends.values:
        heat_ax.axvline(x_pos[line], color="black", linewidth=0.75)

    xtick_pos = np.append([0], x_pos[solution_ends.values][:-1])
    xtick_pos = (xtick_pos + np.roll(xtick_pos, -1)) / 2
    xtick_pos[-1] += x_pos[-1] / 2
    heat_ax.set_xticks(xtick_pos)
    heat_ax.set_xticklabels(
        [x[3:] for x in solution_ends.index],
        ha="center",
        rotation=0,
        va="bottom",
    )
    heat_ax.tick_params(width=0)
    heat_ax.xaxis.set_tick_params(labelbottom=False, labeltop=True, bottom=False)
    heat_ax.set_yticks([])
    heat_ax.set_ylim(len(sample_labels) + 0.5, 0.5)

    logger.debug("Finished plot_multiple_clusters")
    return fig
