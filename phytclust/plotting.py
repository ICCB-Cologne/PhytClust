import os
import logging
import matplotlib.collections as mpcollections
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import logging
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

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
LINEWIDTHS = {
    "copy_numbers": 2,
    "chr_boundary": 1,
    "segment_boundary": 0.5,
}
ALPHAS = {
    "patches": 0.15,
    "patches_wgd": 0.3,
    "clonal": 0.3,
}
SIZES = {
    "tree_marker": 40,
    "ylabel_font": 8,
    "ylabel_tick": 6,
    "xlabel_font": 10,
    "xlabel_tick": 8,
    "chr_label": 8,
}


def plot_tree(
    input_tree: Any,
    label_func: Optional[Callable[[Any], str]] = None,
    title: str = "",
    ax: Optional[plt.Axes] = None,
    output_name: Optional[str] = None,
    outgroup: Optional[str] = None,
    width_scale: float = 1,
    height_scale: float = 1,
    show_terminal_labels: bool = False,
    show_branch_lengths: bool = True,
    show_branch_support: bool = False,
    show_events: bool = False,
    branch_labels: Optional[Union[Dict[Any, str], Callable[[Any], str]]] = None,
    label_colors: Optional[Union[Dict[str, str], Callable[[str], str]]] = None,
    hide_internal_nodes: bool = False,
    marker_size: Optional[int] = None,
    line_width: Optional[float] = None,
    **kwargs: Any,
) -> plt.Figure:
    """Plot the given tree using matplotlib (or pylab).
    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.
    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).
    Example using the pyplot options 'axhspan' and 'axvline'::
        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})
    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).
    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_support : bool
            Whether to display confidence values, if present on the tree.
        ax : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_support`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.
    """

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MEDICCPlotError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None

    import matplotlib.collections as mpcollections

    if ax is None:
        nsamp = len(list(input_tree.find_clades()))
        plot_height = height_scale * nsamp * 0.25
        max_leaf_to_root_distances = np.max(
            [
                np.sum([x.branch_length for x in input_tree.get_path(leaf)])
                for leaf in input_tree.get_terminals()
            ]
        )
        plot_width = 5 + np.max(
            [0, width_scale * np.log10(max_leaf_to_root_distances / 100) * 5]
        )

        # maximum figure size is 250x250 inches
        fig, ax = plt.subplots(figsize=(min(250, plot_width), min(250, plot_height)))

    label_func = (
        label_func
        if label_func is not None
        else lambda x: x.name if hasattr(x, "name") else x
    )

    # options for displaying label colors.
    if label_colors is not None:
        get_label_color = (
            label_colors
            if callable(label_colors)
            else lambda label: label_colors.get(label, "black")
        )
    else:
        clade_colors = {}
        for sample in [x.name for x in list(input_tree.find_clades(""))]:
            is_terminal = input_tree.find_clades(sample).__next__().is_terminal()
            clade_colors[sample] = (
                COLORS["marker_terminal"] if is_terminal else COLORS["marker_internal"]
            )
            if sample == outgroup:
                clade_colors[sample] = COLORS["marker_normal"]

        def get_label_color(label):
            return clade_colors.get(label, "black")

    marker_size = marker_size if marker_size is not None else SIZES["tree_marker"]
    marker_func = lambda x: (
        (marker_size, get_label_color(x.name)) if x.name is not None else None
    )

    ax.axes.get_yaxis().set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
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

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def value_to_str(value):
        if value is None or value == 0:
            return None
        return str(int(value)) if int(value) == value else str(value)

    if not branch_labels:
        if show_branch_lengths:

            def format_branch_label(x):
                return (
                    value_to_str(np.round(x.branch_length, 1))
                    if x.name != "root" and x.name is not None
                    else None
                )

        else:

            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):

        def format_branch_label(clade):
            return branch_labels.get(clade)

    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )

        def format_branch_label(clade):
            return value_to_str(branch_labels(clade))

    if show_branch_support:

        def format_support_value(clade):
            if clade.name == "root" or clade.name is None:
                return None
            try:
                confidences = clade.confidences
            except AttributeError:
                pass
            else:
                return "/".join(value_to_str(cnf.value) for cnf in confidences)
            return (
                value_to_str(clade.confidence) if clade.confidence is not None else None
            )

    def draw_clade_lines(
        use_linecollection: bool = False,
        orientation: str = "horizontal",
        y_here: float = 0,
        x_start: float = 0,
        x_here: float = 0,
        y_bot: float = 0,
        y_top: float = 0,
        color: str = "black",
        lw: float = 0.1,
    ) -> None:
        """Create a line with or without a line collection object.
        """
        if not use_linecollection and orientation == "horizontal":
            ax.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif not use_linecollection and orientation == "vertical":
            ax.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )

    def draw_clade(clade: Any, x_start: float, color: str, lw: float) -> None:
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
            # Add the label handling code here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        if marker_func is not None:
            marker = marker_func(clade)
            if (
                marker is not None
                and clade is not None
                and not (hide_internal_nodes and not clade.is_terminal())
            ):
                marker_size, marker_col = marker_func(clade)
                ax.scatter(x_here, y_here, s=marker_size, color=marker_col, zorder=3)
        label = label_func(str(clade.name))
        ax_scale = ax.get_xlim()[1] - ax.get_xlim()[0]
        if label not in (None, clade.__class__.__name__) and not (
            hide_internal_nodes and not clade.is_terminal()
        ):
            ax.text(
                x_here + min(0.02 * ax_scale, 1),
                y_here,
                f" {label}",
                verticalalignment="center",
                color=get_label_color(label),
            )
        if (
            clade.name is not None
            and (clade.is_terminal() and show_terminal_labels)
            or not hide_internal_nodes
        ):
            label = label_func(str(clade.name))
            ax_scale = ax.get_xlim()[1] - ax.get_xlim()[0]
            ax.text(
                x_here + min(0.02 * ax_scale, 1),
                y_here,
                " %s" % label,  # Display the label
                verticalalignment="center",
                color=get_label_color(label),
            )
        conf_label = format_branch_label(clade)
        if conf_label:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )
        if show_branch_support:
            support_value = format_support_value(clade)
            if support_value:
                ax.text(
                    0.5 * (x_start + x_here),
                    y_here + 0.25,
                    support_value + "%",
                    fontsize="small",
                    color="grey",
                    horizontalalignment="center",
                )
        if show_events and clade.events is not None:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                clade.events,
                fontsize="small",
                color=COLORS["marker_normal"],
                horizontalalignment="center",
            )
        if clade.clades:
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
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

    line_width = (
        line_width if line_width is not None else plt.rcParams["lines.linewidth"]
    )
    draw_clade(input_tree.root, 0, "k", line_width)

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        ax.add_collection(i)
    for i in vertical_linecollections:
        ax.add_collection(i)

    ax.set_xlabel("branch length")
    ax.set_ylabel("taxa")

    # Add margins around the `tree` to prevent overlapping the ax
    xmax = max(x_posns.values())
    # ax.set_xlim(-0.05 * xmax, 1.25 * xmax)
    ax.set_xlim(-0.05 * xmax, 1.05 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    top_margin = 0.5
    ax.set_ylim(max(y_posns.values()) + top_margin, -0.5)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                f'Keyword argument "{key}={value}" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) "
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if output_name is not None:
        plt.savefig(output_name + ".png", bbox_inches="tight")

    return plt.gcf()


def _get_x_positions(tree: Any) -> Dict[Any, float]:
    """Create a mapping of each clade to its horizontal position.
    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths


def _get_y_positions(tree: Any, adjust: bool = False, outgroup: str = "outgroup") -> Dict[Any, float]:
    """Create a mapping of each clade to its vertical position.
        Dict of {clade: y-coord}.
    Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    heights = {
        tip: maxheight - 1 - i
        for i, tip in enumerate(
            reversed(
                [
                    x
                    for x in tree.get_terminals()
                    if outgroup is None or x.name != outgroup
                ]
            )
        )
    }
    if outgroup is not None:
        normal_clades = list(tree.find_clades(outgroup))
        if not normal_clades:
            raise PlotError(f"Normal clade {outgroup} not found in tree")
        heights[normal_clades[0]] = maxheight

    # Internal nodes: place at midpoint of children
    def calc_row(clade: Any) -> None:
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)

    if adjust:
        pos = pd.DataFrame(
            [(clade, val) for clade, val in heights.items()], columns=["clade", "pos"]
        ).sort_values("pos")
        pos["newpos"] = 0
        count = 0
        for i in pos.index:
            if pos.loc[i, "clade"] != tree.root:
                count += 1
            pos.loc[i, "newpos"] = count

        pos.set_index("clade", inplace=True)
        heights = pos.to_dict()["newpos"]

    return heights


def plot_peaks(
    scores_subset: List[float],
    peaks: List[int],
    k_start: int,
    k_end: Optional[int] = None,
) -> None:
    """
    Plots the peaks of a given score subset.

    Args:
        scores_subset (List[float]): Subset of scores to plot.
        peaks (List[int]): Indices of the peaks in the scores_subset.
        k_start (int): Starting value of k.
        k_end (Optional[int], optional): Ending value of k. Defaults to None.

    Returns:
        None
    """
    plt.figure()
    plt.plot(np.arange(len(scores_subset)) + 1 + k_start, scores_subset)
    plt.plot(
        peaks + k_start + 1,
        scores_subset[peaks],
        "x",
        markersize=10,
        label="Top Peaks",
        color="red",
    )
    plt.legend()
    plt.title("Top Peaks for the given Score Range")
    plt.xlabel("k")
    plt.ylabel("Score")

    if k_end is not None:
        plt.xlim(k_start, k_end + 1)
        plt.title(f"Top Peaks for k between {k_start} and {k_end}")
    plt.show()


def plot_cluster(
    cluster: Dict[Any, int],
    tree: Any,
    cmap: Any,
    save: bool = False,
    filename: Optional[str] = None,
    outlier: bool = False,
    hide_internal_nodes: bool = True,
    show_terminal_labels: bool = False,
    width_scale: float = 2,
    height_scale: float = 0.4,
    label_func: Optional[callable] = None,
    show_branch_lengths: bool = False,
    marker_size: int = 50,
    outgroup: Optional[str] = None,
    results_dir: Optional[str] = None,
    **kwargs: Any,
) -> plt.Figure:
    """
    Plots a cluster on a phylogenetic tree.

    Args:
        cluster (dict): Dictionary mapping clades to cluster IDs.
        cluster_number (int): The number of the cluster being plotted.
        tree (Tree): The phylogenetic tree to plot.
        cmap (Colormap): Colormap for the clusters from matplotlib (e.g. tab20).
        save (bool, optional): Whether to save the plot. Defaults to False.
        filename (str, optional): Filename to save the plot. Defaults to None.
        outlier (bool, optional): Whether to highlight outliers. Defaults to False.
        hide_internal_nodes (bool, optional): Whether to hide internal nodes. Defaults to True.
        show_terminal_labels (bool, optional): Whether to show terminal labels. Defaults to False.
        width_scale (float, optional): Scale factor for the width of the plot. Defaults to 2.
        height_scale (float, optional): Scale factor for the height of the plot. Defaults to 0.4.
        label_func (function, optional): Function to generate labels. Defaults to lambda x: None.
        show_branch_lengths (bool, optional): Whether to show branch lengths. Defaults to False.
        marker_size (int, optional): Size of the markers. Defaults to 50.
        outgroup (str, optional): Name of the outgroup. Defaults to None.
        results_dir (str, optional): Directory to save the results. Defaults to None.
        **kwargs: Additional keyword arguments.

    Returns:
        Figure: The generated plot.
    """
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    colors = cmap.colors if hasattr(cmap, 'colors') else cmap(np.linspace(0, 1, cmap.N))

    cluster_sizes = defaultdict(int)
    for cluster_id in cluster.values():
        cluster_sizes[cluster_id] += 1

    clumap = {}
    for clade, cluster_id in cluster.items():
        clade_name = clade.name if hasattr(clade, "name") else None
        if clade_name == outgroup:
            clumap[clade_name] = "grey"
        elif outlier and cluster_sizes[cluster_id] == 1:
            clumap[clade_name] = "black"
        else:
            color_index = cluster_id % len(colors)
            clumap[clade_name] = colors[color_index]

    cluster_number = len(set(cluster.values()))

    fig = plot_tree(
        tree,
        title=(
            f"No. of clusters: {cluster_number}, Score: {kwargs.get('scores')[cluster_number - 1]:.4f}"
            if kwargs.get("scores")
            else ""
        ),
        label_colors=clumap,
        hide_internal_nodes=hide_internal_nodes,
        show_terminal_labels=show_terminal_labels,
        width_scale=width_scale,
        height_scale=height_scale,
        label_func=label_func,
        show_branch_lengths=show_branch_lengths,
        marker_size=marker_size,
        **kwargs,
    )
    if save:
        filename = filename or f"num_clusters_{cluster_number}.png"
        full_path = os.path.join(results_dir, filename)
        plt.savefig(full_path)

    return fig


def plot_multiple_clusters(
    input_df: pd.DataFrame,
    final_tree: Optional[Any] = None,
    y_posns: Optional[Dict[str, int]] = None,
    cmax: Optional[int] = None,
    tree_width_ratio: float = 1,
    cbar_width_ratio: float = 0.05,
    figsize: Tuple[int, int] = (20, 10),
    tree_marker_size: int = 0,
    show_internal_nodes: bool = False,
    title: str = "",
    tree_label_func: Optional[callable] = None,
    cmap: str = "tab20b",
    outgroup: str = "diploid",
    fixed_x_range: Tuple[int, int] = (10000, 50000),
) -> plt.Figure:
    """
    Plots multiple clusters on a phylogenetic tree.

    Args:
        input_df (pd.DataFrame): DataFrame containing the input data.
        final_tree (Optional[Any], optional): The phylogenetic tree to plot. Defaults to None.
        y_posns (Optional[Dict[str, int]], optional): Y positions for the tree. Defaults to None.
        cmax (Optional[int], optional): Maximum value for color normalization. Defaults to None.
        tree_width_ratio (float, optional): Width ratio for the tree plot. Defaults to 1.
        cbar_width_ratio (float, optional): Width ratio for the color bar. Defaults to 0.05.
        figsize (Tuple[int, int], optional): Figure size. Defaults to (20, 10).
        tree_marker_size (int, optional): Size of the tree markers. Defaults to 0.
        show_internal_nodes (bool, optional): Whether to show internal nodes. Defaults to False.
        title (str, optional): Title of the plot. Defaults to "".
        tree_label_func (Optional[callable], optional): Function to generate tree labels. Defaults to None.
        cmap (str, optional): Colormap for the clusters. Defaults to "tab20b".
        outgroup (str, optional): Name of the outgroup. Defaults to "diploid".
        fixed_x_range (Tuple[int, int], optional): Fixed x range for the plot. Defaults to (10000, 50000).

    Returns:
        plt.Figure: The generated plot.
    """
    if cmax is None:
        cmax = np.max(input_df.values.astype(int))

    sample_labels = input_df.index.get_level_values("leaf_name").unique()

    if final_tree is None:
        fig, ax = plt.subplots(
            figsize=figsize,
            ncols=1,
            sharey=False,
            gridspec_kw={"width_ratios": [1]},
        )
        if not show_internal_nodes:
            logger.warning(
                'No tree provided, so "show_internal_nodes=False" is ignored'
            )
        if y_posns is None:
            y_posns = {s: i for i, s in enumerate(sample_labels)}
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
            figsize=figsize,
            ncols=2,
            sharey=False,
            gridspec_kw={"width_ratios": [tree_width_ratio, 1]},
        )
        tree_ax = axs[0]
        ax = axs[1]

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
            label_func=tree_label_func if tree_label_func is not None else lambda x: "",
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

    ind = [y_posns.get(x, -1) for x in sample_labels]
    sample_labels = sample_labels[np.argsort(ind)]
    color_norm = mcolors.Normalize(vmin=1, vmax=cmax)

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
    im = ax.pcolormesh(
        x_pos,
        y_pos,
        input_df.loc[sample_labels]
        .astype(int)
        .unstack("leaf_name")
        .loc[:, "cluster_ID"]
        .loc[:, sample_labels]
        .values.T,
        cmap=cmap,
        norm=color_norm,
    )

    for line in solution_ends.values:
        ax.axvline(x_pos[line], color="black", linewidth=0.75)

    xtick_pos = np.append([0], x_pos[solution_ends.values][:-1])
    xtick_pos = (xtick_pos + np.roll(xtick_pos, -1)) / 2
    xtick_pos[-1] += x_pos[-1] / 2
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(
        [x[3:] for x in solution_ends.index], ha="center", rotation=0, va="bottom"
    )
    ax.tick_params(width=0)
    ax.xaxis.set_tick_params(labelbottom=False, labeltop=True, bottom=False)
    ax.set_yticks([])

    ax.set_ylim(len(sample_labels) + 0.5, 0.5)

    logger.debug("Finished plot_cn_heatmap function")
    return fig


class PlotError(Exception):
    pass
