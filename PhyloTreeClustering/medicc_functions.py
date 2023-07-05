import Bio
from Bio import Phylo # I added that 

import logging

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

COL_ALLELE_A = mpl.colors.to_rgba('orange')
COL_ALLELE_B = mpl.colors.to_rgba('teal')
COL_CLONAL = mpl.colors.to_rgba('lightgrey')
COL_NORMAL = mpl.colors.to_rgba('dimgray')
COL_GAIN = mpl.colors.to_rgba('red')
COL_WGD = mpl.colors.to_rgba('green')
COL_LOSS = mpl.colors.to_rgba('blue')
COL_CHR_LABEL = mpl.colors.to_rgba('grey')
COL_VLINES = '#1f77b4'
COL_MARKER_INTERNAL = COL_VLINES
COL_MARKER_TERMINAL = 'black'
COL_MARKER_NORMAL = 'green'
COL_SUMMARY_LABEL = 'grey'
COL_BACKGROUND = 'white'
COL_BACKGROUND_HATCH = 'lightgray'
COL_PATCH_BACKGROUND = 'white'
LINEWIDTH_COPY_NUMBERS = 2
LINEWIDTH_CHR_BOUNDARY = 1
LINEWIDTH_SEGMENT_BOUNDARY = 0.5
ALPHA_PATCHES = 0.15
ALPHA_PATCHES_WGD = 0.3
ALPHA_CLONAL = 0.3
BACKGROUND_HATCH_MARKER = '/////'
TREE_MARKER_SIZE = 40
YLABEL_FONT_SIZE = 8
YLABEL_TICK_SIZE = 6
XLABEL_FONT_SIZE = 10
XLABEL_TICK_SIZE = 8
CHR_LABEL_SIZE = 8
SMALL_SEGMENTS_LIMIT = 1e7



def _get_x_positions(tree):
    """Create a mapping of each clade to its horizontal position.
    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths

def _get_y_positions(tree, adjust=False, normal_name='diploid'):
    """Create a mapping of each clade to its vertical position.
    Dict of {clade: y-coord}.
Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    heights = {tip: maxheight -1 -i for i,
            tip in enumerate(reversed([x for x in tree.get_terminals() if x.name != normal_name]))}
    heights.update({list(tree.find_clades(normal_name))[0]: maxheight})

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
        
    if adjust:
        pos = pd.DataFrame([(clade, val) for clade, val in heights.items()], columns=['clade','pos']).sort_values('pos')
        pos['newpos'] = 0
        count = 0
        for i in pos.index:
            if pos.loc[i,'clade'].name is not None and pos.loc[i,'clade'].name != 'root':
                count = count+1
            pos.loc[i, 'newpos'] = count

        pos.set_index('clade', inplace=True)
        heights = pos.to_dict()['newpos']

    return heights


def plot_tree(input_tree,
              label_func=None,
              title='',
              ax=None,
              output_name=None,
              normal_name='diploid',
              width_scale=1,
              height_scale=1,
              show_branch_lengths=True,
              show_branch_support=False,
              show_events=False,
              branch_labels=None,
              label_colors=None,
              hide_internal_nodes=False,
              marker_size=None,
              line_width=None,
              **kwargs):
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
        max_leaf_to_root_distances = np.max([np.sum([x.branch_length for x in input_tree.get_path(leaf)])
                            for leaf in input_tree.get_terminals()])
        plot_width = 5 + np.max([0, width_scale * np.log10(max_leaf_to_root_distances / 100) * 5])

        # maximum figure size is 250x250 inches
        fig, ax = plt.subplots(figsize=(min(250, plot_width), min(250, plot_height)))

    label_func=label_func if label_func is not None else lambda x: x

    # options for displaying label colors.
    if label_colors is not None:
        if callable(label_colors):
            def get_label_color(label):
                return label_colors(label)
        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")
    else:
        clade_colors = {}
        for sample in [x.name for x in list(input_tree.find_clades(''))]:
            ## determine if sample is terminal
            is_terminal = True
            matches = list(input_tree.find_clades(sample))
            if len(matches) > 0:
                clade = matches[0]
                is_terminal = clade.is_terminal()
            ## determine if sample is normal
            clade_colors[sample] = COL_MARKER_TERMINAL
            if not is_terminal:
                clade_colors[sample] = COL_MARKER_INTERNAL
            if sample == normal_name:
                clade_colors[sample] = COL_MARKER_NORMAL
        
        def get_label_color(label):
            return clade_colors.get(label, "black")

    if marker_size is None:
        marker_size = TREE_MARKER_SIZE
    marker_func=lambda x: (marker_size, get_label_color(x.name)) if x.name is not None else None

    ax.axes.get_yaxis().set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune=None))
    ax.xaxis.set_tick_params(labelsize=XLABEL_TICK_SIZE)
    ax.xaxis.label.set_size(XLABEL_FONT_SIZE)
    ax.set_title(title, x=0.01, y=1.0, ha='left', va='bottom',
                fontweight='bold', fontsize=16, zorder=10)
    x_posns = _get_x_positions(input_tree)
    y_posns = _get_y_positions(input_tree, adjust=not hide_internal_nodes, normal_name=normal_name)

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def value_to_str(value):
        if value is None or value == 0:
            return None
        elif int(value) == value:
            return str(int(value))
        else:
            return str(value)

    if not branch_labels:
        if show_branch_lengths:
            def format_branch_label(x): 
                return value_to_str(np.round(x.branch_length, 1)) if x.name != 'root' and x.name is not None else None
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
            if clade.name == 'root' or clade.name is None:
                return None
            try:
                confidences = clade.confidences
            # phyloXML supports multiple confidences
            except AttributeError:
                pass
            else:
                return "/".join(value_to_str(cnf.value) for cnf in confidences)
            if clade.confidence is not None:
                return value_to_str(clade.confidence)
            return None


    def draw_clade_lines(
        use_linecollection=False,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
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

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add node marker
        if marker_func is not None:
            marker = marker_func(clade)
            if marker is not None and clade is not None and not(hide_internal_nodes and not clade.is_terminal()):
                marker_size, marker_col = marker_func(clade)
                ax.scatter(x_here, y_here, s=marker_size, c=marker_col, zorder=3)

        # Add node/taxon labels
        label = label_func(str(clade.name))
        ax_scale = ax.get_xlim()[1] - ax.get_xlim()[0]

        if label not in (None, clade.__class__.__name__) and \
                not (hide_internal_nodes and not clade.is_terminal()):
            ax.text(
                x_here + min(0.02*ax_scale, 1),
                y_here,
                " %s" % label,
                verticalalignment="center",
                color=get_label_color(label),
            )
        # Add label above the branch
        conf_label = format_branch_label(clade)
        if conf_label:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )
        # Add support below the branch
        if show_branch_support:
            support_value = format_support_value(clade)
            if support_value:
                ax.text(
                    0.5 * (x_start + x_here),
                    y_here + 0.25,
                    support_value + '%',
                    fontsize="small",
                    color='grey',
                    horizontalalignment="center",
                )
        # Add Events list
        if show_events and clade.events is not None:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                clade.events,
                fontsize="small",
                color=COL_MARKER_NORMAL,
                horizontalalignment="center",
            )
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

    if line_width is None:
        line_width = plt.rcParams["lines.linewidth"]
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
    #ax.set_xlim(-0.05 * xmax, 1.25 * xmax)
    ax.set_xlim(-0.05 * xmax, 1.05 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    ax.set_ylim(max(y_posns.values()) + 0.5, 0.5)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if output_name is not None:
        plt.savefig(output_name + ".png", bbox_inches='tight')

    return plt.gcf()


def plot_cn_heatmap(input_df, final_tree=None, y_posns=None, cmax=8, total_copy_numbers=False,
                    alleles=('cn_a', 'cn_b'), tree_width_ratio=1, cbar_width_ratio=0.05, figsize=(20, 10),
                    tree_line_width=0.5, tree_marker_size=0, show_internal_nodes=False, title='',
                    tree_label_colors=None, tree_label_func=None, cmap='coolwarm', normal_name='diploid',
                    ignore_segment_lengths=False):

    input_df = input_df[alleles].copy()
    # if len(np.intersect1d(cur_sample_labels, input_df.index.get_level_values('sample_id').unique())) != len(cur_sample_labels):
    #     raise MEDICCPlotError("tree nodes and labels in dataframe are not the same")

    if not isinstance(alleles, list) and not isinstance(alleles, tuple):
        alleles = [alleles]
    nr_alleles = len(alleles)

    cmax = min(cmax, np.max(input_df[alleles].values.astype(int)))

    if final_tree is None:
        fig, axs = plt.subplots(figsize=figsize, ncols=1+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': nr_alleles*[1] + [cbar_width_ratio]})

        cur_sample_labels = (input_df.index.get_level_values('sample_id').unique())
        if y_posns is None:
            y_posns = {s: i for i, s in enumerate(cur_sample_labels)}

        cn_axes = axs[:-1]
    else:
        fig, axs = plt.subplots(figsize=figsize, ncols=2+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': [tree_width_ratio] + nr_alleles*[1] + [cbar_width_ratio]})
        tree_ax = axs[0]
        cn_axes = axs[1:-1]

        if show_internal_nodes:
            cur_sample_labels = np.array([x.name for x in list(final_tree.find_clades()) if x.name is not None])
        else:
            cur_sample_labels = np.array([x.name for x in final_tree.get_terminals()])

        y_posns = {k.name:v for k, v in _get_y_positions(final_tree, adjust=show_internal_nodes, normal_name=normal_name).items()}
        
        _ = plot_tree(final_tree, ax=tree_ax, normal_name=normal_name,
                      label_func=tree_label_func if tree_label_func is not None else lambda x: '',
                      hide_internal_nodes=(not show_internal_nodes), show_branch_lengths=False, show_events=False,
                      line_width=tree_line_width, marker_size=tree_marker_size,
                      title=title, label_colors=tree_label_colors)
        tree_ax.set_axis_off()
        tree_ax.set_axis_off()
        fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0.0, wspace=100)

    cax = axs[-1]

    ind = [y_posns.get(x, -1) for x in cur_sample_labels]
    cur_sample_labels = cur_sample_labels[np.argsort(ind)]
    color_norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=(2 if total_copy_numbers else 1), vmax=cmax)

    chr_ends = input_df.loc[cur_sample_labels[0]].copy()
    chr_ends['end_pos'] = np.cumsum([1]*len(chr_ends))
    chr_ends = chr_ends.reset_index().groupby('chrom').max()['end_pos']
    chr_ends = chr_ends.dropna()
    chr_ends = chr_ends.astype(int)


    if ignore_segment_lengths:
        x_pos = np.arange(len(input_df.loc[cur_sample_labels].astype(int).unstack('sample_id'))+1)
    else:
        x_pos = np.append([0], np.cumsum(input_df.loc[cur_sample_labels].astype(int).unstack(
            'sample_id').loc[:, (alleles[0])].loc[:, cur_sample_labels].eval('end-start').values))
    y_pos = np.arange(len(cur_sample_labels)+1)+0.5

    for ax, allele in zip(cn_axes, alleles):
        im = ax.pcolormesh(x_pos, y_pos,
                        input_df.loc[cur_sample_labels].astype(int).unstack(
                            'sample_id').loc[:, (allele)].loc[:, cur_sample_labels].values.T,
                        cmap=cmap,
                        norm=color_norm)

        for _, line in chr_ends.items():
            ax.axvline(x_pos[line], color='black', linewidth=0.75)
        
        xtick_pos = np.append([0], x_pos[chr_ends.values][:-1])
        xtick_pos = (xtick_pos + np.roll(xtick_pos, -1))/2
        xtick_pos[-1] += x_pos[-1]/2
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels([x[3:] for x in chr_ends.index], ha='center', rotation=90, va='bottom')
        ax.tick_params(width=0)
        ax.xaxis.set_tick_params(labelbottom=False, labeltop=True, bottom=False)
        ax.set_yticks([])

    cax.pcolormesh([0, 1],
                   np.arange(0, cmax+2),
                   np.arange(0, cmax+1)[:, np.newaxis],
                   cmap=cmap,
                   norm=color_norm)

    cax.set_xticks([])
    cax.set_yticks(np.arange(0, cmax+1)+0.5)
    if np.max(input_df.values.astype(int)) > cmax:
        cax.set_yticklabels([str(x) + '+' if x == cmax else str(x)
                             for x in np.arange(0, cmax+1)], ha='left')
    else:
        cax.set_yticklabels(np.arange(0, cmax+1), ha='left')
    cax.yaxis.set_tick_params(left=False, labelleft=False, labelright=True)

    for ax in axs[:-1]:
        ax.set_ylim(len(cur_sample_labels)+0.5, 0.5)

    return fig


class MEDICCPlotError(Exception):
    pass



def import_tree(tree_file, normal_name='diploid', file_format='newick'):
    """Loads a phylogenetic tree in the given format and roots it at the normal sample. """
    tree = Bio.Phylo.read(tree_file, file_format)
    input_tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)
    tmpsearch = [c for c in input_tree.find_clades(name = normal_name)]
    normal_name = tmpsearch[0]
    root_path = input_tree.get_path(normal_name)[::-1]

    if len(root_path) > 1:
        new_root = root_path[1]
        input_tree.root_with_outgroup(new_root)
    else:
        pass

    # check that internal node names are unique
    node_names = [c.name for c in input_tree.find_clades() if c.name is not None]
    if len(node_names) != len(np.unique(node_names)):
        logger.warning("Internal node names are not unique. This will cause problems in MEDICC2.")

    return input_tree