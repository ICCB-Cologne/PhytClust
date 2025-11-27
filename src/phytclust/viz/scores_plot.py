from typing import Optional, List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatter, MaxNLocator


def plot_scores(
    pc,
    scores_subset: Optional[np.ndarray] = None,
    k_start: int = 1,
    k_end: Optional[int] = None,
    peaks: Optional[List[int]] = None,
    resolution_on: bool = False,
    num_bins: int = 3,
    fig_width: int = 18,
    fig_height: int = 10,
    log_scale_y: bool = False,
    x_axis_mode: str = "log",
    log_base: Optional[float] = None,
) -> plt.Figure:
    """Split-out of your _plot_scores; logic unchanged except pulling from pc where needed."""
    if scores_subset is None:
        scores_subset = pc.scores
    if len(scores_subset) == 0:
        raise ValueError("Scores are empty, please calculate scores first.")

    title_fontsize = 50
    axis_label_fontsize = 35
    tick_labelsize = 30
    peak_labelsize = 50
    bin_labelsize = 30

    if k_end is None:
        k_end = k_start + len(scores_subset)
    start_idx = k_start - 1
    end_idx = k_end - 1 if k_end is not None else None

    scores_slice = scores_subset[start_idx:end_idx]
    if scores_slice.size == 0:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.text(0.5, 0.5, "No scores to plot", ha="center", va="center")
        ax.set_axis_off()
        plt.tight_layout()
        return fig
    x_indices = np.arange(k_start, k_start + len(scores_slice))
    data_min, data_max = np.nanmin(scores_slice), np.nanmax(scores_slice)
    data_range = data_max - data_min or 1.0

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.plot(x_indices, scores_slice, "o-", label="Scores", markersize=2, linewidth=7)

    if peaks:
        valid_peaks, valid_scores = [], []
        for p in peaks:
            idx = p - k_start
            if 0 <= idx < len(scores_slice):
                valid_peaks.append(p)
                valid_scores.append(scores_slice[idx])

        if valid_peaks:
            ax.plot(valid_peaks, valid_scores, "rx", label="Peaks", markersize=8)
            offset_y = 0.02 * data_range if data_range else 0.02
            for px, py in zip(valid_peaks, valid_scores):
                ax.text(
                    px,
                    py + offset_y,
                    str(px),
                    fontsize=peak_labelsize,
                    color="red",
                    ha="center",
                    va="bottom",
                )

        if resolution_on and num_bins > 0:
            bin_ranges = getattr(pc, "bin_ranges_current", None)
            if bin_ranges is None:
                from ..algo.bins import define_bins

                bin_ranges = define_bins(pc, num_bins, k_lo=k_start, k_hi=k_end)

            boundaries = [
                (bin_ranges[i][1] + bin_ranges[i + 1][0]) / 2
                for i in range(len(bin_ranges) - 1)
            ]
            cont_bins = []
            cont_bins.append((bin_ranges[0][0], boundaries[0]))
            for i in range(1, len(bin_ranges) - 1):
                cont_bins.append((boundaries[i - 1], boundaries[i]))
            cont_bins.append((boundaries[-1], bin_ranges[-1][1]))

            cb_palette = [
                "#0072B2",
                "#009E73",
                "#D55E00",
                "#56B4E9",
                "#E69F00",
                "#CC79A7",
                "#F0E442",
                "#000000",
            ]
            colours = (cb_palette * (num_bins // len(cb_palette) + 1))[:num_bins]

            mid_pts, xtick_labels = [], []
            y_min, y_max = ax.get_ylim()
            y_offset = 0.04 * (y_max - y_min)
            use_log_mid = x_axis_mode == "log"

            for i, ((lo_nom, hi_nom), (lo, hi), colour) in enumerate(
                zip(bin_ranges, cont_bins, colours), 1
            ):
                ax.axvspan(lo, hi, color=colour, alpha=0.20, zorder=-1)
                mid = (lo_nom * hi_nom) ** 0.5 if use_log_mid else (lo_nom + hi_nom) / 2
                mid_pts.append(mid)
                xtick_labels.append(f"[{lo_nom}-{hi_nom}]")
                ax.text(
                    mid,
                    y_max - y_offset,
                    f"CL{i}",
                    ha="center",
                    va="bottom",
                    fontsize=20,
                    weight="bold",
                    color=colour,
                )

            ax.set_xticks(mid_pts)
            # ax.set_xticklabels(xtick_labels, fontsize=tick_labelsize)
            import matplotlib.ticker as mticker

            ax.xaxis.set_major_locator(mticker.FixedLocator(mid_pts))
            ax.xaxis.set_minor_locator(mticker.NullLocator())

    x_min_plot = max(2, x_indices[0])
    x_max_plot = x_indices[-1] + 0.05 * (x_indices[-1] - x_min_plot)
    ax.set_xlim(x_min_plot, x_max_plot)

    if x_axis_mode == "log":
        ax.set_xscale("log")
        import matplotlib.ticker as mticker

        if resolution_on:
            ax.xaxis.set_major_locator(mticker.FixedLocator(mid_pts))
            ax.xaxis.set_minor_locator(mticker.NullLocator())
            # ax.set_xticklabels(xtick_labels, fontsize=tick_labelsize)
        else:
            ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
            ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs="auto", numticks=10))
            ax.xaxis.set_major_formatter(LogFormatter(base=10.0, labelOnlyBase=False))
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha("right")
                # label.set_fontsize(tick_labelsize)
        x_label_str = "No. of Clusters (log)"

    elif x_axis_mode == "linear":
        ax.set_xscale("linear")
        if resolution_on:
            ax.set_xticks(mid_pts)
            # ax.set_xticklabels(xtick_labels, fontsize=tick_labelsize)
        else:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            # for label in ax.get_xticklabels():
            #     label.set_fontsize(tick_labelsize)
        x_label_str = "No. of Clusters"

    else:
        raise ValueError("x_axis_mode must be either 'log' or 'linear'.")

    if log_scale_y:
        ax.set_yscale("log")
        y_label_str = "Scores (log)"
    else:
        y_label_str = "Scores"
    # for label in ax.get_yticklabels():
    #     label.set_fontsize(tick_labelsize)

    ax.tick_params(axis="both", which="both", labelsize=tick_labelsize)
    ax.set_xlabel(x_label_str, fontsize=axis_label_fontsize)
    ax.set_ylabel(y_label_str, fontsize=axis_label_fontsize)
    ax.set_title("Scores", fontsize=title_fontsize, pad=40)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    y_min = np.min(scores_slice)
    y_max = np.max(scores_slice)
    padding = (y_max - y_min) * 0.2
    ax.set_ylim(y_min - padding, y_max + padding)

    plt.tight_layout()
    return fig
