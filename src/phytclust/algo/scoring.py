from typing import Optional, Dict, Any, List
import numpy as np

from ..viz.scores_plot import plot_scores as _plot_scores
from .bins import define_bins as _define_bins


def _single_cluster_score(
    pc, clusters: Optional[Dict[Any, Any]] = None, k: Optional[int] = None
):
    if not pc.max_k or pc.max_k <= 0:
        raise ValueError("max_k must be set and positive to compute cluster scores.")

    active_tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    root = active_tree.root
    root_id = pc.node_to_id[root]
    dp_row = pc.dp_table[root_id]

    if dp_row is None:
        raise ValueError("Root DP row missing (dp_table[root] is None).")

    dp_row = np.asarray(dp_row, dtype=float)
    pc.beta_1 = dp_row[0]
    num_terminals = pc.num_terminals

    if clusters is not None:
        num_clusters = len(set(clusters.values()))
        if num_clusters < 1 or num_clusters > pc.max_k:
            return (float("inf"), float("inf"), float("inf"))
        if num_clusters - 1 >= len(dp_row):
            return (float("inf"), float("inf"), float("inf"))
        beta = dp_row[num_clusters - 1]

    elif k is not None:
        if k < 1 or k > pc.max_k or k - 1 >= len(dp_row):
            return (float("inf"), float("inf"), float("inf"))
        num_clusters = k
        beta = dp_row[k - 1]

    else:
        raise ValueError(
            "Either 'clusters' or 'k' must be provided to compute the score."
        )

    if np.isinf(beta):
        return (beta, float("inf"), 0.0)

    if beta == 0:
        return (beta, float("inf"), 0.0)

    beta_ratios = (pc.beta_1 - beta) / beta
    norm_ratios = (
        (num_terminals - num_clusters) / float(num_clusters)
        if num_clusters
        else float("inf")
    )

    if not np.isfinite(beta_ratios) or not np.isfinite(norm_ratios):
        score = float("inf")
    else:
        score = beta_ratios * norm_ratios

    return (beta, beta_ratios, score)


def calculate_scores(pc, plot: bool = False) -> None:
    """Unified score calculation for (single k | preset clusters | full range)."""
    results = []

    if pc.k is not None:
        from .dp import cluster_map

        cmap = cluster_map(pc, pc.k)
        if cmap is None:
            pc.get_clusters(pc.k)
            from .dp import cluster_map as _cm

            cmap = _cm(pc, pc.k)
        results.append(_single_cluster_score(pc, clusters=cmap))

    elif pc.clusters:
        results.extend(
            _single_cluster_score(pc, clusters=cm) for cm in pc.clusters.values()
        )

    else:
        if not pc.max_k or pc.max_k <= 0:
            raise ValueError(
                "max_k must be set and positive to compute DP-based scores."
            )
        results.extend(
            _single_cluster_score(pc, k=k_val) for k_val in range(1, pc.max_k + 1)
        )

    beta_values, den_list, scores = map(np.array, zip(*results))

    scores[scores < 0] = 0
    beta_values[beta_values < 0] = 0
    beta_values = np.nan_to_num(beta_values, nan=0.0, posinf=0.0, neginf=0.0)

    elbow_scores = [
        (
            (beta_values[i - 1] - beta_values[i])
            / (beta_values[i] - beta_values[i + 1])
            if i > 0 and (beta_values[i] - beta_values[i + 1]) != 0
            else 0.0
        )
        for i in range(len(beta_values) - 1)
    ]
    elbow_scores.append(0.0)
    elbow_scores = np.array(elbow_scores, dtype=float)

    invalid_mask = (
        np.isnan(scores)
        | np.isinf(scores)
        | np.isnan(elbow_scores)
        | np.isinf(elbow_scores)
    )
    valid_mask = ~invalid_mask

    scores_valid = scores[valid_mask]
    beta_valid = beta_values[valid_mask]
    den_valid = den_list[valid_mask]
    elbow_valid = elbow_scores[valid_mask]

    if len(scores_valid) == 0:
        pc.scores = np.array([], dtype=float)
        pc.beta_values = np.array([], dtype=float)
        pc.norm_ratios = np.array([], dtype=float)
        return

    combined_scores = np.nan_to_num(
        elbow_valid * scores_valid, nan=0.0, posinf=0.0, neginf=0.0
    )

    eps = 1e-12
    nonzero_idx = np.where(np.abs(combined_scores) > eps)[0]
    if nonzero_idx.size > 0:
        last_useful = nonzero_idx[-1] + 1
    else:
        last_useful = len(combined_scores)

    combined_scores = combined_scores[:last_useful]
    beta_valid = beta_valid[:last_useful]
    den_valid = den_valid[:last_useful]

    pc.scores = combined_scores
    pc.beta_values = beta_valid
    pc.norm_ratios = den_valid

    if plot:
        pass


def find_score_peaks(
    pc,
    scores: Optional[np.ndarray] = None,
    global_peaks: int = 3,
    peaks_per_bin: int = 1,
    resolution_on: bool = False,
    num_bins: int = 3,
    min_k: int = 2,
    k_start: Optional[int] = None,
    k_end: Optional[int] = None,
    plot: bool = True,
    smooth_window_size: int = 1,
    min_prominence: float = 1e-3,
    ranking_mode: str = "adjusted",
    alpha=0.7,
) -> List[int]:
    """Direct split of your original _find_score_peaks (unchanged logic)."""
    import numpy as np
    from scipy.signal import find_peaks

    if scores is None:
        scores = pc.scores
    if scores is None or len(scores) == 0:
        pc.peaks_by_rank = []
        pc.resolution_info = None
        pc.peaks_by_resolution = None
        if plot:
            # optional: draw a minimal placeholder or skip entirely
            pc.plot_of_scores = _plot_scores(
                pc,
                scores_subset=np.array([], dtype=float),
                peaks=[],
                k_start=1,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )
        return []

    k_start = k_start if k_start is not None else 1
    k_end = k_end if k_end is not None else len(scores)
    if not (0 <= k_start < len(scores)):
        raise ValueError(f"k_start must be between 0 and {len(scores)-1}.")
    if not (k_start < k_end <= len(scores)):
        raise ValueError(f"k_end must be <= {len(scores)} and > k_start.")

    scores_subset = scores[k_start:k_end].astype(float)
    scores_subset = np.where(np.isinf(scores_subset), np.nan, scores_subset)

    if smooth_window_size > 1:
        kernel = np.ones(smooth_window_size) / smooth_window_size
        scores_subset_smoothed = np.convolve(scores_subset, kernel, mode="same")
    else:
        scores_subset_smoothed = scores_subset

    if len(scores_subset_smoothed) > 1 and scores_subset_smoothed[-1] <= 0:
        log_scores = np.log(scores_subset_smoothed[:-1] + 1e-10)
    else:
        log_scores = np.log(scores_subset_smoothed + 1e-10)

    peaks, props = find_peaks(log_scores, prominence=min_prominence)
    prominences = props["prominences"]

    peak_data = []
    for i, peak_idx in enumerate(peaks):
        pk = peak_idx + k_start + 1
        if pk == 2:
            continue
        if pk >= min_k:
            orig_score = scores[pk - 1]
            prom = prominences[i]
            peak_data.append((pk, prom, orig_score))

    if len(peak_data) == 0:
        pc.peaks_by_rank = []
        pc.resolution_info = None
        pc.peaks_by_resolution = None
        if plot:
            pc.plot_of_scores = _plot_scores(
                pc,
                scores_subset=scores_subset,
                peaks=[],
                k_start=k_start,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )
        return pc.peaks_by_rank

    if not resolution_on:
        pc.resolution_info = None
        pc.peaks_by_resolution = None

        if ranking_mode == "raw":
            peak_data.sort(key=lambda x: x[1], reverse=True)
        elif ranking_mode == "adjusted":
            all_prom = [x[1] for x in peak_data]
            all_sc = [x[2] for x in peak_data]
            prom_min, prom_max = min(all_prom), max(all_prom)
            score_min, score_max = min(all_sc), max(all_sc)
            adjusted_data = []
            for pk, prom, sc in peak_data:
                prom_norm = (
                    ((prom - prom_min) / (prom_max - prom_min))
                    if prom_max > prom_min
                    else 1.0
                )
                score_norm = (
                    ((sc - score_min) / (score_max - score_min))
                    if score_max > score_min
                    else 1.0
                )
                combined_metric = alpha * prom_norm + (1 - alpha) * score_norm
                adjusted_data.append((pk, prom, sc, combined_metric))
            adjusted_data.sort(key=lambda x: x[3], reverse=True)
            peak_data = adjusted_data
        else:
            raise ValueError("ranking_mode must be either 'raw' or 'adjusted'.")

        chosen = peak_data[:global_peaks]
        # chosen_sorted = sorted(chosen, key=lambda x: x[0])
        final_peaks = [int(x[0]) for x in chosen]
        pc.peaks_by_rank = final_peaks

    else:
        bin_ranges = _define_bins(pc, num_bins=num_bins)
        pc.bin_ranges_current = bin_ranges
        pc.resolution_info = {}
        pc.peaks_by_resolution = {}

        all_prom = [x[1] for x in peak_data]
        all_sc = [x[2] for x in peak_data]
        prom_min, prom_max = min(all_prom), max(all_prom)
        score_min, score_max = min(all_sc), max(all_sc)

        normed_data = []
        for pk, prom, sc in peak_data:
            prom_norm = (
                ((prom - prom_min) / (prom_max - prom_min))
                if prom_max > prom_min
                else 1.0
            )
            score_norm = (
                ((sc - score_min) / (score_max - score_min))
                if score_max > score_min
                else 1.0
            )
            combined_metric = (alpha) * prom_norm + (1 - alpha) * score_norm
            normed_data.append((pk, prom, sc, combined_metric))

        normed_data.sort(key=lambda x: x[3], reverse=True)

        all_picked_peaks = []
        for i, (start_k, end_k) in enumerate(bin_ranges, start=1):
            bin_label = f"Bin {i}: {start_k}-{end_k}"
            candidates = [
                (pk, pm, sc, combo)
                for (pk, pm, sc, combo) in normed_data
                if start_k <= pk <= end_k
            ]
            chosen = candidates[:peaks_per_bin]
            chosen_kvals = [x[0] for x in chosen]
            pc.resolution_info[bin_label] = chosen
            pc.peaks_by_resolution[bin_label] = chosen_kvals
            all_picked_peaks.extend(chosen_kvals)

        final_peaks = sorted({int(pk) for pk in all_picked_peaks})
        pc.peaks_by_rank = final_peaks

    if plot:
        pc.plot_of_scores = _plot_scores(
            pc,
            peaks=pc.peaks_by_rank,
            k_start=k_start,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )

    return pc.peaks_by_rank
