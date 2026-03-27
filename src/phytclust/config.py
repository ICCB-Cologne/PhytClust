"""Configuration for PhytClust visualization and algorithm parameters."""

# Visualization parameters
VISUALIZATION_CONFIG = {
    # Plotting parameters
    "plot": {
        "default_dpi": 150,
        "default_figsize": (10, 8),
    },
    # Score plot parameters
    "scores_plot": {
        "title_fontsize": 50,
        "axis_label_fontsize": 35,
        "tick_labelsize": 30,
        "peak_labelsize": 50,
        "bin_labelsize": 30,
        "colorblind_palette": [
            "#0072B2",
            "#009E73",
            "#D55E00",
            "#56B4E9",
            "#E69F00",
            "#CC79A7",
            "#F0E442",
            "#000000",
        ],
    },
    # Cluster plot parameters
    "cluster_plot": {
        "default_width_scale": 2.0,
        "default_height_scale": 0.1,
        "marker_size": 40,
        "default_cmap": "tab20",
    },
    # Tree plot parameters
    "tree_plot": {
        "marker_terminal": "green",
        "marker_internal": "blue",
        "marker_normal": "red",
        "marker_size": 8,
        "segment_linewidth": 1.0,
        "default_title_fontsize": 16,
        "default_label_fontsize": 14,
        "xlabel_fontsize": 12,
        "xlabel_tick_fontsize": 10,
    },
}

# Algorithm parameters
ALGORITHM_CONFIG = {
    # DP parameters
    "dp": {
        "dtype_threshold": 80_000,
        "int16_max": 32767,
    },
    # Scoring parameters
    "scoring": {
        "epsilon": 1e-12,  # threshold for zero detection
    },
    # Core clustering parameters
    "core": {
        "default_min_cluster_size": 1,
        "default_max_k_limit": 0.9,
        "default_support_weight": 1.0,
        "default_support_min": 0.05,
    },
}


def get_visualization_param(section: str, param: str, default=None):
    """
    Get a visualization parameter.

    Parameters
    ----------
    section : str
        Configuration section (e.g., 'scores_plot', 'cluster_plot').
    param : str
        Parameter name.
    default : optional
        Default value if parameter not found.

    Returns
    -------
    Value of the parameter or default if not found.
    """
    return VISUALIZATION_CONFIG.get(section, {}).get(param, default)


def get_algorithm_param(section: str, param: str, default=None):
    """
    Get an algorithm parameter.

    Parameters
    ----------
    section : str
        Configuration section (e.g., 'dp', 'scoring').
    param : str
        Parameter name.
    default : optional
        Default value if parameter not found.

    Returns
    -------
    Value of the parameter or default if not found.
    """
    return ALGORITHM_CONFIG.get(section, {}).get(param, default)
