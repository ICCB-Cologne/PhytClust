import numpy as np
from scipy.ndimage import uniform_filter1d


def find_plateau_edges(F, smoothing=4):
    """
    Finds points where score increases and then remains constant or decreases.

    Parameters
    ----------
    F : Signal.
    smoothing: Size of uniform filter 1D applied to F and its derivatives.

    Returns
    -------
    edges: array of indices where the conditions are met, sorted by relative plateau size.
    plateau_sizes: array of relative plateau sizes corresponding to the edges.
    dF: (smoothed) derivative of F
    d2F: (smoothed) Second Derivative of F
    """
    F = np.nan_to_num(F)

    # calculate smooth gradients
    smoothF = uniform_filter1d(F, size=smoothing)
    dF = uniform_filter1d(np.gradient(smoothF), size=smoothing)
    d2F = uniform_filter1d(np.gradient(dF), size=smoothing)

    # Find points where score increases and then remains constant or decreases
    edges = []
    plateau_sizes = []
    for k in range(1, len(F) - 1):
        if (
            np.isfinite(dF[k])
            and np.isfinite(d2F[k])
            and F[k - 1] < F[k]
            and F[k] >= F[k + 1]
        ):
            plateau_size = (abs(F[k] - F[k - 1]) + abs(F[k] - F[k + 1])) / F[k]
            edges.append(k)
            plateau_sizes.append(plateau_size)

    # Sort edges and plateau_sizes by relative plateau size
    sorted_indices = np.argsort(plateau_sizes)[::-1]
    edges = np.array(edges)[sorted_indices]
    plateau_sizes = np.array(plateau_sizes)[sorted_indices]

    return edges, plateau_sizes, dF, d2F


def select_representative_edges(edges, plateau_sizes):
    """
    Selects representative edges with the largest plateau size from groups of close edges.

    Parameters
    ----------
    edges : List of edge indices.
    plateau_sizes: List of plateau sizes corresponding to each edge.

    Returns
    -------
    representatives: List of representative edge indices.
    representative_plateau_sizes: List of plateau sizes corresponding to the representatives.
    """
    # Sort edges and plateau_sizes by edge
    sorted_indices = np.argsort(edges)
    edges = np.array(edges)[sorted_indices]
    plateau_sizes = np.array(plateau_sizes)[sorted_indices]

    # Initialize groups and representatives
    groups = [[edges[0]]]
    group_plateau_sizes = [[plateau_sizes[0]]]
    representatives = []
    representative_plateau_sizes = []

    # Group close edges
    for edge, plateau_size in zip(edges[1:], plateau_sizes[1:]):
        threshold = max(1, edge * 0.01)  # Threshold is 1% of the edge index
        if edge - groups[-1][-1] <= threshold:
            groups[-1].append(edge)
            group_plateau_sizes[-1].append(plateau_size)
        else:
            groups.append([edge])
            group_plateau_sizes.append([plateau_size])

    # Select representative from each group
    for group, group_plateau_size in zip(groups, group_plateau_sizes):
        # Here we select the edge with the largest plateau size as the representative
        representative_index = np.argmax(group_plateau_size)
        representatives.append(group[representative_index])
        representative_plateau_sizes.append(group_plateau_size[representative_index])

    # Pair each representative with its plateau size
    pairs = list(zip(representatives, representative_plateau_sizes))

    # Sort the pairs by plateau size
    pairs.sort(key=lambda pair: pair[1], reverse=True)

    # Unzip the pairs back into two lists
    representatives, representative_plateau_sizes = zip(*pairs)

    return list(representatives), list(representative_plateau_sizes)