import numpy as np
from scipy.ndimage import uniform_filter1d
import numpy.matlib

def calculate_prominence(scores):
    prominences = []
    for i in range(len(scores)):
        if i == 0:
            left_min = scores[i]
        else:
            left_min = min(scores[:i])

        if i == len(scores) - 1:
            right_min = scores[i]
        else:
            right_min = min(scores[i + 1 :])

        prominence = scores[i] - max(left_min, right_min)
        prominences.append(prominence)
    return prominences


def find_plateau_edges(scores_subset, k_start, smoothing=4):
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
    F = scores_subset
    F = np.nan_to_num(F, nan=0, posinf=0, neginf=0)  # replace nan and inf values with 0

    smoothF = uniform_filter1d(F, size=smoothing)
    dF = uniform_filter1d(np.gradient(smoothF), size=smoothing)#changed smoothF to F
    d2F = uniform_filter1d(np.gradient(dF), size=smoothing)
    # dF = np.gradient(F)
    # d2F = np.gradient(dF)
    last_valid_index = len(scores_subset) - np.isnan(scores_subset[::-1]).argmax() - 2

    prominences = calculate_prominence(F)
    # Find points where score increases and then remains constant or decreases
    edges = []
    prominence_edges = []

    for k in range(1, len(F) - 1):
        if (
            np.isfinite(dF[k])
            and np.isfinite(d2F[k])
            and F[k - 1] < F[k]
            and F[k] >= F[k + 1]
        ):
            if k == 1 or k == last_valid_index:
                prominence_edge = 0
            else:
                prominence_edge = prominences[k]
            edges.append(k)
            prominence_edges.append(prominence_edge)

    sorted_indices = np.argsort(prominence_edges)[::-1]
    edges = np.array(edges)[sorted_indices]
    plateau_sizes = np.array(prominence_edges)[sorted_indices]

    return edges, prominence_edges, prominences


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
    sorted_indices = np.argsort(edges)
    edges = np.array(edges)[sorted_indices]
    plateau_sizes = np.array(plateau_sizes)[sorted_indices]

    groups = [[edges[0]]]
    group_plateau_sizes = [[plateau_sizes[0]]]
    representatives = []
    representative_plateau_sizes = []

    for edge, plateau_size in zip(edges[1:], plateau_sizes[1:]):
        threshold = max(1, edge * 0.1)  # Threshold is 1% of the edge index
        if edge - groups[-1][-1] <= threshold:
            groups[-1].append(edge)
            group_plateau_sizes[-1].append(plateau_size)
        else:
            groups.append([edge])
            group_plateau_sizes.append([plateau_size])

    for group, group_plateau_size in zip(groups, group_plateau_sizes):
        representative_index = np.argmax(group_plateau_size)
        representatives.append(group[representative_index])
        representative_plateau_sizes.append(group_plateau_size[representative_index])

    pairs = list(zip(representatives, representative_plateau_sizes))
    pairs.sort(key=lambda pair: pair[1], reverse=True)
    representatives, representative_plateau_sizes = zip(*pairs)

    return list(representatives), list(representative_plateau_sizes)


def normalize(data):
    min_val = min(data)
    max_val = max(data)
    return [(x - min_val) / (max_val - min_val) for x in data]


def elbow_point(data):
    curve = data
    nPoints = len(curve)
    allCoord = np.vstack((range(nPoints), curve)).T
    np.array([range(nPoints), curve])

    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(
        vecFromFirst * np.matlib.repmat(lineVecNorm, nPoints, 1), axis=1
    )
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine**2, axis=1))
    idxOfBestPoint = np.argmax(distToLine)
    return idxOfBestPoint
