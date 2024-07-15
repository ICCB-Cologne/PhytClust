import numpy as np
from scipy.ndimage import uniform_filter1d
import numpy.matlib


def find_peaks_ext(
    scores=None,
    num_terminals=None,
    n=3,
    k_start=0,
    k_end=None,
):
    """
    Finds and plots the top n maxima in the scores within a specified range of k.

    Parameters:
    - n: int, the number of top peaks to find.
    - plot: bool, whether to plot the score graph with the peaks highlighted.
    - method: str, 'default' or 'alt' for different peak finding methods.
    - prominence: tuple or None, the prominence parameter for peak finding, applicable in 'alt' method.
    - k_start: int or None, the starting index for the range of k values to consider.
    - k_end: int or None, the ending index for the range of k values to consider.
    """
    if len(scores) == 0:
        print("Please calculate scores first")
        return []

    # Adjust k_start for 0-based indexing within Python
    k_start = max(k_start - 1, 0)
    k_end = k_end or len(scores)

    if k_end > len(scores) or k_end < k_start:
        raise ValueError(f"Invalid k_end value. Max allowed value is {len(scores)}")

    scores_subset = scores[k_start:k_end]
    scores_subset = np.where(np.isinf(scores_subset), np.nan, scores_subset)

    #####
    smoothing = max(1, int(num_terminals * 0.1))
    edges, prominence_edges, prominences = find_plateau_edges(
        scores_subset, k_start, smoothing)

    if len(edges) == 1:
        top_peaks = list(peak + k_start + 1 for peak in edges)
        print(f"Found only {len(top_peaks)} peak(s)")
        # if plot:
        #     plot_scores(scores_subset, list(edges), k_start)

        return prominences, top_peaks

    representatives, representative_plateau_sizes = select_representative_edges(
        edges, prominence_edges
    )

    # Sort representatives and their plateau sizes by plateau size
    # sorted_indices = np.argsort(representative_plateau_sizes)[::-1]
    representatives = np.array(representatives)  # [sorted_indices]
    representative_plateau_sizes = np.array(
        representative_plateau_sizes
    )  # [sorted_indices]

    top_peaks = list(representatives[: min(n, len(representatives))])
    top_peak_plateau_sizes = list(
        representative_plateau_sizes[: min(n, len(representative_plateau_sizes))]
    )
    top_peak_plateau_sizes = top_peak_plateau_sizes

    # Adjust top_peaks to the original indexing context
    top_peaks = list(peak + k_start + 1 for peak in top_peaks)

    if len(top_peaks) < n:
        print(f"Found only {len(top_peaks)} peak(s)")

    return [peak for peak in top_peaks]

#tbc
def rank_peaks_ext():
    normalized_prominences = normalize(peak_prominence)
    normalized_prominences = np.array(normalized_prominences)
    x = np.arange(len(self.den)).reshape(-1, 1)
    y = self.den
    model = LinearRegression()
    model.fit(x, y)
    predictions = model.predict(x)
    r_squared = r2_score(y, predictions)

    indices = [peak - 1 for peak in self.top_peaks]
    max_product_index = elbow_point(y)
    distances = [abs(i - max_product_index) for i in indices]

    normalized_distances = distances / np.max(distances)
    reciprocal_distances = 1 / (1 + normalized_distances)

    scores = ((1 - r_squared) * (reciprocal_distances)) + (
        (r_squared) * normalized_prominences[indices]
    )

    sorted_indices = np.argsort(scores)[::-1]
    indices = np.array(indices)
    # Use the sorted indices to sort the scores and indices
    scores = scores[sorted_indices]
    indices = indices[sorted_indices]

    # Calculate the rankings
    rankings = np.arange(1, len(scores) + 1)

    # Create a list of tuples (ranking, index, score)
    ranked_data = [
        (rank, (idx + 1), score)
        for rank, idx, score in zip(rankings, indices, scores)
    ]

    self.ranked_peaks = ranked_data


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

    # calculate smooth gradients
    smoothF = uniform_filter1d(F, size=smoothing)
    dF = uniform_filter1d(np.gradient(smoothF), size=smoothing)
    d2F = uniform_filter1d(np.gradient(dF), size=smoothing)
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

    # Sort edges and plateau_sizes by relative plateau size
    sorted_indices = np.argsort(prominence_edges)[::-1]
    edges = np.array(edges)[sorted_indices]
    plateau_sizes = np.array(prominence_edges)[sorted_indices]
    # print(edges, plateau_sizes)

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
        threshold = max(1, edge * 0.1)  # Threshold is 1% of the edge index
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
