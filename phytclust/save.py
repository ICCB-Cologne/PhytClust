import numpy as np
import pandas as pd
import os
import logging
from typing import Any, Dict, List, Optional, Tuple, Union
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def save_clusters(
    scores: Union[List[float], np.ndarray],
    clusters: Dict[int, Dict[Any, int]],
    results_dir: str,
    top_n: int = 1,
    filename: str = "phyclust_results.csv",
    outlier: bool = True,
    selected_peaks: Optional[List[int]] = None,
    n: Optional[int] = None,
    output_all: bool = False,
) -> None:
    """
    Save clustering results to a CSV file.

    Parameters:
    scores (list or np.ndarray): The scores used to identify peaks.
    clusters (dict): A dictionary where keys are cluster indices and values are dictionaries of node-to-cluster mappings.
    results_dir (str): The directory where the results file will be saved.
    top_n (int, optional): The number of top peaks to consider. Defaults to 1.
    filename (str, optional): The name of the output CSV file. Defaults to "phyclust_results.csv".
    outlier (bool, optional): Whether to mark outliers with a cluster ID of -1. Defaults to True.
    selected_peaks (list, optional): Pre-selected peaks to use instead of finding peaks. Defaults to None.
    n (int, optional): Specific cluster index to save. Defaults to None.
    output_all (bool, optional): Whether to output all solutions. Defaults to False.

    Returns:
    None
    """
    if scores is None or (isinstance(scores, np.ndarray) and scores.size == 0):
        logger.error("Please calculate scores first.")
        return

    top_peaks = selected_peaks
    if top_peaks is None:
        logger.error("No peaks found.")
        return

    data = []
    if output_all:
        for k in range(1, len(clusters) + 1):
            cluster = clusters[k - 1]
            for node, cluster_id in cluster.items():
                data.append(
                    {
                        "Node Name": node.name,
                        "Number of Clusters": k,
                        "Cluster Number": cluster_id,
                    }
                )
    else:
        selected_clusters = (
            {n: clusters[n]}
            if n is not None and n in clusters
            else [clusters[i] for i in top_peaks]
        )
        for i, cluster in enumerate(selected_clusters, start=1):
            for node, cluster_id in cluster.items():
                data.append(
                    {
                        "Node Name": node.name,
                        "Number of Clusters": top_peaks[i - 1] if selected_peaks else n,
                        "Cluster Number": cluster_id,
                    }
                )

    df = pd.DataFrame(data)
    df_pivot = df.pivot(
        index="Node Name", columns="Number of Clusters", values="Cluster Number"
    )
    df_pivot.columns = [f"{col}_clusters" for col in df_pivot.columns]
    df_pivot.reset_index(inplace=True)
    df_pivot.to_csv(os.path.join(results_dir, filename), index=False)


def identify_clusters(cluster_ids: List[int], outlier: bool) -> Tuple[set, set]:
    """
    Identify outlier and non-outlier clusters.

    Parameters:
    cluster_ids (List[int]): List of cluster IDs.
    outlier (bool): Whether to identify outliers.

    Returns:
    Tuple[set, set]: A tuple containing sets of outlier and non-outlier clusters.
    """
    if outlier:
        outlier_clusters = {
            cluster_id
            for cluster_id in cluster_ids
            if cluster_ids.count(cluster_id) == 1
        }
    else:
        outlier_clusters = set()
    non_outlier_clusters = (
        set(cluster_ids) - outlier_clusters if outlier else set(cluster_ids)
    )
    return outlier_clusters, non_outlier_clusters


def save_to_csv(
    data: List[Dict[str, Union[str, int]]], results_dir: str, filename: str
) -> None:
    """
    Save data to a CSV file.

    Parameters:
    data (List[Dict[str, Union[str, int]]]): The data to be saved.
    results_dir (str): The directory where the results file will be saved.
    filename (str): The name of the output CSV file.

    Returns:
    None
    """
    df = pd.DataFrame(data)
    full_path = os.path.join(results_dir, filename)
    df.to_csv(full_path, index=False)
    logger.info(f"Cluster data saved to {filename}")
