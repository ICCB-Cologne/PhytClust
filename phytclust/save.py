import numpy as np
import pandas as pd
import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def save_clusters(
    scores, clusters, results_dir, top_n=1, filename="phyclust_results.csv", outlier=True, selected_peaks=None, n=None
):
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

    selected_clusters = {n: clusters[n]} if n is not None and n in clusters else [clusters[i] for i in top_peaks]

    data = []
    all_cluster_ids = set()
    for cluster in selected_clusters:
        cluster_ids = list(cluster.values())
        outlier_clusters, non_outlier_clusters = identify_clusters(cluster_ids, outlier)

        new_cluster_id_map = {old_id: new_id for new_id, old_id in enumerate(sorted(non_outlier_clusters), start=1)}
        data.extend(
            {
                "Node Name": node.name,
                "Cluster Number": (
                    -1
                    if outlier and cluster_id in outlier_clusters
                    else new_cluster_id_map.get(cluster_id, -1)
                ),
            }
            for node, cluster_id in cluster.items()
        )
        # for node, cluster_id in cluster.items():
        #     mapped_cluster_id = -1 if outlier and cluster_id in outlier_clusters else new_cluster_id_map.get(cluster_id, -1)
        #     all_cluster_ids.add(mapped_cluster_id)
        #     data.append({"Node Name": node.name, "Cluster Number": mapped_cluster_id})
    save_to_csv(data, results_dir, filename)


def identify_clusters(cluster_ids, outlier):
    """Identify outlier and non-outlier clusters."""
    if outlier:
        outlier_clusters = {cluster_id for cluster_id in cluster_ids if cluster_ids.count(cluster_id) == 1}
    else:
        outlier_clusters = set()
    non_outlier_clusters = set(cluster_ids) - outlier_clusters if outlier else set(cluster_ids)
    return outlier_clusters, non_outlier_clusters


def save_to_csv(data, results_dir, filename):
    """Save data to a CSV file."""
    df = pd.DataFrame(data)
    full_path = os.path.join(results_dir, filename)
    df.to_csv(full_path, index=False)
    logger.info(f"Cluster data saved to {filename}")














    # if selected_peaks is not None:
    #     top_peaks = selected_peaks
    # else:
    #     peaks, _ = find_peaks(scores, distance=1)
    #     if peaks.size < 1:
    #         print("No peaks found.")
    #         return
    #     top_peaks = peaks[np.argsort(-scores[peaks])][:top_n]

    # if n is not None:
    #     selected_clusters = {n: clusters[n]} if n in clusters else {}
    # else:
    #     selected_clusters = [clusters[i] for i in top_peaks]

    # data = []
    # all_cluster_ids = set()
    # for cluster in selected_clusters:
    #     cluster_ids = list(cluster.values())

    #     if outlier:
    #         outlier_clusters = {
    #             cluster_id
    #             for cluster_id in cluster_ids
    #             if cluster_ids.count(cluster_id) == 1
    #         }
    #     else:
    #         outlier_clusters = set()

    #     non_outlier_clusters = (
    #         set(cluster_ids) - outlier_clusters if outlier else set(cluster_ids)
    #     )

    #     new_cluster_id_map = {
    #         old_id: new_id
    #         for new_id, old_id in enumerate(sorted(non_outlier_clusters), start=1)
    #     }

    #     for node, cluster_id in cluster.items():
    #         if outlier and cluster_id in outlier_clusters:
    #             mapped_cluster_id = -1
    #         else:
    #             mapped_cluster_id = new_cluster_id_map.get(
    #                 cluster_id, -1
    #             )
    #         all_cluster_ids.add(mapped_cluster_id)
    #         data.append({"Node Name": node.name, "Cluster Number": mapped_cluster_id})

    # df = pd.DataFrame(data)
    # full_path = os.path.join(results_dir, filename)
    # df.to_csv(full_path, index=False)
    # print(f"Cluster data saved to {filename}")
