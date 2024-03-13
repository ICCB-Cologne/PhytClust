from scipy.signal import find_peaks
import numpy as np
import pandas as pd
import os


def save_clusters(
    scores, clusters, results_dir, top_n=1, filename="phyclust_results.csv", outliers=True
):
    if scores is None or len(scores) == 0:
        print("Please calculate scores first.")
        return

    peaks, _ = find_peaks(scores, distance=1)
    if peaks.size < 1:
        print("No peaks found.")
        return

    top_peaks = peaks[np.argsort(-scores[peaks])][:top_n]

    selected_clusters = [clusters[i] for i in top_peaks]

    data = []
    all_cluster_ids = set()
    for cluster in selected_clusters:
        cluster_ids = list(cluster.values())
        outlier_clusters = {
            cluster_id
            for cluster_id in cluster_ids
            if cluster_ids.count(cluster_id) == 1
        }
        non_outlier_clusters = set(cluster_ids) - outlier_clusters

        new_cluster_id_map = {
            old_id: new_id
            for new_id, old_id in enumerate(sorted(non_outlier_clusters))
        }

        for node, cluster_id in cluster.items():
            # Nodes that are alone in a cluster are marked as outliers
            if outliers and cluster_id in outlier_clusters:
                cluster_id = -1
            else:
                cluster_id = new_cluster_id_map[cluster_id]

            all_cluster_ids.add(cluster_id)
            data.append({"Node Name": node, "Cluster Number": cluster_id})

    df = pd.DataFrame(data)
    full_path = os.path.join(results_dir, filename)
    df.to_csv(full_path, index=False)
    print(f"Cluster data saved to {filename}")
