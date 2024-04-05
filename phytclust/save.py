from scipy.signal import find_peaks
import numpy as np
import pandas as pd
import os


def save_clusters(
    scores, clusters, results_dir, top_n=1, filename="phyclust_results.csv", outlier=True, selected_peaks=None, n=None
):
    if scores is None or len(scores) == 0:
        print("Please calculate scores first.")
        return

    if selected_peaks is not None:
        top_peaks = selected_peaks
    else:
        peaks, _ = find_peaks(scores, distance=1)
        if peaks.size < 1:
            print("No peaks found.")
            return
        top_peaks = peaks[np.argsort(-scores[peaks])][:top_n]

    if n is not None:
        selected_clusters = {n: clusters[n]} if n in clusters else {}
    else:
        selected_clusters = [clusters[i] for i in top_peaks]

    data = []
    all_cluster_ids = set()
    for cluster in selected_clusters:
        cluster_ids = list(cluster.values())

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

        new_cluster_id_map = {
            old_id: new_id
            for new_id, old_id in enumerate(sorted(non_outlier_clusters), start=1)
        }

        for node, cluster_id in cluster.items():
            if outlier and cluster_id in outlier_clusters:
                mapped_cluster_id = -1
            else:
                mapped_cluster_id = new_cluster_id_map.get(
                    cluster_id, -1
                )
            all_cluster_ids.add(mapped_cluster_id)
            data.append({"Node Name": node.name, "Cluster Number": mapped_cluster_id})

    df = pd.DataFrame(data)
    full_path = os.path.join(results_dir, filename)
    df.to_csv(full_path, index=False)
    print(f"Cluster data saved to {filename}")
