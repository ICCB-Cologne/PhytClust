import os
from typing import Any, Optional
import pandas as pd

from ..algo.dp import cluster_map


def save_clusters(
    pc,
    results_dir: str,
    top_n: int = 1,
    filename: str = "phyclust_results.csv",
    outlier: bool = True,
    n: Optional[int] = None,
    output_all: bool = False,
) -> None:
    """
    Write out cluster assignments:
    - Save score plot if available.
    - Save specific k(s) depending on flags.
    """
    os.makedirs(results_dir, exist_ok=True)

    if getattr(pc, "plot_of_scores", None) is not None:
        pc.plot_of_scores.savefig(os.path.join(results_dir, "scores.png"))

    if pc.k is not None:
        ks = [pc.k]
    elif output_all:
        ks = list(range(1, pc.max_k + 1))
    elif n is not None:
        ks = [n]
    else:
        ks = (pc.peaks_by_rank or [])[:top_n]

    records: list[dict[str, Any]] = []

    for k_val in ks:
        cmap = cluster_map(pc, k_val)
        if cmap is None:
            continue
        counts: dict[int, int] = {}
        for cid in cmap.values():
            counts[cid] = counts.get(cid, 0) + 1

        for node, cid in cmap.items():
            mark = -1 if (outlier and counts[cid] == 1) else cid
            records.append({"Node Name": node.name, "k": k_val, "Cluster ID": mark})

    if not records:
        print("No clusters to save.")
        return

    df = pd.DataFrame(records)
    pivot = df.pivot_table(
        index="Node Name", columns="k", values="Cluster ID", aggfunc="first"
    )
    pivot.columns = [f"clusters_k{col}" for col in pivot.columns]
    pivot.reset_index(inplace=True)
    pivot.to_csv(os.path.join(results_dir, filename), index=False, sep="\t")
    print(f"Wrote clusters to {filename}")

    if getattr(pc, "peaks_by_rank", None):
        with open(os.path.join(results_dir, "peaks_by_rank.txt"), "w") as fh:
            for rank, k_val in enumerate(pc.peaks_by_rank, 1):
                fh.write(f"Rank {rank}: {k_val} clusters\n")
        print("Wrote peaks_by_rank.txt")
