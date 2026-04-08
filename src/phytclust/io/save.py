import os
import logging
from collections import Counter
from typing import Any, Optional

import pandas as pd

from ..algo.dp import cluster_map

logger = logging.getLogger(__name__)


def save_clusters(
    pc,
    results_dir: str,
    top_n: int = 1,
    filename: str = "phytclust_results.tsv",
    outlier: bool = True,
    n: Optional[int] = None,
    output_all: bool = False,
) -> Optional[str]:
    """Write cluster assignments to a TSV file.

    Returns the path to the written file, or None if nothing was saved.
    """
    os.makedirs(results_dir, exist_ok=True)

    if pc.k is not None:
        ks = [pc.k]
    elif output_all:
        ks = list(range(1, pc.max_k + 1))
    elif n is not None:
        ks = [n]
    else:
        ks = (pc.peaks_by_rank or [])[:top_n]

    outlier_thresh = pc.outlier.size_threshold
    records: list[dict[str, Any]] = []

    for k_val in ks:
        cmap = cluster_map(pc, k_val)
        if cmap is None:
            continue
        counts = Counter(cmap.values())

        for node, cid in cmap.items():
            if outlier_thresh is not None:
                mark = -1 if counts[cid] < outlier_thresh else cid
            elif outlier:
                mark = -1 if counts[cid] == 1 else cid
            else:
                mark = cid
            records.append({
                "Node Name": node.name,
                "k": k_val,
                "Cluster ID": mark,
            })

    if not records:
        logger.warning("No clusters to save.")
        return None

    df = pd.DataFrame(records)
    pivot = df.pivot_table(
        index="Node Name", columns="k",
        values="Cluster ID", aggfunc="first",
    )
    pivot.columns = [f"clusters_k{col}" for col in pivot.columns]
    pivot.reset_index(inplace=True)

    out_path = os.path.join(results_dir, filename)
    pivot.to_csv(out_path, index=False, sep="\t")
    logger.info("Wrote clusters to %s", out_path)

    if getattr(pc, "peaks_by_rank", None):
        peaks_path = os.path.join(results_dir, "peaks_by_rank.txt")
        with open(peaks_path, "w") as fh:
            for rank, k_val in enumerate(pc.peaks_by_rank, 1):
                fh.write(f"Rank {rank}: {k_val} clusters\n")
        logger.info("Wrote peaks_by_rank.txt")

    return out_path
