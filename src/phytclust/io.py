import os
import logging
from typing import Any, Optional
import pandas as pd
import json
from datetime import datetime

from .algo.dp import cluster_map

logger = logging.getLogger("Save results")
logger.setLevel(logging.INFO)


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
        logger.info("No clusters to save.")
        return

    df = pd.DataFrame(records)
    pivot = df.pivot_table(
        index="Node Name", columns="k", values="Cluster ID", aggfunc="first"
    )
    pivot.columns = [f"clusters_k{col}" for col in pivot.columns]
    pivot.reset_index(inplace=True)
    pivot.to_csv(os.path.join(results_dir, filename), index=False, sep="\t")
    logger.info(f"Wrote clusters to {filename}")

    if getattr(pc, "peaks_by_rank", None):
        with open(os.path.join(results_dir, "peaks_by_rank.txt"), "w") as fh:
            for rank, k_val in enumerate(pc.peaks_by_rank, 1):
                fh.write(f"Rank {rank}: {k_val} clusters\n")
        logger.info("Wrote peaks_by_rank.txt")


def save_full(
    pc,
    results_dir: str,
    filename: str = "phyclust_results.csv",
) -> None:
    
    if pc._last_result is None:
        logger.info("No last result to save.")
        return
    
    results = {
        "k": pc._last_result.get("k"),
        "ks": pc._last_result.get("ks"),
        "outgroup": pc._last_result.get("outgroup"),
        "newick": pc.tree.format("newick"),
        "scores": list(pc.scores),
        "clusters": list(pc.clusters),
        "date": datetime.now().isoformat(),
    }
    
    out_path = os.path.join(results_dir, filename)
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=4)
    logger.info(f"Wrote full results to {out_path}")
    


def load_full(
    pc,
    results_dir: str,
    filename: str = "phyclust_results.csv",
) -> None:
    from phytclust.algo.core import PhytClust
    
    in_path = os.path.join(results_dir, filename)
    with open(in_path, "r") as fh:
        results = json.load(fh)
    pc = PhytClust()
    pc._last_result = results
    logger.info(f"Wrote full results to {in_path}")

    return pc
