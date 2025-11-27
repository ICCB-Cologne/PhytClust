# bootstrap_stability.py
from collections import defaultdict
from typing import List, Dict, Any, Optional

import numpy as np

from .core import PhytClust  # your class


def _taxon_order(trees):
    """
    Get a sorted, consistent list of taxa present in *all* trees (intersection).
    You can choose union instead if you want, but intersection is safer.
    """
    sets = []
    for t in trees:
        leaves = [term.name for term in t.get_terminals()]
        sets.append(set(leaves))

    common = set.intersection(*sets)
    if not common:
        raise ValueError("No common taxa across bootstrap trees.")

    return sorted(common)


def _labels_from_cmap(cmap, taxa, missing_label: int = -1):
    """
    Convert {leaf_obj -> cluster_id} into a label vector aligned with `taxa` list.
    Assumes cmap keys are leaf objects with .name.
    """
    name_to_cluster = {leaf.name: cid for leaf, cid in cmap.items()}
    labels = np.full(len(taxa), missing_label, dtype=int)
    for i, name in enumerate(taxa):
        if name in name_to_cluster:
            labels[i] = name_to_cluster[name]
    return labels


def compute_coassoc_for_k(
    trees,
    k: int,
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[Dict[str, Any]] = None,
):
    """
    For a given k, run PhytClust on each bootstrap tree and compute:

    - labels: (B, N) int array of cluster IDs (or -1 for missing)
    - coassoc: (N, N) float array of co-association frequencies
    - taxa: list of taxon names in order
    """
    if pc_kwargs is None:
        pc_kwargs = {}

    taxa = _taxon_order(trees)
    B = len(trees)
    N = len(taxa)
    labels = np.full((B, N), -1, dtype=int)

    for b, tree in enumerate(trees):
        pc = PhytClust(
            tree=tree, outgroup=outgroup, min_cluster_size=min_cluster_size, **pc_kwargs
        )
        res = pc.run(k=k)
        cmap = res["clusters"]  # {leaf_obj -> cluster_id}
        labels[b, :] = _labels_from_cmap(cmap, taxa)

    coassoc = np.zeros((N, N), dtype=float)
    counts = np.zeros((N, N), dtype=float)

    for b in range(B):
        row = labels[b]
        for i in range(N):
            if row[i] < 0:
                continue
            for j in range(i + 1, N):
                if row[j] < 0:
                    continue
                counts[i, j] += 1
                counts[j, i] += 1
                if row[i] == row[j]:
                    coassoc[i, j] += 1
                    coassoc[j, i] += 1

    mask = counts > 0
    coassoc[mask] = coassoc[mask] / counts[mask]
    np.fill_diagonal(coassoc, 1.0)

    return taxa, labels, coassoc


def stability_for_k(coassoc: np.ndarray) -> float:
    """
    Simple example: average co-association over all pairs.
    You can design something more clever later.
    """
    N = coassoc.shape[0]
    # exclude diagonal
    triu_idx = np.triu_indices(N, k=1)
    vals = coassoc[triu_idx]
    if vals.size == 0:
        return 0.0
    return float(np.mean(vals))


def choose_k_by_stability(
    trees,
    k_values: List[int],
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[Dict[str, Any]] = None,
):
    """
    For each k in k_values, compute co-association and a stability score.
    Return k* with max stability and all metadata.

    Returns dict:
    {
      "best_k": k_star,
      "scores": {k: score_k, ...},
      "coassoc": {k: coassoc_matrix, ...},
      "taxa": [...],
    }
    """
    scores = {}
    coassoc_by_k = {}
    taxa_ref = None

    for k in k_values:
        taxa, labels, coassoc = compute_coassoc_for_k(
            trees,
            k,
            outgroup=outgroup,
            min_cluster_size=min_cluster_size,
            pc_kwargs=pc_kwargs,
        )
        if taxa_ref is None:
            taxa_ref = taxa
        else:
            if taxa != taxa_ref:
                raise ValueError(
                    "Taxon order mismatch across k; this should not happen."
                )

        score_k = stability_for_k(coassoc)
        scores[k] = score_k
        coassoc_by_k[k] = coassoc

    best_k = max(scores, key=scores.get)
    return {
        "best_k": best_k,
        "scores": scores,
        "coassoc": coassoc_by_k,
        "taxa": taxa_ref,
    }

