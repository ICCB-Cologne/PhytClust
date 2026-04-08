"""Bootstrap co-association stability analysis."""

from typing import Any, Optional

import numpy as np

from .core import PhytClust
from ..exceptions import DataError


def _taxon_order(trees):
    """Sorted list of taxa present in *all* trees (intersection)."""
    sets = [
        {term.name for term in t.get_terminals()}
        for t in trees
    ]
    common = set.intersection(*sets)
    if not common:
        raise DataError("No common taxa across bootstrap trees.")
    return sorted(common)


def _labels_from_cmap(cmap, taxa, missing_label: int = -1):
    """Convert {leaf_obj -> cluster_id} into a label vector aligned with `taxa`."""
    name_to_cluster = {leaf.name: cid for leaf, cid in cmap.items()}
    labels = np.full(len(taxa), missing_label, dtype=int)
    for i, name in enumerate(taxa):
        if name in name_to_cluster:
            labels[i] = name_to_cluster[name]
    return labels


def _coassoc_from_labels(labels: np.ndarray) -> np.ndarray:
    """Compute co-association matrix from a (B, N) label array.

    For each pair (i, j), the co-association is the fraction of bootstrap
    replicates where both taxa were assigned to the same cluster (ignoring
    replicates where either taxon has label < 0).
    """
    B, N = labels.shape
    coassoc = np.zeros((N, N), dtype=float)
    counts = np.zeros((N, N), dtype=float)

    for b in range(B):
        row = labels[b]
        valid = row >= 0
        # Only consider pairs where both are valid
        valid_pair = np.outer(valid, valid)
        counts += valid_pair
        # Same-cluster pairs
        same = np.equal.outer(row, row) & valid_pair
        coassoc += same

    mask = counts > 0
    coassoc[mask] /= counts[mask]
    np.fill_diagonal(coassoc, 1.0)
    return coassoc


def compute_coassoc_for_k(
    trees,
    k: int,
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[dict[str, Any]] = None,
):
    """For a given k, run PhytClust on each bootstrap tree and compute
    co-association.

    Returns
    -------
    taxa : list[str]
        Taxon names in order.
    labels : ndarray of shape (B, N)
        Cluster IDs per replicate (or -1 for missing).
    coassoc : ndarray of shape (N, N)
        Co-association frequencies.
    """
    if pc_kwargs is None:
        pc_kwargs = {}

    taxa = _taxon_order(trees)
    B = len(trees)
    N = len(taxa)
    labels = np.full((B, N), -1, dtype=int)

    for b, tree in enumerate(trees):
        pc = PhytClust(
            tree=tree,
            outgroup=outgroup,
            min_cluster_size=min_cluster_size,
            **pc_kwargs,
        )
        res = pc.run(k=k)
        cmap = res["clusters"][0]
        labels[b, :] = _labels_from_cmap(cmap, taxa)

    coassoc = _coassoc_from_labels(labels)
    return taxa, labels, coassoc


def stability_for_k(coassoc: np.ndarray) -> float:
    """Average co-association over all off-diagonal pairs."""
    N = coassoc.shape[0]
    triu_idx = np.triu_indices(N, k=1)
    vals = coassoc[triu_idx]
    if vals.size == 0:
        return 0.0
    return float(np.mean(vals))


def choose_k_by_stability(
    trees,
    k_values: list[int],
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[dict[str, Any]] = None,
):
    """For each k in k_values, compute co-association and stability.

    Returns
    -------
    dict with keys:
        best_k : int
        scores : dict[int, float]
        coassoc : dict[int, ndarray]
        taxa : list[str]
    """
    scores = {}
    coassoc_by_k = {}
    taxa_ref = None

    for k in k_values:
        taxa, _labels, coassoc = compute_coassoc_for_k(
            trees,
            k,
            outgroup=outgroup,
            min_cluster_size=min_cluster_size,
            pc_kwargs=pc_kwargs,
        )
        if taxa_ref is None:
            taxa_ref = taxa
        elif taxa != taxa_ref:
            raise DataError(
                "Taxon order mismatch across k; this should not happen."
            )

        scores[k] = stability_for_k(coassoc)
        coassoc_by_k[k] = coassoc

    best_k = max(scores, key=scores.get)
    return {
        "best_k": best_k,
        "scores": scores,
        "coassoc": coassoc_by_k,
        "taxa": taxa_ref,
    }
