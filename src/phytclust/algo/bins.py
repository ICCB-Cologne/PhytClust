from typing import Optional, List, Tuple
import numpy as np


def define_bins(
    pc,
    num_bins: int = 3,
    *,
    k_lo: int = 1,
    k_hi: Optional[int] = None,
) -> List[Tuple[int, int]]:
    """Return num_bins log-spaced (inclusive) ranges covering [k_lo … k_hi]."""
    if k_hi is None:
        k_hi = pc.num_terminals
    if k_hi <= k_lo:
        raise ValueError("k_hi must be > k_lo")

    raw = np.geomspace(k_lo, k_hi, num_bins + 1)
    edges = np.unique(np.round(raw).astype(int))
    edges[0], edges[-1] = k_lo, k_hi
    edges = np.maximum.accumulate(edges)

    bins = []
    for i in range(len(edges) - 1):
        lo = edges[i] if i == 0 else bins[-1][1] + 1
        hi = edges[i + 1]
        if lo > hi:
            continue
        bins.append((lo, hi))
    return bins
