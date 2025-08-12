from typing import Optional, List, Tuple
import numpy as np


def define_bins(
    pc,
    num_bins: int = 3,
    *,
    k_lo: int = 1,
    k_hi: Optional[int] = None,
    log_base: Optional[float] = None,
) -> List[Tuple[int, int]]:
    """Return num_bins log-spaced (inclusive) ranges covering [k_lo … k_hi]."""
    if k_hi is None:
        k_hi = pc.num_terminals
    if k_hi <= k_lo:
        raise ValueError("k_hi must be > k_lo")

    if log_base is None:
        log_base = (k_hi / k_lo) ** (1 / num_bins)

    edges = np.geomspace(k_lo, k_hi, num_bins + 1)
    edges = np.round(edges).astype(int)
    edges[0], edges[-1] = k_lo, k_hi

    bin_ranges = []
    for i in range(len(edges) - 1):
        lo = int(edges[i])
        hi = int(edges[i + 1])
        if i:
            lo = bin_ranges[-1][1] + 1
        bin_ranges.append((lo, hi))

    return bin_ranges
