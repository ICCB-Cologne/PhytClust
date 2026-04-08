import pytest
from phytclust.algo.bins import define_bins
from phytclust.exceptions import InvalidKError


class _PC:
    def __init__(self, n):
        self.num_terminals = n


def test_define_bins_basic_properties():
    pc = _PC(100)
    bins = define_bins(pc, num_bins=3, k_lo=1, k_hi=100)

    assert len(bins) >= 1
    assert bins[0][0] == 1
    assert bins[-1][1] == 100

    for (lo1, hi1), (lo2, hi2) in zip(bins, bins[1:]):
        assert lo1 <= hi1
        assert lo2 <= hi2
        assert hi1 + 1 == lo2

    assert all(1 <= lo <= hi <= 100 for lo, hi in bins)


def test_define_bins_rejects_invalid_range():
    pc = _PC(10)
    with pytest.raises(InvalidKError, match="k_hi must be > k_lo"):
        define_bins(pc, num_bins=3, k_lo=5, k_hi=5)
