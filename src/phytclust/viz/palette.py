"""PhytClust brand palette — single source of truth for all colors.

The 8 base colors are derived from the PhytClust logo. When more colors
are needed, lighter/more transparent variants are generated automatically.
"""

from __future__ import annotations

import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap

# 8 logo-derived base colors (same as GUI BASE_COLORS)
BASE_HEX: list[str] = [
    "#b84b4b",  # red
    "#849060",  # olive
    "#3d7c74",  # teal
    "#6e3f8a",  # purple
    "#ceb94b",  # gold
    "#3f648a",  # steel blue
    "#3f408a",  # indigo
    "#da63aa",  # pink
]

# Colorblind-safe accent palette (for score-plot peaks / bin labels)
ACCENT_HEX: list[str] = [
    "#0072B2",
    "#009E73",
    "#D55E00",
    "#56B4E9",
    "#E69F00",
    "#CC79A7",
    "#F0E442",
    "#000000",
]

BASE_RGBA: list[tuple[float, ...]] = [
    mcolors.to_rgba(h) for h in BASE_HEX
]


def expand_palette(
    n: int,
    base: list[str] | None = None,
) -> list[tuple[float, ...]]:
    """Return *n* RGBA colors by cycling base colors with decreasing alpha.

    For n <= len(base), returns the first n base colors.
    For larger n, lightens and reduces alpha in successive rounds,
    matching the GUI's ``generateClusterColors`` behavior.
    """
    base = base or BASE_HEX
    base_rgba = [mcolors.to_rgba(h) for h in base]
    if n <= len(base_rgba):
        return base_rgba[:n]

    colors: list[tuple[float, ...]] = []
    repeats = int(np.ceil(n / len(base_rgba)))
    for r in range(repeats):
        factor = (r - 1) * 0.18 if r > 0 else 0.0
        alpha = max(1.0 - r * 0.15, 0.3)
        for rgba in base_rgba:
            adjusted = _adjust_lightness(rgba, factor)
            colors.append((*adjusted[:3], alpha))
    return colors[:n]


def get_cmap(n: int = 20) -> ListedColormap:
    """Return a PhytClust ``ListedColormap`` with *n* colors."""
    return ListedColormap(expand_palette(n), name="phytclust")


def _adjust_lightness(
    rgba: tuple[float, ...], factor: float
) -> tuple[float, ...]:
    """Shift lightness of an RGBA color by *factor* (-1..1)."""
    r, g, b = rgba[:3]
    if factor > 0:
        r = r + (1 - r) * factor
        g = g + (1 - g) * factor
        b = b + (1 - b) * factor
    elif factor < 0:
        r = r * (1 + factor)
        g = g * (1 + factor)
        b = b * (1 + factor)
    return (
        max(0.0, min(1.0, r)),
        max(0.0, min(1.0, g)),
        max(0.0, min(1.0, b)),
        rgba[3],
    )
