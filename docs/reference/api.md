# Python API Reference

This page documents the public Python interface. For a guided walkthrough, see the [Python Tutorial](../tutorials/python-end-to-end.md).

---

## PhytClust

The main class. Import it from the package root:

```python
from phytclust import PhytClust
```

### Constructor

```python
pc = PhytClust(
    tree,                          # Bio.Phylo tree, Newick string, or file path
    min_cluster_size=None,         # hard lower bound on cluster size
    outlier_size_threshold=None,   # clusters below this → outlier (-1)
    prefer_fewer_outliers=False,   # prefer solutions with fewer outlier groups
    optimize_polytomies=True,      # use native polytomy DP (vs. legacy dummy nodes)
    polytomy_mode="hard",          # "hard" or "soft"
    soft_polytomy_max_degree=18,   # safety limit for soft mode
    no_split_zero_length=False,    # block splitting at zero-length edges
    runtime_config=None,           # RuntimeConfig for plot/save defaults
    peak_config=None,              # default PeakConfig (can be overridden per run)
)
```

### `run(**kwargs)` → dict

The unified entry point. What it computes depends on which arguments you pass:

| Call pattern | Mode | What you get |
|-------------|------|-------------|
| `pc.run(k=5)` | exact-*k* | One cluster map for exactly *k* groups |
| `pc.run(top_n=3, max_k=120)` | global | Top *n* peaks from the score curve |
| `pc.run(by_resolution=True, num_bins=4)` | resolution | One peak per log bin |

You can pass `peak_config=PeakConfig(...)` to override peak detection settings for a single run.

**Return value** — a dictionary with these keys:

| Key | Type | Description |
|-----|------|-------------|
| `mode` | str | `"k"`, `"global"`, or `"resolution"` |
| `k_values` | list[int] | All selected *k* values |
| `selected_k` | int or None | The top-ranked *k* (first element of `k_values`), or `None` if no peaks found |
| `clusters` | dict or list[dict] | Leaf → cluster ID mapping. A single dict for exact-*k*, a list for multi-*k* modes |
| `scores` | numpy array or None | Score vector over all *k* values (None in exact-*k* mode) |

!!! note "Backward-compatible aliases"
    `result["ks"]` and `result["peaks"]` are aliases of `result["k_values"]`. `result["k"]` is present in exact-*k* mode. Prefer `k_values` and `selected_k` in new code.

### `get_clusters(k)` → dict

Compute (or return cached) cluster assignments for a specific *k*. Returns a dict mapping leaf names to cluster IDs.

```python
clusters = pc.get_clusters(k=7)
# {"leaf_A": 0, "leaf_B": 0, "leaf_C": 1, ...}
```

### `plot(**kwargs)`

Generate score curve and/or tree plots. Key arguments:

```python
pc.plot(
    results_dir="results/",  # where to save PNGs
    save=True,               # write to disk (vs. interactive display)
    dpi=150,                 # PNG resolution
)
```

Uses the `RuntimeConfig` set on the object for styling.

### `save(**kwargs)`

Write cluster assignments and metadata to disk:

```python
pc.save(
    results_dir="results/",
    filename="phytclust_results.tsv",
)
```

---

## Configuration classes

All importable from the package root:

```python
from phytclust import (
    PeakConfig,
    RuntimeConfig,
    PlotConfig,
    ScorePlotConfig,
    ClusterPlotConfig,
)
```

See the [Configuration Reference](configuration.md) for all fields and their defaults.

---

## Visualization helpers

### `plot_multiple_k`

Compare several *k* values side by side:

```python
from phytclust.viz.cluster import plot_multiple_k

plot_multiple_k(
    pc,                           # PhytClust object (must have run at least once)
    k_values=[3, 5, 8],          # which k values to plot
    results_dir="results/multi",  # output directory
    save=True,
)
```

---

## Exceptions

PhytClust defines a hierarchy of exceptions for clear error handling:

```python
from phytclust.exceptions import (
    PhytClustError,          # base — catch this for "any PhytClust problem"
    ValidationError,         # bad input (tree format, missing leaves, etc.)
    ConfigurationError,      # invalid config values
    InvalidKError,           # k is out of range or infeasible
    ComputationError,        # DP or scoring failure
    InvalidClusteringError,  # result violates constraints
    MissingDPTableError,     # get_clusters() called before run()
    DataError,               # general data issues
    InvalidTreeError,        # tree can't be parsed or isn't rooted
)
```

In practice, you'll mostly encounter `InvalidKError` (e.g. asking for *k*=100 on a 50-leaf tree) and `InvalidTreeError` (e.g. passing an unrooted tree without `--root-taxon`).

---

## Utilities

These are available for advanced use but aren't part of the main workflow.

### Tree metrics

```python
from phytclust.metrics.indices import (
    colless_index_calc,              # tree balance (Colless index)
    normalized_colless,              # normalized to [0, 1]
    calculate_internal_terminal_ratio,  # branch-length stemmy/tippy ratio
    calculate_variance_branch_length,   # variance of branch lengths
)
```

### Phylogenetic diversity selection

```python
from phytclust.selection.representatives import (
    maximize_pd,           # greedy selection maximizing phylogenetic diversity
    rank_terminal_nodes,   # rank leaves by distance metrics
)
```
