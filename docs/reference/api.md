# Python API Overview

## Main Class

```python
from phytclust import PhytClust
```

Constructor accepts:

- Tree input (`Bio.Phylo` tree, Newick string, or path-like)
- Constraint and behavior flags (`min_cluster_size`, `optimize_polytomies`, `polytomy_mode`, `soft_polytomy_max_degree`, `no_split_zero_length`)
- Optional `runtime_config` and default `peak_config`

## Unified Run Interface

```python
result = pc.run(...)
```

### Exact-k mode

```python
result = pc.run(k=5)
```

Returned keys:

- `mode`: `"k"`
- `k`: selected integer k
- `clusters`: dict mapping clade to cluster id
- `scores`: `None`
- `peaks`: `[k]`

### Global peak mode

```python
result = pc.run(top_n=3, max_k=120)
```

Returned keys:

- `mode`: `"global"`
- `ks`: selected k values by rank
- `clusters`: list of cluster maps (one per selected k)
- `scores`: score vector (numpy array copy or `None`)
- `peaks`: same values as `ks`

### Resolution mode

```python
result = pc.run(by_resolution=True, num_bins=4, max_k=120)
```

Returned keys:

- `mode`: `"resolution"`
- `ks`: one representative k per bin
- `clusters`: list of cluster maps
- `scores`: score vector copy
- `peaks`: same values as `ks`

## Utility Methods

### `get_clusters(k)`

Returns exact cluster mapping for k and caches it.

### `plot(...)`

Dispatches to visualization backend and can save plots to `results_dir`.

### `save(...)`

Writes TSV assignments and optional metadata files.

## Typical Script Skeleton

```python
from Bio import Phylo
from phytclust import PhytClust, PeakConfig

tree = Phylo.read("tree.nwk", "newick")
pc = PhytClust(tree, min_cluster_size=2, outlier_size_threshold=3)

peak_cfg = PeakConfig(lambda_weight=0.5, min_prominence=5.0)
result = pc.run(top_n=3, max_k=120, peak_config=peak_cfg)

pc.plot(results_dir="results", save=True)
pc.save(results_dir="results", filename="phytclust_results.tsv")
```

## Exceptions

Package exports domain exceptions from `phytclust.exceptions`, including:

- `PhytClustError`
- `ValidationError`
- `ConfigurationError`
- `ComputationError`
- `DataError`
- `InvalidKError`
- `InvalidTreeError`
- `MissingDPTableError`
- `InvalidClusteringError`
