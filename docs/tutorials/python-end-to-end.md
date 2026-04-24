# Python Tutorial: Scripting with the PhytClust API

This tutorial covers the same ground as the [CLI tutorial](cli-end-to-end.md), but from Python. Use this when you need to embed clustering in a notebook, loop over many trees, or do custom post-processing on the results.

---

## 1. Load a tree and create a PhytClust object

PhytClust wraps a [Biopython](https://biopython.org/) `Phylo` tree. You can load one from a Newick file, pass a Newick string directly, or hand in an already-parsed tree object:

```python
from Bio import Phylo
from phytclust import PhytClust

tree = Phylo.read("examples/sample_tree.nwk", "newick")

pc = PhytClust(
    tree,
    min_cluster_size=2,
    outlier_size_threshold=3,
    prefer_fewer_outliers=True,
    optimize_polytomies=True,
)
```

All the constraint and behavior flags from the CLI have Python equivalents here. Set them once on the object, and they apply to every `run()` call.

## 2. Fixed-*k* clustering

The simplest case — you know how many clusters you want:

```python
result = pc.run(k=5)

print(result["mode"])       # "k"
print(result["k"])          # 5
print(result["clusters"])   # dict: leaf name -> cluster id
```

The returned `clusters` is a dictionary mapping each leaf to its cluster ID. You can turn it into a DataFrame, merge it with metadata, or pass it to whatever comes next in your pipeline.

## 3. Global peak search

When you don't know the right *k*, ask PhytClust to find it. This computes scores for every *k* from 2 up to `max_k`, then picks the most prominent peaks:

```python
result = pc.run(top_n=3, max_k=120)

print(result["mode"])       # "global"
print(result["k_values"])   # e.g. [7, 23, 45] — the best k values, ranked
print(len(result["clusters"]))  # 3 — one cluster map per selected k
```

`result["scores"]` gives you the raw score vector if you want to plot it yourself or do further analysis. `result["k_values"]` is ordered by prominence — the first entry is the strongest signal.

## 4. Multi-resolution mode

Get one representative *k* per logarithmic bin for a coarse-to-fine view:

```python
result = pc.run(by_resolution=True, num_bins=4, max_k=120)

print(result["mode"])       # "resolution"
print(result["k_values"])   # one k per bin, e.g. [3, 12, 38, 95]
```

This is useful for generating a panel of tree colorings at different scales, or for systematically comparing how clusters merge as you zoom out.

## 5. Fine-tune peak selection

The default peak detection works well for most trees, but you can adjust it with a `PeakConfig`:

```python
from phytclust import PeakConfig

peak_cfg = PeakConfig(
    lambda_weight=0.5,       # balance between raw and adjusted ranking (0-1)
    ranking_mode="adjusted", # "raw" or "adjusted" (accounts for outliers)
    min_prominence=5.0,      # ignore small bumps in the score curve
    min_k=2,                 # don't consider k=1 as a peak
    resolution_fallback_mode="max_score",  # if a bin has no peak, use max score
)

result = pc.run(top_n=3, max_k=120, peak_config=peak_cfg)
print(result["k_values"])
```

See the [Configuration Reference](../reference/configuration.md) for all `PeakConfig` fields and what they control.

## 6. Customize plots and output

Visualization defaults are controlled by `RuntimeConfig`. You set it once on the `PhytClust` object:

```python
from phytclust import (
    RuntimeConfig, PlotConfig,
    ScorePlotConfig, ClusterPlotConfig,
)

runtime_cfg = RuntimeConfig(
    plot=PlotConfig(
        scores=ScorePlotConfig(
            log_scale_y=False,
            fig_width=22,
            fig_height=11,
        ),
        cluster=ClusterPlotConfig(
            height_scale=0.25,
            marker_size=60,
        ),
    )
)

pc = PhytClust(tree, runtime_config=runtime_cfg)
result = pc.run(top_n=3, max_k=100)

# save plots and assignments
pc.plot(results_dir="results/python", save=True)
pc.save(results_dir="results/python", filename="phytclust_results.tsv")
```

## 7. Compare multiple *k* side by side

If you want to visually compare a handful of *k* values — whether they came from peak detection or your own list — use `plot_multiple_k`:

```python
from phytclust.viz.cluster import plot_multiple_k

plot_multiple_k(
    pc,
    k_values=[3, 5, 8],
    results_dir="results/python_multi",
    save=True,
)
```

This generates one tree image per *k*, making it easy to see how clusters split as granularity increases.

## 8. Batch processing: looping over trees

A common pattern in pipelines — run the same analysis across many trees and collect results:

```python
from pathlib import Path
from Bio import Phylo
from phytclust import PhytClust, PeakConfig
import pandas as pd

peak_cfg = PeakConfig(lambda_weight=0.5, min_prominence=5.0)

rows = []
for tree_file in Path("data/trees/").glob("*.nwk"):
    tree = Phylo.read(str(tree_file), "newick")
    pc = PhytClust(tree, min_cluster_size=2, outlier_size_threshold=3)
    result = pc.run(top_n=1, max_k=100, peak_config=peak_cfg)

    rows.append({
        "tree": tree_file.stem,
        "best_k": result["selected_k"],
        "n_leaves": len(list(tree.get_terminals())),
    })

    # optionally save per-tree artifacts
    pc.save(results_dir=f"results/{tree_file.stem}")

summary = pd.DataFrame(rows)
print(summary)
```

The key design choice: create a new `PhytClust` object per tree, but reuse the same `PeakConfig`. This keeps the analysis consistent while letting each tree have its own DP state.

## 9. Working with results programmatically

The result dictionary is designed to be easy to work with:

```python
result = pc.run(top_n=3, max_k=100)

# iterate over solutions
for k, clusters in zip(result["k_values"], result["clusters"]):
    print(f"k={k}: {len(set(clusters.values()))} clusters")

    # convert to DataFrame
    df = pd.DataFrame([
        {"leaf": leaf, "cluster": cid}
        for leaf, cid in clusters.items()
    ])

    # filter out outliers (cluster id -1)
    real_clusters = df[df["cluster"] != -1]
    outliers = df[df["cluster"] == -1]
    print(f"  {len(real_clusters)} leaves in real clusters, {len(outliers)} outliers")
```

!!! tip "Result key aliases"
    For backwards compatibility, `result["ks"]` and `result["peaks"]` are aliases of `result["k_values"]`. Prefer `k_values` in new code.
