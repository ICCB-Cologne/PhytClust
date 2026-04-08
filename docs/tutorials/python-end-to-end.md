# Python End-to-End Tutorial

This tutorial mirrors CLI behavior while giving you scripting control for iterative analyses.

## 1. Load Tree and Construct PhytClust

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

## 2. Fixed-k Mode

```python
result_k = pc.run(k=5)
print(result_k["mode"])      # k
print(result_k["k"])         # 5
print(len(result_k["clusters"]))
```

## 3. Global Peak Mode

```python
result_global = pc.run(top_n=3, max_k=120)
print(result_global["mode"])   # global
print(result_global["ks"])      # selected k values
print(result_global["peaks"])   # same ranking information
```

## 4. Resolution Mode

```python
result_res = pc.run(by_resolution=True, num_bins=4, max_k=120)
print(result_res["mode"])   # resolution
print(result_res["ks"])      # one representative k per bin
```

## 5. Tune Peak Selection with PeakConfig

```python
from phytclust import PeakConfig

peak_cfg = PeakConfig(
    lambda_weight=0.5,
    ranking_mode="adjusted",
    min_prominence=5.0,
    min_k=2,
    resolution_fallback_mode="max_score",
)

result_tuned = pc.run(top_n=3, max_k=120, peak_config=peak_cfg)
print(result_tuned["ks"])
```

## 6. Set Plot and Save Defaults with RuntimeConfig

```python
from phytclust import RuntimeConfig, PlotConfig, ScorePlotConfig, ClusterPlotConfig

runtime_cfg = RuntimeConfig(
    plot=PlotConfig(
        scores=ScorePlotConfig(log_scale_y=False, fig_width=22, fig_height=11),
        cluster=ClusterPlotConfig(height_scale=0.25, marker_size=60),
    )
)

pc_cfg = PhytClust(tree, runtime_config=runtime_cfg)
res = pc_cfg.run(top_n=3, max_k=100)
pc_cfg.plot(results_dir="results/python", save=True)
pc_cfg.save(results_dir="results/python", filename="phytclust_results.tsv")
```

## 7. Side-by-Side k Comparison

```python
from phytclust.viz.cluster import plot_multiple_k

plot_multiple_k(
    pc_cfg,
    k_values=[3, 5, 8],
    results_dir="results/python_multi",
    save=True,
)
```

## 8. Integration Pattern for Pipelines

A robust batch pattern:

1. Build a `PeakConfig` template
2. Iterate trees
3. Run `pc.run(top_n=...)`
4. Persist `result["ks"]`, per-k cluster maps, and score vectors
5. Re-run only selected trees with stricter constraints for sensitivity checks
