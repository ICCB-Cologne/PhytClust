# CLI End-to-End Tutorial

This tutorial is a practical CLI workflow from first run to constrained peak search.

## 1. Run a Baseline Exact-k Clustering

```bash
phytclust examples/sample_tree.nwk --k 5 --save-fig --out-dir results/cli_baseline
```

What this gives you:

- Single fixed-k solution (`k=5`)
- Cluster assignments in TSV
- Colored tree visualization

Use this mode when `k` is known from prior biology or downstream constraints.

## 2. Let PhytClust Propose Good k Values

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --max-k 120 \
  --save-fig \
  --out-dir results/cli_global
```

What to inspect:

- `peaks_by_rank.txt`: selected k values in rank order
- `scores.png`: whether top peaks are well-separated
- `phytclust_results.tsv`: cluster assignments for selected peaks

## 3. Multi-resolution Summary (One k per Bin)

```bash
phytclust examples/sample_tree.nwk \
  --resolution \
  --bins 4 \
  --save-fig \
  --out-dir results/cli_resolution
```

Use this when you need a compact panel of representative granularities rather than many neighboring k values.

## 4. Enforce Hard and Soft Constraints

### Hard minimum cluster size

```bash
phytclust examples/sample_tree.nwk --top-n 3 --min-cluster-size 3 --out-dir results/cli_min_size
```

If a k is infeasible under this constraint, it is not returned.

### Soft outlier handling

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --outlier-size-threshold 3 \
  --prefer-fewer-outliers \
  --out-dir results/cli_outlier
```

Interpretation:

- Clusters smaller than threshold are considered outliers in output handling
- With `--prefer-fewer-outliers`, DP prefers fewer outlier clusters before raw-cost tie-breaking

## 5. Handle Polytomies and Zero-Length Edges Explicitly

### Polytomies

```bash
# default: optimized handling
phytclust tree_with_polytomies.nwk --k 5

# legacy behavior
phytclust tree_with_polytomies.nwk --k 5 --no-optimize-polytomies
```

### Zero-length split guard

```bash
phytclust tree.nwk --k 5 --no-split-zero-length
```

Useful for trees with many collapsed or near-identical leaves.

## 6. Reproducible Run via YAML Config

Create config file:

```yaml
# examples/phytclust.config.yaml
peak:
  lambda_weight: 0.6
  ranking_mode: adjusted
  min_prominence: 5.0
  resolution_fallback_mode: max_score

runtime:
  plot:
    scores:
      log_scale_y: false
      fig_width: 20
    cluster:
      height_scale: 0.25
      marker_size: 55
  save:
    csv_name: phytclust_results.tsv
```

Run with config:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --max-k 120 \
  --config examples/phytclust.config.yaml \
  --save-fig \
  --out-dir results/cli_config
```

## 7. Artifact Quality Checklist

Before downstream analysis, verify:

- Peak separation is visible in `scores.png`
- Returned k values are biologically plausible and not all near `max-k`
- Outlier-marked clusters (`-1`) match expected tiny clades
- Tree plots are legible at your reporting scale (adjust `height_scale` and `marker_size` if needed)
