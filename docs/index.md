# PhytClust Documentation

PhytClust performs dynamic-programming clustering of rooted phylogenetic trees, with a hard monophyly guarantee for every cluster.

It is designed for workflows where you need:

- Exact clustering at a fixed k
- Data-driven selection of good k values using score peaks
- Multiple representative resolutions of the same tree
- Practical controls for outliers, minimum cluster size, zero-length edges, and polytomies

## Start Here

If this is your first time, follow pages in this order:

1. [Getting Started](getting-started.md)
2. [CLI End-to-End](tutorials/cli-end-to-end.md)
3. [Python End-to-End](tutorials/python-end-to-end.md)

## What Makes PhytClust Distinct

- **Monophyly by construction**: final clusters are clades, not arbitrary partitions.
- **Tree-aware optimization**: objective and constraints are solved directly on the phylogeny with DP.
- **Multi-mode execution**: exact-k, global peak search, and resolution bins are all first-class.
- **Constraint-aware behavior**: supports hard constraints (`min_cluster_size`) and soft behavior controls (`outlier_size_threshold`, `prefer_fewer_outliers`).

## Output Artifacts

A typical run can produce:

- `phytclust_results.tsv`: node-to-cluster assignments across one or multiple k values
- `scores.png`: score curve over k with selected peaks
- `tree_k{K}.png`: colored tree rendering for selected K values
- `peaks_by_rank.txt`: selected k values in rank order

## Documentation Scope

- Tutorials are practical and executable.
- Reference pages map directly to the current code behavior in `src/phytclust`.
- Configuration examples are aligned with the runtime and peak config dataclasses.
