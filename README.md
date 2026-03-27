# PhytClust <img src="src/phytclust/phytclust_logo_colour.png" width="120" align="right" alt="PhytClust logo">

Monophyletic, dynamic-programming **clustering of phylogenetic trees**.

PhytClust finds clusterings of the leaves of a phylogenetic tree such that **each cluster corresponds to a monophyletic clade**.

It supports:
- **Exact k-way clustering** (`run(k=...)`)
- **Global peak search** over k using internal validation scores (Calinski–Harabasz + elbow-style index)
- **Multi-resolution clustering**: one representative k per log-spaced resolution bin
- Polytomies, minimum cluster size constraints, support-aware branch lengths, outlier penalties, and more

---

## Requirements

- Python 3.10+
- Biopython, NumPy (installed automatically via pip)

Input trees are expected in **Newick** format. If the tree is unrooted, provide an outgroup (CLI option) or root it beforehand.

---

## Installation

### 1) Recommended: clean conda environment + install from PyPI

```bash
conda create -n phyt_env python=3.10
conda activate phyt_env
pip install phytclust
```

### 2) Install from source (development)

```bash
git clone https://bitbucket.org/schwarzlab/phytclust.git
cd phytclust
pip install -e ".[dev]"
```

---

## Command-line usage

### Exact k clusters

Compute an exact k-way clustering, plot it, and save PNG + CSV under `./results`:

```bash
phytclust tree.nwk --k 5 --save-fig --out-dir results
```

### Global clustering solution

Search for the top 3 peaks up to `k = 200`, save everything in `./out`:

```bash
phytclust tree.nwk --top-n 3 \
  --max-k 200 \
  --save-fig \
  --out-dir out
```

### Multi-resolution clustering

Pick one representative peak per 4 log-spaced bins of k, don’t show plots interactively,
save all k-specific CSVs and plots:

```bash
phytclust tree.nwk --bins 4 \
  --no-plot \
  --save-all-k \
  --save-fig \
  --out-dir out
```

Run `phytclust --help` for the complete CLI reference.

---

## Outlier Penalties and Cluster Size Constraints

PhytClust supports three key parameters to control clustering behavior:

### 1. Minimum Cluster Size (Hard Constraint)

The `--min-cluster-size` option enforces that no final cluster can have fewer leaves than the specified threshold. This is a **hard constraint** — infeasible clusterings are blocked entirely.

```bash
phytclust tree.nwk -k 5 --min-cluster-size 3
```

The DP algorithm prevents any cluster from having fewer than 3 leaves. If it's impossible to partition the tree into 5 clusters with all clusters ≥ 3 leaves, the command will fail with an appropriate error.

### 2. Outlier Penalty (Soft Constraint)

The `--outlier-penalty` and `--outlier-size-threshold` options add a **soft penalty** to the DP objective function when isolating small clusters. This steers the clustering toward grouping small clades with larger ones.

```bash
phytclust tree.nwk --top-n 3 \
  --outlier-penalty 10 \
  --outlier-size-threshold 2
```

Any cluster with fewer than 2 leaves incurs a penalty of 10 in the DP cost. This makes the DP prefer to merge small singletons into larger clusters, but doesn't forbid it entirely. The `--outlier-size-threshold` output flag marks clusters below the threshold with `-1` in the output CSV for easy post-processing.

### 3. Polytomy Optimization (Default On)

The `--no-optimize-polytomies` option controls how polytomous (non-binary) nodes are handled:

- **Default (`optimize_polytomies=True`)**: Polytomous nodes are kept and handled directly in the DP via sequential convolution. This automatically finds the globally optimal split across all possible binary orderings of the polytomy children — equivalent to trying all `(2k−3)!!` resolutions but in polynomial time.

- **With `--no-optimize-polytomies`**: Polytomies are first binarized using 0-length dummy nodes (left-to-right greedy pairing). The DP then treats dummy nodes as unsplittable units. This mode preserves the behavior from earlier versions.

```bash
# Optimal polytomy handling (default)
phytclust tree_with_polytomies.nwk -k 5

# Greedy dummy-node polytomy handling (legacy)
phytclust tree_with_polytomies.nwk -k 5 --no-optimize-polytomies
```

For most use cases, the default convolution-based approach is preferred because it considers all biologically plausible splits within a polytomy, not just the arbitrary grouping imposed by left-to-right pairing.

### Combined Example

```bash
phytclust tree.nwk --top-n 3 \
  --max-k 100 \
  --min-cluster-size 2 \
  --outlier-penalty 5 \
  --outlier-size-threshold 3 \
  --save-fig --out-dir results
```

This finds the top 3 global peaks up to k=100, enforces all clusters have ≥2 leaves, and penalizes clusters with <3 leaves.

---

## Rooting Unrooted Trees

PhytClust requires a rooted tree for dynamic programming. If your tree is unrooted, use the `--root-taxon` option:

### Root at a Specific Taxon

Specify any leaf name as the root. The tree will be rerooted such that this taxon is on the outgroup edge:

```bash
phytclust tree.nwk -k 5 --root-taxon "species_A"
```

### Midpoint Rooting

For unbiased rooting, use `"midpoint"`:

```bash
phytclust tree.nwk -k 5 --root-taxon "midpoint"
```

This roots the tree at the edge equidistant from the two most-distant leaves, a common phylogenetic practice for rooting in the absence of an outgroup.

### Combining with Outgroup

You can specify both a root and an outgroup:

```bash
phytclust tree.nwk -k 5 --root-taxon "species_A" --outgroup "outgroup_species"
```

---

## Handling Zero-Length Clusters

By default, PhytClust allows splitting at any edge. If you want to prevent the DP from splitting nodes where the resulting child clusters would have **only zero-length internal edges** (no evolutionary signal), use `--no-split-zero-length`:

```bash
phytclust tree.nwk -k 5 --no-split-zero-length
```

This treats zero-length clusters like outliers — they can exist as clusters, but splitting them is blocked because there's no phylogenetic signal to justify separation.

This is useful when:
- Your tree contains many zero-length edges (identical sequences, collapsed nodes, redundant copies)
- You want to avoid arbitrary partitioning of phylogenetically indistinguishable taxa
- You need clusters with meaningful evolutionary distance/support

When enabled, any node where both child edges are zero-length cannot be split — the children must stay together.

---

## Visualization: Comparing Multiple Cluster Solutions

### Using the Python API

Compare multiple cluster solutions (different k values) side-by-side:

```python
from phytclust import PhytClust
from phytclust.viz.cluster_plot import plot_multiple_k

pc = PhytClust(tree="tree.nwk")
pc.run(top_n=5)  # Compute top 5 peaks

# Plot the top 3 peak solutions
plot_multiple_k(pc, top_n=3, save=True, results_dir="./results")

# Or specify exact k values
plot_multiple_k(pc, k_values=[3, 5, 8], save=True, results_dir="./results")
```

This generates individual tree plots for each k, saved as `tree_k{k}.png`, allowing you to visually compare how cluster assignments change across different solutions.

### Using the Command-Line

Save cluster plots for the top 3 peaks:

```bash
phytclust tree.nwk --top-n 3 --save-fig --out-dir results
```

This will save `tree_k{k}.png` for each selected k value, enabling side-by-side comparison.

---

## Development

Run tests:

```bash
pytest -q
```

---


## Please cite

> K. Ganesan, E. Billard, T.L. Kaufmann, C. B Strange, M. C. Cwikla, A. Altenhoff, C. Dessimoz, R.F. Schwarz, *PhytClust* (2025), repository: https://bitbucket.org/schwarzlab/phytclust/
