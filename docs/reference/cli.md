# CLI Reference

Complete flag reference for the `phytclust` command. For a guided walkthrough, see the [CLI Tutorial](../tutorials/cli-end-to-end.md).

---

## Positional argument

| Argument | Description |
|----------|-------------|
| `tree` | Path to a Newick file, or `-` to read from stdin |

## Clustering mode

These flags control *what* PhytClust computes. If `--k` is given, it takes priority over peak-based modes.

| Flag | Default | Description |
|------|---------|-------------|
| `-k, --k` | — | Cluster into exactly this many groups (must be >= 2) |
| `--top-n` | `1` | Return this many peaks from the global score curve |
| `--resolution` | off | Enable multi-resolution mode: one peak per log bin |
| `--bins` | `3` | Number of bins for resolution mode |

## Input / output

| Flag | Default | Description |
|------|---------|-------------|
| `-o, --out-dir` | `results/` | Where to write all output files |
| `--save-fig` | off | Write score and tree PNG files |
| `--save-tree` | off | Write tree PNGs only (no score plot) |
| `--save-all-k` | off | Write cluster assignments for every *k* from 1 to max_k |
| `--tsv-name` | `phytclust_results.tsv` | Filename for the cluster assignment table |
| `--no-tsv` | off | Skip writing the TSV entirely |
| `--dpi` | `150` | Resolution for PNG output |

## Tree handling

| Flag | Default | Description |
|------|---------|-------------|
| `--outgroup` | — | Taxon name to exclude from clustering (pruned before DP) |
| `--root-taxon` | — | Re-root the tree on this taxon, or use `midpoint` for midpoint rooting |

## Constraints

| Flag | Default | Description |
|------|---------|-------------|
| `--min-cluster-size` | — | Hard lower bound: no cluster may have fewer leaves than this |
| `--outlier-size-threshold` | — | Clusters below this size are marked as outliers (`-1`) in output |
| `--prefer-fewer-outliers` | off | When ranking solutions, prefer those with fewer outlier clusters |

## Polytomies and zero-length edges

| Flag | Default | Description |
|------|---------|-------------|
| `--no-optimize-polytomies` | off | Fall back to legacy dummy-node insertion instead of native polytomy DP |
| `--polytomy-mode` | `hard` | `hard`: each child goes entirely into one cluster. `soft`: allows partial merges across children |
| `--soft-polytomy-max-degree` | `18` | Safety limit for soft mode — nodes above this degree use hard mode instead |
| `--no-split-zero-length` | off | Prevent the DP from splitting at edges with zero branch length |

## Peak selection

These control how peaks are detected in the score curve. For detailed explanations, see [Configuration > PeakConfig](configuration.md#peakconfig-fields).

| Flag | Default | Description |
|------|---------|-------------|
| `--max-k` | — | Explicit upper bound on *k* to search |
| `--max-k-limit` | `0.9` | If `--max-k` is omitted, search up to `ceil(max_k_limit * n_leaves)` |
| `--lambda-weight` | `0.5` | Blending between raw and adjusted peak ranking (0 = raw only, 1 = adjusted only) |

## Configuration files

| Flag | Default | Description |
|------|---------|-------------|
| `--config` | — | Path to a YAML or JSON config file |
| `--plot / --no-plot` | `--plot` | Show or suppress interactive matplotlib windows |

Config files use a nested structure — see the [Configuration Reference](configuration.md) for the full schema. CLI flags override config file values when both are present.

## Logging and diagnostics

| Flag | Default | Description |
|------|---------|-------------|
| `-v, --verbose` | — | Increase log verbosity (use `-vv` for DEBUG level) |
| `-q, --quiet` | — | Decrease log verbosity (use `-qq` for CRITICAL only) |
| `--no-color` | off | Disable colored terminal output |
| `--log-format` | — | Custom Python logging format string |
| `--time` | off | Print timing for each phase (DP, scoring, plotting, etc.) |
| `--progress` | off | Show a spinner during long computations |
| `--version` | — | Print version and exit |

---

## Examples

```bash
# exact k, with figures
phytclust tree.nwk --k 5 --save-fig --out-dir out

# global top-3 peaks, scanning up to k=120
phytclust tree.nwk --top-n 3 --max-k 120 --save-fig --out-dir out

# resolution mode with 4 bins
phytclust tree.nwk --resolution --bins 4 --max-k 120 --save-fig --out-dir out

# constrained run with config file
phytclust tree.nwk --top-n 3 --min-cluster-size 3 --config settings.yaml --out-dir out

# pipe from stdin
cat tree.nwk | phytclust - --k 5

# verbose timing
phytclust tree.nwk --top-n 1 --time -vv
```
