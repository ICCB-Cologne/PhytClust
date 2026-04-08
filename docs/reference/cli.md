# CLI Flags Reference

This page reflects the current argparse interface in `src/phytclust/cli/_cli.py`.

## Core Mode Selection

- `-k, --k`: exact k-way clustering (must be >=2)
- `--top-n`: number of global peaks (default `1`)
- `--resolution`: enable one-peak-per-log-bin mode
- `--bins`: number of bins for resolution mode (default `3`)

Rules:

- If `--k` is given, it overrides global/resolution peak modes.
- `--top-n` is ignored for exact-k mode.

## Input and Output

- Positional `tree`: Newick file path or `-` for stdin
- `-o, --out-dir`: output directory (default `results`)
- `--save-fig`: write score and tree PNG artifacts
- `--save-tree`: explicitly write tree PNGs
- `--save-all-k`: write assignments for every k from 1..max_k
- `--tsv-name`: output TSV filename (default `phytclust_results.tsv`)
- `--no-tsv`: skip writing TSV
- `--dpi`: PNG resolution (default `150`)

## Tree Handling

- `--outgroup`: taxon excluded from all clusters
- `--root-taxon`: root on taxon name, or `midpoint`

## Constraints and Behavior

- `--min-cluster-size`: hard lower bound on cluster size
- `--outlier-size-threshold`: clusters below threshold considered outliers in output handling
- `--prefer-fewer-outliers`: prioritize solutions with fewer outlier clusters
- `--no-optimize-polytomies`: disable optimal polytomy handling and use legacy dummy-node behavior
- `--polytomy-mode {hard,soft}`: choose strict hard multifurcation handling or subset-merge soft handling
- `--soft-polytomy-max-degree`: complexity guardrail for soft mode (soft mode is exponential in node degree)
- `--no-split-zero-length`: block splitting clades with all-zero edge structure

## Peak Selection / Scoring Controls

- `--max-k`: explicit maximum k for peak search
- `--max-k-limit`: if `--max-k` omitted, use `ceil(max_k_limit * leaves)`
- `--lambda-weight`: peak prominence blending parameter passed into `PeakConfig`

## Runtime and Config Files

- `--config`: optional YAML/JSON config file
- `--plot/--no-plot`: show or suppress interactive plots

Configuration mapping supports:

- `peak`: maps into `PeakConfig`
- `algorithm.polytomy_mode`: maps into `PhytClust(polytomy_mode=...)` unless CLI flag overrides
- `algorithm.soft_polytomy_max_degree`: maps into `PhytClust(soft_polytomy_max_degree=...)` unless CLI flag overrides
- `runtime`: nested runtime dataclass overrides
- Backward-compatible blocks: `plot`, `scores_plot`, `save`

## Logging and Diagnostics

- `-v, --verbose`: increase log verbosity (`-vv` for debug)
- `-q, --quiet`: reduce log verbosity (`-qq` for critical)
- `--no-color`: disable rich color logs
- `--log-format`: logging formatter string
- `--time`: report phase timings
- `--progress`: show spinner UI when supported
- `--version`: print package version

## Example Commands

```bash
# Exact k
phytclust tree.nwk --k 5 --save-fig --out-dir out

# Global top-3 peaks
phytclust tree.nwk --top-n 3 --max-k 120 --save-fig --out-dir out

# Resolution mode
phytclust tree.nwk --resolution --bins 4 --max-k 120 --save-fig --out-dir out
```
