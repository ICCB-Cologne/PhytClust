# Getting Started

## Requirements

- Python 3.10+
- Newick tree input

Install from PyPI:

```bash
pip install phytclust
```

Or install from source in this repository:

```bash
git clone https://bitbucket.org/schwarzlab/phytclust.git
cd phytclust
pip install -e ".[dev]"
```

## Verify Installation

```bash
phytclust --version
phytclust --help
```

## Your First Run (Exact k)

```bash
phytclust examples/sample_tree.nwk --k 5 --save-fig --out-dir results/quickstart
```

This performs a fixed-k clustering and writes:

- `results/quickstart/phytclust_results.tsv`
- `results/quickstart/tree_k5.png`
- `results/quickstart/scores.png` (if score plotting is available for that run mode)

## Common Modes

### Global peak search

```bash
phytclust examples/sample_tree.nwk --top-n 3 --max-k 100 --save-fig --out-dir results/global
```

### Multi-resolution mode

```bash
phytclust examples/sample_tree.nwk --resolution --bins 4 --save-fig --out-dir results/resolution
```

## Rooting Guidance

PhytClust expects rooted trees. For unrooted trees:

```bash
# Root by named taxon
phytclust tree.nwk --k 5 --root-taxon "species_A"

# Midpoint root
phytclust tree.nwk --k 5 --root-taxon midpoint
```

## Next Steps

- For command-focused usage: [CLI End-to-End](tutorials/cli-end-to-end.md)
- For scripting workflows: [Python End-to-End](tutorials/python-end-to-end.md)
