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

## Development

Run tests:

```bash
pytest -q
```

---


## Please cite

> K. Ganesan, E. Billard, T.L. Kaufmann, C. B Strange, M. C. Cwikla, A. Altenhoff, C. Dessimoz, R.F. Schwarz, *PhytClust* (2025), repository: https://bitbucket.org/schwarzlab/phytclust/
