# PhytClust <img src="https://bitbucket.org/schwarzlab/phytclust/raw/HEAD/src/phytclust/phytclust_logo_colour.png" width="120">

Monophyletic, dynamic-programming **clustering of phylogenetic trees**.

PhytClust finds clusterings of the leaves of a rooted tree such that **every cluster is a monophyletic clade**. It supports:

- **Exact k-way clustering** (`run(k=...)`)
- **Global peak search** in k using a Calinski–Harabasz + Elbow 1 Index for validating cluster quality
- **Multi-resolution clustering**: one representative k per log-spaced resolution bin
- Polytomies, minimum cluster size constraints, support-aware branch lengths, outlier penalties and more!

---

## Installation

### 1. Recommended: clean conda environment and install with PyPI

```bash
conda create -n phyt_env python=3.10
conda activate phyt_env
pip install phytclust
```

### 2. Install from source

```bash
git clone https://bitbucket.org/schwarzlab/phytclust.git
cd phytclust
pip install -e .[dev]
```

## Command-line usage

### Exact k clusters

Compute an exact k-way clustering, plot it, and save PNG + CSV under ./results:

```bash
phytclust tree.nwk --k 5 --save-fig --out-dir results
```

### Global clustering solution

Search for the top 3 Calinski–Harabasz + Elbow index peaks up to k = 200, save everything in ./out:

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

(See phytclust --help for the complete CLI.)

## Please cite

Please cite this repository if you use the algorithm in your work:

> K. Ganesan, E. Billard, T.L. Kaufmann, C. B Strange, M. C. Cwikla, A. Altenhoff, C. Dessimoz, R.F. Schwarz, PhytClust, (2025), Bitbucket repository, https://bitbucket.org/schwarzlab/phytclust/

```

```
