# PhytClust &nbsp;<img alt="PyPI" src="https://img.shields.io/pypi/v/phytclust?color=brightgreen"> <img alt="Tests" src="https://img.shields.io/github/actions/workflow/status/schwarzlab/phytclust/ci.yml?label=CI&logo=github">

Monophyletic, dynamic-programming **clustering of phylogenetic trees** –
pick an _exact_ `k`, search for the _global_ Calinski-Harabasz/Bowker
(“CalBow”) peaks, or obtain one representative peak **per log-resolution
bin**.

_This is a modernised fork of the original “PhyloTreeClustering”
algorithm by Schwarz lab._

---

## Installation

### 1. Recommended: clean conda environment + PyPI

```bash
conda create -n phyt_env python=3.10
conda activate phyt_env
pip install phytclust

### 2. Install from source

git clone https://bitbucket.org/schwarzlab/phytclust.git
cd phytclust
pip install -e .[dev]         # editable install + dev extras (black, pytest, nbdev…)


## Output Files


## API Reference

## References



## Please cite

Please cite this repository if you use the algorithm in your work:

> K. Ganesan, E. Billard, T.L. Kaufmann, R.F. Schwarz, PhytClust, (2024), Bitbucket repository, https://bitbucket.org/schwarzlab/phytclust/

```


from phytclust.metrics import normalized_colless, calculate_total_length_variation
from phytclust.selection import select_representative_species

nc = normalized_colless(tree)
rt_tip_sd = calculate_total_length_variation(tree)

reps = select_representative_species(tree, clusters, mode="maximize", distance_ref="mrca")


# exact k=5, show plot, save PNG+CSV under ./results
phytclust tree.nwk k --k 5 --save-fig

# peak search: top 3, cap k at 200, save everything in ./out
phytclust tree.nwk global --top-n 3 --max-k 200 --save-fig --out-dir out

# one peak per 4 log-bins, don’t show plots interactively, save all k CSVs
phytclust tree.nwk resolution --bins 4 --no-plot --save-all-k --save-fig
