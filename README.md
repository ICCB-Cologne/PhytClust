# PhytClust 
Monophyletic Clustering Algorithm for Phylogenetic Trees

### Intallation via pip (recommended)
It is best to use a dedicated conda environment for your PhytClust installation with `conda create -n phyt_env`.
After activating the environment with `conda activate phylo_env` you can install PhyloTreeCLustering via `pip install phyclust`

### Installation from source
Clone the PhytClust repository using `git clone -b Phytclust https://bitbucket.org/schwarzlab/phylotreeclustering`

Then, inside the PhytClust folder, run `pip install . `to install PhytClust to your environment.

## Usage
An example notebook and input trees are provided. As an alternative, a newick string can also be used as input

PhytClust creates the following output files:
-   tree.png: image of the tree with nodes coloured by cluster number. 
-   results.csv : list of all terminals nodes and the cluster they were assigned to. -1 refers to outliers

## References

**For the examples :**

- Original data from Gao et al. 2016 
Gao, R., Davis, A., McDonald, T. et al. 
Punctuated copy number evolution and clonal stasis in triple-negative breast cancer. 
Nat Genet 48, 1119–1130 (2016). https://doi.org/10.1038/ng.3641

- Original data from Minussi et al. 2021
Minussi, D.C., Nicholson, M.D., Ye, H. et al. 
Breast tumours maintain a reservoir of subclonal diversity during expansion. 
Nature 592, 302–308 (2021). https://doi.org/10.1038/s41586-021-03357-x

## Please cite

Please cite this repository if you use the algorithm in your work:

> E. Billard, T.L. Kaufmann, R.F. Schwarz, PhyloTreeClus, (2023), Bitbucket repository, https://bitbucket.org/schwarzlab/phylotreeclustering


