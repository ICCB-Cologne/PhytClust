# PhyloTreeClustering 
Clustering algorithm for phylogenetic cancer evolution trees.

## Installation
Install PhyloTreeClustering via pip or from source.

### Intallation via pip (recommended)
It is best to use a dedicated conda environment for your PhyloTreeClustering installation with `conda create -n phylo_env`.
After activating the environment with `conda activate phylo_env` you can install PhyloTreeCLustering via `pip install phylotreeclus`

### Installation from source
Clone the PhyloTreeClustering repository using `git clone https://bitbucket.org/schwarzlab/phylotreeclustering`. 

Then, inside the PhyloTreeClustering folder, run `pip install . `to install PhyloTreeClustering to your environment.

## Usage
After installing PhyloTreeClustering, you can use PhyloTreeClustering functions in python scripts (through `import PhyloTreeClustering`) and from the command line. General usage from the command line is `phylotreeclus path/to/input/file path/to/output/folder`. Run `phylotreeclus --help` for information on optional arguments.

## Command line Flags

-   `tree_file_path`: tree you want to cluster, in Newick format.

-   `output_path`: path to save figures and files produced by algorithm.

-   `tree`: directly the tree you want to cluster.

-   `dist_matrix`: panda dataframe, matrix of the distances between all leaves of the tree, index are the names of the leaves 
                The default None means the distance matrix will be calculated by the function `calculate_pdm_from_tree`.

-   `threshold`: manually choose the cutting threshold (not recommended).

-   `calinski_harabasz_score`: if True, the algorithm chooses the best threshold based on the clustering with the highest score,
                                True by default.

-   `n_diff_th`: number of different threshold that will be tested.
                The default None means this number will be calculated according to the size of the tree in `init_n_diff_th`.
            
-   `CH_matrix`: distance matrix that will be used to evaluate the clusering with the CH score,
                if no matrix is given, the normal `dist_matrix` will be used.
        
-   `min_size_clus`: minimum size of a cluster. 
                If no value is given, it is calculated in `init_min_size_clus` according to the size of the tree.
        
-   `multi_level_clus`: if True, a second level of clustering will be performed for the bigger clusters.
                The default is False, so only one level of clustering will be performed.
        
-   `dist_type`: distance calculation used to compute the new dendrogram distances.
                The default and recommended distance is `half_max`. The others are `mean`, `median` and `max`.

## Output plot

PhyloTreeClustering creates the following output files:
-   summary_plot.png : dendrogram and tree next to each other, with the clustering shown with colors (outliers in black), the red dotted lines are the min and max threshold     tested, the yellow line is the chosen threshold.
-   sample_labels.txt : list of all the samples and the cluster they were assigned to.
-   sample_labels_2nd_layer.txt : list of all the samples and the cluster they were assigned to if the option multi_layer_clustering is chosen.

## Usage example

For first time users we recommend to have a look at the trees (Newick format) in /examples on Bitbucket to get an idea of how input data should look like. Then run `phylotreeclus examples/Gao_2016_12_final_tree.new path/to/output/folder` as an example of a standard PhyloTreeClustering run. Finally, the notebook `notebooks/simple_example.py` shows how the individual functions in the workflow are used.

![Clustering example](https://bitbucket.org/schwarzlab/phylotreeclustering/raw/master/example_data/example_summary_plot_readme.png)

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

**For the file medicc_functions, copied functions from :**
Kaufmann, T.L., Petkovic, M., Watkins, T.B.K. et al.
MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution.
Genome Biol 23, 241 (2022). https://doi.org/10.1186/s13059-022-02794-9

**For the file calinski_harabasz, copied and modified function from :**
Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.

