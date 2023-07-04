# README #

# PhyloTreeClustering 
Clustering method for phylogenetic trees.

## Installation
Install PhyloTreeClustering via conda.

## Usage
After installing PhyloTreeClustering, you can use PhyloTreeClustering functions in python scripts (through import PhyloTreeClustering) and from the command line. General usage from the command line is PhyloTreeClustering path/to/input/file path/to/output/folder. Run PhyloTreeClustering --help for information on optional arguments.

## Command line Flags

-   tree_file_path: tree to cluster, in new format

-   output_path: path to save figure produced by algorithm

-   tree: directly tree to cluster

-   dist_matrix: panda dataframe, matrix of the distances between all leaves of the tree, index are the names of the leaves 
                The default None means the distance matrix will be calculated by the function calculate_pdm_from_tree.

-   threshold: manually choose the cutting threshold (not recommended)

-   calinski_harabasz_score: if True, the algorithm chooses the best threshold based on the clustering with the highest score

-   n_diff_th: number of repetition of scoring, meaning the number of different threshold cuts will be tested.
                The default None means this number will be calculated according to the size of the tree in init_n_diff_th.
            
-   CH_matrix: distance matrix that will be used to evaluate the clusering with the CH score
                If no matrix is given, the normal dist_matrix will be used
        
-   min_size_clus: minimum size of a cluster. 
                If no value is given, it is calculated in init_min_size_clus according to the size of the tree
        
-   multi_level_clus: if True, make a second level of clusters.
                The default is False, so only one level of clustering will be performed
        
-   dist_type: distance calculation used to compute the new dendrogram distances.
                The default and recommended distance is 'half_max'. The others are 'mean', 'median' and 'max'

## Output plot

PhyloTreeClustering creates the following output files:
-   summary_plot.png : dendrogram and medicc tree next to each other, with the clustering
-   

## Usage example
For first time users we recommend to have a look at examples/simple_example to get an idea of how input data should look like. Then run PhyloTreeClustering examples/simple_example/simple_example.tsv path/to/output/folder as an example of a standard PhyloTreeClustering run. Finally, the notebook notebooks/example_workflows.py shows how the individual functions in the workflow are used.