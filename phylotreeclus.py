#!/usr/bin/env python3

from argparse import ArgumentParser
from PhyloTreeClustering.PhyloTreeClustering import PhyloTreeClustering

if __name__ == "__main__":
    parser = ArgumentParser(
        # name of program
        prog="Dendrogram Clustering",
        # description
        description="Takes a phylogeny tree, transforms it in a dendrogram, and clusters the leaves on the dendrogram",
        # bottom text
        epilog="By Elisa Billard"
    )
    # positional arguments  
    parser.add_argument(
        dest= "tree_file_path",
        help="The input file : the copy number data"
    )
    parser.add_argument(
        dest= "output_path",
        help="The output path to save the figure"
    )

    # optional arguments
    parser.add_argument(
        "--th", dest= "threshold",
        help="The threshold where to cut the dendrogram"
    )
    parser.add_argument(
        "-c", "--calinski_h", action="store_false", dest = "calinski_harabasz_score",
        help="Whether to test different threshold and choose the final one based on CH score"
    )
    parser.add_argument(
        "-n", "--n_diff_th", dest = "n_diff_th",
        help="The number of different thresholds that will be tested in the tuning"
    )
    parser.add_argument(
        "-s", "--min_size_clus", dest = "min_size_clus",
        help="The minimal size of a cluster, in number of leaves"
    )
    parser.add_argument(
        "--dist_type", dest = "dist_type",
        help="The distance type used to calculate the dendrogram distances"
    )
    parser.add_argument(
        "-m", "--multi_level_clus", action= "store_true", # flag, meaning it is False by default
        dest = "multi_level_clus",
        help="Whether to to a 2nd level of clustering on the big clusters"
    )
    
    
    args = parser.parse_args()

    clus = PhyloTreeClustering(**vars(args))
    clus.run()
    clus.plot_summary()
    clus.sample_labels()