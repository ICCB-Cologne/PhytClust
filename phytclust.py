#!/usr/bin/env python3
import argparse
import os
import sys
from Bio import Phylo
from phytclust import PhytClust


def main():
    parser = argparse.ArgumentParser(
        description="Phylogenetic clustering using dynamic programming and resolution-based peak selection."
    )
    parser.add_argument(
        "--tree",
        "-t",
        type=str,
        required=True,
        help="Path to the input Newick tree file.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=None,
        help="Number of clusters to use. If not provided, the algorithm will try 1..max_k.",
    )
    parser.add_argument(
        "--max_k",
        type=int,
        default=None,
        help="Maximum number of clusters to consider. If not provided and k is not set, defaults to 90%% of terminal nodes.",
    )
    parser.add_argument(
        "--resolution_on",
        action="store_true",
        help="Turn on resolution binning when selecting peaks.",
    )

    parser.add_argument(
        "--plot_scores",
        type=bool,
        default=True,
        help="Whether scores should be plotted (default: True).",
    )

    parser.add_argument(
        "--num_bins",
        type=int,
        default=3,
        help="Number of resolution bins to use (default: 3).",
    )
    parser.add_argument(
        "--global_peaks",
        type=int,
        default=3,
        help="Number of global peaks to select if resolution is off (default: 3).",
    )
    parser.add_argument(
        "--peaks_per_bin",
        type=int,
        default=1,
        help="Number of peaks to select per bin when resolution is on (default: 1).",
    )
    parser.add_argument(
        "--min_k",
        type=int,
        default=3,
        help="Minimum k value to consider a peak as a solution (default: 3).",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="results",
        help="Directory where results will be saved (default: results).",
    )
    parser.add_argument("--no_plot", action="store_true", help="Do not display plots.")
    args = parser.parse_args()

    try:
        tree = Phylo.read(args.tree, "newick")
    except Exception as e:
        print(f"Error reading tree file {args.tree}: {e}", file=sys.stderr)
        sys.exit(1)

    clustering = PhytClust(
        tree=tree,
        k=args.k,
        should_plot_scores=args.plot_scores,
        max_k=args.max_k,
        resolution_on=args.resolution_on,
        num_bins=args.num_bins,
    )

    clustering.save(results_dir=args.output_dir)


if __name__ == "__main__":
    main()
