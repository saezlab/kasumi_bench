from utag import utag
import scanpy as sc
import squidpy as sq
import anndata
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse
from argparse import ArgumentParser


def get_arguments():
    # Start with the description
    description = "Run utag on a h5ad file."

    # Create the parser
    parser = ArgumentParser(description=description)

    # Add a group of arguments for the input
    input_group = parser.add_argument_group("Input arguments", "Paths for the input data")
    input_group.add_argument("-i", "--input", required=True,
                             help="Path to the input h5ad file")

    # Add a group of arguments for the tool
    tool_group = parser.add_argument_group("Tool parameters", "Parameters for utag")
    tool_group.add_argument("-d", "--max_dist", type=int, default=20,
                            help="Maximum distance for neighbors [default: 20]")
    tool_group.add_argument("-n", "--normalization", default="l1_norm",
                            help="Normalization method [default: l1_norm]")
    tool_group.add_argument("-c", "--clustering", default="leiden",
                            help="Clustering method [default: leiden]")
    tool_group.add_argument("-r", "--resolutions", nargs='+', type=float,
                            default=[0.5, 0.6, 0.7, 0.8, 0.9],
                            help="Resolutions for clustering [default: 0.5, 0.6, 0.7, 0.8, 0.9]")

    # Add a group of arguments for the output
    output_group = parser.add_argument_group("Output arguments", "Paths for the output data")
    output_group.add_argument("-o", "--output", required=True,
                              help="Path to the output folder")

    # Parse the arguments
    args = parser.parse_args()

    # Standardize the paths
    args.input = Path(args.input).resolve()
    args.output = Path(args.output).resolve()

    return args


def get_relative_spatial_clusters_per_sample(utag_results, sample_key="sample", resolution=0.5):
    # Get the number of cells in each spatial cluster for each sample
    output = pd.DataFrame(
        utag_results.obs.groupby("sample", observed=True).apply(
            lambda x: x[f"UTAG Label_leiden_{resolution}"].value_counts() / len(x), include_groups=False
        )
    )

    # Convert output to a wide format using spatial_cluster as columns
    output = output.unstack(level=1)
    output.columns = output.columns.get_level_values(1)
    output = output.rename(columns={c: f"{c}" for c in output.columns.values})

    return output

def main():
    # Get arguments
    arg  = get_arguments()

    # Read in data
    print("Reading in data")
    adata = anndata.read_h5ad(arg.input)

    # Run UTAG on provided data
    print("Running UTAG")
    utag_results = utag(
        adata,
        slide_key="sample",
        max_dist=arg.max_dist,
        normalization_mode=arg.normalization,
        apply_clustering=True,
        clustering_method=arg.clustering,
        resolutions=arg.resolutions
    )

    # Save results
    print("Saving results")
    for resolution in arg.resolutions:
        output_path = arg.output / f"{arg.input.stem}_{resolution}.csv"
        get_relative_spatial_clusters_per_sample(utag_results, resolution=resolution).to_csv(output_path)


if __name__ == "__main__":
    main()
