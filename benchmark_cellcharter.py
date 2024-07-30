# import the necessary packages
import squidpy as sq
import cellcharter as cc
import scanpy as sc
import anndata
import pandas as pd
from pathlib import Path
from argparse import ArgumentParser


def get_argumnets():
    # Start with the description
    description = "Run cellcharter on a h5ad file"

    # Create the parser
    parser = ArgumentParser(description=description)

    # Add a group of arguments for the input
    input_group = parser.add_argument_group("Input arguments", "Paths for the input data")
    input_group.add_argument("-i", "--input", required=True, help="Path to the input h5ad file")

    # Add a group of arguments for the tool
    tool_group = parser.add_argument_group("Tool parameters", "Parameters for CellCharter")
    tool_group.add_argument("-s", "--steps", type=int, default=1, help="Number of steps to aggregate neighbors")
    tool_group.add_argument("-k", "--k_range", nargs=2, type=int, default=(2, 10), help="Range of k values to test")
    tool_group.add_argument("-m", "--max_runs", type=int, default=10, help="Maximum number of runs for each k")
    tool_group.add_argument("-a", "--accelerator", default="cpu", help="Accelerator to use for training")

    # Add a group of arguments for the output
    output_group = parser.add_argument_group("Output arguments", "Paths for the output data")
    output_group.add_argument("-o", "--output", required=True, help="Path to the output directory")

    # Parse the arguments
    args = parser.parse_args()

    # Standardize the paths
    args.k_range = tuple(args.k_range)
    args.input = Path(args.input).resolve()
    args.output = Path(args.output).resolve()

    # Return the arguments
    return args


def main():
    # Get arguments
    arg = get_argumnets()

    # Set the paths
    output_model_path = arg.output / f"{arg.input.stem}_l{arg.steps}_k{arg.k_range[0]}-{arg.k_range[1]}_model.pkl"
    print(f"Using steps = {arg.steps}")
    print(f"K range     = {arg.k_range}")
    print(f"Max runs    = {arg.max_runs}")
    print(f"Input data  = {arg.input}")
    print(f"Output model= {output_model_path}")


    # Read in data
    print("Reading in data")
    adata = anndata.read_h5ad(arg.input)

    # Add spatial to uns based on the sample name
    adata.uns["spatial"]={sample:{} for sample in adata.obs["sample"].cat.categories}

    # Normalize each sample independently
    print("Normalizing data")
    adata.raw = adata.copy()
    for sample in adata.obs['sample'].cat.categories:
        adata.X[adata.obs['sample'] == sample, :] = sc.pp.scale(adata[adata.obs['sample'] == sample], copy=True).X

    # Compute the delunay graph using the spatial coordinates
    print("Computing delaunay graph")
    sq.gr.spatial_neighbors(adata, library_key="sample", coord_type='generic', delaunay=True)

    # Cut the long links from the delaunay graph
    print("Removing long links")
    cc.gr.remove_long_links(adata)

    # Aggregate neighbors for each layer
    print("Aggregating neighbors")
    cc.gr.aggregate_neighbors(adata, n_layers=arg.steps)

    # Set the model parameters
    print("Running auto K")
    model_params = {
            'random_state': 42,
            'trainer_params': {
                'accelerator': arg.accelerator,
                'enable_progress_bar': False
            },
        }

    # Run the cluster auto K
    autok = cc.tl.ClusterAutoK(n_clusters=arg.k_range,
                               model_class=cc.tl.GaussianMixture,
                               model_params=model_params,
                               max_runs=arg.max_runs)
    autok.fit(adata, use_rep='X_cellcharter')

    # Save the autok model
    print("Saving auto K model")
    autok.save(output_model_path)

    # Predict using the best k
    print(f"Predicting spatial clusters with best k = {autok.best_k}")
    adata.obs['spatial_cluster'] = autok.predict(adata, use_rep='X_cellcharter', k=autok.best_k)

    # Get the number of cells in each spatial cluster for each sample
    print("Computing spatial cluster frequencies")
    output = pd.DataFrame(
        adata.obs.groupby("sample", observed=True).apply(
            lambda x: x["spatial_cluster"].value_counts()/len(x), include_groups=False
        )
    )

    # Convert output to a wide format using spatial_cluster as columns
    output = output.unstack(level=-1)

    # Save the output to a csv file
    print("Saving spatial cluster frequencies")
    output.columns = output.columns.get_level_values(1)
    output = output.rename(columns={c: f"{output.columns.name}_{c}" for c in output.columns.values})
    out_sc_path = arg.output / f"{arg.input.stem}_l{arg.steps}_k{autok.best_k}.csv"
    output.to_csv(out_sc_path)

    # save the .h5ad file with the spatial clusters
    print("Saving h5ad file")
    out_adata_path = arg.output / f"{arg.input.stem}_l{arg.steps}_k{autok.best_k}.h5ad"
    output.to_csv(out_adata_path)


if __name__ == "__main__":
    main()
