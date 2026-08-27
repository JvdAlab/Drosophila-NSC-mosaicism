# Script to perform clustering using the scVI trained model

# %% Import libraries
import logging
from pathlib import Path
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scvi
import pyclustree

# %% Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("04_cluster_data.log"),
    ],
)
logger = logging.getLogger(__name__)

# %% Set up seeds
# Set random seeds for reproducibility
sc.settings.seed = 42
scvi.settings.seed = 42
np.random.seed(42)

# %% Set up paths
obj_dir = Path("results") / "scRNA"
obj_dir.mkdir(parents=True, exist_ok=True)
model_dir = obj_dir / "scVI_model_trained"
figures_dir = obj_dir / "figures_clustering"
figures_dir.mkdir(parents=True, exist_ok=True)

# %% Load data
scvi_latent_path = obj_dir / "scRNA_data_QC_filtered_w_scVI_latent.h5ad"
if not scvi_latent_path.exists():
    raise FileNotFoundError(
        f"{scvi_latent_path} not found. Run 03_train_scVI_model.py first."
    )
scVI_adata = ad.read_h5ad(scvi_latent_path)

# %% Setup scVI model
# Setup anndata
scvi.model.SCVI.setup_anndata(
    scVI_adata,
    batch_key="batch",
    categorical_covariate_keys=["replicate"],
    continuous_covariate_keys=["pct_counts_mt"],
)

# Load scVI model
model = scvi.model.SCVI.load(
    dir_path=str(model_dir),
    adata=scVI_adata,
)

# %% Compute neighbors
sc.pp.neighbors(scVI_adata, use_rep="X_scVI", n_neighbors=15)

# %% Compute UMAP
sc.tl.umap(scVI_adata)
sc.tl.tsne(scVI_adata, use_rep="X_scVI", perplexity=30)

# %%
# Plot UMAP
# UMAP and t-SNE plots
umap_plot_path = figures_dir / "umap_scVI.png"
tsne_plot_path = figures_dir / "tsne_scVI.png"

umap_fig = sc.pl.umap(
    scVI_adata,
    color=["batch", "condition", "doublet_prediction"],
    show=False,
    return_fig=True,
)
umap_fig.savefig(umap_plot_path, bbox_inches="tight", dpi=300)
plt.close(umap_fig)
logger.info(f"UMAP plot saved to {umap_plot_path}")

tsne_fig = sc.pl.tsne(
    scVI_adata,
    color=["batch", "condition", "doublet_prediction"],
    show=False,
    return_fig=True,
)
tsne_fig.savefig(tsne_plot_path, bbox_inches="tight", dpi=300)
plt.close(tsne_fig)
logger.info(f"t-SNE plot saved to {tsne_plot_path}")

# %% Cluster at different resolutions
res_list = [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2]

logger.info(f"Testing different clustering resolutions: {res_list}")

for res in res_list:
    cluster_key = f"leiden_r{res}".replace(".", "p")
    logger.info(f"Clustering with resolution {res}")
    sc.tl.leiden(
        scVI_adata,
        resolution=res,
        key_added=cluster_key,
        neighbors_key="neighbors",
        flavor="leidenalg",
    )
    logger.info(f"Number of clusters found: {scVI_adata.obs[cluster_key].nunique()}")
    plt.figure(figsize=(8, 7))
    sc.pl.umap(
        scVI_adata,
        color=cluster_key,
        legend_loc="on data",
        legend_fontsize=8,
        title=f"Leiden clustering at resolution {res}",
        show=False,
    )
    umap_plot_path = figures_dir / f"umap_leiden_{cluster_key}.png"
    plt.savefig(umap_plot_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"Saved UMAP plot to {umap_plot_path}")

# %% Run clustree
clustree_figure = pyclustree.clustree(
    scVI_adata,
    [f"leiden_r{res}".replace(".", "p") for res in res_list],
    title="Clustree of Leiden clustering",
)
clustree_figure_path = figures_dir / "clustree_leiden_clustering.png"
clustree_figure.savefig(clustree_figure_path, bbox_inches="tight", dpi=300)
logger.info(f"Saved clustree plot to {clustree_figure_path}")

# %% Check marker genes
categorized_marker_genes = {
    "Neural_Progenitors": {
        "Neuroblasts": [
            "N",
            "CycE",
            "dpn",
            "mira",
            "wor",
            "ase",
            "insc",
            "stg",
        ],
        "Ganglion_Mother_Cells": [
            "pros",
            "nerfin-1",
            "Hey",
            "jim",
            "ase",
            "tap",
            "dap",
            "insb",
            "insc",
            "spdo",
        ],
    },
    "Neurons": {
        "Newborn_Immature_Neurons": [
            "Hey",
            "E(spl)m6-BFM",
        ],
        "Cholinergic_Neurons": ["ChAT"],
        "GABAergic_Neurons": ["Gad1"],
        "Glutamatergic_Neurons": ["VGlut"],
        "Dopaminergic_Neurons": ["DAT"],
    },
    "Glia": {
        "Astrocytes": [
            "alrm",
            "Eaat1",
            "Gat",
            "Gs2",
            "wun2",
        ],
        "Cortex_Glia": ["wrapper", "zyd", "hoe1"],
        "Surface_Glia": [
            "Mdr65",
            "CG6126",
            "Indy",
            "moody",
            "AdamTS-A",
        ],
    },
}

flat_categorized_markers_for_plotting = {}
all_unique_markers_to_plot = []

# Iterate in the defined order of categorized_marker_genes to preserve it for plotting
for main_cat_label, sub_dict_or_list_or_genes in categorized_marker_genes.items():
    if isinstance(
        sub_dict_or_list_or_genes, dict
    ):  # Main categories like "Neural_Progenitors", "Neurons", "Glia"
        for sub_cat_label, genes_or_deeper_dict in sub_dict_or_list_or_genes.items():
            if isinstance(genes_or_deeper_dict, list):
                present_genes = [
                    g for g in genes_or_deeper_dict if g in scVI_adata.var_names
                ]
                if present_genes:
                    plot_group_key = sub_cat_label
                    if plot_group_key in flat_categorized_markers_for_plotting:
                        plot_group_key = (
                            f"{main_cat_label}_{sub_cat_label}"  # Make unique if needed
                        )
                    flat_categorized_markers_for_plotting[plot_group_key] = (
                        present_genes
                    )
                    all_unique_markers_to_plot.extend(present_genes)
            elif isinstance(genes_or_deeper_dict, dict):
                for (
                    deeper_sub_cat_label,
                    gene_list,
                ) in genes_or_deeper_dict.items():
                    present_genes = [g for g in gene_list if g in scVI_adata.var_names]
                    if present_genes:
                        plot_group_key = deeper_sub_cat_label
                        if plot_group_key in flat_categorized_markers_for_plotting:
                            plot_group_key = f"{main_cat_label}_{sub_cat_label}_{deeper_sub_cat_label}"
                        flat_categorized_markers_for_plotting[plot_group_key] = (
                            present_genes
                        )
                        all_unique_markers_to_plot.extend(present_genes)
    elif isinstance(sub_dict_or_list_or_genes, list):
        present_genes = [
            g for g in sub_dict_or_list_or_genes if g in scVI_adata.var_names
        ]
        if present_genes:
            flat_categorized_markers_for_plotting[main_cat_label] = present_genes
            all_unique_markers_to_plot.extend(present_genes)
all_unique_markers_to_plot = sorted(list(set(all_unique_markers_to_plot)))


# %% Get normalized counts
plotting_layer_key = "log1p_norm"

if plotting_layer_key not in scVI_adata.layers:
    if "counts" in scVI_adata.layers:
        temp_adata = ad.AnnData(
            X=scVI_adata.layers["counts"].copy(),
            obs=scVI_adata.obs.copy(),
            var=scVI_adata.var.copy(),
        )

        sc.pp.normalize_total(
            temp_adata,
            target_sum=1e4,
            inplace=True,
        )

        sc.pp.log1p(temp_adata)

        scVI_adata.layers[plotting_layer_key] = temp_adata.X.copy()
        del temp_adata

        logger.info(f"Created '{plotting_layer_key}' from 'counts' layer.")
    else:
        logger.warning(f"No 'counts' layer found. Using .X for normalization.")

        scVI_adata.layers[plotting_layer_key] = sc.pp.normalize_total(
            scVI_adata, target_sum=1e4, inplace=False
        )["X"]

        scVI_adata.layers[plotting_layer_key] = sc.pp.log1p(
            scVI_adata.layers[plotting_layer_key]
        )

        logger.info(f"Created '{plotting_layer_key}' from .X.")

# %% Plot marker genes (Dotplot)
chosen_cluster_key = "leiden_r0p4"

logger.info(f"Plotting marker genes for cluster key: {chosen_cluster_key}")

plt.figure()
sc.pl.dotplot(
    scVI_adata,
    flat_categorized_markers_for_plotting,
    groupby=chosen_cluster_key,
    layer=plotting_layer_key,
    standard_scale="var",
    show=False,
)

dotplot_path = figures_dir / f"dotplot_marker_genes_{chosen_cluster_key}.png"
num_clusters = scVI_adata.obs[chosen_cluster_key].nunique()
num_genes_in_plot = len(all_unique_markers_to_plot)

plt.gcf().set_size_inches(
    max(10, num_genes_in_plot * 0.35), max(6, num_clusters * 0.35 + 2)
)
plt.xticks(rotation=90)
plt.savefig(dotplot_path, bbox_inches="tight", dpi=300)
plt.close()
logger.info(f"Categorized dot plot saved to {dotplot_path}")

# %% Plot marker genes (Stacked Violin)
plt.figure()
sc.pl.stacked_violin(
    scVI_adata,
    flat_categorized_markers_for_plotting,
    groupby=chosen_cluster_key,
    layer=plotting_layer_key,
    standard_scale="var",
    show=False,
)

stacked_violin_path = (
    figures_dir / f"stacked_violin_marker_genes_{chosen_cluster_key}.png"
)
num_clusters = scVI_adata.obs[chosen_cluster_key].nunique()
num_genes_in_plot = len(all_unique_markers_to_plot)

plt.gcf().set_size_inches(
    max(10, num_genes_in_plot * 0.35), max(6, num_clusters * 0.35 + 2)
)
plt.xticks(rotation=90)
plt.savefig(stacked_violin_path, bbox_inches="tight", dpi=150)
plt.close()
logger.info(f"Stacked violin plot saved to {stacked_violin_path}")


# %% Plot marker genes (UMAP)
cluster_centroids = {}
if chosen_cluster_key in scVI_adata.obs.columns and "X_umap" in scVI_adata.obsm:
    for cluster_id in scVI_adata.obs[chosen_cluster_key].cat.categories:
        cluster_cells = scVI_adata[scVI_adata.obs[chosen_cluster_key] == cluster_id, :]
        if cluster_cells.n_obs > 0:  # Ensure cluster is not empty
            umap_coords = cluster_cells.obsm["X_umap"]
            median_coord = np.median(umap_coords, axis=0)
            cluster_centroids[cluster_id] = median_coord

umap_marker_dir = figures_dir / "umap_marker_genes"
umap_marker_dir.mkdir(parents=True, exist_ok=True)

for gene in all_unique_markers_to_plot:
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(
        scVI_adata,
        color=gene,
        layer=plotting_layer_key,
        show=False,
        title=f"UMAP Expression of {gene}",
        cmap=sns.blend_palette(["lightgray", sns.xkcd_rgb["blood"]], as_cmap=True),
        ax=ax,
        size=20,
    )
    if cluster_centroids:
        for cluster_id, coord in cluster_centroids.items():
            ax.text(
                coord[0],
                coord[1],
                str(cluster_id),
                fontsize=8,
                fontweight="bold",
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.6, ec="none"),
            )

    umap_gene_path = umap_marker_dir / f"umap_{gene}.png"
    plt.savefig(umap_gene_path, bbox_inches="tight", dpi=300)
    plt.close()

    logger.info(f"UMAP plot for {gene} saved to {umap_gene_path}")


# %% Plot marker genes (Heatmap)
# Calculate the number of cells in each cluster
cluster_sizes = scVI_adata.obs[chosen_cluster_key].value_counts()

base_height = 7
extra_height_per_small_cluster = 0.5
small_cluster_threshold = 50

# Calculate additional height for small clusters
extra_height = sum(
    extra_height_per_small_cluster
    for size in cluster_sizes
    if size < small_cluster_threshold
)

# Set the figure size dynamically
figsize = (
    max(8, len(all_unique_markers_to_plot) * 0.25),
    base_height + extra_height,
)

# Plot the heatmap
plt.figure()
sc.pl.heatmap(
    scVI_adata,
    flat_categorized_markers_for_plotting,
    groupby=chosen_cluster_key,
    layer=plotting_layer_key,
    standard_scale="var",
    show=False,
    figsize=figsize,
)
heatmap_path = figures_dir / f"heatmap_marker_genes_{chosen_cluster_key}.png"
plt.savefig(heatmap_path, bbox_inches="tight", dpi=150)
plt.close()
logger.info(f"Heatmap saved to {heatmap_path}")

# %% Annotate clusters
manual_annotation_key = "cell_type_manual"

annotation_mapping = {
    "0": "Immature Neurons",
    "1": "Cholinergic Neurons",
    "2": "Newborn Neurons",
    "3": "GABAergic Neurons",
    "4": "Ganglion Mother Cells (Early)",
    "5": "Immature Neurons",
    "6": "Dopaminergic Neurons",
    "7": "Astrocytes",
    "8": "Intermediate Neural Progenitors",
    "9": "Immature Neurons",
    "10": "Surface Glia",
    "11": "Glutamatergic Neurons",
    "12": "Cholinergic Neurons",
    "13": "Neuroblasts",
    "14": "Ganglion Mother Cells (Late)",
    "15": "Ganglion Mother Cells (Late)",
    "16": "Cortex Glia",
    "17": "Cholinergic Neurons",
    "18": "Contamination (Muscle Cells)",
}

logger.info(f"Annotating clusters with mapping: {annotation_mapping}")

# Ensure the chosen_cluster_key column is string type for mapping
scVI_adata.obs[chosen_cluster_key] = scVI_adata.obs[chosen_cluster_key].astype(str)

# Map the annotations
scVI_adata.obs[manual_annotation_key] = scVI_adata.obs[chosen_cluster_key].map(
    annotation_mapping
)

defined_categories = list(pd.Series(list(annotation_mapping.values())).unique())

scVI_adata.obs[manual_annotation_key] = pd.Categorical(
    scVI_adata.obs[manual_annotation_key],
    categories=defined_categories,
    ordered=False,
)

logger.info(
    f"Manual cell type annotations added to 'scVI_adata.obs[\"{manual_annotation_key}\"]'."
)
logger.info("Verifying annotations (head):")
print(scVI_adata.obs[[chosen_cluster_key, manual_annotation_key]].head())
logger.info("\nValue counts for manual cell types:")
print(scVI_adata.obs[manual_annotation_key].value_counts(dropna=False))

# %% Plot UMAP with manual annotations
logger.info("Plotting UMAP with manual annotations")
plt.figure(figsize=(12, 10))
sc.pl.umap(
    scVI_adata,
    color="cell_type_manual",
    legend_loc="on data",
    legend_fontsize=8,
    legend_fontoutline=2,
    title=f"Manually Annotated Cell Types ({chosen_cluster_key})",
    show=False,
)
umap_manual_annot_path = (
    figures_dir / f"umap_manual_annotation_{chosen_cluster_key}.png"
)
plt.savefig(umap_manual_annot_path, bbox_inches="tight", dpi=150)
plt.close()
logger.info(f"UMAP with manual annotations saved to {umap_manual_annot_path}")

# %% Marker genes for annotated clusters
logger.info(f"Using layer: {plotting_layer_key} for marker gene plots")

marker_genes_dir = figures_dir / "marker_genes"
marker_genes_dir.mkdir(parents=True, exist_ok=True)

sc.tl.rank_genes_groups(
    scVI_adata,
    groupby=manual_annotation_key,
    method="wilcoxon",
    use_raw=False,
    layer=plotting_layer_key,
    pts=True,
    key_added=f"rank_genes_{manual_annotation_key}",
)

marker_genes = scVI_adata.uns[f"rank_genes_{manual_annotation_key}"]
groups = marker_genes["names"].dtype.names

for group in groups:
    marker_data_to_save = []
    for i in range(len(marker_genes["names"][group])):
        marker_data_to_save.append(
            {
                "cell_type": group,
                "gene": marker_genes["names"][group][i],
                "logfoldchanges": marker_genes["logfoldchanges"][group][i],
                "pvals": marker_genes["pvals"][group][i],
                "pvals_adj": marker_genes["pvals_adj"][group][i],
                "scores": marker_genes["scores"][group][i],
                "pts": marker_genes["pts"][group][i],
            }
        )

    marker_genes_df = pd.DataFrame(marker_data_to_save)
    marker_genes_output_file = marker_genes_dir / f"marker_genes_{group}.tsv"
    marker_genes_df.to_csv(
        marker_genes_output_file,
        sep="\t",
        index=False,
    )

    logger.info(f"Marker genes for {group} saved to {marker_genes_output_file}.")

# %% Remove contaminated cells
print(f"Cells before filtering: {scVI_adata.n_obs}")
scVI_adata = scVI_adata[
    scVI_adata.obs["cell_type_manual"] != "Contamination (Muscle Cells)"
].copy()
print(f"Cells after filtering: {scVI_adata.n_obs}")

# %% Save annotated data
adata_annotated_save_path = (
    obj_dir / "scRNA_data_QC_filtered_w_scVI_latent_annotated.h5ad"
)
scVI_adata.write(
    adata_annotated_save_path,
    compression="gzip",
)
logger.info(f"Annotated data saved to { adata_annotated_save_path}")
