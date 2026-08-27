# Script to create pseudobulks for DEG analysis using edgeR

# %% Import libraries
import logging
from pathlib import Path
import anndata as ad
import pandas as pd

# %% Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("00_create_pseudobulk.log"),
    ],
)
logger = logging.getLogger(__name__)

# %% Set up paths
# input_obj_dir: produced by scRNA_pipeline/04_cluster_data.py (see scRNA_pipeline/README.md)
input_obj_dir = Path("results") / "scRNA"
pseudobulk_dir = Path("results") / "edgeR" / "pseudobulk_all_genes"

pseudobulk_dir.parent.mkdir(parents=True, exist_ok=True)
pseudobulk_dir.mkdir(parents=True, exist_ok=True)


# %% Load the data
annotated_path = input_obj_dir / "scRNA_data_QC_filtered_w_scVI_latent_annotated.h5ad"
if not annotated_path.exists():
    raise FileNotFoundError(
        f"{annotated_path} not found. Run scRNA_pipeline/04_cluster_data.py first."
    )
scVI_adata = ad.read_h5ad(annotated_path)

# %% Create pseudobulk data
# Sanitize cell type names
scVI_adata.obs["cell_type_manual"] = scVI_adata.obs["cell_type_manual"].str.replace(
    " ", "_", regex=False
)

# Convert to categorical variables
scVI_adata.obs["sample_id"] = scVI_adata.obs["sample_id"].astype("category")
scVI_adata.obs["batch"] = scVI_adata.obs["batch"].astype("category")
scVI_adata.obs["cell_type_manual"] = scVI_adata.obs["cell_type_manual"].astype(
    "category"
)

# Convert counts to integers
if scVI_adata.layers["counts"].dtype != "int32":
    scVI_adata.layers["counts"] = scVI_adata.layers["counts"].astype("int32")

# %%
unique_cell_types = scVI_adata.obs["cell_type_manual"].cat.categories.tolist()
# Filter out "Unknown" cell types
unique_cell_types = [ct for ct in unique_cell_types if ct != "Unknown"]

for cell_type in unique_cell_types:
    print(f"Processing cell type: {cell_type}")
    outdir = pseudobulk_dir / cell_type
    outdir.mkdir(parents=True, exist_ok=True)

    # Filter the data for the current cell type
    adata_ct = scVI_adata[scVI_adata.obs["cell_type_manual"] == cell_type].copy()

    # Aggregate counts
    pseudobulk_profiles = {}
    unique_sample_ids = adata_ct.obs["sample_id"].unique()

    for s_id in unique_sample_ids:
        sample_mask_ct = adata_ct.obs["sample_id"] == s_id
        adata_ct_sample = adata_ct[sample_mask_ct].copy()

        summed_counts_vector = adata_ct_sample.layers["counts"].sum(axis=0).A1

        pseudobulk_profiles[s_id] = summed_counts_vector

    pseudobulk_counts_df = pd.DataFrame(pseudobulk_profiles, index=scVI_adata.var.index)
    pseudobulk_counts_df = pseudobulk_counts_df.astype(int)

    # Export count matrix
    pseudobulk_counts_df.to_csv(
        outdir / f"{cell_type}_pseudobulk_counts.tsv",
        sep="\t",
        index=True,
        index_label="gene",
        header=True,
    )
    logger.info(
        f"Exported pseudobulk counts for {cell_type} to {outdir / f'{cell_type}_pseudobulk_counts.tsv'}"
    )

    # Create and export metadata for R
    sample_info_list = []
    for s_id in pseudobulk_counts_df.columns:
        filter_cell_for_sample = scVI_adata.obs[
            scVI_adata.obs["sample_id"] == s_id
        ].iloc[0]
        sample_info_list.append(
            {
                "sample_id": s_id,
                "condition": filter_cell_for_sample["condition"],
                "batch": filter_cell_for_sample["batch"],
            }
        )

    sample_info_df = pd.DataFrame(sample_info_list).set_index("sample_id")

    sample_info_df.to_csv(
        outdir / f"sample_info_{cell_type}.tsv",
        sep="\t",
        index=True,
        index_label="sample_id",
        header=True,
    )

    logger.info(
        f"Exported sample info for {cell_type} to {outdir / f'sample_info_{cell_type}.tsv'}"
    )

    # Export gene metadata
    gene_info_df = scVI_adata.var.copy()
    gene_info_df.to_csv(
        outdir / f"gene_info_{cell_type}.tsv",
        sep="\t",
        index=True,
        index_label="gene_id",
        header=True,
    )

    logger.info(
        f"Exported gene info for {cell_type} to {outdir / f'gene_info_{cell_type}.tsv'}"
    )

    logger.info(
        f"Completed processing for cell type: {cell_type}. Pseudobulk data saved to {outdir}"
    )

logger.info("All cell types processed successfully.")

# %%
logger.info("\n" + "=" * 60)
logger.info("Creating combined glial pseudobulk...")
logger.info("=" * 60)

glial_cell_types = ["Astrocytes", "Cortex_Glia", "Surface_Glia"]
combined_name = "Combined_Glia"

logger.info(f"Combining glial cell types: {', '.join(glial_cell_types)}")

# Check which glial cell types are present
available_glia = [
    ct
    for ct in glial_cell_types
    if ct in scVI_adata.obs["cell_type_manual"].cat.categories
]
missing_glia = [
    ct
    for ct in glial_cell_types
    if ct not in scVI_adata.obs["cell_type_manual"].cat.categories
]

if missing_glia:
    logger.warning(
        f"Warning: The following glial cell types were not found: {', '.join(missing_glia)}"
    )

logger.info(f"Found the following glial cell types: {', '.join(available_glia)}")

if len(available_glia) > 0:
    # Filter for glial cells
    adata_glia = scVI_adata[
        scVI_adata.obs["cell_type_manual"].isin(available_glia)
    ].copy()

    logger.info(f"Total glial cells: {adata_glia.n_obs}")
    for ct in available_glia:
        n_cells = (adata_glia.obs["cell_type_manual"] == ct).sum()
        logger.info(f"  {ct}: {n_cells} cells")

    # Create output directory
    outdir = pseudobulk_dir / combined_name
    outdir.mkdir(parents=True, exist_ok=True)

    # Aggregate counts by sample
    pseudobulk_profiles = {}
    unique_sample_ids = adata_glia.obs["sample_id"].unique()

    for s_id in unique_sample_ids:
        sample_mask = adata_glia.obs["sample_id"] == s_id
        adata_sample = adata_glia[sample_mask].copy()
        summed_counts_vector = adata_sample.layers["counts"].sum(axis=0).A1
        pseudobulk_profiles[s_id] = summed_counts_vector

        n_cells = sample_mask.sum()
        logger.info(f"  Sample {s_id}: {n_cells} glial cells aggregated")

    # Create pseudobulk count matrix
    pseudobulk_counts_df = pd.DataFrame(pseudobulk_profiles, index=scVI_adata.var.index)
    pseudobulk_counts_df = pseudobulk_counts_df.astype(int)

    # Export count matrix
    pseudobulk_counts_df.to_csv(
        outdir / f"{combined_name}_pseudobulk_counts.tsv",
        sep="\t",
        index=True,
        index_label="gene",
        header=True,
    )
    logger.info(
        f"Exported pseudobulk counts to {outdir / f'{combined_name}_pseudobulk_counts.tsv'}"
    )

    # Create and export sample metadata
    sample_info_list = []
    for s_id in pseudobulk_counts_df.columns:
        filter_cell_for_sample = scVI_adata.obs[
            scVI_adata.obs["sample_id"] == s_id
        ].iloc[0]
        n_glia_cells = (adata_glia.obs["sample_id"] == s_id).sum()

        sample_info_list.append(
            {
                "sample_id": s_id,
                "condition": filter_cell_for_sample["condition"],
                "batch": filter_cell_for_sample["batch"],
                "n_glial_cells": n_glia_cells,
            }
        )

    sample_info_df = pd.DataFrame(sample_info_list).set_index("sample_id")
    sample_info_df.to_csv(
        outdir / f"sample_info_{combined_name}.tsv",
        sep="\t",
        index=True,
        index_label="sample_id",
        header=True,
    )
    logger.info(
        f"Exported sample info to {outdir / f'sample_info_{combined_name}.tsv'}"
    )

    # Export gene metadata
    gene_info_df = scVI_adata.var.copy()
    gene_info_df.to_csv(
        outdir / f"gene_info_{combined_name}.tsv",
        sep="\t",
        index=True,
        index_label="gene_id",
        header=True,
    )
    logger.info(f"Exported gene info to {outdir / f'gene_info_{combined_name}.tsv'}")

    # Create cell type composition summary
    composition_list = []
    for s_id in unique_sample_ids:
        sample_glia = adata_glia[adata_glia.obs["sample_id"] == s_id]
        total_cells = len(sample_glia)

        comp_dict = {"sample_id": s_id, "total_glial_cells": total_cells}

        for ct in available_glia:
            n_cells = (sample_glia.obs["cell_type_manual"] == ct).sum()
            prop = n_cells / total_cells if total_cells > 0 else 0
            comp_dict[f"{ct}_n"] = n_cells
            comp_dict[f"{ct}_prop"] = prop

        composition_list.append(comp_dict)

    composition_df = pd.DataFrame(composition_list).set_index("sample_id")
    composition_df.to_csv(
        outdir / f"cell_type_composition_{combined_name}.tsv",
        sep="\t",
        index=True,
        index_label="sample_id",
        header=True,
    )
    logger.info(
        f"Exported cell type composition to {outdir / f'cell_type_composition_{combined_name}.tsv'}"
    )

    logger.info(f"Successfully created combined glial pseudobulk in {outdir}")
else:
    logger.warning(
        "No glial cell types found in the data. Skipping combined glial analysis."
    )

logger.info("=" * 60 + "\n")

# %%
logger.info("\n" + "=" * 60)
logger.info("Creating combined GMC pseudobulk...")
logger.info("=" * 60)

# Define the GMC subsets to combine
gmc_cell_types = ["Ganglion_Mother_Cells_(Early)", "Ganglion_Mother_Cells_(Late)"]
combined_gmc_name = "Combined_GMCs"

logger.info(f"Combining GMC subsets: {', '.join(gmc_cell_types)}")

# Check availability in the dataset
available_gmcs = [
    ct
    for ct in gmc_cell_types
    if ct in scVI_adata.obs["cell_type_manual"].cat.categories
]
missing_gmcs = [
    ct
    for ct in gmc_cell_types
    if ct not in scVI_adata.obs["cell_type_manual"].cat.categories
]

if missing_gmcs:
    logger.warning(
        f"Warning: The following GMC subsets were not found: {', '.join(missing_gmcs)}"
    )

if len(available_gmcs) > 0:
    # Filter for GMCs
    adata_gmcs = scVI_adata[
        scVI_adata.obs["cell_type_manual"].isin(available_gmcs)
    ].copy()

    logger.info(f"Total GMC cells: {adata_gmcs.n_obs}")
    for ct in available_gmcs:
        n_cells = (adata_gmcs.obs["cell_type_manual"] == ct).sum()
        logger.info(f"  {ct}: {n_cells} cells")

    # Create output directory
    outdir_gmc = pseudobulk_dir / combined_gmc_name
    outdir_gmc.mkdir(parents=True, exist_ok=True)

    # Aggregate counts by sample
    pseudobulk_profiles_gmc = {}
    unique_sample_ids = adata_gmcs.obs["sample_id"].unique()

    for s_id in unique_sample_ids:
        sample_mask = adata_gmcs.obs["sample_id"] == s_id
        adata_sample = adata_gmcs[sample_mask].copy()
        summed_counts_vector = adata_sample.layers["counts"].sum(axis=0).A1
        pseudobulk_profiles_gmc[s_id] = summed_counts_vector

        n_cells = sample_mask.sum()
        logger.info(f"  Sample {s_id}: {n_cells} GMCs aggregated")

    # Create and export count matrix
    pseudobulk_counts_gmc_df = pd.DataFrame(
        pseudobulk_profiles_gmc, index=scVI_adata.var.index
    ).astype(int)

    pseudobulk_counts_gmc_df.to_csv(
        outdir_gmc / f"{combined_gmc_name}_pseudobulk_counts.tsv",
        sep="\t",
        index=True,
        index_label="gene",
    )

    # Export sample metadata with GMC-specific counts
    sample_info_gmc = []
    for s_id in pseudobulk_counts_gmc_df.columns:
        meta = scVI_adata.obs[scVI_adata.obs["sample_id"] == s_id].iloc[0]
        n_cells = (adata_gmcs.obs["sample_id"] == s_id).sum()
        sample_info_gmc.append(
            {
                "sample_id": s_id,
                "condition": meta["condition"],
                "batch": meta["batch"],
                "n_gmc_cells": n_cells,
            }
        )

    pd.DataFrame(sample_info_gmc).set_index("sample_id").to_csv(
        outdir_gmc / f"sample_info_{combined_gmc_name}.tsv", sep="\t"
    )

    # Export gene info for consistency
    scVI_adata.var.to_csv(
        outdir_gmc / f"gene_info_{combined_gmc_name}.tsv",
        sep="\t",
        index_label="gene_id",
    )

    logger.info(f"Successfully created combined GMC pseudobulk in {outdir_gmc}")
else:
    logger.warning("No GMC cell types found. Skipping.")

logger.info("=" * 60 + "\n")
# %%
