# Script to count, per gene, how many cells express it in each condition
# within each pseudobulk group. Feeds the is_artifact flag in 02_edgeR_DGE.R.

# %% Import libraries
import logging
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

# %% Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("01_create_pct_expressing.log"),
    ],
)
logger = logging.getLogger(__name__)

# %% Set up paths
# input_obj_dir: produced by scRNA_pipeline/04_cluster_data.py (see scRNA_pipeline/README.md)
input_obj_dir = Path("results") / "scRNA"
pseudobulk_dir = Path("results") / "edgeR" / "pseudobulk_all_genes"

# %% Load the data
annotated_path = input_obj_dir / "scRNA_data_QC_filtered_w_scVI_latent_annotated.h5ad"
if not annotated_path.exists():
    raise FileNotFoundError(
        f"{annotated_path} not found. Run scRNA_pipeline/04_cluster_data.py first."
    )
if not pseudobulk_dir.exists():
    raise FileNotFoundError(
        f"{pseudobulk_dir} not found. Run edgeR/00_create_pseudobulk.py first."
    )

scVI_adata = ad.read_h5ad(annotated_path)

# %% Match the cell-type naming used by 00_create_pseudobulk.py
scVI_adata.obs["cell_type_manual"] = scVI_adata.obs["cell_type_manual"].str.replace(
    " ", "_", regex=False
)

# Composite groups, mirroring 00_create_pseudobulk.py
combined_mapping = {
    "Combined_Glia": ["Astrocytes", "Cortex_Glia", "Surface_Glia"],
    "Combined_GMCs": [
        "Ganglion_Mother_Cells_(Early)",
        "Ganglion_Mother_Cells_(Late)",
    ],
}

# Condition columns expected by 02_edgeR_DGE.R
condition_levels = ["Control", "ND75-KD"]
observed_conditions = set(scVI_adata.obs["condition"].astype(str).unique())
unexpected = observed_conditions - set(condition_levels)
if unexpected:
    raise ValueError(
        f"Unexpected condition value(s): {sorted(unexpected)}; expected {condition_levels}."
    )


def detection_counts(adata_group):
    """Per gene, the number of cells with raw count > 0 in each condition."""
    detected = adata_group.layers["counts"] > 0  # sparse or dense boolean
    condition = adata_group.obs["condition"].astype(str).to_numpy()

    out = {"gene": adata_group.var_names.to_numpy()}
    for level in condition_levels:
        mask = condition == level
        col = f"n_cells_{level.replace('-', '_')}"
        if mask.any():
            out[col] = np.asarray(detected[mask].sum(axis=0)).ravel().astype(int)
        else:
            out[col] = np.zeros(adata_group.n_vars, dtype=int)
    return pd.DataFrame(out)


# %% One TSV per pseudobulk folder created by 00_create_pseudobulk.py
for folder in sorted(pseudobulk_dir.iterdir()):
    if not folder.is_dir() or folder.name == "edgeR_results":
        continue

    name = folder.name
    if name in combined_mapping:
        members = combined_mapping[name]
        group = scVI_adata[scVI_adata.obs["cell_type_manual"].isin(members)].copy()
    else:
        group = scVI_adata[scVI_adata.obs["cell_type_manual"] == name].copy()

    if group.n_obs == 0:
        logger.warning(f"No cells found for '{name}'; skipping.")
        continue

    df = detection_counts(group)
    out_file = folder / f"{name}_percent_expressing.tsv"
    df.to_csv(out_file, sep="\t", index=False)
    logger.info(f"{name}: {group.n_obs} cells -> {out_file}")

logger.info("Per-condition detection counts written for all pseudobulk groups.")
