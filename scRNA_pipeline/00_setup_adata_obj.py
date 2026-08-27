# Script to build per-sample and combined raw AnnData objects from CellRanger output

# %% Import libraries
import logging
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd

# %% Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("00_setup_adata_obj.log"),
    ],
)

logger = logging.getLogger(__name__)
logger.info("Starting setup adata object script.")

# %% Set up paths
raw_input_dir = Path("scRNA_pipeline") / "data" / "raw"
output_dir = Path("results") / "scRNA"
output_dir.mkdir(parents=True, exist_ok=True)

# Path to the samplesheet and CellRanger outputs
samplesheet_file = raw_input_dir / "samplesheet.csv"
cellranger_dir = raw_input_dir / "CellRanger_Outputs"

# %% Load samplesheet
scRNA_data_dict = {}
sample_ids = []

samplesheet = pd.read_csv(samplesheet_file)
required_samplesheet_columns = {
    "SampleID",
    "Condition",
    "Batch",
    "Replicate",
    "CellRangerSample",
}
missing_columns = required_samplesheet_columns.difference(samplesheet.columns)
if missing_columns:
    raise ValueError(
        f"Samplesheet is missing required columns: {sorted(missing_columns)}"
    )
if samplesheet["SampleID"].duplicated().any():
    raise ValueError("Samplesheet contains duplicate SampleID values.")

samplesheet_metadata = (
    samplesheet.set_index("SampleID")
    .rename(
        columns={
            "Condition": "condition",
            "Batch": "batch",
            "Replicate": "replicate",
            "CellRangerSample": "cellranger_sample",
        }
    )
    .to_dict(orient="index")
)

for index, row in samplesheet.iterrows():
    sample_id = row["SampleID"]
    logging.info(f"Processing sample: {sample_id}")
    sample_path = (
        Path(cellranger_dir) / sample_id / "outs" / "filtered_feature_bc_matrix"
    )

    if sample_path.is_dir():
        scRNA_data_sample = sc.read_10x_mtx(
            sample_path, var_names="gene_symbols", cache=True
        )

        scRNA_data_sample.obs.index = scRNA_data_sample.obs.index.astype(str)
        scRNA_data_sample.var.index = scRNA_data_sample.var.index.astype(str)

        scRNA_data_sample.var_names_make_unique()

        scRNA_data_sample.obs["sample_id"] = sample_id

        sample_metadata = samplesheet_metadata.get(sample_id)
        if sample_metadata is None:
            raise ValueError(f"No metadata found for SampleID {sample_id!r}.")

        # Add samplesheet metadata to the AnnData object.
        for metadata_name, metadata_value in sample_metadata.items():
            scRNA_data_sample.obs[metadata_name] = metadata_value

        # Make them categorical
        scRNA_data_sample.obs["condition"] = pd.Categorical(
            scRNA_data_sample.obs["condition"]
        )
        scRNA_data_sample.obs["batch"] = pd.Categorical(scRNA_data_sample.obs["batch"])
        scRNA_data_sample.obs["replicate"] = pd.Categorical(
            scRNA_data_sample.obs["replicate"]
        )

        # Identify mitochondrial genes
        scRNA_data_sample.var["mt"] = scRNA_data_sample.var_names.str.startswith("mt:")

        # Calculate the percentage of mitochondrial genes
        sc.pp.calculate_qc_metrics(
            scRNA_data_sample,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        logging.info("\nAnnData object created successfully.")
        logging.info(scRNA_data_sample.obs.head())

        logging.info("\nValue Counts:")
        logging.info(
            scRNA_data_sample.obs.value_counts(["batch", "condition", "replicate"])
        )

        # Output path for the AnnData object
        scRNA_obj_file = Path(output_dir) / f"{sample_id}_scRNA_data.h5ad"

        # Save the AnnData object
        scRNA_data_sample.write(scRNA_obj_file)

        logging.info(f"AnnData object saved to {scRNA_obj_file}")

        # Append Anndata object to dict
        scRNA_data_dict[sample_id] = scRNA_data_sample

    else:
        logging.warning(f"Sample path does not exist: {sample_path}")
        continue

# %% Concatenate all AnnData objects into one
scRNA_data = ad.concat(
    scRNA_data_dict, label="sample_id", join="outer", index_unique="-"
)

scRNA_data.obs_names_make_unique()

scRNA_data_file = Path(output_dir) / "scRNA_raw_data.h5ad"
scRNA_data.write(scRNA_data_file)
