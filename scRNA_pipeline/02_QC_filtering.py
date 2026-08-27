# Script to perform QC filtering

# %% Import libraries
import logging
from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation

# %% Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("02_QC_filtering.log"),
    ],
)
logger = logging.getLogger(__name__)
logger.info("Starting QC filtering script.")

# %% Set up paths
obj_dir = Path("results") / "scRNA"

# %% Load adata object
adata_list = list(obj_dir.glob("*_w_doublets.h5ad"))

if len(adata_list) == 0:
    logger.error("No adata object found in the specified directory.")
    raise FileNotFoundError(
        f"No *_w_doublets.h5ad files found in {obj_dir}. "
        "Run 01_doublet_detection.py (per sample) first."
    )
else:
    scRNA_w_doublets_dict = {
        Path(file).stem.split("_")[0]: ad.read_h5ad(file) for file in adata_list
    }


# %% Adaptive filtering function
def adaptive_thresholds(metric_data, n_mads=3, log_transform=True, clip_min=0):
    """Calculate adaptive thresholds based on MAD"""

    metric_data = pd.to_numeric(metric_data, errors="coerce")
    metric_data = metric_data.dropna()

    if log_transform:
        metric_data_log = metric_data[metric_data >= 0]
        log_metric = np.log1p(metric_data_log.astype(np.float64))
        median_val = np.median(log_metric)
        mad_val = median_abs_deviation(log_metric, scale="normal", nan_policy="omit")

        lower = np.expm1(median_val - n_mads * mad_val)
        upper = np.expm1(median_val + n_mads * mad_val)

    else:
        median_val = np.median(metric_data)
        mad_val = median_abs_deviation(
            metric_data.astype(np.float64), scale="normal", nan_policy="omit"
        )

        lower = median_val - n_mads * mad_val
        upper = median_val + n_mads * mad_val

    lower = clip_min if np.isnan(lower) or lower < clip_min else round(lower)
    upper = np.inf if np.isnan(upper) or np.isinf(upper) else round(upper)

    # Sanity check: Ensure lower <= upper
    if lower > upper:
        print(
            f"  Warning: Lower threshold ({lower}) > Upper threshold ({upper}). Using broad thresholds."
        )
        return {"lower": clip_min, "upper": np.inf}

    return {"lower": lower, "upper": upper}


# %% Filter out doublets from each sample
# Store filtered singlets
scRNA_singlets_dict = {}

logger.info("Filtering out doublets from each sample...")
# Loop through each sample
for sample_id, adata in scRNA_w_doublets_dict.items():
    adata.obs["sample_id"] = sample_id
    adata = adata[adata.obs["doublet_prediction"] == "singlet"].copy()
    scRNA_singlets_dict[sample_id] = adata

# %% Pool samples by batch
# Concatenate samples by batch
batch_groups = {}
for sample_id, adata in scRNA_singlets_dict.items():
    batch = adata.obs["batch"].unique()[0]
    if batch not in batch_groups:
        batch_groups[batch] = {}
    batch_groups[batch][sample_id] = adata

batch_adata_concat = {}
for batch, samples_dict in batch_groups.items():
    batch_adata_concat[batch] = ad.concat(
        samples_dict, label="sample_id", join="outer", index_unique="-"
    )
    batch_adata_concat[batch].obs_names_make_unique()

# %% Calculate adaptive thresholds per batch
batch_thresholds = {}
for batch, adata in batch_adata_concat.items():
    batch_thresholds[batch] = {
        "n_genes_by_counts": adaptive_thresholds(
            adata.obs["n_genes_by_counts"], n_mads=3, log_transform=True
        ),
        "total_counts": adaptive_thresholds(
            adata.obs["total_counts"], n_mads=3, log_transform=True
        ),
        "pct_counts_mt": adaptive_thresholds(
            adata.obs["pct_counts_mt"], n_mads=3, log_transform=False
        ),
    }
    logger.info(f"Adaptive thresholds for batch {batch}: {batch_thresholds[batch]}")

# Save batch thresholds to file
batch_thresholds_list = []
for batch, metrics in batch_thresholds.items():
    for metric, thresholds in metrics.items():
        # Explicitly extract 'lower' and 'upper' values
        lower = thresholds["lower"] if isinstance(thresholds, dict) else thresholds[0]
        upper = thresholds["upper"] if isinstance(thresholds, dict) else thresholds[1]
        batch_thresholds_list.append(
            {
                "Batch": batch,
                "Metric": metric,
                "Lower": lower,
                "Upper": upper,
            }
        )

batch_thresholds_df = pd.DataFrame(batch_thresholds_list)

# Save batch thresholds to file
batch_thresholds_df.to_csv(obj_dir / "batch_thresholds.csv", index=False, header=True)
logger.info(f"Batch thresholds saved to {obj_dir / 'batch_thresholds.csv'}")

# %% Apply adaptive thresholds to filter cells
# Store QC filtered data
scRNA_data_QC_filtered = {}

for batch, scRNA_data in batch_adata_concat.items():
    scRNA_data = scRNA_data[
        (
            scRNA_data.obs["n_genes_by_counts"]
            >= batch_thresholds[batch]["n_genes_by_counts"]["lower"]
        )
        & (
            scRNA_data.obs["n_genes_by_counts"]
            <= batch_thresholds[batch]["n_genes_by_counts"]["upper"]
        )
        & (
            scRNA_data.obs["total_counts"]
            >= batch_thresholds[batch]["total_counts"]["lower"]
        )
        & (
            scRNA_data.obs["total_counts"]
            <= batch_thresholds[batch]["total_counts"]["upper"]
        )
        & (
            scRNA_data.obs["pct_counts_mt"]
            <= batch_thresholds[batch]["pct_counts_mt"]["upper"]
        )
    ].copy()
    scRNA_data_QC_filtered[batch] = scRNA_data
    logger.info(
        f"Filtered cells for batch {batch}: {scRNA_data.shape[0]} cells remaining."
    )

# %% Save filtered data
# concatenate all batches
scRNA_data_filtered = ad.concat(
    list(scRNA_data_QC_filtered.values()),
    join="outer",
    index_unique=None,
)

scRNA_data_filtered.obs_names_make_unique()
scRNA_data_filtered.write_h5ad(obj_dir / f"scRNA_data_QC_filtered.h5ad")
