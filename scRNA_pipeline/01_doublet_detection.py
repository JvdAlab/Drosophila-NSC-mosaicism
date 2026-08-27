# Script to detect doublets in a single sample via scVI + Solo

# %% Import libraries
import logging
import argparse
from pathlib import Path
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import scvi
import torch

# %% Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("01_doublet_detection.log"),
    ],
)
logger = logging.getLogger(__name__)
logger.info("Starting doublet detection script.")

# %% Set up seeds
# Set random seeds for reproducibility
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)

sc.settings.seed = 42
scvi.settings.seed = 42

# %% Set up paths
obj_dir = Path("results") / "scRNA"
obj_dir.mkdir(parents=True, exist_ok=True)

# %% Check GPU
gpu_available = torch.cuda.is_available()
accelerator = "gpu" if gpu_available else "cpu"
devices = 1 if gpu_available else "auto"

# Log GPU availability
if gpu_available:
    gpu_name = torch.cuda.get_device_name(0)
    logger.info(f"Using GPU: {gpu_name}")
else:
    logger.info("No GPU available, using CPU.")

# %% Parse sampleIDs from command line
parser = argparse.ArgumentParser(description="Doublet detection script for a sample.")
parser.add_argument(
    "--sample_id",
    type=str,
    required=True,
    help="Sample ID for doublet detection.",
)

args = parser.parse_args()
sample_id = args.sample_id

# %% Load adata object
torch.set_num_threads(8)  # Set number of threads for PyTorch

input_path = obj_dir / f"{sample_id}_scRNA_data.h5ad"
if not input_path.exists():
    raise FileNotFoundError(f"{input_path} not found. Run 00_setup_adata_obj.py first.")
scRNA_data_sample = ad.read_h5ad(input_path)

# Use sparse matrix for large datasets
scRNA_data_sample.X = csr_matrix(scRNA_data_sample.X)

scRNA_data_sample.layers["counts"] = scRNA_data_sample.X.copy()

# Basic filtering necessary for scVI/Solo
sc.pp.filter_genes(scRNA_data_sample, min_counts=3)

# Initialize the scVI model
scvi.model.SCVI.setup_anndata(scRNA_data_sample, layer="counts")

logger.info(f"Training scVI model for sample {sample_id}...")

vae_sample = scvi.model.SCVI(
    scRNA_data_sample, n_latent=20, n_layers=1, gene_likelihood="zinb"
)

vae_sample.train(
    max_epochs=200,
    accelerator=accelerator,
    devices=devices,
    train_size=0.9,
    validation_size=0.1,
    early_stopping=True,
    plan_kwargs={"lr": 1e-3},
    check_val_every_n_epoch=10,
)

# Save the trained model
vae_sample.save(obj_dir / "vae_models" / sample_id, overwrite=True)

logger.info(f"Training Solo model for sample {sample_id}...")

# Initialize the Solo model
solo_model = scvi.external.SOLO.from_scvi_model(vae_sample)

# Train the Solo model
solo_model.train(
    max_epochs=50,
    accelerator=accelerator,
    devices=devices,
    train_size=0.9,
    validation_size=0.1,
    early_stopping=True,
    early_stopping_patience=10,
    check_val_every_n_epoch=5,
)

# Save the trained Solo model
solo_model.save(obj_dir / "solo_models" / sample_id, overwrite=True)

logger.info(f"Predicting doublets for sample {sample_id}...")
preds_sample = solo_model.predict(soft=True)

scRNA_data_sample.obs["doublet_prob"] = preds_sample["doublet"]
scRNA_data_sample.obs["singlet_prob"] = preds_sample["singlet"]

# Reset index to turn barcodes into a column for saving with index=False
preds_for_csv = preds_sample.reset_index()

preds_for_csv = preds_for_csv.rename(
    columns={
        preds_for_csv.columns[0]: "barcode",
        "doublet": "doublet_probability",
        "singlet": "singlet_probability",
    }
)

# Add sample_id to the predictions DataFrame
preds_for_csv["sample_id"] = sample_id

preds_for_csv = preds_for_csv[
    ["barcode", "sample_id", "doublet_probability", "singlet_probability"]
]

logger.info("Saving probability predictions to CSV...")
preds_for_csv.to_csv(
    obj_dir / f"{sample_id}_doublet_probs.csv",
    index=False,
)

doublet_thresh = 0.5

logger.info(f"Classifying cells with threshold {doublet_thresh}...")

scRNA_data_sample.obs["doublet_prediction"] = (
    scRNA_data_sample.obs["doublet_prob"] > doublet_thresh
).map({True: "doublet", False: "singlet"})

output_file = obj_dir / f"{sample_id}_scRNA_data_w_doublets.h5ad"
logger.info(f"Saving final AnnData object to {output_file}...")
scRNA_data_sample.write(output_file, compression="gzip")

logger.info(f"Doublet detection completed for sample {sample_id}.")
logger.info(f"Results saved to {output_file}.")
