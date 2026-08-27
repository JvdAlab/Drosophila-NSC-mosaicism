# Script to train scVI model

# %% Import libraries
import logging
from pathlib import Path
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import scvi
import torch

# %% Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("03_train_scVI_model.log"),
    ],
)
logger = logging.getLogger(__name__)
logger.info("Starting scVI training script.")

# %% Set up seeds
# Set random seeds for reproducibility
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)

sc.settings.seed = 42
scvi.settings.seed = 42

# %% Set up paths
obj_dir = Path("results") / "scRNA"
figures_dir = obj_dir / "figures"
figures_dir.mkdir(parents=True, exist_ok=True)

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

# %% Load data
qc_filtered_path = obj_dir / "scRNA_data_QC_filtered.h5ad"
if not qc_filtered_path.exists():
    raise FileNotFoundError(
        f"{qc_filtered_path} not found. Run 02_QC_filtering.py first."
    )
scRNA_data_QC = ad.read_h5ad(qc_filtered_path)

logger.info("Data loaded successfully.")
logger.info(f"AnnData shape: {scRNA_data_QC.shape}")
logger.info(f"First 5 rows of obs:\n{scRNA_data_QC.obs.head()}")

# %% Setup scVI
# Check if the data is sparse and convert to sparse if necessary
if not isinstance(scRNA_data_QC.X, csr_matrix):
    scRNA_data_QC.X = csr_matrix(scRNA_data_QC.X)
    logger.info("Converted data to sparse matrix format.")
else:
    logger.info("Data is already in sparse matrix format.")

# Register data with scVI with batch and categorical covariates (replicate & pct_counts_mt)
scvi.model.SCVI.setup_anndata(
    scRNA_data_QC,
    batch_key="batch",
    categorical_covariate_keys=["replicate"],
    continuous_covariate_keys=["pct_counts_mt"],
)

# %% Train scVI model
# --- Model Architecture Parameters ---
n_latent_dim = 20  # Dimensionality of the latent space. Default is 10.
n_hidden_nodes = 128  # Number of nodes per hidden layer. Default is 128.
n_layers_scvi = 1  # Number of hidden layers. Default is 1.
gene_likelihood_scvi = (
    "zinb"  # Likelihood function. 'zinb' (Zero-Inflated Negative Binomial)
)

# --- Training Parameters ---
max_train_epochs = 400  # Maximum number of epochs for training. Default is 400.
# Early stopping (if enabled) can stop training sooner.
learning_rate_scvi = 1e-3  # Learning rate. Default is 0.001.
use_early_stopping = True  # Highly recommended to prevent overfitting and save time.
early_stopping_monitor_metric = (
    "elbo_validation"  # Metric to monitor for early stopping.
)
early_stopping_patience_epochs = (
    25  # Number of epochs with no improvement after which training stops.
)
# scVI default is 15 if early_stopping=True and this isn't set.
early_stopping_min_delta = (
    0.00  # Minimum change in monitored quantity to qualify as improvement.
)

# Batch size for training. Default is 128.
# Adjust based on dataset size and GPU memory. Larger might speed up training.
train_batch_size = 256
if scRNA_data_QC.shape[0] < 5000:  # Heuristic for smaller datasets
    train_batch_size = 128
if scRNA_data_QC.shape[0] < 1000:
    train_batch_size = 64

logger.info(f"Initializing scVI model with parameters:")
logger.info(f"  n_latent: {n_latent_dim}")
logger.info(f"  n_hidden: {n_hidden_nodes}")
logger.info(f"  n_layers: {n_layers_scvi}")
logger.info(f"  gene_likelihood: {gene_likelihood_scvi}")

scvi_model = scvi.model.SCVI(
    scRNA_data_QC,
    n_hidden=n_hidden_nodes,
    n_latent=n_latent_dim,
    n_layers=n_layers_scvi,
    gene_likelihood=gene_likelihood_scvi,
    # dropout_rate: Default is 0.1.
    # dispersion: Default is 'gene'. Options: 'gene-batch', 'gene-label', 'gene-cell'.
)

logger.info("scVI model initialized.")

# %% Train scVI model
logger.info(f"Starting scVI model training with the following parameters:")
logger.info(f"  max_epochs: {max_train_epochs}")
logger.info(f"  accelerator: {accelerator}, devices: {devices}")
logger.info(f"  batch_size: {train_batch_size}")
logger.info(f"  learning_rate: {learning_rate_scvi}")
logger.info(f"  early_stopping: {use_early_stopping}")
if use_early_stopping:
    logger.info(f"    early_stopping_monitor: {early_stopping_monitor_metric}")
    logger.info(f"    early_stopping_patience: {early_stopping_patience_epochs}")
    logger.info(f"    early_stopping_min_delta: {early_stopping_min_delta}")

scvi_model.train(
    max_epochs=max_train_epochs,
    accelerator=accelerator,
    devices=devices,
    batch_size=train_batch_size,
    plan_kwargs={"lr": learning_rate_scvi},
    early_stopping=use_early_stopping,
    early_stopping_monitor=early_stopping_monitor_metric,
    early_stopping_patience=early_stopping_patience_epochs,
    early_stopping_min_delta=early_stopping_min_delta,
    # check_val_every_n_epoch: Default is 1 (evaluates validation set every epoch)
)

logger.info("scVI model training complete.")

# %% Plot training history (ELBO)
# ELBO (Evidence Lower BOund) is the objective function scVI maximizes.
try:
    train_elbo = scvi_model.history[
        "elbo_train"
    ].dropna()  # elbo_train might have NaNs if not computed every epoch
    val_elbo = scvi_model.history["elbo_validation"]

    plt.figure(figsize=(8, 6))
    if not train_elbo.empty:
        plt.plot(train_elbo.index, train_elbo.values, label="Train ELBO")
    plt.plot(val_elbo.index, val_elbo.values, label="Validation ELBO")
    plt.title("scVI Model Training ELBO")
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.legend()
    plt.grid(True)
    elbo_plot_path = figures_dir / "scvi_training_elbo_plot.png"
    plt.savefig(elbo_plot_path)
    plt.close()  # Close the plot to free memory
    logger.info(f"ELBO training history plot saved to {elbo_plot_path}")
    logger.info(
        f"Final training ELBO: {train_elbo.iloc[-1] if not train_elbo.empty else 'N/A'}"
    )
    logger.info(f"Final validation ELBO: {val_elbo.iloc[-1]}")
except KeyError as e:
    logger.warning(
        f"Could not retrieve ELBO history for plotting (KeyError: {e}). Ensure validation step ran and history keys are correct."
    )
except Exception as e:
    logger.error(f"An error occurred during ELBO plot generation: {e}")

# %% Save the trained model
model_save_path = obj_dir / "scVI_model_trained"
try:
    scvi_model.save(
        str(model_save_path), overwrite=True, save_anndata=False
    )  # save_anndata=False is often preferred to save separately
    logger.info(f"Trained scVI model saved to {model_save_path}")
except Exception as e:
    logger.error(f"Error saving scVI model: {e}")

# %% Add latent representation to AnnData and save AnnData
adata_latent_save_path = obj_dir / "scRNA_data_QC_filtered_w_scVI_latent.h5ad"
try:
    logger.info("Getting latent representation from the model.")
    scRNA_data_QC.obsm["X_scVI"] = scvi_model.get_latent_representation()
    logger.info(f"Writing AnnData with scVI latent space to {adata_latent_save_path}")
    scRNA_data_QC.write_h5ad(adata_latent_save_path)
    logger.info(f"AnnData with scVI latent space saved successfully.")
except Exception as e:
    logger.error(f"Error getting latent representation or saving AnnData: {e}")

logger.info("Script finished.")
