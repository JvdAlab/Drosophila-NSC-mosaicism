# Script to perform the GRN construction step of the pySCENIC pipeline

# %% Import libraries
import argparse
import random
import yaml
import logging
import torch
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import regdiffusion as rd
from pathlib import Path


# %% Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


# %% Set up config
def load_config(config_path):
    """Load pipeline configuration from a YAML file."""
    logger.info(f"Loading configuration from: {config_path}")
    if not Path(config_path).exists():
        raise FileNotFoundError(f"Config file not found at {config_path}")
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# %% GRN construction
def run_grn_inference(scRNA_adata, valid_tfs, top_percentile=85, num_workers=16):
    """Train RegDiffusion on the expression matrix and return the valid-TF edge
    list (columns TF, target, importance) above the top-percentile importance cut,
    sorted by importance."""

    logger.info("Starting GRN inference with RegDiffusion.")

    if "log1p_norm" in scRNA_adata.layers:
        logger.info("Using 'log1p_norm' layer for GRN inference.")
        input_matrix = scRNA_adata.layers["log1p_norm"]
    else:
        logger.info("Using .X (checking if normalization is needed) for GRN inference.")
        if np.max(scRNA_adata.X) > 20:
            sc.pp.normalize_total(scRNA_adata, target_sum=1e4)
            sc.pp.log1p(scRNA_adata)
        input_matrix = scRNA_adata.X

    # Ensure dense matrix (RegDiffusion requirement). Densify the matrix selected
    # above (the log1p_norm layer when present) rather than re-reading raw .X.
    input_matrix = (
        input_matrix.toarray() if hasattr(input_matrix, "toarray") else np.asarray(input_matrix)
    ).astype(np.float32)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info(f"Training on device: {device}")

    rd_trainer = rd.RegDiffusionTrainer(
        input_matrix, device=device, batch_size=8 if device == "cuda" else None
    )

    logger.info("Training RegDiffusion model...")
    rd_trainer.train()

    logger.info("Extracting edges...")
    grn = rd_trainer.get_grn(scRNA_adata.var_names, top_gene_percentile=top_percentile)

    adjacencies = grn.extract_edgelist(k=-1, workers=num_workers)
    adjacencies.columns = ["TF", "target", "importance"]

    # Force float32 to prevent 'object' dtype conversion during pickling
    adjacencies["importance"] = adjacencies["importance"].astype(np.float32)

    # Force string types to avoid mixed-type object columns
    adjacencies["TF"] = adjacencies["TF"].astype(str)
    adjacencies["target"] = adjacencies["target"].astype(str)

    # Filter for valid TFs
    adjacencies = adjacencies[adjacencies["TF"].isin(valid_tfs)]
    logger.info(f"Number of edges after TF filtering: {adjacencies.shape[0]}")

    # Remove self-loops
    adjacencies = adjacencies[adjacencies["TF"] != adjacencies["target"]]

    # Sort by importance
    adjacencies = adjacencies.sort_values(by="importance", ascending=False)

    logger.info(
        f"Final number of edges after self-loop removal: {adjacencies.shape[0]}"
    )

    return adjacencies


# %% Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="pyscenic_config.yaml")
    args = parser.parse_args()

    # Load config
    cfg = load_config(args.config)

    p = cfg["paths"]
    shared = cfg["shared_params"]
    run_cfg = cfg["run"]

    # Seed all RNGs from the single config seed (shared_params.seed) before
    # RegDiffusion runs. RegDiffusion draws from the torch / numpy global RNGs;
    # residual CUDA-kernel nondeterminism may still remain across GPUs.
    seed = int(shared["seed"])
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    # Resolve paths
    OUT_DIR = Path(run_cfg["output_dir"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = OUT_DIR / "01_pyscenic_grn_construction.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    logging.getLogger().addHandler(file_handler)

    logger.info(f"Logging to {log_file}")

    FILES = {
        "PROCESSED_H5AD": OUT_DIR / "preprocessed.h5ad",
        "VALID_TFS": OUT_DIR / "valid_tfs.txt",
        "OUTPUT_ADJ": OUT_DIR / "adjacencies.csv",
    }

    logger.info(f"Loading preprocessed data from {FILES['PROCESSED_H5AD']}")
    scRNA_adata = ad.read_h5ad(FILES["PROCESSED_H5AD"])

    valid_tfs = pd.read_csv(FILES["VALID_TFS"], header=None)[0].tolist()

    adjacencies = run_grn_inference(
        scRNA_adata,
        valid_tfs,
        top_percentile=run_cfg["top_gene_percentile"],
        num_workers=shared["num_workers"],
    )

    adjacencies.to_csv(FILES["OUTPUT_ADJ"], index=False)
    adjacencies.to_pickle(FILES["OUTPUT_ADJ"].with_suffix(".pkl"))
    logger.info(f"Adjacency list saved to {FILES['OUTPUT_ADJ']}")

    logger.info("01_pyscenic_grn_construction completed successfully.")
