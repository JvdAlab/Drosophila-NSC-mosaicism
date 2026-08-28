# Script to perform the AUC scoring step of the pySCENIC pipeline

# %% Import libraries
import argparse
import random
import yaml
import logging
import anndata as ad
import pandas as pd
import numpy as np
from pyscenic.prune import df2regulons
from pyscenic.aucell import aucell
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


# %% AUC scoring
def run_auc_scoring(scRNA_adata, regulons, num_workers=16, seed=42):
    """Score per-cell regulon activity with AUCell; return a cells x regulons AUC matrix."""

    logger.info("Starting AUC scoring...")

    # Core AUCell computation
    auc_mtx = aucell(scRNA_adata.to_df(), regulons, num_workers=num_workers, seed=seed)

    logger.info(f"AUC scoring completed. Matrix shape: {auc_mtx.shape}")
    return auc_mtx


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

    # Seed all RNGs from the single config seed (shared_params.seed); the same
    # value is passed to AUCell below.
    seed = int(shared["seed"])
    random.seed(seed)
    np.random.seed(seed)

    # Resolve paths
    OUT_DIR = Path(run_cfg["output_dir"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = OUT_DIR / "03_pyscenic_auc_scoring.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    logging.getLogger().addHandler(file_handler)

    logger.info(f"Logging to {log_file}")

    FILES = {
        "PROCESSED_H5AD": OUT_DIR / "preprocessed.h5ad",
        "MOTIF_PICKLE": OUT_DIR / "ctx_pruning_results.pkl",
        "AUC_RESULTS": OUT_DIR / "auc_matrix.csv",
    }

    logger.info(f"Loading inputs from {OUT_DIR}")
    scRNA_adata = ad.read_h5ad(FILES["PROCESSED_H5AD"])

    # Load motifs from the ctx pruning step
    df_motifs = pd.read_pickle(FILES["MOTIF_PICKLE"])
    regulons = df2regulons(df_motifs)

    logger.info(f"Number of regulons loaded for scoring: {len(regulons)}")

    auc_mtx = run_auc_scoring(
        scRNA_adata,
        regulons,
        num_workers=cfg["shared_params"]["num_workers"],
        seed=seed,
    )

    auc_mtx.to_csv(FILES["AUC_RESULTS"])
    logger.info(
        f"03_pyscenic_auc_scoring complete. Final matrix saved to {FILES['AUC_RESULTS']}"
    )
