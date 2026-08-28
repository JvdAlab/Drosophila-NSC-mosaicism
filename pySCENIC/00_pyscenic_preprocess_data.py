# Script to preprocess data step of the pySCENIC pipeline

# %% Import libraries
import argparse
import random
import yaml
import logging
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import pyarrow.feather as feather
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


# %% Pre-process data
def preprocess_data(h5ad_path, ranking_db_path, tf_list_path, mapping_path=None):
    """Rescue gene symbols, deduplicate, intersect with the ranking DB, and drop
    near-zero-variance genes; return the filtered AnnData and the TFs present in it."""

    logger.info(f"Loading data from {h5ad_path}")
    scRNA_adata = ad.read_h5ad(h5ad_path)
    logger.info(f"Initial data shape: {scRNA_adata.shape}")

    # --- 1. Map gene names ---
    if mapping_path and Path(mapping_path).exists():
        logger.info(f"Applying gene mapping from {mapping_path}")
        mapping_df = pd.read_csv(mapping_path)

        # Create dictionary of renames present in the dataset
        rename_dict = {
            row["Original_Symbol"]: row["SCENIC_DB_Symbol"]
            for _, row in mapping_df.iterrows()
            if row["Original_Symbol"] in scRNA_adata.var_names
        }

        # Apply the rename to the AnnData index
        scRNA_adata.var.rename(index=rename_dict, inplace=True)
        logger.info(f"Renamed {len(rename_dict)} genes based on rescue map.")

    # --- 2. Deduplicate gene names ---
    # This specifically fixes the 'Reindexing' error by ensuring unique var_names
    if scRNA_adata.var_names.duplicated().any():
        duplicate_mask = scRNA_adata.var_names.duplicated(keep="first")
        duplicated_names = scRNA_adata.var_names[duplicate_mask].unique().tolist()

        logger.warning(
            f"DEDUPLICATION ALERT: {len(duplicated_names)} gene names "
            f"produced collisions and will be deduplicated."
        )

        # Log specific culprits for your records
        logger.info(f"Top collisions handled: {duplicated_names[:10]}")

        # Keep only the first occurrence of each gene name
        scRNA_adata = scRNA_adata[:, ~scRNA_adata.var_names.duplicated()].copy()
        logger.info(f"Shape after deduplication: {scRNA_adata.shape}")

    # Final physical check to ensure unique index for downstream matrix operations
    scRNA_adata.var_names_make_unique()

    # --- 3. Database intersection: Keep only genes present in database ---
    try:
        db_genes = feather.read_table(ranking_db_path, columns=[]).column_names
        db_genes_set = set([g for g in db_genes if g != "motifs"])

        valid_genes = [
            gene for gene in scRNA_adata.var_names
            if gene in db_genes_set
        ]

        logger.info(f"Overlap with SCENIC DB: {len(valid_genes)} genes.")

        # Reindex the AnnData to keep only valid genes in their original order.
        scRNA_adata = scRNA_adata[:, valid_genes].copy()
        logger.info(f"Filtered data shape: {scRNA_adata.shape}")

    except Exception as e:
        logger.error(f"Error during DB filtering: {e}")
        exit(1)

    # --- 4. Variance filtering: Filter out zero-variance genes ---
    logger.info("Checking for zero-variance genes...")
    if hasattr(scRNA_adata.X, "toarray"):
        # Stable variance calculation for sparse matrices
        mean = np.array(scRNA_adata.X.mean(axis=0)).flatten()
        sqr = np.array(scRNA_adata.X.power(2).mean(axis=0)).flatten()
        var = sqr - mean**2
    else:
        var = np.var(scRNA_adata.X, axis=0)

    # Filter: Keep only genes with variance > 0 to avoid model failure
    scRNA_adata = scRNA_adata[:, var > 1e-10].copy()
    logger.info(f"Final shape after variance filtering: {scRNA_adata.shape}")

    # --- 5. TF validation: Keep only TFs present in the database ---
    with open(tf_list_path, "r") as tf_file:
        tf_list = [line.strip() for line in tf_file.readlines()]

    # Ensure TF list is unique and limited to genes actually present in final data
    valid_tfs = sorted(list(set([tf for tf in tf_list if tf in scRNA_adata.var_names])))
    logger.info(f"Final count of valid unique TFs: {len(valid_tfs)}")

    return scRNA_adata, valid_tfs


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

    # Seed all RNGs from the single config seed (shared_params.seed) for
    # run-to-run reproducibility.
    seed = int(shared["seed"])
    random.seed(seed)
    np.random.seed(seed)

    # Resolve paths
    OUT_DIR = Path(run_cfg["output_dir"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = OUT_DIR / "00_pyscenic_preprocess_data.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    logging.getLogger().addHandler(file_handler)

    logger.info(f"Logging to {log_file}")

    FILES = {
        "H5AD": Path(p["work_dir"]) / p["h5ad"],
        "RANKING_DB": Path(p["res_dir"]) / p["ranking_db"],
        "TF_LIST": Path(p["res_dir"]) / p["tf_list"],
        "GENE_MAP": Path(p["res_dir"]) / p["gene_map"],
        "PROCESSED_H5AD": OUT_DIR / "preprocessed.h5ad",
    }

    scRNA_adata, valid_tfs = preprocess_data(
        FILES["H5AD"],
        FILES["RANKING_DB"],
        FILES["TF_LIST"],
        mapping_path=FILES["GENE_MAP"],
    )

    scRNA_adata.write_h5ad(FILES["PROCESSED_H5AD"])
    pd.Series(valid_tfs).to_csv(OUT_DIR / "valid_tfs.txt", index=False, header=False)
    logger.info("00_pyscenic_preprocess_data completed successfully.")
