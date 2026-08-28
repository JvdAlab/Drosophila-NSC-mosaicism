# Script to perform the CTX pruning step of the pySCENIC pipeline

# %% Import libraries
import argparse
import random
import yaml
import logging
import anndata as ad
import pandas as pd
import numpy as np
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from ctxcore.rnkdb import FeatherRankingDatabase
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


# %% CTX Pruning
def run_ctx_pruning(
    adjacencies,
    scRNA_adata,
    ranking_db_path,
    motif_annotations_path,
    config,
    num_workers=16,
    out_motifs_path=None,
    out_regulons_path=None,
):
    """Build co-expression modules from the adjacencies and prune them by cisTarget
    motif enrichment; return the regulons, the motif-enrichment table, and a
    human-readable regulon summary."""

    logger.info("Starting CTX pruning...")

    # Force conversion to float32 and ensure it is a dense numpy-backed DataFrame
    # This specifically addresses the Numba 'No matching signature' error
    ex_matrix = scRNA_adata.to_df().astype(np.float32)
    # Remove any columns that are ALL zeros or have NaN
    ex_matrix = ex_matrix.loc[:, (ex_matrix != 0).any(axis=0)].dropna(axis=1)

    # Force C-contiguous memory layout
    # This specifically addresses the Numba 'No matching signature' error
    ex_matrix_dense = pd.DataFrame(
        data=np.ascontiguousarray(ex_matrix.values),
        index=ex_matrix.index,
        columns=ex_matrix.columns,
    )

    logger.info(f"Cleaned matrix for Numba. Shape: {ex_matrix_dense.shape}")

    # 2. Intersection safety and Type casting
    valid_genes = set(ex_matrix_dense.columns)
    adjacencies = adjacencies[
        (adjacencies["TF"].isin(valid_genes))
        & (adjacencies["target"].isin(valid_genes))
    ].copy()

    # Force numeric types to fix the 'nlargest' object error
    adjacencies["importance"] = pd.to_numeric(
        adjacencies["importance"], errors="coerce"
    )
    adjacencies = adjacencies.dropna(subset=["importance"])
    adjacencies["importance"] = adjacencies["importance"].astype(np.float32)

    # 3. Clean strings
    adjacencies["TF"] = adjacencies["TF"].astype(str)
    adjacencies["target"] = adjacencies["target"].astype(str)

    logger.info(f"Final adjacency count for module generation: {len(adjacencies)}")

    modules = modules_from_adjacencies(
        adjacencies, ex_matrix_dense, rho_mask_dropouts=False
    )

    logger.info(f"Co-expression modules generated: {len(modules)}")

    # Fix database initialization naming
    dbs = [FeatherRankingDatabase(ranking_db_path, name="dm6_v10")]

    df_motifs = prune2df(
        dbs,
        modules,
        motif_annotations_path,
        rank_threshold=config["rank_threshold"],
        nes_threshold=config["nes_threshold"],
        motif_similarity_fdr=shared["motif_similarity_fdr"],
        auc_threshold=shared["auc_threshold"],
        num_workers=num_workers,
        weighted_recovery=False,
        filter_for_annotation=False,
    )

    regulons = df2regulons(df_motifs)

    # Build human-readable summary DataFrame
    regulon_data = []
    for r in regulons:
        context = list(r.context)
        if context[0].endswith(".png"):
            context[0], context[1] = context[1], context[0]

        logo_url = f"https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/logos/{context[1]}"
        regulon_data.append(
            [
                r.name,
                r.transcription_factor,
                context[0],
                len(r.gene2weight),
                logo_url,
                r.score,
                ",".join(r.gene2weight),
            ]
        )

    regulon_df = pd.DataFrame(
        regulon_data,
        columns=[
            "Regulon",
            "TF",
            "TFTargetGenesCorrelation",
            "NbMarkers",
            "Motif",
            "NES",
            "Markers",
        ],
    )

    return regulons, df_motifs, regulon_df.sort_values(by="NbMarkers", ascending=False)


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
    log_file = OUT_DIR / "02_pyscenic_ctx_pruning.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    logging.getLogger().addHandler(file_handler)

    logger.info(f"Logging to {log_file}")

    FILES = {
        "PROCESSED_H5AD": OUT_DIR / "preprocessed.h5ad",
        "OUTPUT_ADJ": OUT_DIR / "adjacencies.csv",
        "MOTIF_RESULTS": OUT_DIR / "ctx_pruning_results.csv",
        "REGULONS_RESULTS": OUT_DIR / "final_regulons.csv",
    }

    logger.info(f"Loading inputs from {OUT_DIR}")
    scRNA_adata = ad.read_h5ad(FILES["PROCESSED_H5AD"])
    adjacencies = pd.read_pickle(FILES["OUTPUT_ADJ"].with_suffix(".pkl"))

    regulons, df_motifs, regulon_df = run_ctx_pruning(
        adjacencies,
        scRNA_adata,
        Path(p["res_dir"]) / p["ranking_db"],
        Path(p["res_dir"]) / p["motif_annotations"],
        config=run_cfg,
        num_workers=cfg["shared_params"]["num_workers"],
    )

    df_motifs.to_csv(FILES["MOTIF_RESULTS"], index=True)
    df_motifs.to_pickle(FILES["MOTIF_RESULTS"].with_suffix(".pkl"))

    regulon_df.to_csv(FILES["REGULONS_RESULTS"], index=False)

    logger.info(f"02_pyscenic_ctx_pruning complete. Results saved in {OUT_DIR}")
