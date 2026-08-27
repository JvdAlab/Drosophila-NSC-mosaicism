# Script to export the annotated AnnData object to Seurat-importable flat files

# %% Export h5ad to Seurat-compatible files
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from pathlib import Path

# %% Load data
obj_dir = Path("results") / "scRNA"
annotated_path = obj_dir / "scRNA_data_QC_filtered_w_scVI_latent_annotated.h5ad"
if not annotated_path.exists():
    raise FileNotFoundError(
        f"{annotated_path} not found. Run 04_cluster_data.py first."
    )
scVI_adata = sc.read_h5ad(annotated_path)

# Create output directory
out_dir = obj_dir / "seurat_conversion"
out_dir.mkdir(exist_ok=True)

# %% 1. Export counts matrix (genes x cells for R)
print("Exporting counts matrix...")
counts = scVI_adata.layers["counts"]
if not isinstance(counts, csr_matrix):
    counts = csr_matrix(counts)
mmwrite(out_dir / "counts.mtx", counts.T)

# %% 2. Export normalized matrix
print("Exporting normalized matrix...")
norm_data = scVI_adata.layers["log1p_norm"]
if not isinstance(norm_data, csr_matrix):
    norm_data = csr_matrix(norm_data)
mmwrite(out_dir / "normalized.mtx", norm_data.T)

# %% 3. Export gene names and barcodes
print("Exporting genes and barcodes...")
pd.DataFrame(scVI_adata.var_names).to_csv(
    out_dir / "genes.csv", header=False, index=False
)
pd.DataFrame(scVI_adata.obs_names).to_csv(
    out_dir / "barcodes.csv", header=False, index=False
)

# %% 4. Export cell metadata
print("Exporting metadata...")
# Select columns to export (exclude internal scvi columns)
metadata_cols = [
    "sample_id",
    "condition",
    "batch",
    "replicate",
    "n_genes_by_counts",
    "total_counts",
    "total_counts_mt",
    "pct_counts_mt",
    "doublet_prob",
    "singlet_prob",
    "doublet_prediction",
    "leiden_r0p4",
    "cell_type_manual",
]
scVI_adata.obs[metadata_cols].to_csv(out_dir / "metadata.csv")

# %% 5. Export UMAP coordinates
print("Exporting UMAP...")
umap_df = pd.DataFrame(
    scVI_adata.obsm["X_umap"], index=scVI_adata.obs_names, columns=["UMAP_1", "UMAP_2"]
)
umap_df.to_csv(out_dir / "umap_coords.csv")

# %% 6. Export tSNE coordinates
print("Exporting tSNE...")
tsne_df = pd.DataFrame(
    scVI_adata.obsm["X_tsne"], index=scVI_adata.obs_names, columns=["tSNE_1", "tSNE_2"]
)
tsne_df.to_csv(out_dir / "tsne_coords.csv")

# %% 7. Export scVI latent space
print("Exporting scVI latent space...")
scvi_df = pd.DataFrame(
    scVI_adata.obsm["X_scVI"],
    index=scVI_adata.obs_names,
    columns=[f"scVI_{i+1}" for i in range(scVI_adata.obsm["X_scVI"].shape[1])],
)
scvi_df.to_csv(out_dir / "scvi_latent.csv")

# %% 8. Export neighborhood graph (optional, for clustering reproducibility)
print("Exporting neighborhood graph...")
mmwrite(out_dir / "connectivities.mtx", scVI_adata.obsp["connectivities"])
mmwrite(out_dir / "distances.mtx", scVI_adata.obsp["distances"])

# %% Summary
print(f"\n=== Export Complete ===")
print(f"Output directory: {out_dir}")
print(f"Cells: {scVI_adata.n_obs}")
print(f"Genes: {scVI_adata.n_vars}")
print(f"\nFiles created:")
for f in sorted(out_dir.glob("*")):
    print(f"  {f.name}")
