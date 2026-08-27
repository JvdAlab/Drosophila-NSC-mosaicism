# Script to import the exported flat files into a Seurat v5 object

# %% Import h5ad-exported data into Seurat
library(Seurat)
library(Matrix)
library(ggplot2)
library(here)

# Set paths
data_dir <- here::here("results", "scRNA")
out_dir <- here(data_dir, "seurat_conversion")

if (!file.exists(here(out_dir, "counts.mtx"))) {
    stop(
        "seurat_conversion/counts.mtx not found. Run 05_convert_adata_to_seurat.py first."
    )
}

# %% 1. Load counts matrix
message("Loading counts matrix...")
counts <- readMM(here(out_dir, "counts.mtx"))
counts <- as(counts, "CsparseMatrix")

# %% 2. Load gene names and barcodes
genes <- read.csv(here(out_dir, "genes.csv"), header = FALSE)$V1
barcodes <- read.csv(here(out_dir, "barcodes.csv"), header = FALSE)$V1

# Set dimension names
rownames(counts) <- genes
colnames(counts) <- barcodes

# %% 3. Create Seurat object
message("Creating Seurat object...")
seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = "OxPhos_scRNA",
    min.cells = 0,
    min.features = 0
)

# %% 4. Add normalized data (Seurat v5 syntax)
message("Adding normalized data...")
norm_data <- readMM(here(out_dir, "normalized.mtx"))
norm_data <- as(norm_data, "CsparseMatrix")
rownames(norm_data) <- genes
colnames(norm_data) <- barcodes

seurat_obj <- SetAssayData(
    seurat_obj,
    layer = "data",
    new.data = norm_data
)

# %% 4b. Finalize Seurat v5 assay structure
message("Joining layers (Seurat v5)...")
seurat_obj <- JoinLayers(seurat_obj)

# %% 5. Add metadata
message("Adding metadata...")
metadata <- read.csv(here(out_dir, "metadata.csv"), row.names = 1)
seurat_obj <- AddMetaData(seurat_obj, metadata)

# %% 6. Add UMAP embedding
message("Adding UMAP...")
umap_coords <- read.csv(here(out_dir, "umap_coords.csv"), row.names = 1)
umap_coords <- as.matrix(umap_coords)
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

seurat_obj[["umap"]] <- CreateDimReducObject(
    embeddings = umap_coords,
    key = "UMAP_",
    assay = "RNA"
)

# %% 7. Add tSNE embedding
message("Adding tSNE...")
tsne_coords <- read.csv(here(out_dir, "tsne_coords.csv"), row.names = 1)
tsne_coords <- as.matrix(tsne_coords)
colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")

seurat_obj[["tsne"]] <- CreateDimReducObject(
    embeddings = tsne_coords,
    key = "tSNE_",
    assay = "RNA"
)

# %% 8. Add scVI latent space
message("Adding scVI latent space...")
scvi_latent <- read.csv(here(out_dir, "scvi_latent.csv"), row.names = 1)
scvi_latent <- as.matrix(scvi_latent)

# Fix column names for Seurat (must be key + integer)
colnames(scvi_latent) <- paste0("scVI_", 1:ncol(scvi_latent))

seurat_obj[["scvi"]] <- CreateDimReducObject(
    embeddings = scvi_latent,
    key = "scVI_",
    assay = "RNA"
)

# %% 9. Set default identity and factor levels
message("Setting identities...")
Idents(seurat_obj) <- seurat_obj$leiden_r0p4

# Convert cluster columns to factors with proper ordering
cluster_cols <- grep("^leiden", colnames(seurat_obj@meta.data), value = TRUE)
for (col in cluster_cols) {
    seurat_obj@meta.data[[col]] <- factor(
        seurat_obj@meta.data[[col]],
        levels = as.character(sort(as.numeric(unique(seurat_obj@meta.data[[col]]))))
    )
}

# %% 10. Verify the object
message("\n=== Seurat Object Summary ===")
print(seurat_obj)

message("\nReductions available:")
print(names(seurat_obj@reductions))

message("\nMetadata columns:")
print(colnames(seurat_obj@meta.data))

# %% 11. Quick verification plots
message("\nGenerating verification plots...")

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition") +
    ggtitle("UMAP by Condition")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "leiden_r0p4") +
    ggtitle("UMAP by Cluster")

print(p1 + p2)

# %% 12. Save the Seurat object
message("Saving Seurat object...")
annotated_seurat_rds <- "scRNA_data_QC_filtered_w_scVI_latent_annotated_seurat.rds"
saveRDS(seurat_obj, here(data_dir, annotated_seurat_rds))

message("\n=== Conversion Complete ===")
message(paste("Cells:", ncol(seurat_obj)))
message(paste("Genes:", nrow(seurat_obj)))
message(paste("Output:", here(data_dir, annotated_seurat_rds)))
