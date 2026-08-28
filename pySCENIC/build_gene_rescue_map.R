#!/usr/bin/env Rscript

# Generate the precomputed SCENIC gene rescue map.

library(org.Dm.eg.db)
library(arrow)
library(tidyverse)
library(here)

logger <- function(message, ...) {
    cat(paste0(Sys.time(), " - ", sprintf(message, ...), "\n"))
}

# --- Repository-relative paths ---------------------------------------------
input_rds <- here::here(
  "results", "scRNA",
  "scRNA_data_QC_filtered_w_scVI_latent_annotated_seurat.rds"
)
out_dir <- here::here("pySCENIC", "resources", "pyscenic")
db_path <- file.path(out_dir, "dm6_v10_clust.genes_vs_motifs.rankings.feather")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_rds)) {
    stop(sprintf(
        "%s not found. Run scRNA_pipeline/06_convert_adata_to_seurat.R first.",
        input_rds
    ))
}
if (!file.exists(db_path)) {
    stop(sprintf(
        "%s not found. See pySCENIC/resources/pyscenic/README.md for the SCENIC ranking database.",
        db_path
    ))
}

# 1. Load all dataset genes.
logger("Loading annotated Seurat object from %s", input_rds)
scRNA_data <- readRDS(input_rds)
all_dataset_genes <- rownames(scRNA_data)

# 2. Load SCENIC headers and normalize them for fuzzy matching.
db_headers <- names(arrow::read_feather(db_path, as_data_frame = FALSE))
db_headers <- db_headers[db_headers != "motifs"]
db_norm <- setNames(gsub("[[:punct:]]| ", "", tolower(db_headers)), db_headers)

# 3. Load the global FlyBase symbol/alias map.
logger("Pulling global FlyBase mapping...")
full_fb_map <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = keys(org.Dm.eg.db, keytype = "SYMBOL"),
    columns = c("SYMBOL", "ALIAS"),
    keytype = "SYMBOL"
)

# 4. Audit missing genes and rescue only collision-free matches.
missing_genes <- all_dataset_genes[!(all_dataset_genes %in% db_headers)]
logger("Auditing %d genes missing from SCENIC DB...", length(missing_genes))

rescue_list <- list()

for (m_gene in missing_genes) {
    target_found <- NULL

    # A. Fuzzy case/punctuation match.
    m_norm <- gsub("[[:punct:]]| ", "", tolower(m_gene))
    match_idx <- which(db_norm == m_norm)
    if (length(match_idx) > 0) {
        target_found <- names(db_norm)[match_idx[1]]
    }

    # B. FlyBase alias match if fuzzy matching failed.
    if (is.null(target_found)) {
        associated_names <- full_fb_map %>%
            filter(SYMBOL == m_gene) %>%
            select(SYMBOL, ALIAS) %>%
            pivot_longer(everything()) %>%
            pull(value) %>%
            unique()

        for (alt_name in associated_names) {
            if (!is.na(alt_name) && alt_name %in% db_headers) {
                target_found <- alt_name
                break
            }
        }
    }

    # Only rescue when the target symbol is not already present in the data.
    if (!is.null(target_found)) {
        if (!(target_found %in% all_dataset_genes)) {
            rescue_list[[m_gene]] <- target_found
        } else {
            logger(
                "Collision: skipping rescue of %s -> %s (target already in dataset)",
                m_gene, target_found
            )
        }
    }
}

# 5. Export the cleaned master map.
global_rescue_df <- data.frame(
    Original_Symbol = names(rescue_list),
    SCENIC_DB_Symbol = unlist(rescue_list)
)

global_rescue_df <- global_rescue_df %>%
    filter(!str_detect(Original_Symbol, "Psi:|rRNA|lncRNA|asRNA|sisRNA"))

write_csv(global_rescue_df, file.path(out_dir, "global_pyscenic_rescue_map.csv"))
logger(
    "Rescued %d genes without nomenclature collisions.",
    nrow(global_rescue_df)
)
