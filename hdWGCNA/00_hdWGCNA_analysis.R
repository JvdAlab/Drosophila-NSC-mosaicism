# Script to perform hdWGCNA analysis

# Import libraries
library(here)
library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(tidyverse)
library(patchwork)

# Metacell aggregation (MetacellsByGroups) samples cells, set seed for
# a reproducible network.
set.seed(42)

# data_dir: produced by scRNA_pipeline/06_convert_adata_to_seurat.R (see scRNA_pipeline/README.md)
data_dir <- here::here("results", "scRNA")
wgcna_dir <- here::here("results", "hdWGCNA")
dir.create(wgcna_dir, showWarnings = FALSE, recursive = TRUE)

# Console + file logging
log_file <- "00_hdWGCNA_analysis.log"

logger <- function(message, ...) {
    line <- paste0(Sys.time(), " - ", sprintf(message, ...), "\n")
    cat(line)
    cat(line, file = log_file, append = TRUE)
}

seurat_rds_path <- here(data_dir, "scRNA_data_QC_filtered_w_scVI_latent_annotated_seurat.rds")
if (!file.exists(seurat_rds_path)) {
    stop(sprintf(
        "%s not found. Run scRNA_pipeline/06_convert_adata_to_seurat.R first.",
        seurat_rds_path
    ))
}
scRNA_data <- readRDS(seurat_rds_path)

# --- FEATURE CLEANING FOR hdWGCNA ---
logger("Starting feature hygiene filtering...")

# 1. Identify patterns for non-essential/noisy genes
# - Psi:CR or rRNA:CR (pseudogenes/rRNAs)
# - sisRNA, asRNA, hpRNA (often non-motif-regulated ncRNAs)
# - FBti (transposable elements)
patterns_to_remove <- c("^Psi:", ":CR", "^rRNA", "^sisRNA", "^asRNA", "^hpRNA", "^FBti")

# Combine patterns into a single regex string
regex_pattern <- paste(patterns_to_remove, collapse = "|")

# 2. Get list of all genes
all_genes <- rownames(scRNA_data)

# 3. Identify genes to keep
genes_to_keep <- all_genes[!grepl(regex_pattern, all_genes)]

# Ensure genes of interest are not filtered out
essential_hubs <- c("Ldh", "sima", "crc", "Xrp1", "ImpL3")
genes_to_keep <- unique(c(genes_to_keep, intersect(all_genes, essential_hubs)))

logger(
    "Filtered out %d genes. Remaining: %d genes.",
    length(all_genes) - length(genes_to_keep),
    length(genes_to_keep)
)

# 4. Subset the Seurat object
scRNA_data <- subset(scRNA_data, features = genes_to_keep)

# Initialize hdWGCNA object
scRNA_obj <- SetupForWGCNA(
    scRNA_data,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "Neuroblast_Project"
)

# Construct Metacells
# Group by both cell type and condition
scRNA_obj <- MetacellsByGroups(
    seurat_obj = scRNA_obj,
    group.by = c("cell_type_manual", "condition"),
    k = 15,
    max_shared = 5,
    ident.group = "cell_type_manual",
    reduction = "scvi"
)

#  Normalize the Metacells
scRNA_obj <- NormalizeMetacells(scRNA_obj)

# Focus on Neuroblasts
scRNA_obj <- SetDatExpr(
    scRNA_obj,
    group_name = "Neuroblasts",
    group.by = "cell_type_manual",
    use_metacells = TRUE,
    assay = "RNA"
)

# Pick soft-thresholding power
scRNA_obj <- TestSoftPowers(scRNA_obj, networkType = "signed")

# Plot the power table
# The chosen power (9, below) is read off this plot
soft_power_plot <- wrap_plots(PlotSoftPowers(seurat_obj = scRNA_obj), ncol = 2)
soft_power_plot_path <- here(wgcna_dir, "soft_power_diagnostics.pdf")
ggsave(soft_power_plot_path, soft_power_plot, width = 12, height = 8)
logger("Saved soft-power diagnostics: %s", soft_power_plot_path)

#  Construct the Co-expression Network
scRNA_obj <- ConstructNetwork(
    scRNA_obj,
    soft_power = 9,
    setDatExpr = FALSE,
    tom_name = "NB_Network",
    tom_outdir = wgcna_dir,
    overwrite_tom = TRUE
)

wgcna_name <- scRNA_obj@misc$active_wgcna
scRNA_obj@misc[[wgcna_name]]$wgcna_net$TOMFiles <-
    file.path(wgcna_dir, "NB_Network_TOM.rda")

dendro_plot_path <- here(wgcna_dir, "NB_dendrogram.pdf")
pdf(dendro_plot_path, width = 12, height = 8)
PlotDendrogram(scRNA_obj, main = "Neuroblast Co-expression Dendrogram")
dev.off()
logger("Saved co-expression dendrogram: %s", dendro_plot_path)

# --- Rebuild the RNA assay from counts, then re-derive data + scale.data ---
# After metacell/network construction the working RNA assay no longer has a
# 'counts'/'data' pair with consistent dimensions, which ModuleEigengenes
# rejects ("invalid class"). Rebuilding the assay from the counts matrix and
# re-running NormalizeData / ScaleData restores a consistent assay. This step
# is required with the pinned Seurat/hdWGCNA versions.
counts_mat <- GetAssayData(scRNA_obj, slot = "counts", assay = "RNA")
scRNA_obj[["RNA"]] <- CreateAssayObject(counts = counts_mat)
scRNA_obj <- NormalizeData(scRNA_obj)

# Scale only the network genes that exist in the rebuilt assay.
network_genes <- GetModules(scRNA_obj)$gene_name
valid_genes <- intersect(network_genes, rownames(scRNA_obj))
scRNA_obj <- ScaleData(
    scRNA_obj,
    features = valid_genes,
    assay = "RNA",
    verbose = FALSE
)

scRNA_obj <- ModuleEigengenes(
    scRNA_obj,
    group.by.vars = NULL,
    assay = "RNA"
)

scRNA_obj <- ModuleConnectivity(
    scRNA_obj,
    group.by = "cell_type_manual",
    group_name = "Neuroblasts",
    assay = "RNA"
)

# Rename the modules to a consistent scheme (NB_M1, NB_M2, ...). ResetModuleNames
# renames the module factor together with the kME_/ME_ columns, so it must run
# after ModuleConnectivity has created them (hdWGCNA >= 0.4.09 requires this
# order; it is also the order used in the hdWGCNA tutorial).
scRNA_obj <- ResetModuleNames(scRNA_obj, new_name = "NB_M")

# Inspect module assignments
modules <- GetModules(scRNA_obj)

ldh_info <- modules %>% filter(gene_name == "Ldh")
ldh_module <- ldh_info$module[1]
ldh_kme <- ldh_info[[paste0("kME_", ldh_module)]]

logger("Ldh is in module: %s (kME = %.3f)", ldh_module, ldh_kme)

ldh_module_genes <- modules %>%
    filter(module == ldh_module) %>%
    arrange(desc(.data[[paste0("kME_", ldh_module)]]))

logger("Module %s contains %d genes", ldh_module, nrow(ldh_module_genes))

logger("Top 10 hub genes in %s module:", ldh_module)
top_hubs_ldh <- ldh_module_genes %>%
    head(10) %>%
    select(gene_name, kME = !!paste0("kME_", ldh_module))
print(top_hubs_ldh)

# Check for other genes of interest
key_genes <- c("Ldh", "Xrp1", "Ets65A", "Hex-A", "kmr", "Tspo", "ND-75")
key_in_module <- ldh_module_genes %>%
    filter(gene_name %in% key_genes) %>%
    select(gene_name, kME = !!paste0("kME_", ldh_module))

logger("Key genes in %s module:", ldh_module)
print(key_in_module)

module_hubs_plot_path <- here(wgcna_dir, paste0(ldh_module, "_module_hubs.pdf"))
pdf(module_hubs_plot_path, width = 12, height = 8)
PlotKMEs(
    scRNA_obj,
    ncol = 5,
    n_hubs = 10,
    text_size = 3
)
dev.off()
logger("Saved module hub kME plot: %s", module_hubs_plot_path)

# Encode condition as a binary trait on the Seurat object so
# ModuleTraitCorrelation can reference it by name.
observed_conditions <- unique(as.character(scRNA_obj$condition))
if (!setequal(observed_conditions, c("Control", "ND75-KD"))) {
    stop(sprintf(
        "Unexpected condition value(s); expected Control and ND75-KD, found: %s",
        paste(sort(observed_conditions), collapse = ", ")
    ))
}
scRNA_obj$ND75_KD <- ifelse(scRNA_obj$condition == "ND75-KD", 1, 0)

# Correlate module eigengenes with the ND75_KD trait, per cell type.
scRNA_obj <- ModuleTraitCorrelation(
    scRNA_obj,
    traits = "ND75_KD",
    group.by = "cell_type_manual"
)

# Get correlation results
mt_cor_list <- GetModuleTraitCorrelation(scRNA_obj)

str(mt_cor_list)

# Extract data for Neuroblasts. hdWGCNA returns cor/pval as named vectors but
# fdr as a 1-row matrix (one trait); handle both.
nb_cor_vec <- mt_cor_list$cor$Neuroblasts
nb_pval_vec <- mt_cor_list$pval$Neuroblasts
nb_fdr_raw <- mt_cor_list$fdr$Neuroblasts
nb_fdr_vec <- if (is.matrix(nb_fdr_raw)) nb_fdr_raw[1, ] else nb_fdr_raw

# Convert to data frame
nb_cor <- data.frame(
    module = names(nb_cor_vec),
    cor = as.numeric(nb_cor_vec),
    pval = as.numeric(nb_pval_vec),
    fdr = as.numeric(nb_fdr_vec),
    stringsAsFactors = FALSE
) %>%
    arrange(desc(cor))

logger("Module correlations with ND75-KD (Neuroblasts):")
print(nb_cor %>% select(module, cor, pval, fdr))

# Highlight Ldh's module
ldh_module_cor <- nb_cor %>% filter(module == ldh_module)
logger("\n=== %s MODULE (Ldh) CORRELATION ===", toupper(ldh_module))
logger("Correlation: %.3f", ldh_module_cor$cor)
logger("P-value: %.2e", ldh_module_cor$pval)
logger("FDR: %.2e", ldh_module_cor$fdr)

# Plot 1: Bar plot of module-trait correlations
p_module_trait <- ggplot(nb_cor, aes(x = reorder(module, cor), y = cor)) +
    geom_col(aes(fill = cor), show.legend = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    # Add significance stars
    geom_text(
        aes(
            label = ifelse(fdr < 0.001, "***",
                ifelse(fdr < 0.01, "**",
                    ifelse(fdr < 0.05, "*", "")
                )
            ),
            y = cor + sign(cor) * 0.05
        ),
        size = 5
    ) +
    # Highlight Ldh's module with a circle
    geom_point(
        data = nb_cor %>% filter(module == ldh_module),
        aes(x = module, y = cor),
        color = "black", size = 4, shape = 21, fill = NA, stroke = 2
    ) +
    scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        name = "Correlation"
    ) +
    coord_flip() +
    labs(
        title = "Module-Trait Correlation: ND75-KD vs Control (Neuroblasts)",
        subtitle = sprintf("* p<0.05, ** p<0.01, *** p<0.001 | %s module (Ldh) highlighted with circle", ldh_module),
        x = "Module",
        y = "Correlation with ND75-KD"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9),
        axis.text = element_text(size = 10)
    )

module_trait_plot_path <- here(wgcna_dir, "module_ND75_KD_correlation.pdf")
ggsave(module_trait_plot_path, p_module_trait, width = 10, height = 6)
logger("Saved module-trait correlation plot: %s", module_trait_plot_path)

# Create a significance plot
p_significance <- ggplot(nb_cor, aes(x = reorder(module, cor), y = -log10(pval))) +
    geom_col(aes(fill = -log10(pval))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_point(
        data = nb_cor %>% filter(module == ldh_module),
        aes(x = module, y = -log10(pval)),
        color = "black", size = 4, shape = 21, fill = NA, stroke = 2
    ) +
    scale_fill_gradient(low = "lightblue", high = "darkred", name = "-log10(p)") +
    coord_flip() +
    labs(
        title = "Statistical Significance",
        x = "Module",
        y = "-log10(p-value)"
    ) +
    theme_classic() +
    theme(axis.text.y = element_blank())

# Combine both plots
combined_cor_plot <- p_module_trait | p_significance

combined_cor_plot_path <- here(wgcna_dir, "module_ND75_KD_correlation_combined.pdf")
ggsave(combined_cor_plot_path, combined_cor_plot, width = 14, height = 6)
logger("Saved combined correlation plot: %s", combined_cor_plot_path)

logger("Saving hdWGCNA results for SCENIC integration...")

# 1. Save the complete Seurat object with all hdWGCNA metadata
complete_obj_path <- here(wgcna_dir, "hdWGCNA_neuroblast_complete_analysis_with_condition.rds")
saveRDS(scRNA_obj, complete_obj_path)
logger("Saved complete hdWGCNA Seurat object: %s", complete_obj_path)

# 2. Extract and save module assignments
modules <- GetModules(scRNA_obj)
module_assignments_rds <- here(wgcna_dir, "hdWGCNA_module_assignments_with_condition.rds")
module_assignments_csv <- here(wgcna_dir, "hdWGCNA_module_assignments_with_condition.csv")
saveRDS(modules, module_assignments_rds)
write.csv(modules, module_assignments_csv, row.names = FALSE)
logger("Saved module assignments: %s, %s", module_assignments_rds, module_assignments_csv)

# 3. Save module eigengenes
MEs <- GetMEs(scRNA_obj, harmonized = FALSE)
module_eigengenes_path <- here(wgcna_dir, "hdWGCNA_module_eigengenes_with_condition.rds")
saveRDS(MEs, module_eigengenes_path)
logger("Saved module eigengenes: %s", module_eigengenes_path)

# 4. Save network connectivity (TOM matrix)
TOM <- GetTOM(scRNA_obj)
tom_matrix_path <- here(wgcna_dir, "hdWGCNA_TOM_matrix_with_condition.rds")
saveRDS(TOM, tom_matrix_path)
logger("Saved TOM matrix: %s", tom_matrix_path)

# 5. Create module-specific gene lists for SCENIC
module_gene_lists <- modules %>%
    split(., .$module) %>%
    lapply(function(x) {
        x %>%
            arrange(desc(.data[[paste0("kME_", x$module[1])]])) %>%
            select(gene_name, starts_with("kME_"))
    })

module_gene_lists_path <- here(wgcna_dir, "hdWGCNA_module_gene_lists_with_condition.rds")
saveRDS(module_gene_lists, module_gene_lists_path)
logger("Saved module-specific gene lists: %s", module_gene_lists_path)

# 6. Extract the Ldh module specifically
ldh_module_genes <- modules %>%
    filter(module == ldh_module) %>%
    arrange(desc(.data[[paste0("kME_", ldh_module)]])) %>%
    select(gene_name, kME = !!paste0("kME_", ldh_module), module)

ldh_module_genes_csv <- here(wgcna_dir, paste0(ldh_module, "_module_genes_for_SCENIC_with_condition.csv"))
write.csv(ldh_module_genes, ldh_module_genes_csv, row.names = FALSE)

# Create a simple gene list (just names) for SCENIC filtering
ldh_genes_only <- ldh_module_genes$gene_name
ldh_module_genes_txt <- here(wgcna_dir, paste0(ldh_module, "_module_genes_list_with_condition.txt"))
writeLines(ldh_genes_only, ldh_module_genes_txt)
logger(
    "Saved %s module genes (%d genes): %s, %s",
    ldh_module, length(ldh_genes_only), ldh_module_genes_csv, ldh_module_genes_txt
)

# 7. Save hub genes (top 25 per module by kME)
hub_genes_all <- GetHubGenes(scRNA_obj, n_hubs = 25)
hub_genes_path <- here(wgcna_dir, "hdWGCNA_hub_genes_with_condition.rds")
saveRDS(hub_genes_all, hub_genes_path)

# Extract hub genes specifically from the Ldh module
ldh_hub_genes <- ldh_module_genes %>%
    head(25) %>%
    pull(gene_name)

ldh_hub_genes_txt <- here(wgcna_dir, paste0(ldh_module, "_hub_genes_with_condition.txt"))
writeLines(ldh_hub_genes, ldh_hub_genes_txt)
logger("Saved hub genes: %s, %s", hub_genes_path, ldh_hub_genes_txt)

# 8. Create a SCENIC-ready target list
# This includes ALL genes in modules that correlate with ND75-KD
significant_modules <- nb_cor %>%
    filter(fdr < 0.05) %>%
    pull(module)

scenic_target_genes <- modules %>%
    filter(module %in% significant_modules) %>%
    pull(gene_name) %>%
    unique()

scenic_target_genes_path <- here(wgcna_dir, "ND75_KD_associated_genes_for_SCENIC_with_condition.txt")
writeLines(scenic_target_genes, scenic_target_genes_path)
logger(
    "Saved %d genes from ND75-KD-associated modules: %s",
    length(scenic_target_genes), scenic_target_genes_path
)

# 9. Save correlation results for reference
module_trait_correlations_path <- here(wgcna_dir, "module_trait_correlations_with_condition.rds")
saveRDS(list(neuroblasts = nb_cor), module_trait_correlations_path)
logger("Saved module-trait correlations: %s", module_trait_correlations_path)

# 10. Create a summary file for easy reference
summary_info <- list(
    analysis_date = Sys.Date(),
    n_modules = length(unique(modules$module)),
    ldh_module = ldh_module,
    ldh_module_size = nrow(ldh_module_genes),
    ldh_module_correlation = ldh_module_cor$cor,
    ldh_module_fdr = ldh_module_cor$fdr,
    significant_modules = significant_modules,
    hub_genes_ldh = ldh_hub_genes,
    key_genes_in_ldh = key_in_module$gene_name
)

summary_rds <- here(wgcna_dir, "analysis_summary_with_condition.rds")
saveRDS(summary_info, summary_rds)
logger("Saved analysis summary (rds): %s", summary_rds)

# Save as human-readable YAML-style text
summary_txt <- here(wgcna_dir, "analysis_summary.txt")
sink(summary_txt)
cat("=== hdWGCNA Analysis Summary ===\n\n")
cat(sprintf("Analysis Date: %s\n", summary_info$analysis_date))
cat(sprintf("Total Modules: %d\n", summary_info$n_modules))
cat(sprintf("Ldh Module: %s\n", summary_info$ldh_module))
cat(sprintf("Ldh Module Size: %d genes\n", summary_info$ldh_module_size))
cat(sprintf(
    "Ldh Module Correlation: %.3f (FDR = %.2e)\n\n",
    summary_info$ldh_module_correlation,
    summary_info$ldh_module_fdr
))
cat("Significant Modules (FDR < 0.05):\n")
cat(paste("  -", significant_modules, collapse = "\n"), "\n\n")
cat("Top Hub Genes in", ldh_module, "module:\n")
cat(paste("  -", ldh_hub_genes[1:10], collapse = "\n"), "\n\n")
cat("Key Genes in", ldh_module, "module:\n")
cat(paste("  -", key_in_module$gene_name, collapse = "\n"), "\n")
sink()
logger("Saved analysis summary (text): %s", summary_txt)

logger("=== SAVE COMPLETE ===")
logger("Files saved to: %s", wgcna_dir)
logger("\nKey outputs for SCENIC:")
logger("  1. Complete analysis: %s", complete_obj_path)
logger("  2. Module genes: %s", ldh_module_genes_csv)
logger("  3. Hub genes: %s", ldh_hub_genes_txt)
logger("  4. ND75-KD genes: %s", scenic_target_genes_path)
