# Metabolomics: per-metabolite differential abundance
#
# Larval brain untargeted metabolomics, ND-42 (Complex I) RNAi vs. Control in
# NSCs. Abundances are normalized to the internal standard (D6-Glutaric Acid)
# and then total-sum normalized upstream (the input workbook); they are then
# log2-transformed with a tiny offset and tested per metabolite with limma.

library(here)
library(readxl)
library(tidyverse)
library(limma)

# Console + file logging, mirroring the other pipeline stages.
log_file <- "00_differential_abundance.log"

logger <- function(message, ...) {
    line <- paste0(Sys.time(), " - ", sprintf(message, ...), "\n")
    cat(line)
    cat(line, file = log_file, append = TRUE)
}

input_path <- here::here(
    "metabolomics", "data", "raw", "metabolite_abundances_normalized.xlsx"
)
output_dir <- here::here("results", "metabolomics")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_path)) {
    stop(sprintf(
        "%s not found. See metabolomics/data/raw/README.md.",
        input_path
    ))
}

# ---- Load data ----
metabolomics_data <- read_xlsx(input_path)

meta_info <- metabolomics_data %>%
    select(Molecule, Molecule.List)

count_df <- metabolomics_data %>%
    select(starts_with("mcherry"), starts_with("ND42"))

count_mat <- as.matrix(count_df)
rownames(count_mat) <- meta_info$Molecule

# Log2 transform with a tiny offset to handle very small abundances
log_mat <- log2(count_mat + 1e-20)

# ---- Experimental design ----
#  6 mcherry controls and 6 ND42 knockdowns
groups <- factor(c(rep("Control", 6), rep("Condition", 6)))
groups <- relevel(groups, ref = "Control")

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

# ---- Fit linear model (limma) ----
fit <- lmFit(log_mat, design)

# Contrast: Condition minus Control
contrast_mat <- makeContrasts(Condition - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)

# Empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results (Benjamini-Hochberg adjusted p-values)
res_df <- topTable(fit2, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("gene") %>%
    left_join(meta_info, by = c("gene" = "Molecule")) %>%
    rename(Pathway = Molecule.List)

# Add per-group mean abundance (used by the volcano plot in
# reports/plot_manuscript_figures.qmd to annotate average abundance).
control_cols <- grep("mcherry", colnames(count_mat))
condition_cols <- grep("ND42", colnames(count_mat))

control_means <- rowMeans(count_mat[, control_cols], na.rm = TRUE)
condition_means <- rowMeans(count_mat[, condition_cols], na.rm = TRUE)

res_df <- res_df %>%
    mutate(
        avg_control = control_means[gene],
        avg_condition = condition_means[gene]
    )

# ---- Save ----
write_csv(
    res_df,
    file.path(output_dir, "metabolomics_differential_abundance_results.csv")
)

saveRDS(count_mat, file.path(output_dir, "metabolomics_count_matrix.rds"))
saveRDS(meta_info, file.path(output_dir, "metabolomics_meta_info.rds"))
saveRDS(log_mat, file.path(output_dir, "metabolomics_log2_count_matrix.rds"))

logger("Wrote differential abundance results and matrix objects to %s", output_dir)
logger("Top metabolites by P-value:")
print(res_df %>% arrange(P.Value) %>% head())
