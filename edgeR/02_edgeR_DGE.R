#  Script to perform edgeR DGE analysis

# %%
# Import libraries
library(here)
library(edgeR)
library(readr)
library(tidyverse)
library(statmod)
library(EnhancedVolcano)

pseudobulk_dir <- here::here("results", "edgeR", "pseudobulk_all_genes")
output_dir <- here::here("results", "edgeR", "edgeR_results")
dir.create(output_dir, showWarnings = FALSE)

# Console + file logging, mirroring the Python stages' logging setup.
log_file <- "02_edgeR_DGE.log"

logger <- function(message, ...) {
    line <- paste0(Sys.time(), " - ", sprintf(message, ...), "\n")
    cat(line)
    cat(line, file = log_file, append = TRUE)
}


lfc_threshold_ev <- 1.0
fdr_threshold_ev <- 0.05
min_cells_required <- 5

dirs <- dir(pseudobulk_dir, full.names = TRUE)
dirs <- dirs[basename(dirs) != "edgeR_results"]

for (dir in dirs) {
    logger("Starting analysis for %s", dir)

    counts_file <- list.files(
        dir,
        pattern = "_pseudobulk_counts\\.tsv$",
        full.names = TRUE
    )
    sample_info_file <- list.files(
        dir,
        pattern = "sample_info_",
        full.names = TRUE
    )
    # Load cell count info
    pct_file <- list.files(
        dir,
        pattern = "_percent_expressing\\.tsv$",
        full.names = TRUE
    )

    if (length(pct_file) == 0) {
        logger("WARN: percent-expressing file not found for %s. Run edgeR/01_create_pct_expressing.py first; proceeding with is_artifact = FALSE.", dir)
        cell_counts_df <- NULL
    } else {
        cell_counts_df <- read_tsv(pct_file, show_col_types = FALSE) %>%
            select(gene, n_cells_Control, n_cells_ND75_KD)
    }

    #  Load counts data
    counts_data <- read_tsv(
        counts_file,
        col_types = cols(
            gene = col_character(),
            .default = col_double()
        )
    ) %>%
        column_to_rownames(var = "gene")

    counts_matrix <- as.matrix(counts_data)

    # Load sample info
    sample_info <- read_tsv(
        sample_info_file,
        col_types = cols(sample_id = col_character(), .default = col_character())
    ) %>%
        column_to_rownames(var = "sample_id")

    # Reorder sample_info to match column order of counts_matrix
    sample_info <- sample_info[colnames(counts_matrix), , drop = FALSE]

    # Create DGEList object
    sample_info$condition <- factor(
        sample_info$condition,
        levels = c("Control", "ND75-KD")
    )
    if (any(is.na(sample_info$condition))) {
        stop("Unexpected condition value; expected Control or ND75-KD.")
    }
    sample_info$batch <- factor(
        sample_info$batch,
        levels = c("Batch1", "Batch2")
    )

    y <- DGEList(
        counts = counts_matrix,
        samples = sample_info,
        group = sample_info$condition
    )

    ct_output_dir <- here(output_dir, basename(dir))
    dir.create(ct_output_dir, showWarnings = FALSE)

    # Design matrix. R sanitizes the hyphen in ND75-KD to a period in the
    # model-matrix column name, so select the condition coefficient from the
    # design itself rather than constructing a name from the biological label.
    design <- model.matrix(~ batch + condition, data = y$samples)
    condition_coefficients <- grep("^condition", colnames(design), value = TRUE)
    if (length(condition_coefficients) != 1) {
        stop(sprintf(
            "Expected exactly one condition coefficient; found: %s",
            paste(condition_coefficients, collapse = ", ")
        ))
    }

    if (ncol(design) >= nrow(y$samples)) {
        logger("Design matrix has more columns than rows. Skipping analysis for %s.", dir)
        next
    }

    # Filter lowly expressed genes
    keep <- filterByExpr(y, design = design)
    y_filtered <- y[keep, , keep.lib.sizes = FALSE]

    logger("Filtering: %d genes kept out of %d", sum(keep), nrow(y))
    if (sum(keep) < 2) {
        logger("Not enough genes to perform analysis. Skipping %s.", dir)
        next
    }

    # Normalize counts
    y_norm <- calcNormFactors((y_filtered))

    # Plot MDS (Multidimensional Scaling) plot
    mds_plot_file <- file.path(
        ct_output_dir,
        sprintf("%s_MDS_plot.png", basename(dir))
    )

    png(mds_plot_file, width = 800, height = 600, res = 150)
    plotMDS(y_norm, col = as.numeric(y_norm$samples$condition), pch = 19)
    title(main = paste("MDS Plot for", basename(dir)))
    dev.off()

    # Estimate dispersion
    y_disp <- estimateDisp(y_norm, design = design, robust = TRUE)

    # Plot Biological Coefficient of Variation (BCV)
    bcv_plot_file <- file.path(
        ct_output_dir,
        sprintf("%s_BCV_plot.png", basename(dir))
    )
    png(bcv_plot_file, width = 800, height = 600, res = 150)
    plotBCV(y_disp)
    title(main = paste("BCV Plot for", basename(dir)))
    dev.off()

    # Fit GLM (Quasi-Likelihood Model)
    fit <- glmQLFit(y_disp, design, robust = TRUE)

    # Plot QL dispersion
    ql_disp_plot_file <- file.path(
        ct_output_dir,
        sprintf("%s_QL_dispersion_plot.png", basename(dir))
    )
    png(ql_disp_plot_file, width = 800, height = 600, res = 150)
    plotQLDisp(fit)
    title(main = paste("QL Dispersion Plot for", basename(dir)))
    dev.off()

    # Perform differential expression analysis
    coef_to_test <- condition_coefficients[[1]]

    qlf <- glmQLFTest(fit, coef = coef_to_test)

    # Get results
    top_degs_table <- topTags(
        qlf,
        n = Inf,
        adjust.method = "BH",
        sort.by = "PValue"
    )

    top_degs_table_df <- as_tibble(top_degs_table$table, rownames = "gene")

    if (!is.null(cell_counts_df)) {
        top_degs_table_df <- top_degs_table_df %>%
            left_join(cell_counts_df, by = "gene") %>%
            mutate(
                # Logic: Is this result suspicious?
                # If upregulated (logFC > 0): trust only if n_cells_ND75_KD >= 5
                # If downregulated (logFC < 0): trust only if n_cells_Control >= 5
                is_artifact = case_when(
                    logFC > 0 & n_cells_ND75_KD < min_cells_required ~ TRUE,
                    logFC < 0 & n_cells_Control < min_cells_required ~ TRUE,
                    TRUE ~ FALSE
                )
            )

        # Count how many significant genes are artifacts
        n_artifacts <- sum(top_degs_table_df$FDR < fdr_threshold_ev & top_degs_table_df$is_artifact, na.rm = TRUE)
        logger("Identified %d significant DEGs as potential single-cell artifacts (< %d cells).", n_artifacts, min_cells_required)
    } else {
        top_degs_table_df$is_artifact <- FALSE
        top_degs_table_df$n_cells_Control <- NA
        top_degs_table_df$n_cells_ND75_KD <- NA
    }

    logger("Found %d DEGs with FDR < 0.05 for %s.", sum(top_degs_table_df$FDR < fdr_threshold_ev), basename(dir))

    #  Save results
    results_file <- file.path(
        ct_output_dir,
        sprintf("%s_edgeR_results.tsv", basename(dir))
    )
    write_tsv(
        top_degs_table_df,
        results_file
    )

    logger("Finished analysis for %s", dir)

    # Volcano Plot - Highlighting Artifacts
    plot_title_suffix_ev <- levels(y_disp$samples$condition)[2]
    main_plot_title_ev <- paste(
        basename(dir),
        "-",
        plot_title_suffix_ev,
        "vs",
        levels(y_disp$samples$condition)[1]
    )

    # Filter out artifacts for the plot OR color them differently
    # Here we plot everything but you can filter 'top_degs_table_df' before passing if preferred
    clean_plot_data <- top_degs_table_df %>% filter(!is_artifact)

    ev_plot_object <- EnhancedVolcano(
        toptable = clean_plot_data, # Plotting only clean data
        lab = clean_plot_data$gene,
        x = "logFC",
        y = "FDR",
        title = paste(main_plot_title_ev, "(Artifacts Removed)"),
        subtitle = paste0(
            "LFC cutoff: |", lfc_threshold_ev,
            "|; FDR cutoff: ", fdr_threshold_ev,
            "; Min Cells: ", min_cells_required
        ),
        pCutoff = fdr_threshold_ev,
        FCcutoff = lfc_threshold_ev,
        pointSize = 2.5,
        labSize = 3.0
    )

    # Save volcano plot
    volcano_plot_file <- file.path(
        ct_output_dir,
        sprintf("%s_edgeR_volcano_plot.png", basename(dir))
    )
    png(volcano_plot_file, width = 1600, height = 1600, res = 150)
    print(ev_plot_object)
    dev.off()
}
