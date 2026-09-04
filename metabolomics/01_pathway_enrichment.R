# Metabolomics: pathway-level enrichment (MSEA)
#
# Metabolite-set enrichment on the per-metabolite differential abundance
# ranking, with metabolite classes as the gene sets. Ranking metric is the
# signed -log10(P-value); enrichment is fgsea.

library(here)
library(tidyverse)
library(fgsea)

log_file <- "01_pathway_enrichment.log"

logger <- function(message, ...) {
    line <- paste0(Sys.time(), " - ", sprintf(message, ...), "\n")
    cat(line)
    cat(line, file = log_file, append = TRUE)
}

processed_dir <- here::here("metabolomics", "data", "processed")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
input_path <- file.path(processed_dir, "metabolomics_differential_abundance_results.csv")

if (!file.exists(input_path)) {
    stop(sprintf(
        "%s not found. Run metabolomics/00_differential_abundance.R first.",
        input_path
    ))
}

res_df <- read_csv(input_path, show_col_types = FALSE)

run_msea <- function(res_df) {
    # Create ranked list
    res_df <- res_df %>%
        filter(!is.na(Pathway), Pathway != "Internal Std") %>%
        mutate(rank_metric = sign(logFC) * -log10(P.Value))

    ranks <- setNames(res_df$rank_metric, res_df$gene)
    ranks <- sort(ranks, decreasing = TRUE)

    # Create pathway list (metabolite class -> its metabolites)
    pathways <- split(res_df$gene, res_df$Pathway)

    # Run fgsea (new API - no nPermutations)
    set.seed(42)
    fgsea_res <- fgsea(
        pathways = pathways,
        stats = ranks,
        minSize = 2,
        maxSize = 500
    )

    fgsea_res <- fgsea_res %>%
        arrange(pval) %>%
        mutate(
            direction = ifelse(ES > 0, "Up", "Down"),
            neg_log10_p = -log10(pval)
        )

    return(fgsea_res)
}

msea_results <- run_msea(res_df)

# fgsea returns leadingEdge as a list-column; collapse it for CSV storage.
msea_results <- msea_results %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write_csv(msea_results, file.path(processed_dir, "metabolomics_fgsea_results.csv"))

logger(
    "Wrote %d pathways to %s", nrow(msea_results),
    file.path(processed_dir, "metabolomics_fgsea_results.csv")
)
print(msea_results %>% select(pathway, ES, NES, pval, padj, size))
