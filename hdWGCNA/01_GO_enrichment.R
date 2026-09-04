# Script to perform GO enrichment and functional labelling of the hdWGCNA
# Neuroblast co-expression modules (NB_M1..NB_M11).

# Import libraries
library(here)
library(tidyverse)
library(Seurat)
library(hdWGCNA)
library(org.Dm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(cowplot)
library(grid)
library(scales)

# enrichplot's treeplot term clustering and layout are deterministic, but a
# seed is set for reproducibility and to mirror 00_hdWGCNA_analysis.R.
set.seed(42)

# --- 1. PATHS AND LOGGING ---

# wgcna_dir: produced by hdWGCNA/00_hdWGCNA_analysis.R (see hdWGCNA/README.md)
wgcna_dir  <- here::here("results", "hdWGCNA")
output_dir <- here::here("results", "hdWGCNA", "GO_exploration_top50")

rds_dir <- file.path(output_dir, "GO_objects")
pdf_dir <- file.path(output_dir, "GO_analysis")
csv_dir <- file.path(output_dir, "label_review")

for (d in c(rds_dir, pdf_dir, csv_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

prefix        <- "NB_hdWGCNA"
ontologies    <- c("BP", "MF", "CC")
module_order  <- paste0("NB_M", 1:11)

# Console + file logging
log_file <- "01_GO_enrichment.log"

logger <- function(message, ...) {
    line <- paste0(Sys.time(), " - ", sprintf(message, ...), "\n")
    cat(line)
    cat(line, file = log_file, append = TRUE)
}

# --- 2. LOAD hdWGCNA OBJECT AND DERIVE HUB GENE LISTS ---

hdWGCNA_rds_path <- file.path(
    wgcna_dir, "hdWGCNA_neuroblast_complete_analysis_with_condition.rds"
)
if (!file.exists(hdWGCNA_rds_path)) {
    stop(sprintf(
        "%s not found. Run hdWGCNA/00_hdWGCNA_analysis.R first.",
        hdWGCNA_rds_path
    ))
}

logger("Loading hdWGCNA object: %s", hdWGCNA_rds_path)
hdWGCNA_results <- readRDS(hdWGCNA_rds_path)

# Module names (NB_M1..NB_Mn) are already set by 00_hdWGCNA_analysis.R.
nb_modules  <- GetModules(hdWGCNA_results)
MEs_all     <- GetMEs(hdWGCNA_results)
all_modules <- setdiff(colnames(MEs_all), "grey")

if (!setequal(all_modules, module_order)) {
    stop(sprintf(
        "Module set mismatch: expected %s, found %s",
        paste(module_order, collapse = ", "),
        paste(sort(all_modules), collapse = ", ")
    ))
}
all_modules <- module_order
logger("Modules found: %s", paste(all_modules, collapse = ", "))

# Hub genes: top 50 per module by kME. The lollipop plot uses
# the top 10; the wider set here is for GO power.
HUB_N <- 50

nb_hub_genes <- lapply(all_modules, function(mod) {
    kme_col <- paste0("kME_", mod)
    if (!kme_col %in% colnames(nb_modules)) {
        logger("  Warning: %s not found in module table", kme_col)
        return(NULL)
    }
    nb_modules %>%
        filter(module == mod) %>%
        arrange(desc(.data[[kme_col]])) %>%
        head(HUB_N) %>%
        pull(gene_name)
})
names(nb_hub_genes) <- all_modules

walk2(names(nb_hub_genes), nb_hub_genes, function(mod, genes) {
    logger("  %s: %d hub genes", mod, length(genes))
})

# Background universe: all non-grey module genes.
network_universe <- nb_modules %>%
    filter(module != "grey") %>%
    pull(gene_name) %>%
    unique()
logger("Background universe: %d non-grey module genes", length(network_universe))

# --- 3. GOSemSim SEMANTIC-SIMILARITY DATABASES ---

logger("Building GOSemSim databases (BP, MF, CC)...")
d_sim_list <- lapply(setNames(ontologies, ontologies), function(ont) {
    GOSemSim::godata(
        annoDb    = "org.Dm.eg.db",
        ont       = ont,
        computeIC = TRUE,
        keytype   = "SYMBOL"
    )
})
logger("GOSemSim databases ready.")

# --- 4. GO OVER-REPRESENTATION PER MODULE x ONTOLOGY ---
#
# Ontology-aware thresholds. BP is denser than MF/CC in Drosophila,
# so BP uses stricter cutoffs:
#   BP      : p < 0.05, q < 0.20, minGSSize = 5, Wang simplify cutoff 0.85
#   MF, CC  : p < 0.10, q < 0.30, minGSSize = 3, Wang simplify cutoff 0.90
# maxGSSize = 500, pAdjustMethod = "BH", background universe = network_universe.

TEXT_SIZE <- 14

go_thresholds <- function(ont) {
    if (ont == "BP") {
        list(pvalue = 0.05, qvalue = 0.20, minGSSize = 5, simplify = 0.85)
    } else {
        list(pvalue = 0.10, qvalue = 0.30, minGSSize = 3, simplify = 0.90)
    }
}

logger("=== Running GO enrichment pipeline ===")
all_top_terms <- list()

for (mod_name in names(nb_hub_genes)) {

    hubs <- nb_hub_genes[[mod_name]]

    if (is.null(hubs) || length(hubs) < 5) {
        logger("[%s] Skipping: only %d hub genes", mod_name, length(hubs))
        next
    }

    logger("[%s] Processing %d hub genes...", mod_name, length(hubs))

    pdf_path <- file.path(pdf_dir, paste0(prefix, "_Explore_", mod_name, ".pdf"))
    pdf(pdf_path, width = 22, height = 18)

    tryCatch({

        for (ont in ontologies) {

            th <- go_thresholds(ont)
            logger("  [%s | %s] enrichGO (p<%.2f, q<%.2f, minGSSize=%d)...",
                   mod_name, ont, th$pvalue, th$qvalue, th$minGSSize)

            ego <- tryCatch({
                enrichGO(
                    gene          = hubs,
                    universe      = network_universe,
                    OrgDb         = org.Dm.eg.db,
                    ont           = ont,
                    keyType       = "SYMBOL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = th$pvalue,
                    qvalueCutoff  = th$qvalue,
                    minGSSize     = th$minGSSize,
                    maxGSSize     = 500
                )
            }, error = function(e) {
                logger("    enrichGO error: %s", e$message)
                return(NULL)
            })

            if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
                logger("  [%s | %s] No significant terms.", mod_name, ont)
                next
            }

            n_full <- nrow(as.data.frame(ego))
            logger("  [%s | %s] %d terms before simplification", mod_name, ont, n_full)

            ego_simple <- clusterProfiler::simplify(
                ego,
                cutoff     = th$simplify,
                by         = "p.adjust",
                select_fun = min
            )
            n_simple <- nrow(as.data.frame(ego_simple))
            logger("  [%s | %s] %d terms after simplification", mod_name, ont, n_simple)

            saveRDS(ego,
                    file.path(rds_dir, paste0(prefix, "_", mod_name, "_", ont, "_full.rds")))
            saveRDS(ego_simple,
                    file.path(rds_dir, paste0(prefix, "_", mod_name, "_", ont, "_simplified.rds")))

            # Diagnostic plots + treeplot RDS (reused by the composite panel).
            tryCatch({
                ego_sim        <- pairwise_termsim(ego,        semData = d_sim_list[[ont]])
                ego_simple_sim <- pairwise_termsim(ego_simple, semData = d_sim_list[[ont]])

                suppressWarnings({
                    if (n_simple >= 2) {
                        n_clusters <- min(5, n_simple)

                        p_tree_base <- treeplot(
                            ego_simple_sim,
                            showCategory   = min(15, n_simple),
                            cluster.params = list(n = n_clusters),
                            fontsize       = 5
                        )
                        saveRDS(
                            p_tree_base,
                            file.path(rds_dir,
                                      paste0(prefix, "_", mod_name, "_", ont, "_treeplot.rds"))
                        )

                        p_tree <- p_tree_base +
                            ggtitle(sprintf("%s | %s - Hierarchy (n=%d terms)",
                                            mod_name, ont, n_simple)) +
                            theme(
                                text       = element_text(size = TEXT_SIZE),
                                plot.title = element_text(size = TEXT_SIZE, face = "bold")
                            )
                        print(p_tree)
                    } else {
                        p_dot <- dotplot(ego_simple, showCategory = n_simple) +
                            ggtitle(sprintf("%s | %s - Dotplot (n=%d terms)",
                                            mod_name, ont, n_simple)) +
                            theme_minimal(base_size = TEXT_SIZE)
                        print(p_dot)
                    }

                    if (n_simple >= 1) {
                        print(
                            cnetplot(
                                ego_simple,
                                showCategory       = min(5, n_simple),
                                cex_label_gene     = 1.0,
                                cex_label_category = 1.5
                            ) +
                                ggtitle(sprintf("%s | %s - Gene-Concept Network", mod_name, ont)) +
                                theme(text = element_text(size = TEXT_SIZE))
                        )
                    }

                    if (n_simple >= 2) {
                        print(
                            emapplot(ego_simple_sim, showCategory = min(10, n_simple),
                                     cex_label_category = 1.5) +
                                ggtitle(sprintf("%s | %s - Enrichment Map", mod_name, ont)) +
                                theme(text = element_text(size = TEXT_SIZE))
                        )
                    }

                    p_bar <- barplot(ego_simple, showCategory = min(15, n_simple)) +
                        ggtitle(sprintf("%s | %s - Top Terms", mod_name, ont)) +
                        theme_minimal(base_size = TEXT_SIZE)
                    print(p_bar)
                })

            }, error = function(e) {
                logger("    Plotting error for %s | %s: %s", mod_name, ont, e$message)
            })

            all_top_terms[[paste0(mod_name, "_", ont)]] <- as.data.frame(ego_simple) %>%
                arrange(p.adjust) %>%
                head(5) %>%
                mutate(module = mod_name, ontology = ont, n_hub_genes = length(hubs))

        } # end ontology loop

    }, finally = {
        dev.off()
        logger("  [%s] Diagnostic PDF written: %s", mod_name, pdf_path)
    })

} # end module loop

top_terms_path <- file.path(csv_dir, "GO_top_terms_per_module.csv")
if (length(all_top_terms) > 0) {
    bind_rows(all_top_terms) %>%
        dplyr::select(module, ontology, n_hub_genes, ID, Description,
                      GeneRatio, p.adjust, qvalue, Count) %>%
        write_csv(top_terms_path)
    logger("Top terms per module written: %s", top_terms_path)
} else {
    logger("WARNING: no enriched terms collected; %s not written", top_terms_path)
}

# --- 5. TERM-COUNT AUDIT ---

logger("=== Term count audit ===")
audit <- purrr::map_dfr(all_modules, function(mod) {
    purrr::map_dfr(ontologies, function(ont) {
        fp <- file.path(rds_dir, paste0(prefix, "_", mod, "_", ont, "_simplified.rds"))
        if (!file.exists(fp)) {
            return(tibble(module = mod, ont = ont, n_terms = 0L, status = "NO FILE"))
        }
        n <- nrow(as.data.frame(readRDS(fp)))
        status <- dplyr::case_when(
            n == 0 ~ "EMPTY",
            n == 1 ~ "SINGLE TERM",
            n <  5 ~ "FEW TERMS",
            TRUE   ~ "OK"
        )
        tibble(module = mod, ont = ont, n_terms = n, status = status)
    })
})
print(audit, n = Inf)
write_csv(audit, file.path(csv_dir, "term_count_audit.csv"))
logger("Term count audit written: %s", file.path(csv_dir, "term_count_audit.csv"))

# --- 6. AUTO LABEL DERIVATION (for provenance / review) ---
#
# For each module, take the top non-generic simplified term per ontology and
# pick a primary label (BP preferred, then MF, then CC). These auto labels are
# recorded for provenance; the authoritative module labels are the curated
# manuscript set in section 8. NB_M10 is annotated from hub-gene identity.

logger("=== Deriving auto labels from enrichment results ===")

generic_terms <- c(
    "biological process", "cellular process", "metabolic process",
    "cellular metabolic process", "organic substance metabolic process",
    "macromolecule metabolic process", "primary metabolic process",
    "nitrogen compound metabolic process",
    "cellular nitrogen compound metabolic process", "gene expression",
    "biosynthetic process", "cellular biosynthetic process",
    "regulation of biological process", "regulation of cellular process",
    "positive regulation of biological process",
    "negative regulation of biological process", "molecular function",
    "catalytic activity", "binding", "molecular adaptor activity", "transport",
    "localization", "intracellular", "cell", "intracellular organelle",
    "organelle", "membrane"
)

extract_top_terms <- function(fp, n_top = 5) {
    if (!file.exists(fp)) return(NULL)
    df <- as.data.frame(readRDS(fp))
    if (nrow(df) == 0) return(NULL)
    df %>%
        filter(!tolower(Description) %in% tolower(generic_terms)) %>%
        filter(nchar(Description) > 8) %>%
        arrange(p.adjust, desc(Count)) %>%
        head(n_top) %>%
        dplyr::select(ID, Description, p.adjust, Count, GeneRatio)
}

clean_term <- function(term) {
    term %>%
        str_remove_all("\\bprocess\\b") %>%
        str_remove_all("\\bregulation of\\b") %>%
        str_remove_all("\\bpositive regulation of\\b") %>%
        str_remove_all("\\bnegative regulation of\\b") %>%
        str_remove_all("\\bcellular\\b") %>%
        str_remove_all("\\bbiological\\b") %>%
        str_remove_all("\\bmacromolecule\\b") %>%
        str_replace_all("\\s{2,}", " ") %>%
        str_trim() %>%
        str_to_sentence()
}

derived_labels <- purrr::map_dfr(all_modules, function(mod) {

    if (mod == "NB_M10") {
        return(tibble(
            module          = mod,
            label_auto      = "mtDNA-encoded ETC subunits",
            source          = "Manual - hub gene identity (mtDNA)",
            top_bp          = "mt:CoI/II/III, mt:Cyt-b, mt:ATPase6, mt:ND1/3/4/5",
            top_mf          = NA_character_,
            top_cc          = NA_character_,
            label_term_id   = NA_character_,
            label_term_padj = NA_real_
        ))
    }

    bp_terms <- extract_top_terms(
        file.path(rds_dir, paste0(prefix, "_", mod, "_BP_simplified.rds")))
    mf_terms <- extract_top_terms(
        file.path(rds_dir, paste0(prefix, "_", mod, "_MF_simplified.rds")))
    cc_terms <- extract_top_terms(
        file.path(rds_dir, paste0(prefix, "_", mod, "_CC_simplified.rds")))

    primary <- if (!is.null(bp_terms) && nrow(bp_terms) > 0) {
        bp_terms[1, ]
    } else if (!is.null(mf_terms) && nrow(mf_terms) > 0) {
        mf_terms[1, ]
    } else if (!is.null(cc_terms) && nrow(cc_terms) > 0) {
        cc_terms[1, ]
    } else {
        NULL
    }

    tibble(
        module          = mod,
        label_auto      = if (!is.null(primary)) clean_term(primary$Description)
                          else paste0(mod, " - no significant terms"),
        source          = "GO enrichment (auto)",
        top_bp          = if (!is.null(bp_terms)) paste(head(bp_terms$Description, 3), collapse = "; ") else NA_character_,
        top_mf          = if (!is.null(mf_terms)) paste(head(mf_terms$Description, 3), collapse = "; ") else NA_character_,
        top_cc          = if (!is.null(cc_terms)) paste(head(cc_terms$Description, 3), collapse = "; ") else NA_character_,
        label_term_id   = if (!is.null(primary)) primary$ID       else NA_character_,
        label_term_padj = if (!is.null(primary)) primary$p.adjust else NA_real_
    )
})

# --- 7. CONSENSUS-LABEL RESOLVER ---
#
# The per-cluster tip labels on the Supp Fig 5D treeplots were shortened by hand
# against the BP/MF/CC treeplot hierarchies during the original analysis. That
# curation is frozen in the `consensus_term_overrides` map in section 9 (applied
# by shorten_consensus_label()). This repository records the analysis as it was
# run, so there is no interactive review round-trip: `consensus_label_resolved`
# is left empty and every tip label is resolved from the frozen override map.

consensus_label_resolved <- tibble::tibble(
    module         = character(),
    ont            = character(),
    cluster        = integer(),
    label_resolved = character()
)

# --- 8. CURATED MANUSCRIPT MODULE LABELS ---
#
# These are assigned from the top non-generic enriched terms and
# refined manually against the GO treeplot hierarchies, not recomputed here.
# This vector is the authoritative copy;
# hdWGCNA/resources/hdWGCNA_module_manuscript_labels.csv holds the same
# content in table form (label_source = figure_2E).

module_labels <- c(
    NB_M1  = "Mitochondrial ATP synthesis",
    NB_M2  = "Protein degradation & Chromatin",
    NB_M3  = "Cytoplasmic translation",
    NB_M4  = "Synaptic signalling & Ion transport",
    NB_M5  = "Chromatin remodelling & Gene silencing",
    NB_M6  = "Transcriptional repression",
    NB_M7  = "Stress response & Growth signalling",
    NB_M8  = "Epigenetic silencing & Chromatin",
    NB_M9  = "Mitochondrial biogenesis & NAD/carbohydrate metabolism",
    NB_M10 = "Mitochondrial OxPhos core subunits",
    NB_M11 = "Ribosome biogenesis & Protein folding"
)

if (!setequal(names(module_labels), all_modules)) {
    stop("module_labels does not cover exactly the detected modules.")
}

# Cross-check, if present.
manuscript_labels_csv <- here::here("hdWGCNA", "resources", "hdWGCNA_module_manuscript_labels.csv")
if (file.exists(manuscript_labels_csv)) {
    csv_labels <- read_csv(manuscript_labels_csv, show_col_types = FALSE)
    cmp <- tibble(module = names(module_labels), in_script = unname(module_labels)) %>%
        left_join(csv_labels %>% dplyr::select(module, in_csv = label), by = "module") %>%
        filter(is.na(in_csv) | in_script != in_csv)
    if (nrow(cmp) > 0) {
        for (i in seq_len(nrow(cmp))) {
            logger("WARNING: label mismatch for %s: script='%s' csv='%s'",
                   cmp$module[i], cmp$in_script[i],
                   ifelse(is.na(cmp$in_csv[i]), "<missing>", cmp$in_csv[i]))
        }
    } else {
        logger("Module labels match %s", manuscript_labels_csv)
    }
} else {
    logger("Note: %s not found; skipping label cross-check.", manuscript_labels_csv)
}

# Provenance table: auto-derived label alongside the curated label.
provenance_path <- file.path(csv_dir, "module_label_provenance.csv")
derived_labels %>%
    left_join(tibble(module = names(module_labels),
                     label_manuscript = unname(module_labels)),
              by = "module") %>%
    dplyr::select(module, label_auto, label_manuscript, source,
                  label_term_id, label_term_padj, top_bp, top_mf, top_cc) %>%
    write_csv(provenance_path)
logger("Module label provenance written: %s", provenance_path)

logger("Curated module labels:")
for (mod in all_modules) logger("  %s: %s", mod, module_labels[[mod]])

# --- 9. SUPP FIG 5D COMPOSITE PANEL ---
#
# Supp Fig 5D legend: "Gene Ontology Biological Process enrichment ... Tree plots
# show relationships among enriched GO terms". The composite renders
# the BP treeplot row only (MF/CC enrichment is still computed and saved in
# section 4).

panel_ontologies <- "BP"
n_ont            <- length(panel_ontologies)

# Saved treeplot layer indices (confirmed across tested plots).
CONS_TEXT_IDX <- 3L
CONS_BAR_IDX  <- 4L
HILIGHT_IDX   <- 5L

# Global colour / size scale ranges across the plotted ontologies.
all_padj  <- c()
all_count <- c()
for (mod in module_order) {
    for (ont in panel_ontologies) {
        fp <- file.path(rds_dir, paste0(prefix, "_", mod, "_", ont, "_simplified.rds"))
        if (!file.exists(fp)) next
        df <- as.data.frame(readRDS(fp))
        if (nrow(df) == 0) next
        all_padj  <- c(all_padj, df$p.adjust)
        all_count <- c(all_count, df$Count)
    }
}
if (length(all_padj) == 0 || length(all_count) == 0) {
    stop(sprintf(
        "No simplified %s enrichment objects found under %s; cannot build the composite panel.",
        paste(panel_ontologies, collapse = "/"), rds_dir
    ), call. = FALSE)
}
color_limits <- range(all_padj)
size_limits  <- range(all_count)

# ---- treeplot loading -------------------------------------------------------

build_treeplot_base <- function(ego_simple, ont, d_sim_list) {
    n_simple <- nrow(as.data.frame(ego_simple))
    if (n_simple < 2) return(NULL)
    ego_simple_sim <- tryCatch(
        pairwise_termsim(ego_simple, semData = d_sim_list[[ont]]),
        error = function(e) NULL
    )
    if (is.null(ego_simple_sim)) return(NULL)
    tryCatch(
        treeplot(
            ego_simple_sim,
            showCategory   = min(15, n_simple),
            cluster.params = list(n = min(5, n_simple)),
            fontsize       = 5
        ),
        error = function(e) NULL
    )
}

load_saved_treeplot <- function(mod, ont) {
    tree_fp <- file.path(rds_dir, paste0(prefix, "_", mod, "_", ont, "_treeplot.rds"))
    simp_fp <- file.path(rds_dir, paste0(prefix, "_", mod, "_", ont, "_simplified.rds"))
    if (file.exists(tree_fp)) return(readRDS(tree_fp))
    if (!file.exists(simp_fp)) return(NULL)
    ego_s <- readRDS(simp_fp)
    if (nrow(as.data.frame(ego_s)) < 2) return(NULL)
    logger("  rebuilding treeplot %s | %s (no saved object)", mod, ont)
    build_treeplot_base(ego_s, ont, d_sim_list)
}

logger("Loading saved treeplots for the composite panel...")
treeplot_cache <- setNames(vector("list", length(module_order)), module_order)
for (mod in module_order) {
    treeplot_cache[[mod]] <- setNames(vector("list", n_ont), panel_ontologies)
    for (ont in panel_ontologies) {
        treeplot_cache[[mod]][[ont]] <- load_saved_treeplot(mod, ont)
        logger("  %s | %s -> %s", mod, ont,
               if (is.null(treeplot_cache[[mod]][[ont]])) "missing" else "ok")
    }
}

# ---- consensus label shortening -------------------------------------------

consensus_term_overrides <- c(
    "generation of precursor metabolites and energy" = "Precursor metabolites",
    "ligase activity" = "Ligase",
    "catalytic activity" = "Catalytic",
    "monoatomic ion channel activity" = "Ion channel",
    "inorganic cation transmembrane transporter activity" = "Cation transporter",
    "proton-transporting ATP synthase complex, rotational mechanism" = "ATP synthase",
    "mitochondrion" = "Mitochondrion",
    "proton-transporting ATP synthase complex, coupling factor F(o)" = "ATP synthase F(o)",
    "cytochrome complex" = "Cytochrome complex",
    "mitochondrial membrane" = "Mito. membrane",
    "proton-transporting two-sector ATPase complex" = "Proton ATPase",
    "catabolic process" = "Protein catabolism",
    "protein-containing complex organization" = "Complex organisation",
    "nucleosome assembly" = "Nucleosome assembly",
    "protein-DNA complex organization" = "Protein-DNA complex",
    "proteasome-mediated ubiquitin-dependent protein catabolic process" = "Proteasomal degradation",
    "threonine-type endopeptidase activity" = "Endopeptidase activity",
    "threonine-type peptidase activity" = "Thr-type peptidase",
    "structural constituent of chromatin" = "Chromatin constituent",
    "protein-containing complex binding" = "Complex binding",
    "catalytic complex" = "Catalytic complex",
    "chromosome" = "Chromosome",
    "NURF complex" = "NURF complex",
    "U4/U5 snRNP" = "U4/U5 snRNP",
    "proteasome complex" = "Proteasome core",
    "ribosome biogenesis" = "Ribosome biogenesis",
    "ribonucleoprotein complex biogenesis" = "RNP biogenesis",
    "ribosomal small subunit biogenesis" = "SSU biogenesis",
    "ribosomal large subunit biogenesis" = "LSU biogenesis",
    "cytoplasmic translation" = "Cytoplasmic translation",
    "rRNA binding" = "rRNA binding",
    "large ribosomal subunit rRNA binding" = "LSU rRNA binding",
    "DNA-(apurinic or apyrimidinic site) endonuclease activity" = "AP endonuclease",
    "structural constituent of ribosome" = "Ribosome constituent",
    "RNA binding" = "RNA binding",
    "cytosolic small ribosomal subunit" = "Cytosolic SSU",
    "small ribosomal subunit" = "Small ribosomal subunit",
    "small-subunit processome" = "SSU processome",
    "cytosolic ribosome" = "Cytosolic ribosome",
    "cytosolic large ribosomal subunit" = "Cytosolic LSU",
    "chemical synaptic transmission" = "Synaptic transmission",
    "negative regulation of synaptic transmission" = "Synaptic inhibition",
    "protein folding" = "Protein folding",
    "rRNA processing" = "rRNA processing",
    "ribonucleoprotein complex assembly" = "RNP complex",
    "unfolded protein binding" = "Unfolded protein binding",
    "nucleic acid binding" = "Nucleic acid binding",
    "protein folding chaperone" = "Chaperone",
    "preribosome" = "Preribosome",
    "preribosome, large subunit precursor" = "Pre-LSU",
    "mitochondrial outer membrane translocase complex" = "TOM complex",
    "prefoldin complex" = "Prefoldin complex",
    "organelle lumen" = "Organelle lumen"
)

normalize_override_match <- function(x) {
    x %>%
        str_to_lower() %>%
        str_replace_all("[-/()]", " ") %>%
        str_replace_all("[[:punct:]]", " ") %>%
        str_replace_all("\\b(of|and|the|a|an|to|for|via|by)\\b", " ") %>%
        str_replace_all("\\s{2,}", " ") %>%
        str_trim()
}

tokenise_override_match <- function(x) {
    x %>%
        normalize_override_match() %>%
        str_split("\\s+", simplify = FALSE) %>%
        purrr::pluck(1) %>%
        unique()
}

shorten_consensus_label <- function(labels) {
    cleaned <- labels %>%
        str_replace_all("\n", " ") %>%
        str_replace_all("\\.{3,}|\u2026", " ") %>%
        str_replace_all("\\s{2,}", " ") %>%
        str_trim()

    override_keys <- names(consensus_term_overrides) %>%
        str_replace_all("\n", " ") %>%
        str_replace_all("\\s{2,}", " ") %>%
        str_trim()

    override_norm   <- purrr::map_chr(override_keys, normalize_override_match)
    override_tokens <- purrr::map(override_keys, tokenise_override_match)

    map_chr(cleaned, function(lbl) {
        if (is.na(lbl) || lbl == "") return("")

        lbl_norm <- normalize_override_match(lbl)
        idx <- which(override_norm == lbl_norm)
        if (length(idx) > 0) return(consensus_term_overrides[idx[1]])

        lbl_tokens <- tokenise_override_match(lbl)
        if (length(lbl_tokens) > 1) {
            overlap_scores <- purrr::map_dbl(override_tokens, function(tok) {
                if (length(tok) <= 1) return(0)
                length(intersect(lbl_tokens, tok)) / length(tok)
            })
            best_idx <- which.max(overlap_scores)
            if (length(best_idx) == 1 && is.finite(overlap_scores[best_idx]) &&
                overlap_scores[best_idx] >= 0.75) {
                return(consensus_term_overrides[best_idx])
            }
        }

        s <- lbl %>%
            str_remove_all("\\bprocess\\b") %>%
            str_remove_all("\\bbiological\\b") %>%
            str_remove_all("\\bcellular\\b") %>%
            str_remove_all("\\bactivity\\b") %>%
            str_remove_all("\\bregulation\\b") %>%
            str_remove_all("\\bpositive\\b") %>%
            str_remove_all("\\bnegative\\b") %>%
            str_remove_all("\\borganization\\b") %>%
            str_remove_all("\\borganelle\\b") %>%
            str_replace_all("\\.{3,}|\u2026", " ") %>%
            str_replace_all("\\s{2,}", " ") %>%
            str_trim() %>%
            str_to_sentence()

        if (is.na(s) || s == "") s <- str_sub(lbl, 1, 18)
        str_wrap(s, width = 14)
    })
}

sanitize_consensus_label_field <- function(labels) {
    labels %>%
        str_replace_all("[\r\n\t]", " ") %>%
        str_replace_all("\\s{2,}", " ") %>%
        str_trim()
}

format_consensus_plot_label <- function(labels, width = 11) {
    labels %>%
        sanitize_consensus_label_field() %>%
        str_wrap(width = width)
}

normalize_consensus_label <- function(labels) {
    labels %>%
        str_to_lower() %>%
        str_replace_all("\n", " ") %>%
        str_replace_all("\\.{3,}|\u2026", " ") %>%
        str_replace_all("\\s{2,}", " ") %>%
        str_trim()
}

# ---- treeplot layer helpers ---------------------------------------------

get_consensus_text_data    <- function(p) if (is.null(p)) NULL else { d <- p$layers[[CONS_TEXT_IDX]]$data; if (is.data.frame(d)) d else NULL }
get_consensus_bar_data     <- function(p) if (is.null(p)) NULL else { d <- p$layers[[CONS_BAR_IDX]]$data;  if (is.data.frame(d)) d else NULL }
get_consensus_hilight_data <- function(p) if (is.null(p)) NULL else { d <- p$layers[[HILIGHT_IDX]]$data;   if (is.data.frame(d)) d else NULL }

get_tree_x_range <- function(p) {
    if (is.null(p)) return(tibble(xmin = NA_real_, xmax = NA_real_))
    x_vals <- numeric(0)
    if (is.data.frame(p$data) && "x" %in% names(p$data)) x_vals <- c(x_vals, p$data$x)
    hil <- get_consensus_hilight_data(p)
    if (!is.null(hil)) {
        if ("xmin" %in% names(hil)) x_vals <- c(x_vals, hil$xmin)
        if ("xmax" %in% names(hil)) x_vals <- c(x_vals, hil$xmax)
    }
    x_vals <- x_vals[is.finite(x_vals)]
    if (length(x_vals) == 0) return(tibble(xmin = NA_real_, xmax = NA_real_))
    tibble(xmin = min(x_vals), xmax = max(x_vals))
}

get_tip_xmax <- function(p) {
    if (is.null(p) || !is.data.frame(p$data)) return(NA_real_)
    d <- p$data
    if (!all(c("x", "isTip") %in% names(d))) return(NA_real_)
    tip_x <- d$x[!is.na(d$isTip) & d$isTip]
    tip_x <- tip_x[is.finite(tip_x)]
    if (length(tip_x) == 0) return(NA_real_)
    max(tip_x)
}

shift_x_columns <- function(d, shift) {
    if (!is.data.frame(d) || !is.finite(shift) || shift == 0) return(d)
    for (col in intersect(c("x", "xend", "xmin", "xmax", "xintercept"), names(d))) {
        d[[col]] <- d[[col]] + shift
    }
    d
}

shift_y_columns <- function(d, shift) {
    if (!is.data.frame(d) || !is.finite(shift) || shift == 0) return(d)
    for (col in intersect(c("y", "yend", "ymin", "ymax", "yintercept"), names(d))) {
        d[[col]] <- d[[col]] + shift
    }
    d
}

align_treeplot_tips <- function(p, tip_target_x) {
    if (is.null(p) || !is.finite(tip_target_x)) return(p)
    tip_xmax <- get_tip_xmax(p)
    if (!is.finite(tip_xmax)) return(p)
    shift <- tip_target_x - tip_xmax
    if (!is.finite(shift) || shift == 0) return(p)
    p2 <- p
    p2$data <- shift_x_columns(p2$data, shift)
    for (i in seq_along(p2$layers)) {
        p2$layers[[i]]$data <- shift_x_columns(p2$layers[[i]]$data, shift)
    }
    p2
}

standardize_treeplot_layers <- function(p, tip_anchor_x,
                                        highlight_gap = 0.04, highlight_width = 0.78,
                                        tree_linewidth = 0.42, bar_linewidth = 0.7,
                                        tip_point_size = 1.1, highlight_alpha = 0.28) {
    if (is.null(p)) return(p)
    p2 <- p
    for (i in intersect(c(1L, 2L), seq_along(p2$layers))) {
        p2$layers[[i]]$aes_params$linewidth <- tree_linewidth
        p2$layers[[i]]$computed_geom_params$linewidth <- tree_linewidth
        p2$layers[[i]]$aes_params$lineend <- "round"
        p2$layers[[i]]$computed_geom_params$lineend <- "round"
    }
    if (CONS_BAR_IDX <= length(p2$layers)) {
        p2$layers[[CONS_BAR_IDX]]$aes_params$linewidth <- bar_linewidth
        p2$layers[[CONS_BAR_IDX]]$computed_geom_params$linewidth <- bar_linewidth
        p2$layers[[CONS_BAR_IDX]]$aes_params$lineend <- "round"
        p2$layers[[CONS_BAR_IDX]]$computed_geom_params$lineend <- "round"
    }
    if (HILIGHT_IDX <= length(p2$layers)) {
        hil <- get_consensus_hilight_data(p2)
        if (!is.null(hil) && all(c("xmin", "xmax") %in% names(hil))) {
            hil_right <- tip_anchor_x - highlight_gap
            hil$xmax <- hil_right
            hil$xmin <- hil_right - highlight_width
            p2$layers[[HILIGHT_IDX]]$data <- hil
        }
        p2$layers[[HILIGHT_IDX]]$aes_params$alpha <- highlight_alpha
        p2$layers[[HILIGHT_IDX]]$aes_params$linewidth <- 0
        p2$layers[[HILIGHT_IDX]]$computed_geom_params$alpha <- highlight_alpha
        p2$layers[[HILIGHT_IDX]]$computed_geom_params$linewidth <- 0
        p2$layers[[HILIGHT_IDX]]$aes_params$linejoin <- "round"
        p2$layers[[HILIGHT_IDX]]$computed_geom_params$linejoin <- "round"
    }
    if (6L <= length(p2$layers)) {
        p2$layers[[6]]$aes_params$size <- tip_point_size
        p2$layers[[6]]$aes_params$alpha <- 1
    }
    p2
}

set_consensus_bar_layer <- function(p, d) {
    bar_mapping <- aes(
        x = x, y = y, xend = xend, yend = yend,
        colour_ggnewscale_2 = color.x
    )
    p$layers[[CONS_BAR_IDX]]$data            <- d
    p$layers[[CONS_BAR_IDX]]$mapping         <- bar_mapping
    p$layers[[CONS_BAR_IDX]]$computed_mapping <- bar_mapping
    p$layers[[CONS_BAR_IDX]]$inherit.aes     <- FALSE
    p
}

position_consensus_bars <- function(p, bar_offset = 0.18, merge_bars = TRUE,
                                    anchor_x = NULL, y_shift = -0.12) {
    p2  <- p
    seg <- get_consensus_bar_data(p2)
    hil <- get_consensus_hilight_data(p2)
    if (is.null(seg) || is.null(hil)) return(p2)

    target_x <- if (is.null(anchor_x)) max(hil$xmax, na.rm = TRUE) + bar_offset else anchor_x

    d <- seg %>%
        mutate(
            label_raw = if ("labels" %in% names(.)) labels else NA_character_,
            y_min = pmin(y, yend),
            y_max = pmax(y, yend)
        )

    if (merge_bars && "label_raw" %in% names(d)) {
        d <- d %>%
            filter(!is.na(label_raw), label_raw != "") %>%
            group_by(label_raw) %>%
            summarise(
                x = target_x, xend = target_x,
                y = min(y_min, na.rm = TRUE), yend = max(y_max, na.rm = TRUE),
                angle = angle[1], labels = label_raw[1], cluster = cluster[1],
                color.x = color.x[1], node = node[1], parent = parent[1],
                branch.length = branch.length[1], isTip = isTip[1], branch = branch[1],
                group = group[1], color.y = color.y[1], count = count[1],
                .groups = "drop"
            )
    } else {
        d$x    <- target_x
        d$xend <- target_x
    }

    d <- shift_y_columns(d, y_shift)
    set_consensus_bar_layer(p2, d)
}

build_label_plan <- function(p, mod, ont, label_sot = NULL) {
    seg <- get_consensus_bar_data(p)
    if (is.null(seg)) return(tibble())

    raw_labels <- if ("labels" %in% names(seg)) seg$labels
                  else if ("label" %in% names(seg)) seg$label
                  else NULL
    if (is.null(raw_labels)) return(tibble())

    d <- seg %>%
        mutate(
            label_raw = raw_labels,
            cluster   = if ("cluster" %in% names(.)) as.integer(cluster) else NA_integer_,
            label_key = normalize_consensus_label(label_raw),
            x_bar = pmax(x, xend),
            y_min = pmin(y, yend),
            y_max = pmax(y, yend)
        ) %>%
        filter(!is.na(label_raw), label_raw != "")
    if (nrow(d) == 0) return(tibble())

    combined <- d %>%
        group_by(cluster, label_raw, label_key) %>%
        summarise(
            x_bar = mean(x_bar, na.rm = TRUE),
            y_min = min(y_min, na.rm = TRUE),
            y_max = max(y_max, na.rm = TRUE),
            y_mid = (min(y_min, na.rm = TRUE) + max(y_max, na.rm = TRUE)) / 2,
            .groups = "drop"
        ) %>%
        arrange(desc(y_mid)) %>%
        mutate(module = mod, ont = ont)

    if (!is.null(label_sot) && nrow(label_sot) > 0) {
        combined <- combined %>%
            left_join(
                label_sot %>% dplyr::select(module, ont, cluster, label_resolved),
                by = c("module", "ont", "cluster")
            )
    } else {
        combined$label_resolved <- NA_character_
    }

    combined %>%
        mutate(
            label_final = dplyr::coalesce(label_resolved, shorten_consensus_label(label_raw)),
            label_final = format_consensus_plot_label(label_final),
            is_visible  = label_final != ""
        )
}

apply_label_plan <- function(p, plan,
                             label_gap = 0.03, label_angle = 90,
                             label_size = 1.5, label_colour = "grey15",
                             label_hjust = 0.5, label_vjust = 1) {
    p2 <- p

    if (CONS_TEXT_IDX <= length(p2$layers)) {
        d <- p2$layers[[CONS_TEXT_IDX]]$data
        if (is.data.frame(d)) {
            if ("label"  %in% names(d)) p2$layers[[CONS_TEXT_IDX]]$data$label  <- rep("", nrow(d))
            if ("labels" %in% names(d)) p2$layers[[CONS_TEXT_IDX]]$data$labels <- rep("", nrow(d))
        }
    }

    protected_layers <- c(CONS_TEXT_IDX, CONS_BAR_IDX, HILIGHT_IDX)
    for (i in seq_along(p2$layers)) {
        if (i %in% protected_layers) next
        cls <- class(p2$layers[[i]]$geom)[1]
        d   <- p2$layers[[i]]$data
        if (grepl("Text|Label", cls, ignore.case = TRUE)) {
            if (is.data.frame(d)) {
                if ("label"  %in% names(d)) p2$layers[[i]]$data$label  <- rep("", nrow(d))
                if ("labels" %in% names(d)) p2$layers[[i]]$data$labels <- rep("", nrow(d))
            }
            p2$layers[[i]]$aes_params$colour <- "transparent"
            p2$layers[[i]]$aes_params$size   <- 0
        }
    }

    label_df <- plan %>%
        filter(is_visible, label_final != "") %>%
        transmute(x = x_bar + label_gap, y = y_mid, label = label_final)

    if (nrow(label_df) > 1) {
        label_df <- label_df %>%
            arrange(desc(y)) %>%
            mutate(n_lines = stringr::str_count(label, "\n") + 1L)
        min_gap <- 0.78
        y_adj <- label_df$y
        for (i in 2:length(y_adj)) {
            needed_gap <- min_gap + 0.20 * max(
                label_df$n_lines[i - 1] - 1L, label_df$n_lines[i] - 1L, 0L
            )
            y_adj[i] <- min(y_adj[i], y_adj[i - 1] - needed_gap)
        }
        label_df$y <- y_adj
        label_df <- label_df %>% select(-n_lines)
    }

    if (nrow(label_df) == 0) return(p2)

    p2 +
        geom_text(
            data = label_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, angle = label_angle, size = label_size,
            colour = label_colour, hjust = label_hjust, vjust = label_vjust,
            lineheight = 0.85
        )
}

# ---- global tree x alignment -------------------------------------------

treeplot_tip_xmax <- map_dfr(module_order, function(mod) {
    map_dfr(panel_ontologies, function(ont) {
        tibble(module = mod, ont = ont, tip_xmax = get_tip_xmax(treeplot_cache[[mod]][[ont]]))
    })
}) %>% pull(tip_xmax)
treeplot_tip_xmax <- treeplot_tip_xmax[is.finite(treeplot_tip_xmax)]
global_tip_xmax   <- if (length(treeplot_tip_xmax)) max(treeplot_tip_xmax) else NA_real_

treeplot_aligned_x_ranges <- map_dfr(module_order, function(mod) {
    map_dfr(panel_ontologies, function(ont) {
        p <- treeplot_cache[[mod]][[ont]]
        xr <- get_tree_x_range(p)
        tip_x <- get_tip_xmax(p)
        shift <- if (is.finite(tip_x)) global_tip_xmax - tip_x else NA_real_
        xr %>% mutate(module = mod, ont = ont,
                      xmin_aligned = xmin + shift, xmax_aligned = xmax + shift)
    })
})
global_tree_xmin <- min(treeplot_aligned_x_ranges$xmin_aligned, na.rm = TRUE)
global_tree_xmax <- max(treeplot_aligned_x_ranges$xmax_aligned, na.rm = TRUE)

# ---- final per-cell treeplot builder ---------------------------------

font_scale <- 1.5
fs <- function(size) size * font_scale

make_empty_plot <- function(msg = "\u2014") {
    ggplot() +
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = msg,
                 color = "grey65", size = 2.5, hjust = 0.5, fontface = "italic") +
        xlim(0, 1) + ylim(0, 1) +
        theme(plot.background = element_rect(fill = "#fafafa", color = "grey90", linewidth = 0.2))
}

make_treeplot_final <- function(mod, ont, bar_offset = 0.60, label_gap = 1, bar_y_shift = -0.12) {
    simp_fp <- file.path(rds_dir, paste0(prefix, "_", mod, "_", ont, "_simplified.rds"))
    if (!file.exists(simp_fp)) return(make_empty_plot("No data"))

    ego_s   <- readRDS(simp_fp)
    n_terms <- nrow(as.data.frame(ego_s))
    if (n_terms == 0) return(make_empty_plot("No terms"))
    if (n_terms < 2)  return(make_empty_plot(as.data.frame(ego_s)$Description[1]))

    p_base <- treeplot_cache[[mod]][[ont]]
    if (is.null(p_base)) return(make_empty_plot("Plot failed"))

    p_aligned <- align_treeplot_tips(p_base, global_tip_xmax)
    p_styled  <- standardize_treeplot_layers(p_aligned, tip_anchor_x = global_tip_xmax)

    bar_x   <- global_tip_xmax + bar_offset
    x_upper <- bar_x + label_gap + 0.35

    p_shifted <- position_consensus_bars(
        p_styled, bar_offset = bar_offset, merge_bars = TRUE,
        anchor_x = bar_x, y_shift = bar_y_shift
    )
    plan <- build_label_plan(p_shifted, mod = mod, ont = ont,
                             label_sot = consensus_label_resolved)

    p_clean <- apply_label_plan(
        p_shifted, plan,
        label_gap = label_gap, label_angle = 90,
        label_size = fs(2.5), label_colour = "grey15",
        label_hjust = 0.5, label_vjust = 1
    )

    p_clean +
        scale_color_gradient(
            low = "#E31A1C", high = "#1F78B4",
            limits = color_limits, oob = scales::squish, name = "adj. P"
        ) +
        scale_size_continuous(limits = size_limits, range = c(0.6, 3.0), name = "Gene count") +
        coord_cartesian(xlim = c(global_tree_xmin, x_upper), clip = "off") +
        theme_void() +
        theme(
            legend.position = "none",
            plot.margin     = margin(2, 14, 2, 6),
            plot.background = element_rect(fill = "white", color = NA)
        )
}

logger("Building final treeplots...")
plot_grid_list <- setNames(vector("list", length(module_order)), module_order)
for (mod in module_order) {
    plot_grid_list[[mod]] <- setNames(vector("list", n_ont), panel_ontologies)
    for (ont in panel_ontologies) {
        logger("  %s | %s", mod, ont)
        plot_grid_list[[mod]][[ont]] <- make_treeplot_final(mod, ont)
    }
}

# ---- shared legend -----------------------------------------------------

color_breaks <- pretty(color_limits, n = 4)
color_breaks <- color_breaks[color_breaks >= color_limits[1] & color_breaks <= color_limits[2]]
count_breaks <- pretty(size_limits, n = 4)
count_breaks <- count_breaks[count_breaks >= size_limits[1] & count_breaks <= size_limits[2]]

legend_proxy <- ggplot(
    data.frame(
        padj  = seq(color_limits[1], color_limits[2], length.out = 50),
        count = seq(size_limits[1], size_limits[2], length.out = 50),
        y = 1
    ),
    aes(x = padj, y = y, color = padj, size = count)
) +
    geom_point(alpha = 0) +
    scale_color_gradient(
        low = "#E31A1C", high = "#1F78B4",
        limits = color_limits, breaks = color_breaks,
        labels = function(x) formatC(x, format = "f", digits = 3),
        name = "FDR",
        guide = guide_colorbar(
            direction = "horizontal", title.position = "top", title.hjust = 0.5,
            barwidth = unit(12, "lines"), barheight = unit(0.7, "lines"),
            ticks.colour = "grey30", frame.colour = "grey50"
        )
    ) +
    scale_size_continuous(
        limits = size_limits, range = c(1.5, 5),
        breaks = count_breaks, labels = as.character(count_breaks),
        name = "Gene count",
        guide = guide_legend(
            direction = "horizontal", nrow = 1, title.position = "top",
            title.hjust = 0.5, override.aes = list(color = "grey40", alpha = 1)
        )
    ) +
    theme_void() +
    theme(
        legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.spacing.x = unit(40, "pt"),
        legend.spacing.y = unit(2, "pt"),
        legend.title     = element_text(size = fs(10), face = "bold", color = "grey15", hjust = 0.5),
        legend.text      = element_text(size = fs(9.5), color = "grey20"),
        legend.key.size  = unit(10, "pt")
    )

shared_legend <- cowplot::get_legend(legend_proxy)

# ---- hub-gene kME strip ---------------------------------------------

# Shared module palette (D3 category20, grey tones excluded; NB_M9 = the red
# used for Ldh/glycolysis in the manuscript). Matches fig2e_module_colors in
# reports/gene_network_analysis.qmd.
mod_colors <- c(
    NB_M1  = "#1f77b4", NB_M2  = "#ff7f0e", NB_M3  = "#2ca02c",
    NB_M4  = "#9467bd", NB_M5  = "#8c564b", NB_M6  = "#e377c2",
    NB_M7  = "#bcbd22", NB_M8  = "#17becf", NB_M9  = "#e31a1c",
    NB_M10 = "#aec7e8", NB_M11 = "#98df8a"
)

hub_data <- nb_modules %>%
    dplyr::filter(module != "grey") %>%
    dplyr::rename(gene_symbol = dplyr::any_of(c("gene_name", "gene", "symbol"))) %>%
    dplyr::select(gene_symbol, module, tidyselect::starts_with("kME_")) %>%
    tidyr::pivot_longer(
        cols = tidyselect::starts_with("kME_"),
        names_to = "kME_mod", values_to = "kME"
    ) %>%
    dplyr::mutate(
        kME_mod = gsub("kME_", "", kME_mod),
        module  = factor(module, levels = module_order)
    ) %>%
    dplyr::filter(module == kME_mod) %>%
    dplyr::group_by(module) %>%
    dplyr::arrange(dplyr::desc(kME), .by_group = TRUE) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup()

make_kme_module_plot <- function(module_id, show_y = FALSE) {
    mod_df <- hub_data %>%
        dplyr::filter(module == module_id) %>%
        dplyr::mutate(gene_symbol = factor(
            gene_symbol, levels = gene_symbol[order(kME, decreasing = TRUE)]
        ))

    ggplot(mod_df, aes(x = gene_symbol, y = kME)) +
        geom_segment(aes(xend = gene_symbol, y = 0, yend = kME),
                     color = mod_colors[[module_id]], linewidth = 1.0, alpha = 0.6) +
        geom_point(fill = mod_colors[[module_id]], shape = 21, size = 3.0,
                   color = "white", stroke = 0.7) +
        geom_text(aes(label = sprintf("%.2f", kME)),
                  angle = 90, hjust = -0.45, size = fs(2.8), color = "black") +
        scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.5),
                           expand = c(0, 0), position = "right") +
        coord_cartesian(clip = "off") +
        labs(x = NULL, y = if (show_y) "kME" else NULL) +
        theme_minimal(base_size = fs(11.5)) +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       size = fs(6.2), color = "black", margin = margin(t = 1)),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = if (show_y) element_text(angle = 90, size = fs(10), color = "black") else element_blank(),
            axis.ticks.y = if (show_y) element_line(color = "black", linewidth = 0.4) else element_blank(),
            axis.title.y.right = if (show_y) element_text(angle = 90, vjust = 0.8, size = fs(11)) else element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.y = if (show_y) element_line(color = "black", linewidth = 0.4) else element_blank(),
            axis.line.x = element_line(color = "black", linewidth = 0.4),
            plot.margin = margin(2, 8, 2, 2)
        )
}

kme_axis_labels <- stats::setNames(
    vapply(module_order, function(mod) paste0(mod, "\n", module_labels[[mod]]), character(1)),
    module_order
)

kme_plot_list <- stats::setNames(
    lapply(seq_along(module_order), function(i) {
        make_kme_module_plot(module_order[i], show_y = i == length(module_order))
    }),
    module_order
)

kme_grob_list <- stats::setNames(lapply(kme_plot_list, ggplotGrob), names(kme_plot_list))
kme_common_heights <- do.call(grid::unit.pmax, lapply(kme_grob_list, function(g) g$heights))
kme_grob_list <- stats::setNames(
    lapply(kme_grob_list, function(g) { g$heights <- kme_common_heights; g }),
    names(kme_grob_list)
)

# ---- layout ----------------------------------------------------------

main_content_w <- 24
legend_strip_w <- 2.8
content_w      <- main_content_w + legend_strip_w
n_mods         <- length(module_order)
row_title_w    <- 1.1
tree_col_w     <- (main_content_w - row_title_w) / n_mods
kme_row_h      <- 6.8
tree_row_h     <- 12
content_h      <- kme_row_h + tree_row_h * n_ont

col_widths  <- unit(c(row_title_w, rep(tree_col_w, n_mods)), "inches")
row_heights <- unit(c(rep(tree_row_h, n_ont), kme_row_h), "inches")

ont_row_map <- stats::setNames(seq_along(panel_ontologies), panel_ontologies)
ROW_KME     <- n_ont + 1
ont_full_name <- c(BP = "GO: Biological Process",
                   MF = "GO: Molecular Function",
                   CC = "GO: Cellular Compartment")

draw_panel_content <- function() {
    pushViewport(viewport(
        x = unit(main_content_w / 2, "inches"), y = unit(content_h / 2, "inches"),
        width = unit(main_content_w, "inches"), height = unit(content_h, "inches"),
        just = c("center", "center")
    ))
    pushViewport(viewport(layout = grid.layout(
        nrow = n_ont + 1, ncol = n_mods + 1, widths = col_widths, heights = row_heights
    )))

    for (ont in panel_ontologies) {
        pushViewport(viewport(layout.pos.row = ont_row_map[[ont]], layout.pos.col = 1))
        grid.text(ont_full_name[[ont]], rot = 90,
                  gp = gpar(fontsize = fs(9.5), fontface = "bold", col = "grey35"))
        popViewport()
    }

    for (i in seq_along(module_order)) {
        mod   <- module_order[i]
        col_i <- i + 1

        pushViewport(viewport(layout.pos.row = ROW_KME, layout.pos.col = col_i))
        pushViewport(viewport(x = unit(0.8, "npc"), y = unit(0.79, "npc"),
                              width = unit(1, "npc"), height = unit(0.40, "npc"),
                              just = c("center", "center")))
        grid.draw(kme_grob_list[[mod]])
        popViewport()
        grid.text(kme_axis_labels[[mod]],
                  x = unit(0.74, "npc"), y = unit(0.56, "npc"),
                  just = c("center", "bottom"), hjust = 1, rot = 90,
                  gp = gpar(fontsize = fs(9.2), fontface = "bold", col = "black", lineheight = 0.9))
        popViewport()

        for (ont in panel_ontologies) {
            pushViewport(viewport(layout.pos.row = ont_row_map[[ont]], layout.pos.col = col_i))
            pushViewport(viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                                  width = unit(0.92, "npc"), height = unit(1, "npc"),
                                  just = c("center", "center")))
            p <- plot_grid_list[[mod]][[ont]]
            if (!is.null(p)) {
                if (inherits(p, "grob") || inherits(p, "gTree") || inherits(p, "gtable")) {
                    grid.draw(p)
                } else {
                    grid.draw(ggplotGrob(p))
                }
            }
            popViewport()
            popViewport()
        }
    }

    pushViewport(viewport(
        x = unit(main_content_w + legend_strip_w / 2, "inches"),
        y = unit(kme_row_h + (content_h - kme_row_h) / 2, "inches"),
        width = unit(legend_strip_w, "inches"),
        height = unit(content_h - kme_row_h, "inches"),
        just = c("center", "center"), clip = "off"
    ))
    pushViewport(viewport(x = unit(0.34, "npc"), y = unit(0.5, "npc"),
                          width = unit(5.8, "inches"), height = unit(1.1, "inches"),
                          just = c("center", "center"), clip = "off"))
    pushViewport(viewport(angle = 90))
    grid.draw(shared_legend)
    popViewport(); popViewport(); popViewport()

    popViewport(); popViewport()
}

# ---- render --------------------------------------------------------

png_land <- file.path(output_dir, "GO_panel_supplementary_landscape.png")
png(png_land, width = content_w, height = content_h, units = "in", res = 300, bg = "white")
grid.newpage()
draw_panel_content()
invisible(dev.off())
logger("Landscape PNG written: %s", png_land)

pdf_rot <- file.path(output_dir, "GO_panel_supplementary_rotated.pdf")
pdf(pdf_rot, width = content_h, height = content_w, bg = "white")
grid.newpage()
pushViewport(viewport(y = unit(0.48, "npc"), angle = -90,
                      width = unit(content_w, "inches"), height = unit(content_h, "inches")))
draw_panel_content()
popViewport()
invisible(dev.off())
logger("Rotated PDF written: %s", pdf_rot)

# --- 10. DONE ---

logger("=== GO ENRICHMENT COMPLETE ===")
logger("Outputs written under: %s", output_dir)
logger("  1. Enrichment objects : %s", rds_dir)
logger("  2. Diagnostic PDFs    : %s", pdf_dir)
logger("  3. Label review       : %s", csv_dir)
logger("  4. Supp Fig 5D panel  : %s", png_land)
