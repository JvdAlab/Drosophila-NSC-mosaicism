# hdWGCNA Co-Expression Network Analysis

## Overview

Neuroblast-focused co-expression network analysis (hdWGCNA) comparing ND-75
knockdown ("ND75-KD") vs. Control, built around the module containing *Ldh*
and its hub genes.

**Input:**
- `results/scRNA/scRNA_data_QC_filtered_w_scVI_latent_annotated_seurat.rds`
  (annotated Seurat object, produced by
  [`scRNA_pipeline/06_convert_adata_to_seurat.R`](../scRNA_pipeline/06_convert_adata_to_seurat.R)).
- [`resources/hdWGCNA_module_manuscript_labels.csv`](./resources/hdWGCNA_module_manuscript_labels.csv)
  (committed, hand-curated — used by stage 01 only as an optional cross-check
  against the frozen `module_labels` vector; see Analysis Details).

**Output:**
- The hdWGCNA network object plus derived tables, gene lists, and plots, all
  under `results/hdWGCNA/` (git-ignored, regenerated each run; see
  [Expected Output](#expected-output)).

---

## Workflow context

Run from the repository root.

### Environment

```bash
conda env create -f hdWGCNA/environment.yml
conda activate hdWGCNA_pipeline
```

### Pipeline commands

```bash
Rscript hdWGCNA/00_hdWGCNA_analysis.R
Rscript hdWGCNA/01_GO_enrichment.R
```

Stage 01 runs in a single pass. Module labels are not re-derived at run time:
the curated `NB_M*` → biological-process labels used in Figure 2E were assigned
from the top non-generic enriched terms and refined by hand against the BP/MF/CC
GO treeplot hierarchies, and are frozen in `01_GO_enrichment.R` (the
`module_labels` vector). The short per-cluster tip labels on the Supp Fig 5D
treeplots are likewise frozen, in the `consensus_term_overrides` map.

---

## Pipeline Scripts

All outputs are written under `results/hdWGCNA/`.

| Stage | Script | Input | Outputs |
|:---|:---|:---|:---|
| 00: Network construction | [`00_hdWGCNA_analysis.R`](./00_hdWGCNA_analysis.R) | Annotated Seurat `.rds` | hdWGCNA Seurat object, TOM matrix, module assignments/eigengenes, hub-gene lists, candidate-gene lists, module-trait correlation tables/plots, `analysis_summary.txt` (see [Expected Output](#expected-output)) |
| 01: GO enrichment + module labels | [`01_GO_enrichment.R`](./01_GO_enrichment.R) | `hdWGCNA_neuroblast_complete_analysis_with_condition.rds` (from stage 00) | Per module × ontology (BP/MF/CC): `GO_exploration_top50/GO_objects/*_{full,simplified,treeplot}.rds`, `GO_analysis/*.pdf` diagnostics, `label_review/{term_count_audit,module_label_provenance}.csv`, and the `GO_panel_supplementary_landscape.png` / `_rotated.pdf` composite (Supp Fig 5C–D: hub-gene kME strip + BP treeplot row). Top 50 hub genes/module by kME; ontology-specific cutoffs (BP `minGSSize = 5`, MF/CC `minGSSize = 3`). |

---

## Analysis Details

- **Feature filtering:** removes pseudogene/rRNA/ncRNA/transposable-element gene name patterns (`^Psi:`, `:CR`, `^rRNA`, `^sisRNA`, `^asRNA`, `^hpRNA`, `^FBti`) before network construction
- **Network construction:** gene selection `fraction = 0.05`; metacells `k = 15`, `max_shared = 5`; network type `signed`; `soft_power = 9` (fixed, chosen from `TestSoftPowers` output)
- **Network scope:** metacells grouped by cell type + condition, network constructed on the Neuroblasts group only, assay `RNA`
- **Module of interest:** after `ResetModuleNames`, modules are named `NB_M1`, `NB_M2`, …; the module containing *Ldh* is then identified (`ldh_module`) and drives the Ldh-focused gene lists, hub lists, and plots
- **Condition variable:** binary `ND75_KD` derived from `condition == "ND75-KD"` — an ND-75 knockdown vs. Control comparison
- **Trait correlation:** calculated within Neuroblasts; genes from modules with FDR < 0.05 are written to a combined candidate-gene list (`ND75_KD_associated_genes_for_SCENIC_with_condition.txt`)
- **GO enrichment (stage 01):** top 50 hub genes/module by kME; background universe = all non-grey module genes; `clusterProfiler::enrichGO` with `pAdjustMethod = "BH"`, `maxGSSize = 500`; ontology-specific thresholds — BP `p < 0.05, q < 0.20, minGSSize = 5`, MF/CC `p < 0.10, q < 0.30, minGSSize = 3`; redundant terms collapsed with Wang semantic similarity (`clusterProfiler::simplify` cutoff 0.85 BP, 0.90 MF/CC). BP/MF/CC are all computed and saved (`GO_objects/*_simplified.rds`); the composite panel's treeplot row (Supp Fig 5D) renders **BP** only (Supp Fig 5C).
- **Module labels:** `NB_M10` is labelled from hub-gene identity (mtDNA-encoded); the other ten `NB_M*` labels are frozen in the `module_labels` vector in `01_GO_enrichment.R` (assigned from the top non-generic enriched terms, refined by hand against the GO treeplot hierarchies) and cross-checked against [`resources/hdWGCNA_module_manuscript_labels.csv`](./resources/hdWGCNA_module_manuscript_labels.csv) (`label_source = figure_2E`). Per-cluster treeplot tip labels are shortened by the frozen `consensus_term_overrides` map.

---

## Dependencies

The pipeline dependencies are defined in [`environment.yml`](./environment.yml).

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | 4.4.3 | [`00_hdWGCNA_analysis.R`](./00_hdWGCNA_analysis.R), [`01_GO_enrichment.R`](./01_GO_enrichment.R) |
| [Seurat](https://satijalab.org/seurat/) | 5.4.0 | Single-cell object handling |
| [hdWGCNA](https://smorabit.github.io/hdWGCNA/) | 0.4.11 | Co-expression network construction |
| [WGCNA](https://bioconductor.org/packages/WGCNA/) | 1.74 | Underlying network/module detection |
| [tidyverse](https://www.tidyverse.org/) | 2.0.0 | Data manipulation |
| [ggplot2](https://ggplot2.tidyverse.org/) | 4.0.2 | Diagnostic and correlation plots |
| [patchwork](https://patchwork.data-imaginist.com/) | 1.3.2 | Plot composition |
| [here](https://cran.r-project.org/package=here) | 1.0.2 | Repository-relative paths |
| [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) | 4.14.6 | GO over-representation (`01_GO_enrichment.R`) |
| [enrichplot](https://bioconductor.org/packages/enrichplot/) / [GOSemSim](https://bioconductor.org/packages/GOSemSim/) | 1.26.6 / 2.32.0 | GO tree plots + semantic similarity (`01`) |
| [org.Dm.eg.db](https://bioconductor.org/packages/org.Dm.eg.db/) | 3.20.0 | *D. melanogaster* GO annotation (`01`) |
| [cowplot](https://cran.r-project.org/package=cowplot) | 1.2.0 | Shared-legend extraction for the composite panel (`01`) |


---

## Expected Output

`results/hdWGCNA/`:

**Objects**
- `hdWGCNA_neuroblast_complete_analysis_with_condition.rds` — complete Seurat object
- `NB_Network_TOM.rda` — the TOM written by `ConstructNetwork`
- `hdWGCNA_module_assignments_with_condition.{rds,csv}`
- `hdWGCNA_module_eigengenes_with_condition.rds`
- `hdWGCNA_TOM_matrix_with_condition.rds`
- `hdWGCNA_module_gene_lists_with_condition.rds`
- `hdWGCNA_hub_genes_with_condition.rds`

**Gene lists**
- `<module>_module_genes_for_SCENIC_with_condition.csv`, `<module>_module_genes_list_with_condition.txt`
- `<module>_hub_genes_with_condition.txt`
- `ND75_KD_associated_genes_for_SCENIC_with_condition.txt`

**Correlation tables**
- `module_trait_correlations_with_condition.rds`

**Plots**
- `soft_power_diagnostics.pdf` — `TestSoftPowers` table (evidence for `soft_power = 9`)
- `NB_dendrogram.pdf` — co-expression module dendrogram
- `<module>_module_hubs.pdf`
- `module_ND75_KD_correlation.pdf`, `module_ND75_KD_correlation_combined.pdf`

**Summary**
- `analysis_summary.txt`, `analysis_summary_with_condition.rds`

### `results/hdWGCNA/GO_exploration_top50/` (stage 01)

- `GO_objects/NB_hdWGCNA_<module>_<BP|MF|CC>_{full,simplified,treeplot}.rds` — enrichment objects and cached base treeplots
- `GO_analysis/*.pdf` — per-module GO diagnostic pages
- `label_review/GO_top_terms_per_module.csv` (top 5 simplified terms per module × ontology), `label_review/term_count_audit.csv` (per module × ontology term counts), and `label_review/module_label_provenance.csv` (auto-derived label alongside the curated manuscript label)
- `GO_panel_supplementary_landscape.png`, `GO_panel_supplementary_rotated.pdf` — the kME + BP-treeplot composite (Supp Fig 5C–D)
