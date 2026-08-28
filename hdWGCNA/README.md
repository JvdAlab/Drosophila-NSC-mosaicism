# hdWGCNA Co-Expression Network Analysis

## Overview

Neuroblast-focused co-expression network analysis (hdWGCNA) comparing ND-75
knockdown ("ND75-KD") vs. Control, built around the module containing *Ldh*
and its hub genes.

**Input:**
- `results/scRNA/scRNA_data_QC_filtered_w_scVI_latent_annotated_seurat.rds`
  (annotated Seurat object, produced by
  [`scRNA_pipeline/06_convert_adata_to_seurat.R`](../scRNA_pipeline/06_convert_adata_to_seurat.R)).

**Output:**
- The hdWGCNA network object plus derived tables, gene lists, and plots, all
  under `results/hdWGCNA/` (see
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
```

---

## Pipeline Scripts

All outputs are written under `results/hdWGCNA/`.

| Stage | Script | Input | Outputs |
|:---|:---|:---|:---|
| 00: Network construction | [`00_hdWGCNA_analysis.R`](./00_hdWGCNA_analysis.R) | Annotated Seurat `.rds` | hdWGCNA Seurat object, TOM matrix, module assignments/eigengenes, hub-gene lists, candidate-gene lists, module-trait correlation tables/plots, `analysis_summary.txt` (see [Expected Output](#expected-output)) |

---

## Analysis Details

- **Feature filtering:** removes pseudogene/rRNA/ncRNA/transposable-element gene name patterns (`^Psi:`, `:CR`, `^rRNA`, `^sisRNA`, `^asRNA`, `^hpRNA`, `^FBti`) before network construction
- **Network construction:** gene selection `fraction = 0.05`; metacells `k = 15`, `max_shared = 5`; network type `signed`; `soft_power = 9` (fixed, chosen from `TestSoftPowers` output)
- **Network scope:** metacells grouped by cell type + condition, network constructed on the Neuroblasts group only, assay `RNA`
- **Module of interest:** after `ResetModuleNames`, modules are named `NB_M1`, `NB_M2`, …; the module containing *Ldh* is then identified programmatically (`ldh_module`) and drives the Ldh-focused gene lists, hub lists, and plots
- **Condition variable:** binary `ND75_KD` derived from `condition == "ND75-KD"` — an ND-75 knockdown vs. Control comparison
- **Trait correlation:** calculated within Neuroblasts; genes from modules with FDR < 0.05 are written to a combined candidate-gene list (`ND75_KD_associated_genes_for_SCENIC_with_condition.txt`)

---

## Dependencies

The pipeline dependencies are defined in [`environment.yml`](./environment.yml).

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | 4.4.3 | [`00_hdWGCNA_analysis.R`](./00_hdWGCNA_analysis.R) |
| [Seurat](https://satijalab.org/seurat/) | 5.4.0 | Single-cell object handling |
| [hdWGCNA](https://smorabit.github.io/hdWGCNA/) | 0.4.11 | Co-expression network construction |
| [WGCNA](https://bioconductor.org/packages/WGCNA/) | 1.74 | Underlying network/module detection |
| [tidyverse](https://www.tidyverse.org/) | 2.0.0 | Data manipulation |
| [ggplot2](https://ggplot2.tidyverse.org/) | 4.0.2 | Diagnostic and correlation plots |
| [patchwork](https://patchwork.data-imaginist.com/) | 1.3.2 | Plot composition |
| [here](https://cran.r-project.org/package=here) | 1.0.2 | Repository-relative paths |


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
