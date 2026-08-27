# edgeR Pseudobulk Differential Expression

## Overview

Builds per-cell-type pseudobulk profiles from the annotated scRNA-seq object
and tests ND75-KD vs. Control differential gene expression with edgeR.

**Input:**
- `results/scRNA/scRNA_data_QC_filtered_w_scVI_latent_annotated.h5ad`
  (annotated AnnData object, produced by
  [`scRNA_pipeline/04_cluster_data.py`](../scRNA_pipeline/04_cluster_data.py))

**Output:**
- Pseudobulk count/metadata/gene-info TSVs per cell type
- Per-condition detection counts per cell type (`<cell_type>_percent_expressing.tsv`)
- edgeR differential expression tables, MDS/BCV/QL-dispersion diagnostic plots,
  and volcano plots per cell type

---

## Workflow context

All commands run from the repository root.

### Environment

```bash
conda env create -f edgeR/environment.yml
conda activate edgeR_pipeline
```

### Pipeline commands

```bash
python edgeR/00_create_pseudobulk.py
python edgeR/01_create_pct_expressing.py
Rscript edgeR/02_edgeR_DGE.R
```

---

## Pipeline Scripts

All outputs are written under `results/edgeR/`.

| Stage | Script | Input | Outputs |
|:---|:---|:---|:---|
| 00: Pseudobulk | [`00_create_pseudobulk.py`](./00_create_pseudobulk.py) | Annotated `.h5ad` | Per group in `pseudobulk_all_genes/<group>/`: `<group>_pseudobulk_counts.tsv`<br>`sample_info_<group>.tsv`<br>`gene_info_<group>.tsv` |
| 01: Detection counts | [`01_create_pct_expressing.py`](./01_create_pct_expressing.py) | Annotated `.h5ad` and the stage-00 group directories | `pseudobulk_all_genes/<group>/<group>_percent_expressing.tsv` (`gene`, `n_cells_Control`, `n_cells_ND75_KD`) |
| 02: edgeR DGE | [`02_edgeR_DGE.R`](./02_edgeR_DGE.R) | Stage-00 and stage-01 outputs | Per group in `edgeR_results/<group>/`: `<group>_edgeR_results.tsv` (with `is_artifact`)<br>`<group>_MDS_plot.png`<br>`<group>_BCV_plot.png`<br>`<group>_QL_dispersion_plot.png`<br>`<group>_edgeR_volcano_plot.png` |

The scripts run in numeric order. `00` also drops any unannotated / `"Unknown"`
cell type if present and builds two composite groups:
`Combined_Glia` = Astrocytes + Cortex_Glia + Surface_Glia;
`Combined_GMCs` = Ganglion_Mother_Cells_(Early) + Ganglion_Mother_Cells_(Late).

---

## Analysis Details

- **Condition factor:** `Control`, `ND75-KD`; **batch factor:** `Batch1`, `Batch2`
- **Design matrix:** `~ batch + condition`, tested via `glmQLFTest` on the `condition` coefficient
- **Low-replicate skip:** any group where the design matrix has ≥ as many columns as samples is skipped (too few replicates to fit)
- **Gene filtering:** `filterByExpr` on the design; volcano cutoffs FDR < 0.05, \|logFC\| ≥ 1.0
- **Artifact flag:** a significant DEG is marked `is_artifact` if its direction of change isn't supported by ≥5 cells expressing in the relevant condition (upregulated genes need ≥5 cells in ND75-KD, downregulated need ≥5 in Control); flagged genes are excluded from the volcano plot but retained in the results table

The per-condition detection counts used by the artifact flag come from
[`01_create_pct_expressing.py`](./01_create_pct_expressing.py) (`<group>_percent_expressing.tsv`). If that file
is missing for a group, [`02_edgeR_DGE.R`](./02_edgeR_DGE.R) logs a warning and sets
`is_artifact = FALSE` for every gene in that group.

---

## Dependencies

The pipeline dependencies are defined in [`environment.yml`](./environment.yml).

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [Python](https://www.python.org/downloads/) | 3.10.17 | [`00_create_pseudobulk.py`](./00_create_pseudobulk.py), [`01_create_pct_expressing.py`](./01_create_pct_expressing.py) |
| [anndata](https://anndata.readthedocs.io/) | 0.11.4 | Reading the annotated AnnData object |
| [numpy](https://numpy.org/) | 1.26.4 | Detection-count arrays in [`01_create_pct_expressing.py`](./01_create_pct_expressing.py) |
| [pandas](https://pandas.pydata.org/) | 1.5.3 | Pseudobulk aggregation |
| [R](https://www.r-project.org/) | 4.4.3 | [`02_edgeR_DGE.R`](./02_edgeR_DGE.R) |
| [edgeR](https://bioconductor.org/packages/edgeR/) | 4.4.2 | Differential expression testing |
| [limma](https://bioconductor.org/packages/limma/) | 3.62.2 | edgeR backend |
| [EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano/) | 1.24.0 | Volcano plots |
| [statmod](https://cran.r-project.org/package=statmod) | 1.5.0 | Quasi-likelihood dispersion estimation |
| [tidyverse](https://www.tidyverse.org/) | 2.0.0 | dplyr/tibble verbs in [`02_edgeR_DGE.R`](./02_edgeR_DGE.R) |
| [readr](https://readr.tidyverse.org/) | 2.2.0 | TSV I/O in [`02_edgeR_DGE.R`](./02_edgeR_DGE.R) |
| [here](https://cran.r-project.org/package=here) | 1.0.2 | Repository-relative paths |

Each script specifies its own output paths directly and expects the repository
root as its working directory. [`00_create_pseudobulk.py`](./00_create_pseudobulk.py) and
[`01_create_pct_expressing.py`](./01_create_pct_expressing.py) read `results/scRNA/`, produced by
[`scRNA_pipeline/04_cluster_data.py`](../scRNA_pipeline/04_cluster_data.py)
(see Input above)

---

## Expected Output

**Pseudobulk** (`results/edgeR/pseudobulk_all_genes/<cell_type>/`):
- `<cell_type>_pseudobulk_counts.tsv`, `sample_info_<cell_type>.tsv`, `gene_info_<cell_type>.tsv`
- `<cell_type>_percent_expressing.tsv` — `gene`, `n_cells_Control`, `n_cells_ND75_KD` (from [`01_create_pct_expressing.py`](./01_create_pct_expressing.py))

**edgeR results** (`results/edgeR/edgeR_results/<cell_type>/`):
- `<cell_type>_edgeR_results.tsv` — full DEG table with `is_artifact` flag
- `<cell_type>_MDS_plot.png`, `<cell_type>_BCV_plot.png`, `<cell_type>_QL_dispersion_plot.png`
- `<cell_type>_edgeR_volcano_plot.png`
