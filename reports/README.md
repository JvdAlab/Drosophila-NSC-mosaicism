# Analysis Reports

Quarto reports that consume the derived outputs of the pipeline modules
([`scRNA_pipeline/`](../scRNA_pipeline/), [`edgeR/`](../edgeR/), [`hdWGCNA/`](../hdWGCNA/), [`pySCENIC/`](../pySCENIC/)) and produce the
manuscript's computational figure panels. Where a rendered `.html` sits
alongside a `.qmd`, that HTML is the record of the run (exact values, tables,
figure captions, and `sessionInfo()`).

| Report | Manuscript panels | Key inputs |
|:---|:---|:---|
| [`gene_network_analysis.qmd`](./gene_network_analysis.qmd) | Fig 2E (module–condition correlation), 2F (hub-gene network), 2G (edgeR ∩ hdWGCNA overlap), 2M (SCENIC TF activity vs. Ldh coupling); Supp Fig 5C (hub-gene kME), 6A (regulon activity heatmap), 6B (regulon activity vs. Ldh coupling). Fig 2E module labels come from [`hdWGCNA/01_GO_enrichment.R`](../hdWGCNA/01_GO_enrichment.R). | `results/hdWGCNA/` (from [`hdWGCNA/00_hdWGCNA_analysis.R`](../hdWGCNA/00_hdWGCNA_analysis.R)); the curated [`hdWGCNA/resources/hdWGCNA_module_manuscript_labels.csv`](../hdWGCNA/resources/hdWGCNA_module_manuscript_labels.csv); the Neuroblast-subset Seurat `.rds`; `results/edgeR/edgeR_results/Neuroblasts/`; `results/pyscenic/`; [`pySCENIC/resources/pyscenic/global_pyscenic_rescue_map.csv`](../pySCENIC/resources/pyscenic/global_pyscenic_rescue_map.csv) |
| [`plot_manuscript_figures.qmd`](./plot_manuscript_figures.qmd) | scRNA-seq panel assembly — Fig 2B (cell-type UMAP), 2C (marker dot plot), 2D (Neuroblast edgeR volcano), 2I (*Ldh* violin, edgeR FDR stars); Supp Fig 4A–H (QC composite, condition/batch/cluster UMAPs, marker feature-density, transgene dot plot); Supp Fig 5A–B (GMC + combined-glia edgeR volcanoes); Supp Fig 8E/F (glial transporter dot plots); Fig 4E/F, 4G, Supp Fig 8I (metabolite-set enrichment, metabolomics volcano, *Eaat1/2* violins, pathway-grouped metabolite heatmap). | Annotated Seurat object, per-cell-type `edgeR_results/` tables, QC and doublet CSVs, and metabolomics `.rds` matrices, all staged under [`reports/data/`](#reportsdata--local-render-inputs); the committed metabolomics DA/fgsea tables from [`metabolomics/data/processed/`](../metabolomics/data/processed/); the committed [`transporter_gene_sets.csv`](./resources/transporter_gene_sets.csv) for Supp Fig 8E/F |

The clone/brain/cell-level mixed-model statistics for the confocal/genetics
figure panels (Fig 1B/1G/1P, 3E/3G, 4C/4K/4L, 5E and Supp Fig S1I-J/S7C/S7F/S7J/S9A/S9I)
are a self-contained analysis with committed inputs and are **not** in this
directory — see [`../statistical_analysis/`](../statistical_analysis/).

## Rendering

The hdWGCNA/network reports resolve paths with `here::here("results", …)` from
the repository root. `results/` is not populated in this package, so the
reports are **rendered by the authors in the analysis environment** and the
resulting `.html` is provided here. A reader consults that HTML; it cannot be
regenerated from this repository alone without the archived intermediate data.

### `reports/data/` — local render inputs

[`plot_manuscript_figures.qmd`](./plot_manuscript_figures.qmd) reads its
external rendering inputs from [`reports/data/`](./data/), which is
git-ignored — nothing in it is committed. Except for what's listed below,
everything under it is a source data table or object, not code; they are
pulled in so the report can be rendered locally, and the rendered `.html` in
this directory is the record of the run. (The clone/brain/cell-level
statistics live in their own module,
[`../statistical_analysis/`](../statistical_analysis/), with committed
inputs — nothing for them is staged here.)

| Local path (under `reports/data/`) | Purpose |
|---|---|
| `py_scRNA_objs/scRNA_data_QC_filtered_w_scVI_latent_annotated_seurat_dec_2025.rds` (393 MB) | annotated Seurat object (UMAPs, dot plots, violins, feature plots) |
| `edgeR_results/<cell_type>/<cell_type>_edgeR_results.tsv` (16 cell types) | per-cell-type pseudobulk edgeR tables (Fig 2D / S5A–B volcano panels; Fig 2I edgeR FDR stars). `Combined_GMCs` = GMC early+late, `Combined_Glia` = Astro+Cortex+Surface. |
| `py_scRNA_objs/QC_raw_singlets_for_plotting.csv` | Supp Fig 4A–C QC composite |
| `py_scRNA_objs/complete_thresholds_for_QC_plots.csv` | Supp Fig 4A–C QC composite (threshold lines) |
| `py_scRNA_objs/{Control1,Control3,Dysfunction1,Dysfunction3}_doublet_probs.csv` | Supp Fig 4C doublet histograms |
| `metabolomics/metabolomics_log2_count_matrix.rds` | Supp Fig 8I pathway-grouped metabolite heatmap (log2 abundances) |

The two metabolomics result tables the report reads —
`metabolomics_differential_abundance_results.csv` (Fig 4F) and
`metabolomics_fgsea_results.csv` (Fig 4E) — are **not** staged here; they are
committed under [`metabolomics/data/processed/`](../metabolomics/data/processed/)
instead (produced by `metabolomics/00_differential_abundance.R` and
`01_pathway_enrichment.R`).

### Neuroblast-subset objects

[`gene_network_analysis.qmd`](./gene_network_analysis.qmd) reads a
Neuroblast-only Seurat object,
`results/scRNA/scRNA_objects/scRNA_data_QC_filtered_w_scVI_latent_annotated_nb_subset_seurat.rds`.
No packaged script produces it; it is a one-line subset of the annotated object
from [`scRNA_pipeline/06_convert_adata_to_seurat.R`](../scRNA_pipeline/06_convert_adata_to_seurat.R):

```r
nb <- subset(seurat_obj, subset = cell_type_manual == "Neuroblasts")
saveRDS(nb, "…/scRNA_data_QC_filtered_w_scVI_latent_annotated_nb_subset_seurat.rds")
```

The equivalent AnnData subset feeds pySCENIC and is documented, with its
derivation and checksum, in
[`pySCENIC/data/pyscenic_input/README.md`](../pySCENIC/data/pyscenic_input/README.md).

The Supp Fig 5D module GO tree plots and the curated Fig 2E module labels are
produced by [`hdWGCNA/01_GO_enrichment.R`](../hdWGCNA/01_GO_enrichment.R) (a
pipeline-module script, not a report), which has its own manual label-review
step — see [`hdWGCNA/README.md`](../hdWGCNA/README.md).

## Curated inputs

### [`resources/transporter_gene_sets.csv`](./resources/transporter_gene_sets.csv) (Supp Fig 8E/F)

The amino-acid- and monocarboxylate-transporter gene sets shown in
`plot_manuscript_figures.qmd` Supp Fig 8E/F. Columns: `panel` (S8E / S8F),
`class` (`amino_acid` / `monocarboxylate`), `gene`, `note`.

**Provenance.** These sets were originally compiled from an external transporter
gene list whose source file predates version control and could not be recovered
(its composition — many `CG` identifiers spanning several SLC families, plus
mitochondrial carriers and sideroflexins — points to a FlyBase gene-group or GO
molecular-function export, but the exact query is not known). The primary source
of record is therefore the **published Supplementary Figure 8E,F itself**: the
lists in the CSV were transcribed verbatim from that figure. Genes below the
plot's detection threshold in these four cell types drop off at render time, as
in the published panels.

A small number of monocarboxylate-panel entries have broad annotations and are
not canonical transporters (secreted phospholipases `sPLA2` / `GIIIspla2` /
`GXIVsPLA2`, `Nrx-1`, `Prestin`); they are kept to reproduce the published
figure and flagged in the `note` column. `VGAT` and `Gat` appear in both panels,
matching the figure.

## Environment

```bash
conda env create -f reports/environment.yml
conda activate reports_pipeline
```

Versions are pinned from the `sessionInfo()` recorded in the rendered
`gene_network_analysis.html` (R 4.4.3, Bioconductor 3.20). This environment is
larger than [`../hdWGCNA/environment.yml`](../hdWGCNA/environment.yml) because it
adds the GO-analysis, network-visualisation, single-cell-figure, and
HTML-rendering packages (clusterProfiler, rrvgo, enrichplot, GOSemSim, ggsci,
ggraph, tidygraph, eulerr, ComplexHeatmap, reactable, SCpubr, ggrastr,
knitr/rmarkdown, …). The statistical-analysis module has its own
[`../statistical_analysis/environment.yml`](../statistical_analysis/environment.yml)
(R 4.6.1).
