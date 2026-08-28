# Metabolomics Differential Abundance

## Overview

Untargeted LC-MS metabolomics of whole larval CNS, comparing Complex I
(*ND-42*) RNAi in neural stem cells against Control. Stage 00 runs
per-metabolite differential abundance with `limma`; stage 01 runs
metabolite-set enrichment (`fgsea`) on the stage-00 ranking, with metabolite
classes as the sets.

**Input:**
- [`metabolomics/data/raw/metabolite_abundances_normalized.xlsx`](./data/raw/metabolite_abundances_normalized.xlsx):
  peak areas normalized to the internal standard (D6-Glutaric Acid) and
  total-sum normalized upstream

**Output** (`results/metabolomics/`):
- `metabolomics_differential_abundance_results.csv`: per-metabolite `limma`
  table: log2 fold change, moderated *t*, raw and BH-adjusted *p*, per-group
  mean abundance, and metabolite class
- `metabolomics_fgsea_results.csv`: per-class `fgsea` enrichment: ES, NES,
  *p*, BH `padj`, leading-edge metabolites, and direction
- `metabolomics_count_matrix.rds`, `metabolomics_log2_count_matrix.rds`,
  `metabolomics_meta_info.rds`: intermediate matrices reused by
  [`../reports/plot_manuscript_figures.qmd`](../reports/plot_manuscript_figures.qmd)

---

## Workflow context

All commands run from the repository root.

### Environment

```bash
conda env create -f metabolomics/environment.yml
conda activate metabolomics_pipeline
```

The same R 4.4.3 / Bioconductor 3.20 stack as the other R modules, adding
`limma` and `fgsea` (`readxl` comes with `tidyverse`).

### Pipeline commands

```bash
Rscript metabolomics/00_differential_abundance.R
Rscript metabolomics/01_pathway_enrichment.R
```

---

## Pipeline Scripts

All outputs are written under `results/metabolomics/`.

| Stage | Script | Input | Outputs |
|:---|:---|:---|:---|
| 00: Differential abundance | [`00_differential_abundance.R`](./00_differential_abundance.R) | `data/raw/metabolite_abundances_normalized.xlsx` | `metabolomics_differential_abundance_results.csv`<br>`metabolomics_count_matrix.rds`<br>`metabolomics_log2_count_matrix.rds`<br>`metabolomics_meta_info.rds` |
| 01: Pathway enrichment | [`01_pathway_enrichment.R`](./01_pathway_enrichment.R) | `metabolomics_differential_abundance_results.csv` | `metabolomics_fgsea_results.csv` |

---

## Analysis Details

- **Design:** 6 Control (`mcherry*`) vs. 6 *ND-42* KD (`ND42*`) samples;
  abundances log2-transformed with a `+ 1e-20` offset
- **Differential abundance (00):** `limma` linear model `~ 0 + group` with the
  `Condition − Control` contrast, `eBayes` moderation, and Benjamini-Hochberg
  adjusted *p*-values (`topTable`); per-group mean raw abundance is appended
  for the downstream volcano plot
- **Enrichment (01):** metabolite classes (`Molecule.List`) as sets; ranking
  metric `sign(logFC) × −log10(P.Value)`; `fgsea` with `minSize = 2`,
  `maxSize = 500`; the `Internal Std` class and unclassified metabolites are
  dropped before ranking
- **Non-metabolite rows:** the `D6-Glutaric Acid` internal-standard row and the
  `Sum` total-ion row carry no differential-abundance result
- **Reproducibility:** `fgsea` is seeded (`set.seed(42)`); stage 00 is
  deterministic

---

## Dependencies

Pinned in [`environment.yml`](./environment.yml) (`metabolomics_pipeline`) —
the R 4.4.3 / Bioconductor 3.20 stack shared with the other R modules.

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | 4.4.3 | `00_differential_abundance.R`, `01_pathway_enrichment.R` |
| [tidyverse](https://www.tidyverse.org/) | 2.0.0 | Data manipulation and CSV I/O |
| [readxl](https://readxl.tidyverse.org/) | via tidyverse | Reading the abundance workbook |
| [limma](https://bioconductor.org/packages/limma/) | 3.62.2 | Per-metabolite differential abundance |
| [fgsea](https://bioconductor.org/packages/fgsea/) | 1.32.0 | Metabolite-set enrichment |
| [here](https://cran.r-project.org/package=here) | 1.0.2 | Repository-relative paths |

---

## Expected Output

`results/metabolomics/`:

**Differential abundance (stage 00)**
- `metabolomics_differential_abundance_results.csv` — one row per metabolite
  (`gene`): `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, `Pathway`,
  `avg_control`, `avg_condition`
- `metabolomics_count_matrix.rds`, `metabolomics_log2_count_matrix.rds` —
  metabolite × sample abundance matrices (raw and log2)
- `metabolomics_meta_info.rds` — metabolite-to-class table

**Enrichment (stage 01)**
- `metabolomics_fgsea_results.csv` — one row per metabolite class: `pathway`,
  `pval`, `padj`, `log2err`, `ES`, `NES`, `size`, `leadingEdge` (`;`-joined),
  `direction`, `neg_log10_p`
