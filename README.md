# Drosophila NSC mosaicism

<!-- Add on Zenodo release:
[![DOI](https://zenodo.org/badge/DOICODE.svg)](https://doi.org/10.5281/zenodo.XXXXXXXX)
-->

Analysis code accompanying:

> Petridi, S., Dhawanjewar, A., Dubal, D., Chandrasegaram, R.,
> Dickmänken, H., Bala-Muraly, N., Wilson, B. A., Eve, T., Marzullo, B.,
> Hynes-Allen, A., Jones, S. A., King, M. S., Butler, R., Sciacovelli, M.,
> Kunji, E. R. S., Tennant, D. A. & van den Ameele, J.
> **Tissue-wide metabolic buffering confers resilience to mosaic
> mitochondrial dysfunction.** bioRxiv 2026.08.27.747448 (2026).
> <https://doi.org/10.64898/2026.08.27.747448>

---

## Overview

Computational analyses for a study of how the *Drosophila* larval brain
tolerates mosaic mitochondrial (Complex I) dysfunction in its neural stem
cell (NSC) niche. Single-cell RNA-seq of the NSC niche feeds sample-level
pseudobulk differential expression (edgeR), co-expression network
analysis (hdWGCNA), and gene-regulatory-network inference (pySCENIC). An
independent untargeted metabolomics assay of whole larval brains compares
Complex I / *ND-42* RNAi in NSCs against control.

---

## Analyses

| Analysis | Description | Report | Source |
|----------|-------------|--------|--------|
| Single-cell RNA-seq | Raw object setup, doublet detection, adaptive QC, scVI representation, clustering, annotation, Seurat conversion | [plot_manuscript_figures.html](reports/plot_manuscript_figures.html) | [scRNA_pipeline/](scRNA_pipeline/) |
| Pseudobulk DGE | Sample-level pseudobulk construction and edgeR differential expression | [plot_manuscript_figures.html](reports/plot_manuscript_figures.html) | [edgeR/](edgeR/) |
| Co-expression networks | Neuroblast hdWGCNA co-expression network construction | [gene_network_analysis.html](reports/gene_network_analysis.html) | [hdWGCNA/](hdWGCNA/) |
| Gene regulatory networks | pySCENIC GRN inference, motif pruning, regulon AUC scoring (DVC pipeline) | [gene_network_analysis.html](reports/gene_network_analysis.html) | [pySCENIC/](pySCENIC/) |
| Metabolomics | Larval-brain differential abundance (limma) and metabolite-set enrichment (fgsea) | [plot_manuscript_figures.html](reports/plot_manuscript_figures.html) | [metabolomics/](metabolomics/) |
| Statistical analysis | Clone/brain/cell-level mixed-model statistics for the confocal/genetics figure panels (lme4/lmerTest, Type III ANOVA, emmeans, DHARMa) | [statistical_analysis.html](statistical_analysis/statistical_analysis.html) | [statistical_analysis/](statistical_analysis/) |
| Figure reports | Quarto reports assembling network, SCENIC, GO, and metabolomics figure panels | [reports/](reports/) | [reports/](reports/) |
| Image analysis | ImageJ/Fiji panorama-analysis Java component | — | [Deadpanorama/](Deadpanorama/) |

---

## Repository Structure

```
Drosophila-NSC-mosaicism/
├── scRNA_pipeline/          # scRNA-seq: QC → scVI → clustering → annotation → Seurat export
├── edgeR/                   # Pseudobulk + edgeR differential expression
├── hdWGCNA/                 # Co-expression network construction
├── pySCENIC/                # GRN inference (DVC pipeline)
├── metabolomics/            # Differential abundance (limma) + enrichment (fgsea)
├── statistical_analysis/    # Mixed-model statistics for the confocal/genetics panels
├── reports/                 # Quarto figure-panel reports
└── Deadpanorama/            # ImageJ/Fiji Java component
```

Pipeline outputs are written to a local `results/` directory at run time; it
is git-ignored and not part of this checkout (see Data Availability below).

See each module's `README.md` for full parameters, dependencies, and usage.

---

## Requirements

Each module carries its own conda `environment.yml` and builds a
self-named environment. Create the ones you need, or run
[`./setup_envs.sh`](setup_envs.sh) from the repository root to build them
all. All scripts
expect the repository root as the working directory; the pySCENIC
pipeline is driven by DVC (`dvc repro`) and targets a Python/GPU or HPC
environment.

| Module | Conda env |
|--------|-----------|
| [`scRNA_pipeline/`](scRNA_pipeline/) | `scRNA_pipeline` |
| [`edgeR/`](edgeR/) | `edgeR_pipeline` |
| [`hdWGCNA/`](hdWGCNA/) | `hdWGCNA_pipeline` |
| [`metabolomics/`](metabolomics/) | `metabolomics_pipeline` |
| [`statistical_analysis/`](statistical_analysis/) | `statistical_analysis` |
| [`reports/`](reports/) | `reports_pipeline` |
| [`pySCENIC/`](pySCENIC/) | `pyscenic_pipeline` |

### Recreate the environment

```bash
conda env create -f <module>/environment.yml
conda activate <conda env name>
```

---

## Data Availability

Large inputs, reference databases, and generated results are external
artifacts, deposited or archived separately.

| Dataset | Repository | Accession |
|:--------|:-----------|:----------|
| scRNA-seq (10x, larval CNS neuroblast niche) | GEO | pending |
| Larval-brain untargeted metabolomics | Metabolomics Workbench / MetaboLights | normalized abundance table committed in repo; raw spectra pending |
| SCENIC input AnnData and reference resources | Zenodo / archive | pending (SHA256 checksums in module READMEs) |

Each `data/` or `resources/` directory that expects an external file carries a
README naming exactly what belongs there and, where applicable, its SHA256
checksum; see that module's own README (linked from the Analyses table
above) for what it holds and where it comes from.

---

## Citation

If you use this code, please cite the study:

> Petridi, S., Dhawanjewar, A., Dubal, D., Chandrasegaram, R.,
> Dickmänken, H., Bala-Muraly, N., Wilson, B. A., Eve, T., Marzullo, B.,
> Hynes-Allen, A., Jones, S. A., King, M. S., Butler, R., Sciacovelli, M.,
> Kunji, E. R. S., Tennant, D. A. & van den Ameele, J.
> **Tissue-wide metabolic buffering confers resilience to mosaic
> mitochondrial dysfunction.** bioRxiv 2026.08.27.747448 (2026).
> <https://doi.org/10.64898/2026.08.27.747448>

---

## License

The analysis code in this repository is released under the
[MIT License](LICENSE).

[`Deadpanorama/`](Deadpanorama/) is a standalone ImageJ/Fiji plugin,
released separately under the
[GNU General Public License v3.0 (or later)](Deadpanorama/LICENSE)

---

## Contact

- Corresponding Author: [Jelle van den Ameele](mailto:jv361@cam.ac.uk)
- Maintainer: [Abhilesh Dhawanjewar](mailto:ad2347@cam.ac.uk)
