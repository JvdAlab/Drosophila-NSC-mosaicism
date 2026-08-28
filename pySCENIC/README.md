# pySCENIC Regulon Inference

## Overview

Infers transcription-factor regulons and scores per-cell regulon activity in
the annotated Neuroblast subset. Stages 00–03 are the pySCENIC steps packaged as a
[DVC](https://dvc.org/) pipeline: RegDiffusion gene-regulatory-network
inference, cisTarget motif pruning, and AUCell activity scoring; stage 04
computes per-regulon differential-activity statistics from the AUCell output.

**Input:**
- `pySCENIC/data/pyscenic_input/scRNA_scVI_annotated_nb_subset.h5ad`: annotated
  Neuroblast subset of the scRNA-seq pipeline output (see
  [`data/pyscenic_input/README.md`](./data/pyscenic_input/README.md))
- SCENIC reference resources under `pySCENIC/resources/pyscenic/`: ranking
  database, motif table, FlyBase TF list, and gene rescue map. The ranking
  database and motif table must be downloaded first — see
  [`resources/pyscenic/README.md`](./resources/pyscenic/README.md)

**Output:**
- Filtered expression matrix and validated TF list (`preprocessed.h5ad`, `valid_tfs.txt`)
- Inferred TF→target GRN edge list (`adjacencies.csv` / `.pkl`)
- cisTarget-pruned regulons and their motif-enrichment table (`final_regulons.csv`, `ctx_pruning_results.csv` / `.pkl`)
- Per-cell regulon activity matrix from AUCell (`auc_matrix.csv`)
- Per-regulon differential-activity statistics, a regulon–*Ldh* correlation
  table, and a ranked annotated table (`regulon_metrics.csv`,
  `regulon_Ldh_correlation.csv`, `top_regulons_annotated.csv`)

<br>

> **This arm is a DVC pipeline:** Stages 00–03 here are defined
> in [`dvc.yaml`](./dvc.yaml) and locked in [`dvc.lock`](./dvc.lock), and are
> re-run with `dvc repro` (which skips stages whose code, parameters, and
> inputs are unchanged). Stage 04 is a post-DVC step run directly.

---

## Workflow Diagram

```mermaid
flowchart LR
    INPUT[Annotated Neuroblast h5ad and SCENIC resources]
    PRE[00 Preprocess]
    GRN[01 GRN construction]
    CTX[02 Context pruning]
    AUC[03 AUC scoring]
    MET[04 Regulon metrics]
    OUTPUT[Regulons, AUC matrix, DRA statistics]

    INPUT --> PRE
    PRE --> GRN
    GRN --> CTX
    CTX --> AUC
    AUC --> MET
    MET --> OUTPUT
```

---

## Workflow context

Run from the repository root.

### Gene rescue map (prerequisite)

[`00_pyscenic_preprocess_data.py`](./00_pyscenic_preprocess_data.py) and [`04_compute_regulon_metrics.py`](./04_compute_regulon_metrics.py) both read
[`resources/pyscenic/global_pyscenic_rescue_map.csv`](./resources/pyscenic/global_pyscenic_rescue_map.csv), which maps dataset gene symbols
onto the SCENIC/cisTarget namespace. It is included in the repository, so the packaged run should work. If it is missing, or if the annotated
object or ranking database changes, rebuild it before running the pipeline:

```bash
Rscript pySCENIC/build_gene_rescue_map.R
```

[`build_gene_rescue_map.R`](./build_gene_rescue_map.R) needs the annotated
Seurat object
(`results/scRNA/scRNA_data_QC_filtered_w_scVI_latent_annotated_seurat.rds`,
from [`scRNA_pipeline/06_convert_adata_to_seurat.R`](../scRNA_pipeline/06_convert_adata_to_seurat.R)) and the ranking database.

### Environment

```bash
conda env create -f pySCENIC/environment.yml
conda activate pyscenic_pipeline
```

### Pipeline commands

`dvc` is installed by the conda environment above, however DVC must be initialized once before the first run:

```bash
# Run dvc init the first time and exclude from subsequent runs
# Will create .dvc
dvc init
./pySCENIC/run_pyscenic_dvc.sh
python pySCENIC/04_compute_regulon_metrics.py --config pySCENIC/pyscenic_config.yaml
```

Notes:

- [`run_pyscenic_dvc.sh`](./run_pyscenic_dvc.sh) just calls `dvc repro pySCENIC/dvc.yaml`. Add `--force`
  to rerun every stage; set `DVC_BIN` to point at a specific `dvc` executable.
- On a fresh checkout the `results/pyscenic/` outputs are absent, so the first
  `dvc repro` runs every stage and rewrites `dvc.lock` with that run's hashes.
- The four DVC stages can also be run individually
  (`python pySCENIC/0N_*.py --config pySCENIC/pyscenic_config.yaml`); doing so
  does **not** update `dvc.lock`.

---

## Pipeline Scripts

All outputs are written under `results/pyscenic/`.

| Stage | Script | Input | Outputs |
|:---|:---|:---|:---|
| 00: Preprocess | [`00_pyscenic_preprocess_data.py`](./00_pyscenic_preprocess_data.py) | Annotated `.h5ad`, ranking database, TF list, gene rescue map | `preprocessed.h5ad`<br>`valid_tfs.txt` |
| 01: GRN construction | [`01_pyscenic_grn_construction.py`](./01_pyscenic_grn_construction.py) | `preprocessed.h5ad`, `valid_tfs.txt` | `adjacencies.csv`<br>`adjacencies.pkl` |
| 02: Context pruning | [`02_pyscenic_ctx_pruning.py`](./02_pyscenic_ctx_pruning.py) | `adjacencies.pkl`, `preprocessed.h5ad`, ranking database, motif annotations | `ctx_pruning_results.csv`<br>`ctx_pruning_results.pkl`<br>`final_regulons.csv` |
| 03: AUC scoring | [`03_pyscenic_auc_scoring.py`](./03_pyscenic_auc_scoring.py) | `ctx_pruning_results.pkl`, `preprocessed.h5ad` | `auc_matrix.csv` |
| 04: Regulon metrics | [`04_compute_regulon_metrics.py`](./04_compute_regulon_metrics.py) | `auc_matrix.csv`, `ctx_pruning_results.csv`, input `.h5ad`, gene rescue map | `regulon_metrics.csv`<br>`regulon_Ldh_correlation.csv`<br>`top_regulons_annotated.csv` |

Stages 00–03 are DVC-tracked ([`dvc.yaml`](./dvc.yaml) / [`dvc.lock`](./dvc.lock));
stage 04 is a post-DVC analysis step. The gene rescue map consumed by stages 00
and 04 must exist first — see
[Gene rescue map](#gene-rescue-map-prerequisite) above.

---

## Analysis Details

- **Input scope:** annotated Neuroblast subset, 707 cells (Control 368,
  ND75-KD 339); GRN inference uses the object's `log1p_norm` layer
- **GRN inference (01):** RegDiffusion; edges taken from the top
  `top_gene_percentile = 85` importance genes, filtered to valid TFs, with
  self-loops removed
- **Context pruning (02):** cisTarget against the dm6 v10 clustered ranking
  database; `rank_threshold = 3000`, `nes_threshold = 3.5`,
  `motif_similarity_fdr = 0.001`, `auc_threshold = 0.01`
- **AUC scoring (03):** AUCell on the pruned regulons
- **Regulon metrics (04):** per-regulon differential regulon activity between
  conditions: Mann-Whitney U, rank-biserial r, Cohen's d, Benjamini-Hochberg
  FDR and Spearman correlation between each regulon's AUCell activity and
  *Ldh* expression; TF symbols absent from the dataset fall back to the gene
  rescue map
- **Condition factor:** `Control`, `ND75-KD`
- **Reproducibility:** a single seed (`shared_params.seed = 42`) seeds `random`,
  `numpy`, and `torch` in every stage; RegDiffusion on GPU is not fully
  deterministic across hardware
- **Parameters:** all paths and the values above are defined in
  [`pyscenic_config.yaml`](./pyscenic_config.yaml)

---

## Dependencies

The pipeline dependencies are defined in [`environment.yml`](./environment.yml).

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [Python](https://www.python.org/downloads/) | 3.9.23 | Pipeline scripts |
| [anndata](https://anndata.readthedocs.io/) | 0.8.0 | Single-cell data containers |
| [NumPy](https://numpy.org/) | 1.23.5 | Numerical arrays |
| [pandas](https://pandas.pydata.org/) | 2.2.3 | Tabular data processing |
| [PyArrow](https://arrow.apache.org/docs/python/) | 20.0.0 | Feather ranking-database access |
| [PyYAML](https://pyyaml.org/) | 6.0.2 | Configuration parsing |
| [Scanpy](https://scanpy.readthedocs.io/) | 1.9.3 | AnnData and single-cell utilities |
| [SciPy](https://scipy.org/) | via Scanpy | Mann-Whitney U and Spearman tests in `04` |
| [statsmodels](https://www.statsmodels.org/) | via Scanpy | Benjamini-Hochberg FDR in `04` |
| [pySCENIC](https://pyscenic.readthedocs.io/) | 0.12.1 | cisTarget pruning and AUCell scoring |
| [ctxcore](https://github.com/aertslab/ctxcore) | 0.2.0 | cisTarget ranking-database access |
| [RegDiffusion](https://github.com/aertslab/RegDiffusion) | 0.1.1 | Gene-regulatory-network inference |
| [PyTorch](https://pytorch.org/) | 2.6.0 with CUDA 12.4 | RegDiffusion backend |
| [DVC](https://dvc.org/) | 3.66.1 | Pipeline execution and provenance tracking |

Each script specifies its own output paths directly and expects the repository
root as its working directory. The input `.h5ad` is a manually-prepared
Neuroblast subset of the scRNA-seq pipeline's annotated output ([`data/pyscenic_input/README.md`](./data/pyscenic_input/README.md)).

---

## Expected Output

`results/pyscenic/`:

**GRN and regulons**
- `preprocessed.h5ad`, `valid_tfs.txt` — filtered expression matrix and validated TF list
- `adjacencies.csv` / `.pkl` — inferred TF→target edge list
- `ctx_pruning_results.csv` / `.pkl` — cisTarget motif-enrichment results
- `final_regulons.csv` — human-readable pruned regulons

**Activity**
- `auc_matrix.csv` — cell × regulon AUCell activity matrix

**Differential activity (stage 04)**
- `regulon_metrics.csv` — per-regulon differential-activity statistics
- `regulon_Ldh_correlation.csv` — per-regulon Spearman correlation with *Ldh*
- `top_regulons_annotated.csv` — regulons ranked by \|effect size\|, joined with the *Ldh*-correlation columns

**Logs**
- `NN_*.log` — one per stage
