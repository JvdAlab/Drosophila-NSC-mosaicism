# scRNA-seq raw data

`CellRanger_Outputs/<sample>/outs/filtered_feature_bc_matrix/` holds the
per-sample Cell Ranger count matrices consumed by
[`scRNA_pipeline/00_setup_adata_obj.py`](../../00_setup_adata_obj.py). Four samples were processed across two
sequencing runs; [`samplesheet.csv`](./samplesheet.csv) holds their sample
IDs, biological conditions, batch mapping, replicate labels, and external
Cell Ranger FASTQ sample names.

| Sample | Sequencing run | Condition |
|:---|:---|:---|
| Control-1 | SLX-21964 (Batch1) | mCherry (Control) |
| ND75-KD-1 | SLX-21964 (Batch1) | ND75 (ND75-KD) |
| Control-3 | SLX-22049 (Batch2) | mCherry (Control) |
| ND75-KD-3 | SLX-22049 (Batch2) | ND75 (ND75-KD) |

## Cell Ranger preprocessing

FASTQ files were processed with
[Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest)
v6.1.2 against a `cellranger mkref` index built from Ensembl
`Drosophila_melanogaster.BDGP6.32.105` — assembly BDGP6.32 (dm6; GenBank
GCA_000001215.4), annotation release 105, from the genome FASTA and GTF. Each
sample was quantified with:

```bash
cellranger count --id=<sample> \
  --transcriptome=<path to cellranger_index_updated> \
  --fastqs=<path to sample's raw FASTQ directory> \
  --sample=<ND75|mCherry> \
  --localcores=8 --localmem=24
```

## GEO submission

Both the raw FASTQ files and the processed `filtered_feature_bc_matrix/`
matrices for the four samples are submitted to GEO. Accession numbers: **to
be added once the GEO record is created**.
