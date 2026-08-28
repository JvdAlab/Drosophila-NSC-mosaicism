# pySCENIC reference resources

Reference resources for the pySCENIC workflow, with SHA256 checksums for
verifying downloaded copies.

## Downloaded resources

External reference databases, not produced by this repository:

- `dm6_v10_clust.genes_vs_motifs.rankings.feather`: gene-vs-motif rankings
  from the cisTarget v10 *Drosophila melanogaster* (`dm6`, clustered)
  collection
- `motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl`: the matching
  `v10nr_clust` motif-to-TF annotation (`m0.001-o0.0`)
- [`flybase_dmel_TFs.txt`](./flybase_dmel_TFs.txt): FlyBase
  *Drosophila melanogaster* transcription-factor list

The `.feather` and `.tbl` databases must be downloaded
separately. Get them from the aertslab cisTarget resource collection at
<https://resources.aertslab.org/cistarget/>.
The rankings are under `databases/`,
the annotations are under `motif2tf/`, place both in this directory before
running the pipeline. Verify each file against its SHA256 below after
downloading.

| File | SHA256 |
|:---|:---|
| `dm6_v10_clust.genes_vs_motifs.rankings.feather` | `3353dcf9396bbd91b84ed1d84b3b504271b48e144b91f4ea7e49297057fa4787` |
| `motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl` | `91284e94b0317b764dc2f8d8147d30db707605df757ff7de16fb0953c63fda2a` |
| `flybase_dmel_TFs.txt` | `cff167fe42c7344002bf23023b37b51bccb6876bac8ecf8cd9251005ae47995f` |

## Generated locally

| File | SHA256 |
|:---|:---|
| `global_pyscenic_rescue_map.csv` | `14f7893210281c832f0a67dbe02de0831ebe60fc30c5fc7739d8fcdf9e223139` |

Produced by [`../build_gene_rescue_map.R`](../build_gene_rescue_map.R) from the
annotated Seurat object and the ranking database above; consumed by
[`00_pyscenic_preprocess_data.py`](../../00_pyscenic_preprocess_data.py) and
[`04_compute_regulon_metrics.py`](../../04_compute_regulon_metrics.py).
