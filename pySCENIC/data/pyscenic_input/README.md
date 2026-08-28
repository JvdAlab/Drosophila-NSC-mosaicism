# pySCENIC input data

`scRNA_scVI_annotated_nb_subset.h5ad`: the annotated Neuroblast subset used as
input to the pySCENIC pipeline.

It is the Neuroblast population
(`adata[adata.obs["cell_type_manual"] == "Neuroblasts"]`) of the annotated
AnnData written by
[`scRNA_pipeline/04_cluster_data.py`](../../../scRNA_pipeline/04_cluster_data.py),
carrying the `counts` and `log1p_norm` layers and the `condition` /
`sample_id` values `Control` and `ND75-KD`.

Deposited separately (accession pending).

```text
cf7af5e0a8055a73b2681b06c673dea049402cf0ecce7cb5580f98449312ac93  scRNA_scVI_annotated_nb_subset.h5ad
```
