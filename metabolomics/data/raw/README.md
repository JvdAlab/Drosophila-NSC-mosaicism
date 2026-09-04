# Metabolomics input data

Untargeted LC-MS metabolomics of whole larval CNS (brain + ventral nerve
cord), wandering third-instar larvae, comparing Control (mCherry RNAi) with
Complex I / *ND-42* RNAi in neural stem cells (`UAS-ND-42-RNAi`, BDSC
#32998). Six biological replicates per condition, ~10 mg brains each.

## `metabolite_abundances_normalized.xlsx`

Collaborator-supplied peak areas, normalized to the internal standard
(D6-Glutaric Acid) then total-sum normalized (originally
`Normalised to internal standard and sum.xlsx`). One sheet, 105 rows ×
17 columns; the direct input to
[`metabolomics/00_differential_abundance.R`](../../00_differential_abundance.R).

The four columns the pipeline reads:

| Column(s) | Contents |
|:---|:---|
| `Molecule` | Metabolite name; becomes the `gene` identifier in the results and the row key of the abundance matrices. |
| `Molecule.List` | Metabolite class; carried through as `Pathway` in the results and used as the gene-set definition for enrichment (stage 01). The value `Internal Std` marks the `D6-Glutaric Acid` row and is dropped before enrichment. |
| `mcherry_65b.Area` … `mcherry_90.Area` (6) | Normalized abundance for the six Control replicates; selected by `starts_with("mcherry")`. One half of the analysis matrix. |
| `ND42_51b.Area` … `ND42_84b.Area` (6) | Normalized abundance for the six *ND-42* KD replicates; selected by `starts_with("ND42")`. The other half. |

Rows: 103 metabolites plus a `D6-Glutaric Acid` (internal standard) row and a
`Sum` (total-ion) row; those two carry no differential-abundance result.

