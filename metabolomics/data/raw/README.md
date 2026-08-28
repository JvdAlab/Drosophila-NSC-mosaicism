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

| Column(s) | Contents |
|:---|:---|
| `Molecule` | Metabolite name |
| `Molecule.List` | Metabolite class (used as the `fgsea` gene set) |
| `mcherry_*.Area` (6) | Control sample abundances |
| `ND42_*.Area` (6) | *ND-42* KD sample abundances |
| `Average Normalised …`, `fold change` | Precomputed summaries, unused by the analysis |

The 105 rows are 103 metabolites plus a `D6-Glutaric Acid` (internal
standard) row and a `Sum` (total-ion) row; those two carry no differential
abundance.

LC-MS platform and acquisition parameters: see the manuscript Materials and
Methods.
