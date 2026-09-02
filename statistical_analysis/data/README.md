# Statistical-analysis input tables

Data tables per analysis in [`../statistical_analysis.qmd`](../statistical_analysis.qmd).
Each row is a single clone, brain, or cell; `value` is the response that is
modelled. Render the analysis report with
`quarto render statistical_analysis/statistical_analysis.qmd`.

| file | manuscript panel(s) | replication unit | response (`value`) |
|---|---|---|---|
| `fig_01b_lineage_size.csv` | Fig 1B | one NSC clone | cells per NSC lineage |
| `fig_01g_heteroplasmy.csv` | Fig 1G | one neuroblast | EdU+ progeny per neuroblast |
| `fig_01p_heteroplasmic_clones.csv` | Fig 1P | one NSC clone | EdU+ cells per clone |
| `fig_03e_ucp_rescue.csv` | Fig 3E (UCP1 rows) + Supp Fig 7J (UCP2 rows) | one brain (VNC) | NSC mitotic index (% pH3+Dpn+ / Dpn+) |
| `fig_03g_lbnox_rescue.csv` | Fig 3G | one brain (VNC) | NSC mitotic index |
| `fig_04c_inx2_edu.csv` | Fig 4C | one NSC clone | EdU+ cells per clone |
| `fig_s01i_j_edu.csv` | Fig S1I–J | one NSC clone | EdU+ cells per clone |
| `fig_s07c_pampk.csv` | Fig S7C | one NSC | normalised pAMPK intensity (a.u.) |
| `fig_s07f_sod2_rescue.csv` | Fig S7F | one brain (VNC) | NSC mitotic index |
| `edu_glial_screen.csv` | Fig 4D, 4L, S8D, S9A + shared controls | one NSC clone | EdU+ cells per clone |

`edu_glial_screen.csv` is a combined data table as we fit a single pooled
model across the whole screen (joint Benjamini–Hochberg FDR); the dedicated
rescue experiments (Fig 4K, 5E, S9A, S9I) are subsets of it selected by
genotype.

## Columns

| column | meaning |
|---|---|
| `manuscript_panel` | manuscript figure panel the row belongs to |
| `measurement`, `unit` | what `value` is, and its unit |
| `value` | the measured response modelled in the report |
| `genotype` | full genotype label; the glial screen (`edu_glial_screen.csv`) uses a `<NSC background> x <glial transgene>` convention |
| `rnai` | NSC RNAi condition / background (`Control`, `Complex I`, `Complex V`) |
| `rescue_transgene` | rescue construct (`GFP`/`UCP1`/`UCP2`; `None`/`SOD2`; `nlsLbNox`/`mitoLbNox`; …) |
| `panel_side` | mosaic system: `left` = most NSCs affected, `right` = few NSCs affected |
| `screen_type` | `control` / `transporter` / `enzyme` / … (glial screen only) |
| `in_combined_screen` | `TRUE` if the row enters the pooled joint-FDR model |
| `collection_date` | imaging-session id (`session_01`, …) or `undated`; used only for the shared-session random effect and session counts |
| `mito_type` | `MitoGFP` (control) or `MitoXho1` (Fig 1P) |
| `mtdna_background` | `Homoplasmic WT` or `Heteroplasmic` (Fig 1P) |
| `glial_gene` | glial knockdown target (Fig 4C) |
| `heteroplasmy_pct`, `heteroplasmy_fraction` | in situ % / fraction mutant mtDNA (Fig 1G predictor) |
| `d_mel`, `d_yak` | *D. melanogaster* / *D. yakuba* mtDNA smFISH spot counts (Fig 1G) |
| `dpn_positive`, `ph3_positive` | Dpn+ and pH3+Dpn+ NSC counts behind the mitotic index |
| `brain_label` | brain id `B1`–`B4` (Fig 1G) |
| `replicate_id` | biological-replicate (brain) code `b01`, `b02`, …; the random-effect grouping |
