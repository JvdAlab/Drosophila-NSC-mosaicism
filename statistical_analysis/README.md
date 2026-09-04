# Statistical Analysis

Clone-, brain-, and cell-level mixed-model statistics for the confocal /
genetics figure panels of the manuscript: linear and generalized
least-squares models (`stats::lm`, `nlme::gls`), linear mixed-effects models
(`lme4`, `lmerTest`), Type III ANOVA (`car`, `lmerTest`), estimated marginal
means and multiplicity-corrected contrasts (`emmeans`), and simulation-based
residual diagnostics (`DHARMa`).

Covers Fig 1B/1G/1P, 3E/3G, 4C/4K, 5E and Supp Fig S1I-J/S7C/S7F/S7J/S9A/S9I,
plus the pooled glial RNAi screen behind Fig 4D/4L/S8D/S9A (one joint model,
Benjamini-Hochberg FDR across all knockdown contrasts).

## Contents

| file | what it is |
|:---|:---|
| [`statistical_analysis.qmd`](./statistical_analysis.qmd) | the analysis document |
| [`statistical_analysis.html`](./statistical_analysis.html) | committed HTML render |
| [`data/`](./data/) | data tables for each analysis (one row per clone / brain / cell); see [`data/README.md`](./data/README.md) |
| [`environment.yml`](./environment.yml) | environment specification |

## Rendering

From the repository root:

```bash
conda env create -f statistical_analysis/environment.yml
conda activate statistical_analysis
quarto render statistical_analysis/statistical_analysis.qmd
```

Reads `data/*.csv`; writes result tables to `results/statistical_analysis/`
(git-ignored, regenerated each render) and the self-contained `.html`.

### Result tables

Machine-readable copies of every result table, written at render time by the
`write-all-results` chunk. Values are the raw model outputs (not the
rounded/annotated display versions).

| file | section | contents |
|---|---|---|
| `fig_01b_results.csv` | Fig 1B | lineage-size contrasts (Complex I / Complex V vs Control), per panel side |
| `fig_01g_anova.csv` | Fig 1G | Type III ANOVA (`progeny ~ heteroplasmy + brain`) |
| `fig_01p_results.csv` | Fig 1P | MitoXho1 vs MitoGFP within each mtDNA background |
| `fig_01p_anova.csv` | Fig 1P | Type III ANOVA (2 × 2 factorial) |
| `fig_03e_within_condition_results.csv` | Fig 3E / S7J | UCP1 / UCP2 vs GFP within each RNAi condition |
| `fig_03e_diff_in_diff_results.csv` | Fig 3E / S7J | difference-in-differences (Complex V vs Complex I) |
| `fig_03g_within_condition_results.csv` | Fig 3G | nlsLbNox / mitoLbNox vs Control within each RNAi condition |
| `fig_03g_diff_in_diff_results.csv` | Fig 3G | difference-in-differences (Complex I vs Control) |
| `fig_04c_results.csv` | Fig 4C | the three glial-Inx2 contrasts |
| `combined_screen_results.csv` | Fig 4D & 4L | per-knockdown contrast vs `ND42 x Luciferase` for the pooled glial RNAi screen, joint BH-FDR; keyed by `manuscript_panel` |
| `combined_phenotype_contrast.csv` | Fig 4D & 4L | `ND42 x Luciferase` vs `mCherry x Luciferase` phenotype contrast (uncorrected) |
| `fig_04k_left_results.csv` | Fig 4K | Eaat1 KD vs wild-type (left-panel LME) |
| `fig_04k_right_results.csv` | Fig 4K | CI-background Eaat1 KD + Eaat2-GFP rescue contrasts (right-panel GLS) |
| `fig_05e_ldh_alone.csv` | Fig 5E | Ldh KD vs wild-type |
| `fig_05e_ldh_worsens.csv` | Fig 5E | Ldh KD + CI vs CI alone |
| `fig_05e_rescue_vs_baseline.csv` | Fig 5E | each rescue construct vs `ND42 x Ldh` |
| `fig_05e_rescue_vs_control.csv` | Fig 5E | each rescue construct vs `ND42 x Luciferase` |
| `fig_s01i_j_results.csv` | Fig S1I-J | proliferation-rate contrasts (Complex I / Complex V vs Control), per panel side |
| `fig_s07c_results.csv` | Fig S7C | pairwise pAMPK ratio contrasts |
| `fig_s07f_results.csv` | Fig S7F | SOD2 vs no SOD2 within each RNAi background |
| `fig_s09a_results.csv` | Fig S9A | NSC-background simple effects within each glial knockdown |
| `fig_s09a_vs_control_results.csv` | Fig S9A | glial knockdown vs pooled control within each NSC background |
| `fig_s09i_results.csv` | Fig S9I | the three glial-LbNox contrasts |

## Analysis details

- **Replication unit** is the brain (VNC): where a brain contributes multiple
  clone- or cell-level measurements, brain identity is a random intercept
  (`(1 | brain_uid)`); where each brain yields a single value (per-brain
  mitotic index), an ordinary linear model is used.
- **Diagnostics:** scaled residuals from 1,000 DHARMa simulations are tested
  for uniform distribution (KS), dispersion, and outliers; a transform or a
  heteroscedastic GLS (`nlme::gls` with `varIdent`) is used where the Gaussian
  assumption is violated.
- **Multiplicity:** small pre-specified contrast families use a
  multivariate-*t* adjustment; the large glial RNAi screen uses
  Benjamini-Hochberg FDR.
- **Reproducibility:** `set.seed(42)`, a fixed DHARMa seed, and
  `emmeans::emm_options(rng.seed = 42)` so the multivariate-*t* p-values are
  deterministic across renders.
