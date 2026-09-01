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
