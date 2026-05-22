# lomda <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/heli-math/lomda/workflows/R-CMD-check/badge.svg)](https://github.com/heli-math/lomda/actions)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

**L**ongitudinal **O**mics **M**ultivariate **D**imension-reduction **A**nalysis

`lomda` is an R package for analysing high-dimensional longitudinal omics data
(metabolomics, proteomics, transcriptomics) measured at repeated clinical visits.
It implements a principled two-stage approach:

| Stage | Method | Purpose |
|-------|--------|---------|
| 1 | **PCA** | Compress high-dimensional omics data into a small number of interpretable PCs |
| 2 | **LMM** (random intercept) | Model PC trajectories over time; test for significant time trends |

---

## Installation

```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("heli-math/lomda")
```

---

## Quick Example

```r
library(lomda)

# Simulate a longitudinal metabolomics dataset
dat <- simulate_lomda_data(
  n_subjects  = 50,
  n_visits    = 3,
  n_features  = 50,   # omics features
  n_signal    = 5,    # features with a true time trend
  time_effect = 1,
  seed        = 2026
)

# Fit the two-stage model
fit <- lomda(dat, n_pc = 3)
print(fit)
summary(fit)

# ── Inference ────────────────────────────────
# Likelihood Ratio Test (time effect)
lomda_lrt(fit)

# Wald Test
lomda_wald(fit)

# ── Plots ────────────────────────────────────
plot(fit, type = "scores")       # PC1 vs PC2
plot(fit, type = "loadings")     # Feature loadings biplot
plot(fit, type = "variance")     # Scree plot
plot(fit, type = "trajectory")   # PC score over time
```

---

## Data Format

The input data frame must have this column layout:

| Column position | Name | Description |
|----------------|------|-------------|
| 1 | `ID` | Subject identifier |
| 2 | `visit` | Visit number (integer: 1, 2, 3, …) |
| 3 | `age` | Subject age at that visit |
| 4+ | `M1`, `M2`, … | Omics features (metabolites, proteins, etc.) |

One row per subject × visit.

---

## Statistical Model

### Stage 1 — PCA

All omics measurements (pooled across subjects and visits) are centred and
scaled, then decomposed via singular value decomposition. The top *K* PCs
explain the maximal variance in the data.

### Stage 2 — LMM

For each PC *k* and subject *i* at visit *j*:

$$
t_{ik} = \beta_0 + \beta_1 \cdot k + b_i + g_{ik}
$$

- $\beta_1$ — **time-effect slope** (primary inferential target)
- $b_i \sim N(0,\,\Sigma_b)$ — subject random intercept
- $g_{ik} \sim N(0,\,\Sigma_g)$ — residual error

### Inference

| Test | Function | Description |
|------|----------|-------------|
| Likelihood Ratio Test | `lomda_lrt()` | Full model vs. null (time dropped) |
| Wald Test | `lomda_wald()` | $\hat{\beta}_1 / \text{SE}$, Satterthwaite df via `lmerTest` |

---

## Downstream Plots

| Function | Plot |
|----------|------|
| `plot_scores()` | PC score scatter (any two PCs, coloured by metadata) |
| `plot_loadings()` | Feature loadings bar chart or biplot |
| `plot_trajectory()` | Mean PC score over visits ± CI + individual lines |
| `plot_variance_explained()` | Scree plot with cumulative variance |

All plots return `ggplot2` objects and can be further customised.

---

## Package Structure

```
lomda/
├── R/
│   ├── lomda-package.R      # Package documentation
│   ├── lomda.R              # Main two-stage fit: lomda()
│   ├── stages.R             # lomda_pca(), lomda_lmm()
│   ├── inference.R          # lomda_lrt(), lomda_wald()
│   ├── simulate.R           # simulate_lomda_data()
│   └── plots.R              # plot.lomda(), plot_scores(), …
├── tests/testthat/          # Unit tests
├── vignettes/               # Package vignette
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

---

## Dependencies

- [`lme4`](https://CRAN.R-project.org/package=lme4) — mixed model fitting
- [`lmerTest`](https://CRAN.R-project.org/package=lmerTest) — Satterthwaite df + Wald tests
- [`ggplot2`](https://CRAN.R-project.org/package=ggplot2) — visualisations
- [`ggrepel`](https://CRAN.R-project.org/package=ggrepel) — non-overlapping labels
- [`dplyr`](https://CRAN.R-project.org/package=dplyr) / [`tidyr`](https://CRAN.R-project.org/package=tidyr) — data manipulation

---

## Citation

If you use `lomda` in your research, please cite:

```
He Li, Said el Bouhaddani, Jeanine Houwing-Duistermaat (2026). lomda: Longitudinal Omics Multivariate Dimension-reduction Analysis.
R package version 0.1.0. https://github.com/heli-math/lomda
```

---

## Funding

This work was supported by a STSM Grant from COST Action CA21169, supported by COST (European Cooperation in Science and Technology).

