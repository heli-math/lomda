# lomda <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/heli-math/lomda/workflows/R-CMD-check/badge.svg)](https://github.com/heli-math/lomda/actions)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

**L**ongitudinal **O**mics **M**ultivariate **D**imension-reduction
**A**nalysis

`lomda` is an R package for analysing high-dimensional longitudinal omics data,
including metabolomics, proteomics, transcriptomics, and related repeated
measurement data.

The package uses a two-stage workflow. Stage 1 is shared by all analyses:
Principal Component Analysis (PCA) reduces the high-dimensional omics block to
a small number of PC scores. Stage 2 models those estimated PC scores using the
time scale chosen by the user.

| Stage | Method | Purpose |
|-------|--------|---------|
| 1 | PCA | Compress correlated omics features into PC scores |
| 2a | LMM over `visit` | Model visit-indexed PC trajectories with random intercepts |
| 2b | FDA over `age` | Smooth age-indexed PC trajectories and predict age-specific scores |

## Installation

```r
# install.packages("devtools")
devtools::install_github("heli-math/lomda")
```

## Quick Example

```r
library(lomda)

dat <- simulate_lomda_data(
  n_subjects         = 50,
  n_visits           = 3,
  n_features         = 50,
  n_pc               = 3,
  age_effect         = 0.3,
  missing_visit_prob = 0.3,
  seed               = 2026
)

# Stage 2a: PCA + LMM over visit
fit_visit <- lomda(dat, n_pc = 3, time = "visit")
print(fit_visit)
summary(fit_visit)

lomda_lrt(fit_visit)
lomda_wald(fit_visit)

plot(fit_visit, type = "scores")
plot(fit_visit, type = "trajectory")
plot(fit_visit, type = "all")

# Stage 2b: PCA + FDA over age
fit_age <- lomda(dat, n_pc = 3, time = "age", method_fda = "spline")
summary(fit_age)

predict(fit_age, ages = c(35, 45, 55, 65))
lomda_important(fit_age, n_top = 10)

plot(fit_age, type = "variance")
plot_fda_trajectory(fit_age, pcs = 1:3, n_subjects = 50)
plot_fda_importance(fit_age)
plot(fit_age, type = "all", pcs = 1:3, n_subjects = 50)
```

## Choosing the Stage 2 Model

Use one public front door:

```r
lomda(dat, time = "visit")
lomda(dat, time = "age", method_fda = "spline")
```

For visit-indexed analysis, `lomda()` fits PCA followed by LMM:

```r
fit_visit <- lomda(dat, time = "visit")
```

For age-indexed analysis, `lomda()` fits PCA followed by FDA:

```r
fit_age <- lomda(dat, time = "age", method_fda = "spline")
fit_age <- lomda(dat, time = "age", method_fda = "face")
```

The lower-level `lomda_fda()` function remains available for advanced users and
for backward-compatible scripts.

## Data Format

The input data frame should have one row per subject-visit observation. LOMDA
expects these metadata columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | `ID` | Subject identifier |
| 2 | `visit` | Visit number |
| 3 | `age` | Subject age at that visit |
| 4+ | omics features | Metabolites, proteins, genes, or other omics variables |

One row is one subject-visit observation. Rename your metadata columns to
`ID`, `visit`, and `age` before calling `lomda()`.

## Statistical Model

### Stage 1: PCA

All omics measurements are pooled across subjects and visits, then centred and
optionally scaled. PCA maps the omics block to PC scores. These PC scores are
the estimated low-dimensional molecular profiles used in Stage 2.

### Stage 2a: LMM Over Visit

When `time = "visit"`, each PC score is modelled with a linear mixed model with
a subject-specific random intercept. This is useful for formal testing of
visit-related change.

Use:

```r
lomda_lrt(fit_visit)
lomda_wald(fit_visit)
```

### Stage 2b: FDA Over Age

When `time = "age"`, each PC score is smoothed as a function of age. This is
useful when biological age is the scientific time scale and the goal is to
describe, predict, or visualize smooth molecular trajectories.

Use:

```r
predict(fit_age, ages = c(40, 50, 60))
predict(fit_age, newdata = dat[1:5, ], type = "deviation")
lomda_important(fit_age)
```

## Downstream Plots

| Function | Plot |
|----------|------|
| `plot_scores()` | PC score scatter plot |
| `plot_loadings()` | Feature loading bar chart or biplot |
| `plot_trajectory()` | Mean PC trajectory over visit |
| `plot_variance_explained()` | Scree plot |
| `plot_fda_trajectory()` | Age-smoothed PC score trajectory |
| `plot_fda_importance()` | Important features for an FDA-smoothed PC |

For large age-indexed datasets, limit the number of drawn participants:

```r
plot_fda_trajectory(fit_age, pcs = 1:3, n_subjects = 100)
```

All plot functions return `ggplot2` objects and can be customized.

## Package Structure

```text
lomda/
|-- R/
|   |-- lomda-package.R
|   |-- lomda.R
|   |-- stages.R
|   |-- fda.R
|   |-- inference.R
|   |-- simulate.R
|   `-- plots.R
|-- tests/testthat/
|-- vignettes/
|-- DESCRIPTION
|-- NAMESPACE
`-- README.md
```

## Citation

If you use `lomda` in your research, please cite:

```text
He Li, Said el Bouhaddani, Jeanine Houwing-Duistermaat (2026).
lomda: Longitudinal Omics Multivariate Dimension-reduction Analysis.
R package version 0.1.0. https://github.com/heli-math/lomda
```

## Funding

This work was supported by a STSM Grant from COST Action CA21169, supported by
COST (European Cooperation in Science and Technology).
