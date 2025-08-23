<div align="left">
  <h1 style="display: inline-block;">PERSEO: Model-Aware Differential Expression for Omics Data</h1>
  <img src="assets/logo_perseo.png" alt="Perseo Logo" width="150" align="right">
  <br/>

  <a href="https://www.r-project.org/">
    <img src="https://img.shields.io/badge/R-%3E=4.2.0-blue?style=flat&logo=R" alt="R">
  </a>
  <a href="https://www.gamlss.com/">
    <img src="https://img.shields.io/badge/GAMLSS-supported-lightgrey?logo=R&style=flat" alt="GAMLSS">
  </a>
  <a href="https://www.tidyverse.org/">
    <img src="https://img.shields.io/badge/tidyverse-compatible-brightgreen?style=flat&logo=tidyverse" alt="tidyverse">
  </a>
  <a href="https://github.com/OPSergio/Perseo">
    <img src="https://img.shields.io/badge/version-0.0.1-orange?style=flat&logo=GitHub" alt="version">
  </a>
</div>


---

Automated model selection and differential expression analysis for omics data using **GAMLSS** (Generalized Additive Models for Location, Scale and Shape). It supports overdispersed, skewed, or otherwise non-standard distributions, allowing for better model fit and more accurate inference.

---

## Table of Contents

- [Why PERSEO?](#why-perseo)
- [Installation](#installation)
- [Data Requirements](#data-requirements)
- [Workflow Overview](#workflow-overview)
- [Function Reference](#function-reference)
  - [find_families()](#find_families)
  - [fit_gamlss_models()](#fit_gamlss_models)
  - [Utilities](#utilities)
- [End‑to‑End Example](#end-to-end-example)
- [Interpretation & Tips](#interpretation--tips)
- [Troubleshooting](#troubleshooting)
- [Roadmap](#roadmap)
- [License](#license)

---

## Why PERSEO?

Classical DE pipelines assume a single likelihood family for all features (often Normal or
Negative Binomial). In high‑dimensional omics this is brittle: supports and shapes **vary by
feature** (counts vs positive continuous vs real‑valued; skewness/zero‑inflation/overdispersion).
PERSEO tackles this by:

- Testing **multiple GAMLSS families per feature**.
- Using **strict transforms** + a **common mask** so all families are compared on the **same observations**.
- Applying **Jacobian correction** so ICs (AIC/BIC/GAIC) are comparable on the original scale.
- Selecting the **best model** per feature and extracting **per‑term stats** from `summary(fit, what="mu")`.
- Scaling with **future.apply**.

---

## Installation (coming soon)

```r
# Not published yet
# remotes::install_github("OPSergio/Perseo")
```

---

## Data Requirements

- **counts_matrix**: numeric matrix **features × samples** (rows = features, cols = samples).
  Values must be numeric, typically non‑negative (RNA‑seq counts) or continuous (normalized assays).
- **design_matrix**: model matrix with **nrow = ncol(counts_matrix)**. Build with `model.matrix()`.
  It must encode your covariates of interest (e.g., tissue type, age, batch, sex…).

Example:

```r
metadata$tissue_type <- stats::relevel(factor(metadata$tissue_type), ref = "Normal")

design <- model.matrix(
  ~ tissue_type + age_at_diagnosis + gender + initial_weight +
    laterality + ajcc_pathologic_stage,
  data = metadata
)

# Align counts to design rows
counts_matrix <- as.matrix(counts[, rownames(design)])
mode(counts_matrix) <- "numeric"
```

---

## Functions

### `find_families()`

Identify frequent, plausible GAMLSS families via bootstrap sampling. For each sampled feature,
it fits **intercept‑only** models across candidate families using **strict transforms**, a **common mask**,
and **Jacobian‑corrected** IC.

**Usage**

```r
ff <- find_families(
  counts_matrix,
  n_genes          = 200,
  n_boot           = 10,
  top_n            = 4,
  families         = NULL,     # if NULL, a default panel is used
  criterion        = "BIC",    # "GAIC" or "AIC"
  gaic_k           = NULL,     # if GAIC and NULL -> log(n_valid)
  min_n            = 5,
  seed             = NULL,
  group_by_support = FALSE,    # do NOT hard-filter by inferred support
  binom_bd         = NULL,     # optional denominator for BI/BB
  filter_beta_inflated = TRUE,
  thr_zero         = 0.005,
  thr_one          = 0.005,
  verbose          = TRUE
)
```

**Arguments (key)**

- `counts_matrix`: numeric matrix (features × samples).
- `families`: character vector of families to test; if `NULL`, uses a default panel covering
  unit interval, counts (incl. zero‑inflated), positive continuous, and real‑valued.
- `criterion`: `"BIC"` (recommended), `"GAIC"`, `"AIC"`.
- `group_by_support`: if `TRUE`, restricts candidates by inferred empirical support; if `FALSE`
  (recommended), the support is kept as **metadata** only (no hard gating).
- `min_n`: minimum valid samples after the common mask.

**Returns**

A list with:
- `top_families_overall`: character vector (length = `top_n`).
- `top_families_by_support`: list with top families per support class (metadata).
- `freq_table_overall`, `prop_table_overall`: frequency/proportion tables.
- `freq_by_support`, `prop_by_support`: by support class.
- `sampled_results`: tibble with one row per feature×bootstrap pull, including
  `feature`, `family`, `skipped`, `n_valid`, `support`, `bootstrap`.

**Example**

```r
candidate_pool <- c("PO","NBI","ZIP","ZINBI","GA","GG","IG","LOGNO","NO","TF")

ff <- find_families(
  counts_matrix    = counts_matrix,
  n_genes          = 300,
  n_boot           = 5,
  top_n            = 6,
  families         = candidate_pool,
  criterion        = "BIC",
  min_n            = 20,
  seed             = 123,
  group_by_support = FALSE
)

ff$top_families_overall
#> [1] "GG"    "LOGNO" "NBI"   "GA"    "TF"    "NO"

ff$freq_table_overall
#> GG   LOGNO   NBI   GA   TF  NO
#> 142     91    64   28    7   2

ff$prop_table_overall
#>       GG    LOGNO     NBI      GA      TF      NO
#> 0.473333 0.303333 0.213333 0.093333 0.023333 0.006667

ff$top_families_by_support
#> $count
#> [1] "NBI" "PO"  "ZINBI"
#>
#> $unit
#> [1] "BE" "BEINF"
#>
#> $positive
#> [1] "GG"    "LOGNO" "GA"    "IG"
#>
#> $real
#> [1] "TF" "NO"

head(ff$sampled_results)
# A tibble: 6 × 6
#  bootstrap feature         family skipped n_valid support
#      <int> <chr>           <chr>  <lgl>      <int> <chr>
#1         1 GENE_0001       GG     FALSE        213 positive
#2         1 GENE_0002       NBI    FALSE        208 count
#3         1 GENE_0003       LOGNO  FALSE        214 positive
#4         1 GENE_0004       GA     FALSE        210 positive
#5         1 GENE_0005       TF     FALSE        216 real
#6         1 GENE_0006       GG     FALSE        211 positive

with(ff$sampled_results[!ff$sampled_results$skipped, ], table(support, family))
```

---

### `fit_gamlss_models()`

Fit multiple GAMLSS families per feature on a **common set of observations**, apply
**Jacobian correction**, select the **best family** by IC, and return **per‑term statistics**
from `summary(fit, what="mu")`. This “no‑VCOV mode” avoids fragile covariance extraction and
is robust across families.

**Usage**

```r
options(progressr.enable = TRUE)
if (requireNamespace("cli", quietly = TRUE)) progressr::handlers("cli")

fit <- fit_gamlss_models(
  counts_matrix      = counts_matrix,
  design_matrix      = design,
  candidate_families = candidate_families,  # e.g., ff$top_families_overall
  criterion          = "BIC",
  gaic_k             = NULL,
  min_n              = 20,
  p_adjust           = "BH",
  workers            = max(1, future::availableCores() - 1),
  show_progress      = TRUE,
  progress_label     = "Fitting genes"
)
```

**Arguments**

- `counts_matrix`, `design_matrix`: see [Data Requirements](#data-requirements).
- `candidate_families`: character vector of GAMLSS families to compare.
- `criterion`: `"BIC"` (recommended), `"GAIC"`, `"AIC"`.
- `min_n`: minimum valid samples after common mask.
- `p_adjust`: method passed to `p.adjust()` to compute **FDR per term across features** (default `"BH"`).
- `workers`: number of parallel workers (`future::multisession`).
- `show_progress`, `progress_label`: progressr control.
- `contrast_matrix`: **ignored** in no‑VCOV mode (reserved for future arbitrary contrasts).

**Outputs**

Returns a list with two tibbles:

- `selection` (one row per feature)
  - `feature`: feature ID (rownames of `counts_matrix`).
  - `best_family`: winning family by corrected IC.
  - `n_valid_obs`: number of samples used after the common mask.
  - `ic_value`: **Jacobian‑corrected** IC for the best family.

- `results` (multiple rows per feature; per‑term statistics on `mu`)
  - `feature`: feature ID.
  - `term`: coefficient name from the design (e.g., `tissue_typeTumor`, `age_at_diagnosis`).
  - `effect`: coefficient estimate (on the link of `mu`).
  - `se`: standard error (as reported by `summary()`; may be `NA` if not available).
  - `stat`: Wald statistic (often `t`).
  - `pval`: p‑value of the term.
  - `padj`: FDR (adjusted across features **by term**).

Example output: 

```r

head(fit$selection)
# A tibble: 6 × 4
#  feature            best_family n_valid_obs ic_value
#  <chr>              <chr>             <int>    <dbl>
#1 ENSG00000153002.12 GG                 1176   19590.
#2 ENSG00000210082.2  GG                 1176   32447.
#3 ENSG00000276168.1  GG                 1176   19424.
#4 ENSG00000211896.7  GG                 1176   28663.
#5 ENSG00000198804.2  GG                 1176   32901.
#6 ENSG00000108821.14 GG                 1176   32407.

head(fit$results)
# A tibble: 6 × 7
#  feature            term                 effect      se    stat     pval     padj
#  <chr>              <chr>                 <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
#1 ENSG00000153002.12 X.Intercept.        6.89e+0 6.45e-1 10.7    1.77e-25 1.85e-25
#2 ENSG00000153002.12 tissue_typeTumor   -7.72e-1 4.05e-1 -1.90   5.71e- 2 7.11e- 2
#3 ENSG00000153002.12 age_at_diagnosis    1.05e-6 2.44e-5  0.0429 9.66e- 1 9.79e- 1
#4 ENSG00000153002.12 gendermale          1.03e+0 1.17e+0  0.877  3.81e- 1 7.56e- 1
#5 ENSG00000153002.12 X.Intercept..1      1.39e+0 2.02e-2 69.1    0        0
#6 ENSG00000153002.12 X.Intercept..2      9.02e-2 1.11e-2  8.13   1.08e-15 2.15e-15

```

**Notes**

- `fit_gamlss_models()` removes an explicit `"(Intercept)"` column from `design_matrix` to
  avoid double intercepts when using `y ~ .`.
- The **common mask** combines family domain validity and `complete.cases(design)` to prevent
  optimizer failures (e.g., `sigma` working vector issues).

---

### Utilities

- `transform_for_family_strict(y, fam, eps = 1e-6, allow_eps = TRUE)`  
  Returns list with `y` (transformed), `mask` (valid obs), `logJ_per_obs`, and `meta`.
- `inverse_transform(z, meta)`  
  Invert strict transforms back to original scale.
- `family_groups()`  
  Family names grouped by theoretical support (count / unit / positive / real).

---

## Contributions

Open an issue or pull request with suggestions, improvements, or feedback. This is an early-stage open science tool.

---

## License

The license will be specified soon. We welcome feedback on open-source licensing options.
