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

## Overview

`Perseo` is an R package that automates the process of:

- Inferring the appropriate distribution per feature (gene, protein, metabolite, etc.)
- Fitting GAMLSS models across a variety of families
- Selecting the best model per feature via AIC/BIC
- Performing hypothesis testing via user-defined contrasts
- Returning results in tidy tabular format, with optional p-value adjustment

Tailored for high-dimensional **omics datasets** where modeling assumptions vary across features.

---

## Functions

### `find_families()`

This function helps identify the best-fitting GAMLSS families for a representative subset of features (e.g, genes) using AIC. It is useful to explore the most appropiate distributions to be used in subsequent modeling

```r
top_families <- find_families(
  counts_matrix = counts_filtered,
  n_genes = 100,
  top_n = 4,
  verbose = TRUE
  )
```
Example output 

```r
Running model fitting on 100 genes using parallel execution...

===== Summary Report =====
Genes analyzed: 100
Genes skipped (all 0s): 5
Genes successfully fitted: 95
Most frequent family: TF (89 genes)
 Top 4 families:
fitted
 TF NBI  BI 
 89   5   1 
```
---

### `fit_gamlss_models()`

This function fits multiple GAMLSS models across genes (features) and selects the best-fitting distribution per gene based on user-defined criteria (AIC, BIC, GAIC, or logLik). It is designed to be scalable and memory-efficient by only storing summary statistics and diagnostics, not full models.

```r
X <- model.matrix(~ group, data = sample_metadata)

# Define candidate families (can be from top_families or fixed list)
families <- c("TF", "NBI", "BI", "NO")

# Run GAMLSS fitting
results_tbl <- fit_gamlss_models(
  counts_matrix = counts_filtered,
  X = X,
  families = families,
  criterion = "GAIC",
  timeout = 5,
  verbose = TRUE
)

```

- `counts_matrix`: matrix of expression/abundance values
- `X`: design matrix
- `families`: list of distribution families to test
- `timeout`: time in seconds for each model
- `verbose`: print feedback

Output

```r
| gene  | best\_family | AIC   | BIC   | GAIC3 | logLik | df  | ks\_p | skewness | kurtosis |
| ----- | ------------ | ----- | ----- | ----- | ------ | --- | ----- | -------- | -------- |
| BRCA1 | TF           | 512.3 | 520.1 | 517.9 | -252.1 | 4   | 0.72  | 0.03     | -0.89    |
| TP53  | NBI          | 438.7 | 447.4 | 443.2 | -215.3 | 5   | 0.59  | -0.18    | 1.02     |
| ...   | ...          | ...   | ...   | ...   | ...    | ... | ...   | ...      | ...      |
```

When `verbose = TRUE`, a final report is printed: 

```
── Fitting GAMLSS models across genes ──────────────────────────────────────────────
Models fitted for 95 of 100 genes
5 genes skipped due to all-zero counts
Most frequent family: TF (89 genes)

```

---

### `select_best_model()`

Selects the best-fitting model per feature.

```r
select_best_model(fit_result, criterion = "AIC")
```

- `fit_result`: list of model fits
- `criterion`: `"AIC"` or `"BIC"`

---

### `contrast_test()`

Performs a Wald test for a coefficient in a GAMLSS model.

```r
contrast_test(model, coef_name)
```

- `model`: fitted GAMLSS object
- `coef_name`: the name of the coefficient (must match the design)

---

## Integrated Function

### `differential_expression_gamlss()`

Runs the full pipeline efficiently across all features:

```r
differential_expression_gamlss(
  counts,
  design,
  families = c("NBI", "GA", "NO", "PO"),
  coef_name,
  criterion = "AIC",
  p_adjust_method = "BH",
  timeout = 10,
  verbose = TRUE
)
```

- `counts`: features × samples numeric matrix
- `design`: model matrix
- `coef_name`: name of coefficient to test
- `p_adjust_method`: `"BH"`, `"bonferroni"`, etc.
- `timeout`: per-feature timeout
- `verbose`: show progress

Optimized with `future_lapply()` for multicore parallelism. If needed:

```r
options(future.globals.maxSize = 2 * 1024^3)  # 2GB limit
```

---

## Input Requirements

PERSEO requires:

### 1. Counts or abundance matrix

A numeric matrix-like object (e.g., `matrix` or `data.frame`) where rows are features and columns are samples. Entries must be **non-negative**. It supports transcriptomics, metabolomics, proteomics, etc.

### 2. Design matrix

The experimental design should be encoded as a model matrix (no intercept if all variables are explicitly encoded):

```r
metadata <- data.frame(
  tissue_type = factor(c("Normal", "Tumor", "Tumor")),
  age = c(55, 63, 70),
  gender = factor(c("female", "female", "male"))
)

design <- model.matrix(~ tissue_type + age + gender, data = metadata)
```

### 3. Contrast matrix (optional)

To test a contrast (e.g., Tumor vs Normal):

```r
contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
rownames(contrast_matrix) <- colnames(design)
colnames(contrast_matrix) <- "Tumor_vs_Normal"
contrast_matrix["tissue_typeTumor", 1] <- 1
```

Then pass `coef_name = "tissue_typeTumor"`.

---

## Example

```r
library(PERSEO)

results <- differential_expression_gamlss(
  counts = counts,
  design = design,
  coef_name = "tissue_typeTumor"
)

head(results)
```

---

## Installation (coming soon)

```r
# Not yet published
# devtools::install_github("OPSergio/Perseo")
```

---

## Contributions

Open an issue or pull request with suggestions, improvements, or feedback. This is an early-stage open science tool.

---

## License

The license will be specified soon. We welcome feedback on open-source licensing options.
