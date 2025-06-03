# PERSEO: Model-Aware Differential Expression for Omics Data

[![R](https://img.shields.io/badge/R-%3E=4.2.0-blue?style=flat&logo=R)](https://www.r-project.org/)
[![gamlss](https://img.shields.io/badge/GAMLSS-supported-lightgrey?logo=R&style=flat)](https://www.gamlss.com/)
[![tidyverse](https://img.shields.io/badge/tidyverse-compatible-brightgreen?style=flat&logo=tidyverse)](https://www.tidyverse.org/)
[![version](https://img.shields.io/badge/version-0.0.1-orange?style=flat&logo=GitHub)](https://github.com/your_username/gamlssOmics)

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

### `fit_gamlss_models()`

Fits multiple GAMLSS models to a single feature across different families.

```r
fit_gamlss_models(y, X, families = c("NBI", "GA", "NO", "PO"), timeout = 10, verbose = TRUE)
```

- `y`: numeric vector of expression/abundance values
- `X`: design matrix
- `families`: list of distribution families to test
- `timeout`: time in seconds for each model
- `verbose`: print feedback

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

## Installation (coming soon)

```r
# Not yet published
# devtools::install_github("yourusername/PERSEO")
```

---

## Contributions

Open an issue or pull request with suggestions, improvements, or feedback. This is an early-stage open science tool.

---

## License

The license will be specified soon. We welcome feedback on open-source licensing options.

Automated model selection and differential expression analysis for omics data using **GAMLSS** (Generalized Additive Models for Location, Scale and Shape). It supports overdispersed, skewed, or otherwise non-standard count distributions, allowing for better model fit and more accurate inference.

---

## Overview

`Perseo` is an R package that automates the process of:

- Inferring data type per feature (gene, metabolite, etc.)
- Selecting appropriate families of distributions from the GAMLSS framework
- Fitting multiple GAMLSS models per feature
- Selecting the best model via AIC/BIC
- Performing statistical contrasts across conditions
- Exporting results and diagnostics

The package is tailored for **high-dimensional omics datasets** (transcriptomics, metabolomics, etc.), where modeling assumptions may vary across features.

---


### fit_gamlss_models()

Fits multiple GAMLSS models to a single gene across different families.

Usage:

```r
fit_gamlss_models(y, X, families = c("NBI", "GA", "NO", "PO"), timeout = 10, verbose = TRUE)
```

Arguments:
- y: numeric vector of gene expression values
- X: design matrix
- families: list of distribution families to test
- timeout: time in seconds to allow for each model fit
- verbose: print progress or not

### select_best_model()

Selects the best-fitting model based on AIC or BIC.

Usage:

```r
select_best_model(fit_result, criterion = "AIC")
```

Arguments:
- fit_result: output from fit_gamlss_models
- criterion: "AIC" or "BIC"

### contrast_test()

Tests the specified model coefficient (e.g. tumor vs normal).

Usage:

```r
contrast_test(model, coef_name)
```

Arguments:
- model: a gamlss object
- coef_name: name of the coefficient to test

## Integrated Function: differential_expression_gamlss()

Runs the full pipeline across all genes.

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

Arguments:
- counts: genes x samples count matrix
- design: model matrix
- coef_name: coefficient to test
- p_adjust_method: method for multiple testing correction
- timeout: timeout per model (in seconds)
- verbose: print progress

This function is optimized with `future_lapply()` for parallel execution.

Memory options:

```r
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB limit
```

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

## Requirements

- R >= 4.0
- Packages: gamlss, future, future.apply, tibble, stats

## Installation (coming soon)

```r
devtools::install_github("yourusername/PERSEO")
```

## Contributions

Contributions and feedback are welcome. Please open an issue or pull request.

## Input Data

PERSEO requires two main inputs:

1. **Counts matrix (`counts`)**: A matrix-like object (e.g. `data.frame` or `matrix`) where rows are genes and columns are samples. The values should be raw integer counts, typically from RNA-Seq data.

2. **Design matrix (`design`)**: A model matrix that represents the experimental design. It is generated from a metadata `data.frame` (usually called `colData`) using the `model.matrix()` function in R. It must not contain an intercept if all variables are explicitly specified.

Example:
```r
# Sample metadata
colData <- data.frame(
  tissue_type = factor(c("Normal", "Tumor", "Tumor")),
  age = c(55, 63, 70),
  gender = factor(c("female", "female", "male"))
)

# Build design matrix
design <- model.matrix(~ tissue_type + age + gender, data = colData)
```

3. **Contrast matrix**: You can define specific contrasts (i.e. hypotheses to test) with a contrast matrix. You must ensure the coefficient you want to test exists in your design.

Example to test "Tumor vs Normal":
```r
contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
rownames(contrast_matrix) <- colnames(design)
colnames(contrast_matrix) <- "Tumor_vs_Normal"
contrast_matrix["tissue_typeTumor", 1] <- 1
```

The `coef_name` argument in the PERSEO functions should match the row name of the contrast column (e.g., `"tissue_typeTumor"`).
