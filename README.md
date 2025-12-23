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
- [Quick Start](#quick-start)
- [Function Reference](#function-reference)
  - [run_perseo() - High-Level Orchestration](#run_perseo---high-level-orchestration)
  - [find_families() - Family Selection](#find_families---family-selection)
  - [fit_gamlss_models() - Differential Expression](#fit_gamlss_models---differential-expression)
  - [Utilities](#utilities)
- [Advanced Usage](#advanced-usage)
  - [Transformation Modes: Strict vs Safe](#transformation-modes-strict-vs-safe)
- [Interpretation & Tips](#interpretation--tips)
- [Troubleshooting](#troubleshooting)
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

## Quick Start

**PERSEO** provides a high-level orchestration function `run_perseo()` that handles the complete workflow in a single call:

```r
library(perseo)

# Prepare your data
metadata$condition <- relevel(factor(metadata$condition), ref = "Control")
design <- model.matrix(~ condition + age + batch, data = metadata)
counts_matrix <- as.matrix(counts[, rownames(design)])
mode(counts_matrix) <- "numeric"

# Run complete PERSEO pipeline
results <- run_perseo(
  counts_matrix = counts_matrix,
  design_matrix = design,
  n_genes = 200,           # Features to sample for family selection
  n_boot = 10,             # Bootstrap iterations
  top_n = 4,               # Top families to select
  criterion = "BIC",       # Model selection criterion
  min_n = 20,              # Minimum valid observations
  p_adjust_method = "BH", # Multiple testing correction
  verbose = TRUE,          # Show progress messages
  seed = 123               # Reproducibility
)

# View results summary
print(results)

# Access components
head(results$differential_expression$results)  # DE results with FDR
head(results$differential_expression$selection) # Best family per feature
results$family_selection$top_families_overall   # Selected families
results$summary                                  # Execution summary

# Filter significant results (FDR < 0.05)
sig_results <- results$differential_expression$results %>%
  filter(p_adj < 0.05, grepl("^condition", term))
```

**What `run_perseo()` does:**

1. **Family Selection**: Bootstraps a subset of features to identify the most frequently selected GAMLSS families
2. **Differential Expression**: Fits selected families to all features and picks the best model per feature
3. **Multiple Testing Correction**: Applies global FDR adjustment across all tests

**Returns** a `perseo_results` object with three components:
- `$family_selection`: Family selection bootstrap results
- `$differential_expression`: DE models with adjusted p-values
- `$summary`: Execution metadata and status

---

## Function Reference

### `run_perseo()` - High-Level Orchestration

One-function workflow combining family selection, differential expression, and p-value adjustment.

**Usage**

```r
results <- run_perseo(
  counts_matrix,
  design_matrix,
  n_genes = 200,
  n_boot = 10,
  top_n = 4,
  families = NULL,
  group_by_support = TRUE,
  criterion = "GAIC",
  gaic_k = NULL,
  min_n = 5,
  binom_bd = NULL,
  filter_beta_inflated = TRUE,
  p_adjust_method = "BH",
  verbose = TRUE,
  seed = NULL
)
```

**Key Arguments**

- `counts_matrix`: numeric matrix (features × samples)
- `design_matrix`: model matrix with nrow = ncol(counts_matrix)
- `n_genes`: number of features to sample for family selection
- `n_boot`: bootstrap iterations for family selection
- `top_n`: number of top families to select
- `families`: candidate families (NULL = use defaults)
- `criterion`: "GAIC", "BIC", or "AIC"
- `p_adjust_method`: "BH", "bonferroni", "holm", etc.
- `verbose`: show progress messages
- `seed`: for reproducibility

**Returns**

`perseo_results` S3 object with:
- `family_selection`: output from `find_families()`
- `differential_expression`: output from `fit_gamlss_models()` with global p-value adjustment
- `summary`: list with execution metadata

---

### `find_families()` - Family Selection

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

### `fit_gamlss_models()` - Differential Expression

Fit multiple GAMLSS families per feature on a **common set of observations**, apply
**Jacobian correction**, select the **best family** by IC, and return **per‑term statistics**
from `summary(fit, what="mu")`. This "no‑VCOV mode" avoids fragile covariance extraction and
is robust across families.

> **Note**: When using `run_perseo()`, this function is called automatically with selected families.
> Use directly only for advanced workflows or when you already know which families to test.

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

### Working with Multi-Level Factors

When your design includes a factor with >2 levels (e.g., `Status` with levels Healthy, Mild, Severe), the default model output provides coefficients for each level **relative to the reference level** (treatment coding). To obtain **custom contrasts** between levels (e.g., Mild-Severe, or average of Mild+Severe vs Healthy), pass a **contrast matrix** `makeContrasts`:

**Example: Custom Contrast Matrix**

```r
# Simulate data with 3-level Status factor
metadata <- data.frame(
  Status = factor(rep(c("Healthy", "Mild", "Severe"), each = 20)),
  Age = rnorm(60, 50, 10)
)

# Create design matrix (Healthy is reference by default)
design <- model.matrix(~ Status + Age, data = metadata)
colnames(design)
#> [1] "(Intercept)" "StatusMild"  "StatusSevere" "Age"

# Define contrasts of interest
# Each row is a linear combination of coefficients
C <- rbind(
  "Mild_vs_Severe" = c(0, 1, -1, 0),      # StatusMild - StatusSevere
  "MS_avg_vs_Healthy" = c(0, 0.5, 0.5, 0) # (Mild + Severe)/2 - Healthy
)
colnames(C) <- colnames(design)

# Fit with custom contrasts
fit <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,
  candidate_families = c("NBI", "GG", "LOGNO"),
  contrast_matrix = C,  # Pass contrast matrix directly
  workers = 4
)

# Standard coefficient results (Mild-Healthy, Severe-Healthy, Age effect)
head(fit$results)
# A tibble: 6 × 7
#  feature   term          effect     se  stat    pval    padj
#  <chr>     <chr>          <dbl>  <dbl> <dbl>   <dbl>   <dbl>
#1 GENE_001  StatusMild     0.234  0.089  2.63 0.00856 0.0234
#2 GENE_001  StatusSevere   0.567  0.092  6.16 7.3e-10 2.1e-9
#3 GENE_001  Age            0.012  0.005  2.40 0.0164  0.0352

# Custom contrasts (Mild-Severe, avg of Mild+Severe vs Healthy)
head(fit$contrasts)
# A tibble: 6 × 8
#  feature  family contrast             estimate    se     z p_value  p_adj
#  <chr>    <chr>  <chr>                   <dbl> <dbl> <dbl>   <dbl>  <dbl>
#1 GENE_001 NBI    Mild_vs_Severe         -0.333 0.095 -3.51 0.00045 0.00128
#2 GENE_001 NBI    MS_avg_vs_Healthy       0.401 0.064  6.27 3.6e-10 1.1e-9
#3 GENE_002 GG     Mild_vs_Severe         -0.078 0.088 -0.89 0.375   0.523
#4 GENE_002 GG     MS_avg_vs_Healthy       0.084 0.061  1.38 0.168   0.287
```

**Outputs**

When `contrast_matrix` is provided, the output list includes a third element:

- `contrasts` (multiple rows per feature; one per contrast)
  - `feature`: feature ID
  - `family`: best-fitting family for this feature
  - `contrast`: contrast name (from `rownames(C)`)
  - `estimate`: point estimate of the linear combination
  - `se`: standard error (computed from variance-covariance matrix)
  - `z`: z-statistic
  - `p_value`: two-sided p-value
  - `p_adj`: FDR-adjusted p-value (BH correction **by contrast** across features)

**Building Contrast Matrices**

Each row of the contrast matrix defines a linear combination of model coefficients:

```r
# For design with columns: (Intercept), StatusMild, StatusSevere, Age
# Coefficients represent:
# - (Intercept): baseline (Healthy group at Age=0)
# - StatusMild: Mild - Healthy
# - StatusSevere: Severe - Healthy
# - Age: linear age effect

# Example contrasts:
C <- rbind(
  # Mild vs Severe: (StatusMild) - (StatusSevere)
  "Mild_vs_Severe" = c(0, 1, -1, 0),
  
  # Healthy vs average of diseased: 0 - (StatusMild + StatusSevere)/2
  "Healthy_vs_Disease" = c(0, -0.5, -0.5, 0),
  
  # Severe vs Healthy (same as coefficient StatusSevere)
  "Severe_vs_Healthy" = c(0, 0, 1, 0)
)
colnames(C) <- c("(Intercept)", "StatusMild", "StatusSevere", "Age")
```

**Important Notes**

- Contrasts are computed using `vcov(fit, what = "mu")`. If the covariance matrix extraction fails or contains non-finite values, contrasts will be `NA` for that feature (the function does **not** error).
- The `p_adj` column in `contrasts` is computed **by contrast** across all features (e.g., all "Mild_vs_Severe" comparisons are adjusted together).
- Standard coefficient results in `fit$results` remain unchanged—contrasts are an **additional output**.
- The model is **re-fitted** for the best family when contrasts are requested to extract `vcov`, which adds computational overhead.

---

### Utilities

- `transform_for_family_strict(y, fam, eps = 1e-6, allow_eps = TRUE)`  
  Returns list with `y` (transformed), `mask` (valid obs), `logJ_per_obs`, and `meta`.
- `inverse_transform(z, meta)`  
  Invert strict transforms back to original scale.
- `family_groups()`  
  Family names grouped by theoretical support (count / unit / positive / real).

---

## Advanced Usage

### Transformation Modes: Strict vs Safe

PERSEO supports two transformation modes for adapting data to GAMLSS family requirements:

#### **Strict Mode** (Default with `group_by_support = TRUE`)

**Behavior:**
- Enforces theoretical domain requirements without modifying data
- Invalid observations (e.g., zeros for positive families, non-integers for count families) are **excluded via mask**
- No data repair or clipping
- Conservative, domain-preserving approach

**When to use:**
- You want support-consistent, conservative analysis
- Domain violations should exclude observations rather than transform them
- Default for most analyses when `group_by_support = TRUE`

**Example:**
```r
ff <- find_families(
  counts_matrix = counts_matrix,
  group_by_support = TRUE,
  transform_mode = "strict",  # Can be omitted (default)
  n_genes = 200,
  criterion = "BIC"
)
```

#### **Safe Mode** (Default with `group_by_support = FALSE`)

**Behavior:**
- Applies **global, deterministic, reversible affine transformations**: y* = a·y + b (a > 0)
- No observation-wise clipping or rounding
- All transformations are invertible with Jacobian correction
- Allows flexible family exploration across support boundaries

**Family-specific transformations:**
- **Positive families** (GA, GG, LOGNO, IG): Global shift if min(y) ≤ 0
  ```
  b = -min(y) + eps
  a = 1
  ```
- **Unit interval families** (BE, BEINF, BEO, etc.): Global min-max scaling
  ```
  a = 1 / (max(y) - min(y))
  b = -min(y) * a
  ```
- **Real-valued families** (NO, TF, GU): Z-score standardization (same as strict)
- **Count families** (PO, NBI, ZIP, etc.): Identity (no rounding by default)

**When to use:**
- Exploratory modeling where you're willing to compare models on transformed scale
- `group_by_support = FALSE` (testing families across support boundaries)
- You want to avoid excluding observations due to domain violations

**Example:**
```r
ff <- find_families(
  counts_matrix = counts_matrix,
  group_by_support = FALSE,
  transform_mode = "safe",  # Can be omitted (default with group_by_support = FALSE)
  n_genes = 200,
  criterion = "BIC"
)
```

#### **Comparing Modes on Same Data**

```r
library(perseo)

# Same data, strict mode
strict_result <- find_families(
  counts_matrix = counts_matrix,
  group_by_support = TRUE,
  transform_mode = "strict",
  n_genes = 100,
  n_boot = 5,
  criterion = "BIC",
  seed = 123
)

# Same data, safe mode
safe_result <- find_families(
  counts_matrix = counts_matrix,
  group_by_support = FALSE,
  transform_mode = "safe",
  n_genes = 100,
  n_boot = 5,
  criterion = "BIC",
  seed = 123
)

# Compare selected families
strict_result$top_families_overall
#> [1] "NBI"   "PO"    "ZINBI" "ZIP"

safe_result$top_families_overall
#> [1] "GG"    "LOGNO" "GA"    "TF"

# Check which mode was used
strict_result$transform_mode  #> "strict"
safe_result$transform_mode    #> "safe"
```

#### **Key Differences**

| Aspect | Strict | Safe |
|--------|--------|------|
| **Data modification** | None | Global affine transformation |
| **Invalid observations** | Excluded (via mask) | Included (after transformation) |
| **Jacobian correction** | Yes | Yes |
| **Invertibility** | N/A | Yes (via metadata) |
| **Use case** | Conservative, support-aware | Exploratory, flexible |
| **Observation-wise ops** | No | No (only global) |
| **Default when** | `group_by_support = TRUE` | `group_by_support = FALSE` |

#### **Important Notes**

- Both modes use **Jacobian correction** to ensure ICs are comparable
- Safe mode **never** applies observation-wise clipping or silent rounding
- All safe transformations are **global, affine, and invertible**
- You can explicitly set `transform_mode` regardless of `group_by_support`
- Transformation metadata is included in output for transparency

---

### Step-by-Step Workflow

For fine-grained control, you can call the internal functions directly:

```r
# Step 1: Select families via bootstrap
family_results <- find_families(
  counts_matrix = counts_matrix,
  n_genes = 300,
  n_boot = 10,
  top_n = 6,
  criterion = "BIC",
  min_n = 20,
  seed = 123
)

selected_families <- family_results$top_families_overall
print(selected_families)
#> [1] "GG"    "LOGNO" "NBI"   "GA"    "TF"    "NO"

# Step 2: Fit models with selected families
de_results <- fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = design,
  candidate_families = selected_families,
  criterion = "BIC",
  min_n = 20,
  p_adjust = "BH",  # Per-term FDR within fit_gamlss_models
  workers = 4,
  show_progress = TRUE
)

# Step 3: Manual global p-value adjustment (optional)
de_results$results$p_adj_global <- p.adjust(
  de_results$results$pval,
  method = "BH"
)

# Access results
head(de_results$selection)  # Best family per feature
head(de_results$results)    # Per-term statistics
```

### Parallel Processing

```r
library(future)

# Configure parallel backend
plan(multisession, workers = 8)

# Enable progress bars
options(progressr.enable = TRUE)
if (requireNamespace("cli", quietly = TRUE)) {
  progressr::handlers("cli")
}

# Run with parallelization
results <- run_perseo(
  counts_matrix = counts_matrix,
  design_matrix = design,
  n_genes = 500,
  n_boot = 20,
  verbose = TRUE
)

# Reset to sequential
plan(sequential)
```

### Custom Family Panel

```r
# Define custom families
my_families <- c("PO", "NBI", "ZIP", "ZINBI",  # Count distributions
                 "GA", "GG", "IG", "LOGNO",    # Positive continuous
                 "NO", "TF")                    # Real-valued

results <- run_perseo(
  counts_matrix = counts_matrix,
  design_matrix = design,
  families = my_families,
  n_genes = 200,
  top_n = 5,
  criterion = "BIC"
)
```

---

## Interpretation & Tips

### Understanding Results

**Family Selection**: 
- `top_families_overall`: Most frequently selected families across bootstrap samples
- `freq_table_overall`: How many times each family was selected
- Higher frequency = more robust/versatile for your dataset

**Differential Expression**:
- `results$pval`: Raw p-values from Wald tests
- `results$padj` or `results$p_adj`: FDR-adjusted p-values (depends on workflow)
- `results$effect`: Coefficient on link scale (not directly interpretable as fold-change)
- `selection$best_family`: Winning family per feature after IC comparison

### Recommended Settings

- **Sample size < 50**: Use `criterion = "AIC"` (less penalty)
- **Sample size 50-200**: Use `criterion = "BIC"` (balanced)
- **Sample size > 200**: Use `criterion = "GAIC"` with `gaic_k = log(n)`
- **Bootstrap**: `n_boot = 10-20` usually sufficient; higher for robustness
- **Sampling**: `n_genes = 200-500` captures diversity without excessive runtime

### Common Pitfalls

1. **All features skipped**: Check `min_n` - you may need fewer valid observations
2. **No families selected**: Increase `n_genes` or `n_boot`
3. **NAs in results**: Normal for some features; check `selection$n_valid_obs`
4. **Slow runtime**: Use parallelization with `future::plan(multisession)`

---

## Troubleshooting

**Error: "unused argument (binom_bd = ...)"**
- This parameter is only for `find_families()`, not `fit_gamlss_models()`

**Warning: "No families selected"**
- Data may be too sparse
- Try increasing `n_genes` or reducing `min_n`
- Check for extreme outliers or batch effects

**Slow performance**
- Enable parallel processing: `future::plan(multisession, workers = N)`
- Reduce `n_genes` or `n_boot` for family selection
- Use `verbose = FALSE` to suppress messages

**Memory issues**
- Process features in batches
- Reduce number of workers
- Use `criterion = "BIC"` (faster than GAIC)

---

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) —
see the `LICENSE` file for details.
