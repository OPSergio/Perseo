<div align="left">
  <h1 style="display: inline-block;">PERSEO: Model-Aware Differential Analysis for Omics Data</h1>
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
    <img src="https://img.shields.io/badge/version-1.0.0-orange?style=flat&logo=GitHub" alt="version">
  </a>
</div>

---

**PERSEO** Automated model selection and differential expression analysis for omics data using **GAMLSS** (Generalized Additive Models for Location, Scale and Shape). It supports overdispersed, skewed, or otherwise non-standard distributions, allowing for better model fit and more accurate inference.

---

## Table of Contents

- [Why PERSEO?](#why-perseo)
- [Installation](#installation)
- [Core Concepts](#core-concepts)
- [Quick Start](#quick-start)
  - [One-Function Workflow: run_perseo()](#one-function-workflow-run_perseo)
  - [Step-by-Step Workflow](#step-by-step-workflow)
- [Input Specifications](#input-specifications)
  - [Required Inputs](#required-inputs)
  - [Workflow A: Formula + Automatic Contrasts](#workflow-a-formula--automatic-contrasts)
  - [Workflow B: Design Matrix + Manual Contrasts](#workflow-b-design-matrix--manual-contrasts)
- [Visualizations](#visualizations)
  - [Volcano Plot](#volcano-plot)
  - [MA Plot](#ma-plot)
- [Interactive HTML Report](#interactive-html-report)
- [Advanced Features](#advanced-features)
  - [Hierarchical Testing (Omnibus + Contrasts)](#hierarchical-testing-omnibus--contrasts)
  - [Transformation Modes](#transformation-modes)
  - [Parallel Processing](#parallel-processing)
- [Interpretation Guide](#interpretation-guide)
- [FAQ](#faq)
- [License](#license)

---

## Why PERSEO?

Classical differential expression pipelines assume a single likelihood family for all features:

- **RNA-seq tools**: Negative Binomial (DESeq2, edgeR)
- **Microarray tools**: Normal distribution (limma)
- **Generic tools**: Log-normal, Poisson, etc.

This is **brittle** for high-dimensional omics data where:

- Features vary in **support** (counts vs. positive continuous vs. real-valued)
- Features vary in **shape** (overdispersed, zero-inflated, skewed, heavy-tailed)
- A single family cannot fit all features well

**PERSEO solves this by:**

1. **Testing multiple GAMLSS families per feature** (Negative Binomial, Gamma, Log-Normal, Beta, etc.)
2. **Using strict transformations + common masking** so all families are compared on the same observations
3. **Applying Jacobian correction** so information criteria (BIC/AIC/GAIC) are comparable on the original scale
4. **Selecting the best family per feature** and extracting coefficient statistics
5. **Scaling efficiently** via parallel processing

---

## Installation

```r
# Install from GitHub
remotes::install_github("OPSergio/Perseo")
```

**Dependencies**: R ≥ 4.2.0, `gamlss`, `gamlss.dist`, `dplyr`, `tibble`, `future`, `future.apply`, `progressr`

---

## Core Concepts

### Family Selection

PERSEO uses **bootstrap sampling** to identify commonly selected distribution families across a representative subset of features. This avoids the computational cost of testing all families on all features during exploration.

- **Bootstrap mode** (`bootstrap = TRUE`): Fast, samples `n_genes` features per pull for `n_boot` iterations
- **Full evaluation mode** (`bootstrap = FALSE`): Comprehensive, evaluates all families on all features

### Model Selection

For each feature, PERSEO:

1. Fits all candidate families on the **same observations** (common mask)
2. Applies family-specific transformations with **Jacobian correction**
3. Compares models using **BIC/AIC/GAIC** (on original scale)
4. Selects the best family and extracts **per-term statistics** from `summary(fit, what = "mu")`

### Contrasts

PERSEO supports two approaches for computing linear contrasts:

- **Automatic**: Specify a categorical variable (`contrast_variable`) → all pairwise contrasts generated
- **Manual**: Provide a custom contrast matrix (`contrast_matrix`) → full control over comparisons

### Hierarchical Testing

For multi-level factors (3+ levels), PERSEO offers **omnibus testing**:

1. **Stage 1 (Omnibus)**: Test if factor has *any* effect (multi-df test)
2. **Stage 2 (Contrasts)**: Only compute pairwise contrasts for significant features

This mirrors ANOVA + post-hoc workflow and reduces multiple testing burden.

---

## Quick Start

### One-Function Workflow: `run_perseo()`

The simplest way to use PERSEO is via `run_perseo()`, which handles family selection, differential expression, and contrast computation in one call.

```r
library(PERSEO)
library(dplyr)

# Prepare metadata
metadata <- data.frame(
  sample_id = colnames(counts_matrix),
  condition = factor(c(rep("Control", 50), rep("Treatment", 50))),
  age = rnorm(100, 50, 10),
  batch = factor(rep(1:4, each = 25)),
  row.names = colnames(counts_matrix)
)

# Set reference level
metadata$condition <- relevel(metadata$condition, ref = "Control")

# Run complete pipeline
results <- run_perseo(
  counts_matrix = counts_matrix,
  design_matrix = "~ condition + age + batch",  # formula string
  metadata = metadata,                           # required for formula
  n_genes = 200,                                 # features to sample
  n_boot = 10,                                   # bootstrap iterations
  top_n = 4,                                     # top families to select
  criterion = "BIC",                             # model selection criterion
  contrast_variable = "condition",               # auto-generate contrasts
  p_adjust_method = "BH",                        # FDR correction
  parallel = TRUE,                               # enable parallelization
  workers = 4,                                   # number of cores
  verbose = TRUE,
  seed = 123
)

# View results
print(results)

# Access components
head(results$differential_expression$results)    # Coefficient statistics
head(results$differential_expression$selection)  # Best family per feature
head(results$differential_expression$contrasts)  # Pairwise contrasts
results$family_selection$top_families_overall    # Selected families

# Filter significant results (FDR < 0.05)
sig_de <- results$differential_expression$results %>%
  filter(p_adj < 0.05, grepl("^condition", term))

sig_contrasts <- results$differential_expression$contrasts %>%
  filter(p_adj < 0.05)
```

---

### Step-by-Step Workflow

For fine-grained control, call functions separately. See [`docs/find_families.md`](docs/find_families.md) for detailed family selection documentation.

```r
# Step 1: Family Selection
family_results <- find_families(
  counts_matrix = counts_matrix,
  n_genes = 300,
  n_boot = 10,
  top_n = 6,
  criterion = "BIC",
  seed = 123
)

selected_families <- family_results$top_families_overall
print(selected_families)
#> [1] "GG"    "LOGNO" "NBI"   "GA"    "TF"    "NO"

# Step 2: Differential Expression with selected families
de_results <- fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = "~ condition + age + batch",
  metadata = metadata,
  candidate_families = selected_families,
  contrast_variable = "condition",
  criterion = "BIC",
  parallel = TRUE,
  workers = 4
)

# Access results
head(de_results$results)     # Per-term coefficients
head(de_results$selection)   # Best family per feature
head(de_results$contrasts)   # Pairwise contrasts
```

---

## Input Specifications

### Required Inputs

**`counts_matrix`**: Numeric matrix (features × samples)

- Rows = features (genes, proteins, metabolites, etc.)
- Columns = samples
- Row names = feature identifiers
- Column names = sample identifiers
- Values = numeric (finite)

**`design_matrix`**: Model specification (two options)

- **Option 1**: Formula string (e.g., `"~ condition + batch"`)
- **Option 2**: Pre-built design matrix from `model.matrix()`

**`metadata`**: Data frame with sample metadata (required for formula input or automatic contrasts)

- Row names must match column names of `counts_matrix`
- Must have columns referenced in formula or `contrast_variable`

---

### Workflow A: Formula + Automatic Contrasts

**When to use**: You want PERSEO to build the design matrix and automatically generate all pairwise contrasts for a categorical variable.

**Requirements**:

- `design_matrix` is a **formula string** or **formula object**
- `metadata` is provided (data.frame)
- `contrast_variable` names a column in metadata

**Example**:

```r
# Prepare metadata
metadata <- data.frame(
  sample_id = colnames(counts_matrix),
  tissue_type = factor(c(rep("Normal", 40), rep("Tumor_A", 30), rep("Tumor_B", 30))),
  age = rnorm(100, 60, 10),
  batch = factor(rep(1:4, each = 25)),
  row.names = colnames(counts_matrix)
)

# Set reference level
metadata$tissue_type <- relevel(metadata$tissue_type, ref = "Normal")

# Run with automatic contrast generation
results <- fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = "~ tissue_type + age + batch",  # formula string
  metadata = metadata,                             # REQUIRED
  candidate_families = c("NBI", "GG", "LOGNO"),
  contrast_variable = "tissue_type",               # auto-generate all pairwise
  criterion = "BIC",
  parallel = TRUE,
  workers = 4
)

# Automatically generates contrasts:
# - Tumor_A_vs_Normal
# - Tumor_B_vs_Normal
# - Tumor_B_vs_Tumor_A

head(results$contrasts)
```

**Benefits**:

- Simple and concise
- No need to manually build design matrix
- All pairwise contrasts generated automatically
- Coefficient names match metadata column names

---

### Workflow B: Design Matrix + Manual Contrasts

**When to use**: You need explicit control over the design matrix and want custom contrasts (e.g., average of groups, specific linear combinations).

**Requirements**:

- `design_matrix` is a **numeric matrix** from `model.matrix()`
- `contrast_matrix` is provided (numeric matrix)
- Column names of `contrast_matrix` match coefficient names from `design_matrix`

**Example**:

```r
# Build design matrix manually
metadata$tissue_type <- relevel(factor(metadata$tissue_type), ref = "Normal")
design <- model.matrix(~ tissue_type + age + batch, data = metadata)

# Inspect coefficient names
colnames(design)
#> [1] "(Intercept)" "tissue_typeTumor_A" "tissue_typeTumor_B" "age" "batch2" "batch3" "batch4"

# Define custom contrast matrix
C <- matrix(0, nrow = 3, ncol = ncol(design))
colnames(C) <- colnames(design)
rownames(C) <- c("Tumor_A_vs_Normal", "Tumor_B_vs_Tumor_A", "AvgTumor_vs_Normal")

# Tumor_A vs Normal (just the Tumor_A coefficient)
C["Tumor_A_vs_Normal", "tissue_typeTumor_A"] <- 1

# Tumor_B vs Tumor_A (Tumor_B - Tumor_A)
C["Tumor_B_vs_Tumor_A", "tissue_typeTumor_B"] <- 1
C["Tumor_B_vs_Tumor_A", "tissue_typeTumor_A"] <- -1

# Average of tumors vs Normal: (Tumor_A + Tumor_B) / 2
C["AvgTumor_vs_Normal", "tissue_typeTumor_A"] <- 0.5
C["AvgTumor_vs_Normal", "tissue_typeTumor_B"] <- 0.5

# Run with explicit contrast matrix
results <- fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = design,                          # pre-built matrix
  candidate_families = c("NBI", "GG", "LOGNO"),
  contrast_matrix = C,                             # explicit contrasts
  criterion = "BIC",
  parallel = TRUE,
  workers = 4
)

head(results$contrasts)
```

**Benefits**:

- Full control over design matrix construction
- Custom contrasts (averages, complex combinations)
- Explicit reference level handling
- Useful for non-standard designs

---

### Invalid Combinations

**DO NOT**:

```r
# INVALID: Design matrix + contrast_variable without metadata
design <- model.matrix(~ condition + batch, data = metadata)

fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = design,              # matrix (not formula)
  contrast_variable = "condition",     # ERROR: automatic generation not supported
  candidate_families = c("NBI", "GG")
)
# ERROR: When contrast_variable is provided and contrast_matrix is NULL,
#        design_matrix must be a formula + metadata must be provided.
```

**DO**:

```r
# VALID: Use formula + metadata + contrast_variable
fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = "~ condition + batch",  # formula
  metadata = metadata,                     # required
  contrast_variable = "condition",         # OK
  candidate_families = c("NBI", "GG")
)

# OR: Use design matrix + explicit contrast_matrix
C <- matrix(...)
fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = design,       # matrix
  contrast_matrix = C,           # explicit
  candidate_families = c("NBI", "GG")
)
```

---

## Visualizations

PERSEO provides two ggplot2-based visualization functions for exploring differential expression results. Both accept the output of `run_perseo()` or `fit_gamlss_models()` and return `ggplot` objects that can be further customized or saved.

**Install optional label package for cleaner plots:**

```r
install.packages("ggrepel")  # non-overlapping text labels
```

### Volcano Plot

`plot_volcano()` displays log fold change against -log10 adjusted p-value, with four significance categories colour-coded:

| Category | Criterion |
|----------|-----------|
| **Not sig** | FDR ≥ threshold AND \|LFC\| < threshold |
| **LFC only** | \|LFC\| ≥ threshold, FDR not significant |
| **FDR only** | FDR significant, \|LFC\| < threshold |
| **Sig & LFC** | Both FDR and \|LFC\| thresholds met |

```r
# Basic volcano plot for one contrast
p <- plot_volcano(results, contrast = "Treatment_vs_Control")
print(p)

# Customise thresholds and labels
p <- plot_volcano(
  results,
  contrast       = "TreatmentA_vs_Control",
  fdr_threshold  = 0.05,    # FDR significance cutoff
  lfc_threshold  = 1,       # |LFC| cutoff (log scale)
  label_top      = 15,      # label top N significant features
  point_size     = 2,
  alpha          = 0.6,
  title          = "My Volcano"
)

# Save to file
ggplot2::ggsave("volcano_TreatA_vs_Control.png", p, width = 8, height = 5.5, dpi = 150)

# Loop over all contrasts
contrasts <- unique(results$differential_expression$contrasts$contrast)
for (cname in contrasts) {
  p <- plot_volcano(results, contrast = cname, label_top = 10)
  ggplot2::ggsave(paste0("volcano_", cname, ".png"), p, width = 8, height = 5.5)
}
```

### MA Plot

`plot_ma()` displays mean expression (log2 scale) against log fold change. Points are coloured by significance. Providing the original `counts_matrix` places features on the true mean expression axis; without it, features are ranked by |LFC|.

```r
# MA plot — true mean expression on x-axis
p <- plot_ma(
  results,
  contrast      = "Treatment_vs_Control",
  counts_matrix = counts_matrix      # original counts for x-axis
)
print(p)

# MA plot without counts (uses |LFC| rank on x-axis)
p <- plot_ma(results, contrast = "Treatment_vs_Control")

# Customise
p <- plot_ma(
  results,
  contrast      = "TreatmentB_vs_TreatmentA",
  counts_matrix = counts_matrix,
  fdr_threshold = 0.01,
  label_top     = 20
)

ggplot2::ggsave("ma_TreatB_vs_TreatA.png", p, width = 8, height = 5.5, dpi = 150)
```

**Parameters shared by both functions:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `contrast` | `NULL` | Contrast name to plot; auto-selected if only one exists |
| `fdr_threshold` | `0.05` | Adjusted p-value threshold |
| `label_top` | `10` | Top N significant features to label (0 = none) |
| `point_size` | `1.8` | Base point size |
| `alpha` | `0.7` | Point transparency |
| `title` | `NULL` | Plot title; auto-generated if NULL |

---

## Interactive HTML Report

`report_perseo()` generates a self-contained HTML report from the output of `run_perseo()`. The report includes:

- **Header**: PERSEO logo, version, run timestamp
- **Analysis Parameters**: Full table of all `run_perseo()` inputs (data, families, models, contrasts, multiple testing, parallelisation)
- **KPI Summary**: Total features, significant hits (FDR < 5%), families tested, contrasts computed
- **Family Selection**: Interactive bar chart (ggiraph) of family frequency distribution + summary table
- **Differential Expression**: Searchable, sortable DT table of all coefficient results
- **Contrasts**: Per-contrast tabset with interactive volcano plot, MA plot, and top-hits table
- **Omnibus** *(if enabled)*: Table of omnibus test results

**Install required packages:**

```r
install.packages(c("rmarkdown", "ggiraph", "DT", "htmltools"))
install.packages("ggrepel")  # optional, for labelled plots
```

**Usage:**

```r
# Generate report in current directory
report_perseo(results)

# Custom output location and title
report_perseo(
  x           = results,
  output_file = "my_analysis.html",
  output_dir  = "reports/",
  title       = "Tumour vs Normal — PERSEO Report",
  open        = TRUE    # auto-open in browser
)

# Suppress browser opening (e.g. in scripts)
path <- report_perseo(results, open = FALSE)
cat("Report saved to:", path, "\n")
```

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `x` | — | Output of `run_perseo()` |
| `output_file` | `"perseo_report.html"` | Output filename |
| `output_dir` | `"."` | Output directory (created if absent) |
| `title` | `"PERSEO Analysis Report"` | Report title in header |
| `open` | `TRUE` | Open in browser after rendering |
| `quiet` | `FALSE` | Suppress rmarkdown render messages |

**Full example — 4-group analysis:**

```r
library(PERSEO)

# 1. Build metadata and run pipeline
metadata <- data.frame(
  group = factor(rep(c("Control", "TreatA", "TreatB", "TreatC"), each = 20))
)
results <- run_perseo(
  counts_matrix     = counts_matrix,
  design_matrix     = "~ group",
  metadata          = metadata,
  contrast_variable = "group",   # auto-generates all 6 pairwise contrasts
  n_genes           = 150,
  n_boot            = 8,
  top_n             = 4,
  criterion         = "BIC",
  p_adjust_method   = "BH",
  seed              = 2024
)

# 2. Visualise individual contrasts
plot_volcano(results, "TreatA_vs_Control", label_top = 15)
plot_ma(results, "TreatA_vs_Control", counts_matrix = counts_matrix)

# 3. Generate full interactive report
report_perseo(
  results,
  output_file = "perseo_4group_report.html",
  output_dir  = "reports/",
  title       = "4-Group PERSEO Analysis"
)
```

---

## Advanced Features

### Hierarchical Testing (Omnibus + Contrasts)

For multi-level factors (3+ levels), hierarchical testing reduces multiple testing burden:

1. **Omnibus test**: Does the factor have *any* effect? (multi-df test)
2. **Contrasts**: Only for features passing omnibus → compute pairwise comparisons

**When to use**:

- Factor has 3+ levels
- Want stricter false positive control
- Large datasets (1000+ features)

**When to skip**:

- Factor has only 2 levels (omnibus = single pairwise test)
- Exploratory analysis (want to see all contrasts)
- Pre-specified contrasts of interest

**Example**:

```r
# Multi-level factor
metadata$tissue_type <- factor(rep(c("Normal", "Luminal", "Basal"), length.out = 100))

results <- fit_gamlss_models(
  counts_matrix = counts_matrix,
  design_matrix = "~ tissue_type + age + batch",
  metadata = metadata,
  candidate_families = c("NBI", "GG"),
  contrast_variable = "tissue_type",
  omnibus = TRUE,                 # Enable hierarchical testing
  omnibus_threshold = 0.05,       # Filter threshold
  omnibus_test = "Wald",          # "Wald" (fast) or "LRT" (robust)
  parallel = TRUE,
  workers = 4
)

# View omnibus results
head(results$omnibus)
#> # A tibble: 6 × 7
#>   feature  family test_type statistic    df p_value pass 
#>   <chr>    <chr>  <chr>         <dbl> <int>   <dbl> <lgl>
#> 1 Gene_001 NBI    Wald          12.5      2  0.0019 TRUE 
#> 2 Gene_002 GG     Wald           3.2      2  0.201  FALSE
#> 3 Gene_003 NBI    Wald          18.7      2  0.0001 TRUE 

# Contrasts only for features passing omnibus
head(results$contrasts)

# Filter significant contrasts
sig_contrasts <- results$contrasts %>%
  filter(p_adj < 0.05)
```

**Omnibus test comparison**:

| Test | Speed | Best For | Notes |
|------|-------|----------|-------|
| **Wald** | Fast | Large samples (n > 30/group), standard families | Uses vcov from fitted model; asymptotic approximation |
| **LRT** | ~2× slower | Small samples (n < 30/group), complex families | Refits reduced model; more robust |

---

### Transformation Modes

PERSEO supports two transformation strategies:

#### `transform_mode = "strict"` (Default)

**Behavior**:

- Domain-preserving transformations
- Invalid observations **excluded via masking**
- No data modification
- Conservative, support-consistent approach

**When to use**:

- Default for most analyses
- Domain violations should exclude observations
- Support-aware family selection

**Example**:

```r
results <- find_families(
  counts_matrix = counts_matrix,
  transform_mode = "strict",  # default
  group_by_support = TRUE,
  criterion = "BIC"
)
```

#### `transform_mode = "safe"`

**Behavior**:

- Global affine transformations: y* = a·y + b (a > 0)
- **Invertible** with Jacobian correction
- No observation exclusion
- All observations retained after transformation

**When to use**:

- Want more flexibility
- Testing families across support boundaries
- Want to avoid observation exclusion

**Example**:

```r
results <- find_families(
  counts_matrix = counts_matrix,
  transform_mode = "safe",
  group_by_support = FALSE,  # flexible family exploration
  criterion = "BIC"
)
```

**Key differences**:

| Aspect | Strict | Safe |
|--------|--------|------|
| **Data modification** | None | Global affine transformation |
| **Invalid observations** | Excluded | Included (after transform) |
| **Jacobian correction** | Yes | Yes |
| **Invertibility** | N/A | Yes (via metadata) |
| **Use case** | Conservative | Flexibility |

See [`docs/utils_transformations.md`](docs/utils_transformations.md) for technical details.

---

### Parallel Processing

PERSEO includes built-in parallel processing support via `future` backend.

**Automatic mode** (recommended):

```r
results <- run_perseo(
  counts_matrix = counts_matrix,
  design_matrix = design,
  parallel = TRUE,     # Auto-configure parallel backend
  workers = 8,         # Number of cores (auto-detected if NULL)
  verbose = TRUE
)

# Automatically:
# 1. Sets up future::plan(multisession)
# 2. Runs analysis in parallel
# 3. Resets to sequential on completion
```

**Manual mode** (advanced):

```r
library(future)
plan(multisession, workers = 8)

results <- run_perseo(
  counts_matrix = counts_matrix,
  design_matrix = design
)

plan(sequential)
```

**Recommended settings**:

| Dataset Size | Parallel | Workers |
|--------------|----------|---------|
| < 1000 features | `FALSE` | N/A |
| 1000–10,000 features | `TRUE` | 4–8 |
| > 10,000 features | `TRUE` | 8–16 |

**Memory considerations**:

- Each worker loads a copy of the data
- More workers = more memory usage
- If memory issues occur, reduce `workers`

---

## Interpretation Guide

### Understanding Results

**Family Selection** (`family_selection`):

- `top_families_overall`: Most frequently selected families across bootstrap samples
- `freq_table_overall`: Selection frequency per family
- Higher frequency → more robust/versatile for your dataset

**Differential Expression** (`differential_expression$results`):

- `feature`: Feature identifier
- `term`: Coefficient name (e.g., `conditionTreatment`, `age`, `batch2`)
- `effect`: Coefficient estimate (on link scale)
- `se`: Standard error
- `stat`: Wald statistic (typically t or z)
- `pval`: Raw p-value
- `padj`: **FDR-adjusted p-value (by term across features)**

**Selection** (`differential_expression$selection`):

- `feature`: Feature identifier
- `best_family`: Winning family per feature
- `n_valid_obs`: Number of observations used (after common mask)
- `ic_value`: Jacobian-corrected information criterion value

**Contrasts** (`differential_expression$contrasts`):

- `feature`: Feature identifier
- `family`: Best family for this feature
- `contrast`: Contrast name (e.g., `Treatment_vs_Control`)
- `estimate`: Point estimate of linear combination
- `se`: Standard error
- `z`: z-statistic
- `p_value`: Raw p-value
- `p_adj`: **FDR-adjusted p-value (by contrast across features)**

**Omnibus** (`differential_expression$omnibus`, when enabled):

- `feature`: Feature identifier
- `family`: Best family
- `test_type`: `"Wald"` or `"LRT"`
- `statistic`: Test statistic
- `df`: Degrees of freedom
- `p_value`: Omnibus p-value
- `pass`: Logical; `TRUE` if p_value < omnibus_threshold

### Multiple Testing Strategy

**Standard workflow** (`omnibus = FALSE`):

- Coefficient p-values adjusted **by term** across features
- Contrast p-values adjusted **by contrast** across features
- Each term/contrast treated independently

**Hierarchical workflow** (`omnibus = TRUE`):

1. **Stage 1 (Omnibus)**: Test each feature for factor effect (1 test per feature)
2. **Stage 2 (Contrasts)**: Only for significant features, compute pairwise tests

**Example**: 500 features, 3-level factor

| Approach | Stage 1 Tests | Stage 2 Tests | Total Tests | Power |
|----------|--------------|---------------|-------------|-------|
| Standard | 0 | 500 × 3 = 1500 | 1500 | Lower |
| Hierarchical | 500 | 150 × 3 = 450 | 950 | Higher |

**Note**: Omnibus test acts only as a **gatekeeper**. Contrast p-values are still adjusted by contrast across features (unchanged from standard workflow).

### Common Patterns

**Filter significant DE results**:

```r
sig_results <- results$differential_expression$results %>%
  filter(p_adj < 0.05, grepl("^condition", term)) %>%
  arrange(p_adj)
```

**Filter significant contrasts**:

```r
sig_contrasts <- results$differential_expression$contrasts %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)
```

**Check family distribution**:

```r
table(results$differential_expression$selection$best_family)
#> GG LOGNO   NBI    GA    TF    NO 
#> 142    91    64    28     7     2
```

**Omnibus pass rate**:

```r
sum(results$differential_expression$omnibus$pass, na.rm = TRUE) / 
  nrow(results$differential_expression$omnibus) * 100
#> 32.5% passed omnibus filter
```

---

## FAQ

**Q: When should I use `bootstrap = TRUE` vs. `bootstrap = FALSE`?**

A: Use `bootstrap = TRUE` (default) for fast exploratory analysis. The bootstrap samples a subset of features to identify common families, then applies them to all features. Use `bootstrap = FALSE` for comprehensive analysis where you want to evaluate all families on all features (slower but thorough).

**Q: How do I choose between Workflow A and Workflow B?**

A: Use **Workflow A** (formula + automatic contrasts) for simple designs with standard pairwise comparisons. Use **Workflow B** (design matrix + manual contrasts) when you need custom contrasts (e.g., averages, complex linear combinations) or explicit control over reference levels.

**Q: What's the difference between Wald and LRT omnibus tests?**

A: **Wald** is faster (uses vcov from fitted model) and suitable for large samples (n > 30/group). **LRT** is more robust (refits reduced model) and better for small samples (n < 30/group) or complex families. LRT is ~2× slower.

**Q: Why do I get NA values in contrast results?**

A: NAs occur when:
- Variance-covariance matrix extraction fails
- Vcov contains non-finite values
- Model did not converge properly

This is normal for some features; check `selection$n_valid_obs` to ensure sufficient data.

**Q: How do I interpret `effect` values in results?**

A: `effect` is the coefficient estimate on the **link scale** (not directly interpretable as fold-change). For count families (NBI, PO), the link is typically log. For positive families (GG, LOGNO), also log. For interpretation, focus on:
- Sign (positive/negative effect)
- Statistical significance (padj)
- Contrast estimates for direct comparisons

**Q: Can I use PERSEO with non-count data?**

A: Yes! PERSEO supports:
- **Count data**: RNA-seq counts → NBI, PO, ZIP, ZINBI
- **Positive continuous**: Proteomics intensities → GG, GA, LOGNO, IG
- **Unit interval**: Beta values, proportions → BE, BEINF, BEO
- **Real-valued**: Normalized/transformed data → NO, TF

Use `group_by_support = TRUE` to restrict families by empirical support.

**Q: How do I speed up analysis?**

A: 
1. Enable parallel processing: `parallel = TRUE, workers = 8`
2. Use `criterion = "BIC"` (faster than GAIC)
3. Reduce `n_genes` and `n_boot` for family selection
4. Use `omnibus_test = "Wald"` (faster than LRT)
5. Skip omnibus filtering if not needed (`omnibus = FALSE`)

**Q: What if all features are skipped?**

A: Check `min_n` parameter. You may need fewer valid observations. Also check for:
- Extreme outliers
- Batch effects
- Data quality issues

**Q: How are contrast p-values adjusted?**

A: Contrast p-values are adjusted **by contrast** across all features using the specified FDR method (default: BH). Each contrast is treated as an independent hypothesis family. Example:
- `Treatment_vs_Control`: All features with this contrast are adjusted together
- `High_vs_Low`: All features with this contrast are adjusted together

This is true regardless of `omnibus` setting. The omnibus test serves only as a gatekeeper, not as part of the adjustment procedure.

---

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) — see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use PERSEO in your research, please cite:

Olmos-Piñero, S. & Fernández-Lanza Val, F.  
*PERSEO: Flexible differential expression analysis using GAMLSS*.  
Manuscript in preparation.

---

## Contributing

Contributions are welcome! Please open an issue or pull request on [GitHub](https://github.com/OPSergio/Perseo).

---

## Contact

For questions or support:

- **Issues**: [GitHub Issues](https://github.com/OPSergio/Perseo/issues)
- **Email**: solmos97@gmail.com

---

**Note**: This package is under active development. API and functionality may change in future versions.
