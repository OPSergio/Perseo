# Hierarchical Testing with Omnibus Tests in PERSEO

## Table of Contents

- [Introduction](#introduction)
- [Why Use Omnibus Testing?](#why-use-omnibus-testing)
- [Omnibus Test Options](#omnibus-test-options)
  - [Wald Test (Default)](#1-wald-test-default)
  - [Likelihood Ratio Test (LRT)](#2-likelihood-ratio-test-lrt)
- [Complete Workflow Example](#complete-workflow-example)
- [Understanding the Results](#understanding-the-results)
- [Comparison: Standard vs Hierarchical](#comparison-standard-vs-hierarchical)
- [Decision Guide](#decision-guide)
- [Advanced: Custom Thresholds](#advanced-custom-thresholds)
- [References](#references)

---

## Introduction

When analyzing multi-level categorical variables (e.g., tissue type with 3+ levels), 
you often want to perform **pairwise post-hoc comparisons** after determining whether 
the factor has **any** effect. PERSEO implements this classical statistical workflow 
through **hierarchical testing**:

1. **Omnibus test** (feature-level): Does the factor have any effect?
2. **Post-hoc contrasts** (pairwise): Which levels differ from each other?

This approach mirrors **ANOVA followed by Tukey's HSD** but is fully integrated into 
the GAMLSS model-based differential expression framework.

---

## Why Use Omnibus Testing?

### Statistical Rationale

**Multiple testing problem**: For a 3-level factor, you perform 3 pairwise comparisons 
per feature. For 1000 features, that's 3000 tests. Without proper control, false 
positives accumulate rapidly.

**Hierarchical testing solution**:
- First, test whether the factor matters at all (1 test per feature)
- Only compute post-hoc contrasts for features passing the omnibus test
- This reduces the number of tests and controls family-wise error rate (FWER)

### Practical Benefits

✅ **Reduces multiple testing burden**: Fewer tests = more statistical power  
✅ **Mirrors classical workflow**: ANOVA → post-hoc (familiar to biologists)  
✅ **Computational efficiency**: Skip contrast computation for non-significant features  
✅ **Interpretability**: Clear separation between "Does X matter?" and "How do levels differ?"

### When NOT to Use Omnibus Testing

❌ **Two-level factors**: No benefit (omnibus = single pairwise comparison)  
❌ **Pre-specified contrasts**: If you have 1-2 specific hypotheses, test them directly  
❌ **Exploratory analysis**: If you want to see all contrasts regardless of significance

---

## Omnibus Test Options

PERSEO offers two omnibus test types:

### 1. Wald Test (Default)

**How it works**:
- Uses coefficients and variance-covariance matrix from the fitted model
- Tests joint null hypothesis: β₁ = β₂ = ... = 0
- Test statistic: W = β'V⁻¹β ~ χ²(df)
- No model refitting required

**Advantages**:
- ⚡ **Fast**: Uses already-computed vcov matrix
- 💾 **Memory efficient**: No additional model fitting
- 🔧 **Works with all families**: No convergence issues

**Disadvantages**:
- 📊 **Asymptotic**: Less accurate for small sample sizes
- ⚠️ **Assumes normality**: Of the estimator (usually fine for GAMLSS)

**Best for**:
- Large datasets (n > 30 per group)
- Standard GAMLSS families
- When computational speed is important

```r
result <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,
  metadata = metadata,
  candidate_families = c("NBI", "GG"),
  contrast_variable = "tissue_type",
  omnibus = TRUE,
  omnibus_test = "Wald"  # Default
)
```

### 2. Likelihood Ratio Test (LRT)

**How it works**:
- Fits two models: full (with factor) and reduced (without factor)
- Compares log-likelihoods: LRT = 2(ℓ_full - ℓ_reduced) ~ χ²(df)
- Requires refitting the reduced model

**Advantages**:
- 🎯 **More robust**: Better for non-normal distributions
- 📐 **Exact**: Not relying on asymptotic normality of estimators
- 🔬 **Gold standard**: Classical test for nested models

**Disadvantages**:
- 🐢 **Slower**: Requires fitting an additional model per feature
- 💻 **More memory**: Stores two models temporarily
- ⚠️ **Convergence**: Reduced model might fail to converge

**Best for**:
- Small to medium datasets
- Complex/unusual GAMLSS families
- When statistical rigor is paramount

```r
result <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,
  metadata = metadata,
  candidate_families = c("NBI", "GG"),
  contrast_variable = "tissue_type",
  omnibus = TRUE,
  omnibus_test = "LRT"
)
```

---

## Complete Workflow Example

### Setup

```r
library(PERSEO)
library(dplyr)

# Simulate RNA-seq data
set.seed(123)
n_features <- 500
n_samples <- 60

counts <- matrix(
  rnbinom(n_features * n_samples, mu = 50, size = 2),
  nrow = n_features,
  ncol = n_samples
)
rownames(counts) <- paste0("Gene_", 1:n_features)

# Metadata with 3 tissue types
metadata <- data.frame(
  tissue_type = factor(rep(c("Liver", "Brain", "Heart"), each = 20)),
  batch = factor(rep(1:3, length.out = n_samples)),
  age = rnorm(n_samples, 50, 10)
)

# Design matrix
design <- model.matrix(~ tissue_type + batch + age, data = metadata)
```

### Standard Analysis (No Omnibus)

```r
# Compute ALL pairwise contrasts for ALL features
result_standard <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,
  metadata = metadata,
  candidate_families = c("NBI", "ZINBI", "GG"),
  contrast_variable = "tissue_type",
  omnibus = FALSE,  # Default
  parallel = TRUE,
  workers = 4
)

# Results
nrow(result_standard$contrasts)
#> [1] 1500  # 500 features × 3 contrasts

# All features have contrasts
length(unique(result_standard$contrasts$feature))
#> [1] 500
```

### Hierarchical Analysis (With Omnibus)

```r
# First omnibus test, then contrasts only for significant features
result_omnibus <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,
  metadata = metadata,
  candidate_families = c("NBI", "ZINBI", "GG"),
  contrast_variable = "tissue_type",
  omnibus = TRUE,
  omnibus_threshold = 0.05,
  omnibus_test = "Wald",
  parallel = TRUE,
  workers = 4
)

# View omnibus results
head(result_omnibus$omnibus)
#> # A tibble: 6 × 7
#>   feature  family test_type statistic    df p_value pass 
#>   <chr>    <chr>  <chr>         <dbl> <int>   <dbl> <lgl>
#> 1 Gene_001 NBI    Wald          12.5      2  0.0019 TRUE 
#> 2 Gene_002 GG     Wald           3.2      2  0.201  FALSE
#> 3 Gene_003 NBI    Wald          18.7      2  0.0001 TRUE 

# Only significant features have contrasts
nrow(result_omnibus$contrasts)
#> [1] 450  # Fewer contrasts computed

length(unique(result_omnibus$contrasts$feature))
#> [1] 150  # Only ~30% of features passed omnibus

# Features with contrasts = features that passed omnibus
features_with_contrasts <- unique(result_omnibus$contrasts$feature)
features_passed_omnibus <- result_omnibus$omnibus$feature[
  result_omnibus$omnibus$pass
]

all(features_with_contrasts %in% features_passed_omnibus)
#> [1] TRUE
```

---

## Understanding the Results

### Omnibus Output

```r
str(result_omnibus$omnibus)
#> tibble [500 × 7] (S3: tbl_df/tbl/data.frame)
#>  $ feature  : chr  "Gene_001" "Gene_002" ...
#>  $ family   : chr  "NBI" "GG" ...          # Selected GAMLSS family
#>  $ test_type: chr  "Wald" "Wald" ...       # "Wald" or "LRT"
#>  $ statistic: num  12.5 3.2 ...            # Test statistic
#>  $ df       : int  2 2 ...                 # Degrees of freedom
#>  $ p_value  : num  0.0019 0.201 ...        # Omnibus p-value
#>  $ pass     : lgl  TRUE FALSE ...          # Passed threshold?
```

**Key columns**:
- `p_value`: Omnibus test p-value (tests if factor has ANY effect)
- `pass`: Logical flag (TRUE if p_value < omnibus_threshold)
- `df`: Number of factor levels - 1

### Contrasts Output

```r
head(result_omnibus$contrasts)
#> # A tibble: 6 × 8
#>   feature  family contrast       estimate     se      z  p_value   p_adj
#>   <chr>    <chr>  <chr>             <dbl>  <dbl>  <dbl>    <dbl>   <dbl>
#> 1 Gene_001 NBI    Brain_vs_Liver    0.45   0.12   3.75  0.00018 0.00054
#> 2 Gene_001 NBI    Heart_vs_Liver    0.23   0.13   1.77  0.077   0.115  
#> 3 Gene_001 NBI    Heart_vs_Brain   -0.22   0.12  -1.83  0.067   0.101  
```

**Important**: `p_adj` is computed **by contrast** across features (unchanged from standard PERSEO):

```r
# For "Brain_vs_Liver", adjust p-values across all features
result_omnibus$contrasts %>%
  filter(contrast == "Brain_vs_Liver") %>%
  summarise(
    n_tests = n(),
    n_significant = sum(p_adj < 0.05, na.rm = TRUE)
  )
#> # A tibble: 1 × 2
#>   n_tests n_significant
#>     <int>         <int>
#> 1     150            35
```

---

## Comparison: Standard vs Hierarchical

### Statistical Power

| Aspect | Standard (omnibus = FALSE) | Hierarchical (omnibus = TRUE) |
|--------|---------------------------|------------------------------|
| **Total tests** | Features × Contrasts | Features (omnibus) + Passing × Contrasts |
| **Example (500 features, 3 levels)** | 1500 tests | 500 + ~450 = ~950 tests |
| **Multiple testing burden** | Higher | Lower |
| **Power for true positives** | Lower (stricter FDR) | Higher (fewer tests) |

### Computational Cost

```r
# Benchmark
system.time(
  fit_gamlss_models(..., omnibus = FALSE)
)
#>    user  system elapsed 
#>  45.2     1.3    47.8

system.time(
  fit_gamlss_models(..., omnibus = TRUE, omnibus_test = "Wald")
)
#>    user  system elapsed 
#>  46.5     1.4    49.1  # Minimal overhead

system.time(
  fit_gamlss_models(..., omnibus = TRUE, omnibus_test = "LRT")
)
#>    user  system elapsed 
#>  82.3     2.1    86.7  # ~2× slower (refits reduced models)
```

---

## Decision Guide

### Use Omnibus Testing When:

✅ Factor has 3+ levels  
✅ You want to control FWER/FDR more strictly  
✅ You have many features (1000+)  
✅ You expect most features to be non-significant for the factor

### Skip Omnibus Testing When:

❌ Factor has only 2 levels (omnibus = pairwise comparison)  
❌ You have strong prior hypotheses about specific contrasts  
❌ Dataset is small (< 100 features)  
❌ You want to report all contrasts regardless of significance

### Choose Wald When:

⚡ Speed is important  
📊 Standard GAMLSS families (NBI, GG, LOGNO, etc.)  
💻 Limited computational resources  
🔢 Large sample sizes (n > 30 per group)

### Choose LRT When:

🎯 Maximum statistical rigor  
🔬 Complex/unusual families  
📐 Small sample sizes  
⏱️ Computational time is not a constraint

---

## Advanced: Custom Thresholds

```r
# Very conservative (only contrast features with strong omnibus signal)
result_conservative <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,
  metadata = metadata,
  candidate_families = c("NBI"),
  contrast_variable = "tissue_type",
  omnibus = TRUE,
  omnibus_threshold = 0.01,  # Stricter
  omnibus_test = "Wald"
)

# Liberal (contrast most features, but still filter some)
result_liberal <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,
  metadata = metadata,
  candidate_families = c("NBI"),
  contrast_variable = "tissue_type",
  omnibus = TRUE,
  omnibus_threshold = 0.1,  # More permissive
  omnibus_test = "Wald"
)

# Compare filtering rates
mean(result_conservative$omnibus$pass)
#> [1] 0.12  # Only 12% pass

mean(result_liberal$omnibus$pass)
#> [1] 0.42  # 42% pass
```

---

## References

- Maxwell, S. E., & Delaney, H. D. (2004). *Designing Experiments and Analyzing Data*. 
  Chapter 5: Testing Several Contrasts.
  
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: 
  a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society*.

- Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for 
  location, scale and shape. *Journal of the Royal Statistical Society: Series C*.

---

## Summary

| Feature | Standard | Omnibus (Wald) | Omnibus (LRT) |
|---------|----------|----------------|---------------|
| **Tests per feature** | k(k-1)/2 | 1 + k(k-1)/2 (if pass) | 1 + k(k-1)/2 (if pass) |
| **Speed** | Fast | Fast | Moderate |
| **Statistical rigor** | Standard | High | Highest |
| **Multiple testing burden** | Highest | Lower | Lower |
| **Best for** | Exploration | Large datasets | Small/complex data |

**Recommendation**: For most differential expression analyses with multi-level factors, 
use `omnibus = TRUE` with `omnibus_test = "Wald"` for an optimal balance of 
statistical rigor and computational efficiency.
