# `find_families`: Family Selection via Bootstrap Sampling

## Overview

`find_families()` performs automated distribution family selection for high-dimensional omics data using bootstrap sampling and information criterion comparison. The function evaluates **intercept-only** GAMLSS models across candidate families, applying Jacobian-corrected information criteria and common masking to ensure valid cross-family comparisons.

## Workflow Context

`find_families()` is **Step 1** in the PERSEO pipeline:

1. **`find_families()`**: Identify frequently selected families on a representative subset of features
2. **`fit_gamlss_models()`**: Apply selected families to all features with full regression models (covariates + contrasts)

> **Important**: `find_families()` uses **intercept-only models** (no covariates) to identify robust families. The selected families are then used in `fit_gamlss_models()` with the full design matrix. See the main [README](../README.md) for complete workflows.

**For differential expression analysis**, see:
- [README - Workflow A: Formula + Automatic Contrasts](../README.md#workflow-a-formula--automatic-contrasts)
- [README - Workflow B: Design Matrix + Manual Contrasts](../README.md#workflow-b-design-matrix--manual-contrasts)

## Statistical Framework

### Model Selection Problem

For each feature `i` with observations `yᵢ = (yᵢ₁, ..., yᵢₙ)`, select the family `f*` from candidates `ℱ` that minimizes:

```
f* = argmin_{f ∈ ℱ} IC(Mf)
```

where `IC(Mf)` is the information criterion for model `Mf` fitted with family `f`.

### Bootstrap Sampling Approach

1. Draw `B` bootstrap samples of `m` features each (without replacement within sample)
2. For each sampled feature, fit intercept-only models with all candidate families
3. Select best family per feature based on IC
4. Aggregate family frequencies across bootstrap samples
5. Return top `k` most frequent families

### Transformation and Jacobian Correction

For transformation `z = g(y)`, the log-likelihood on the original scale is:

```
log L(θ; y) = log L(θ; z) + Σⱼ log|∂zⱼ/∂yⱼ|
```

Information criterion becomes:

```
IC(Mf) = -2[log L(θ̂; z) + Σⱼ log|Jⱼ|] + penalty(df, n)
```

### Common Mask Rationale

All candidate families are evaluated on identical observations via intersection of validity masks:

```
mask_common = ⋂_{f ∈ ℱ} mask_f
```

This ensures:
- Consistent effective sample size `n_valid` across families
- Valid IC comparisons (same data likelihood base)
- Unbiased selection (no confounding between fit quality and data availability)

---

## Function Signature

```r
find_families(
  counts_matrix,
  n_genes = 200,
  n_boot = 10,
  top_n = 4,
  families = NULL,
  verbose = TRUE,
  min_n = 5,
  seed = NULL,
  group_by_support = TRUE,
  binom_bd = NULL,
  criterion = c("GAIC", "BIC", "AIC"),
  gaic_k = NULL,
  filter_beta_inflated = TRUE,
  thr_zero = 0.005,
  thr_one = 0.005,
  transform_mode = NULL
)
```

---

## Arguments

### Data Arguments

**`counts_matrix`** (required)
- Type: Numeric matrix (features × samples)
- Structure: Rows = features, Columns = samples
- Requirements: Rownames should contain feature identifiers
- Values: Numeric (finite values required for analysis)

**`design_matrix`**
- Note: Not used in `find_families()` (intercept-only models)
- See `fit_gamlss_models()` for covariate inclusion

### Sampling Arguments

**`n_genes`** (default: `200`)
- Type: Integer > 0
- Description: Number of features sampled per bootstrap pull
- Constraint: Must be ≤ `nrow(counts_matrix)`
- Effect: Larger values increase computational cost but improve coverage

**`n_boot`** (default: `10`)
- Type: Integer > 0
- Description: Number of bootstrap pulls
- Effect: More pulls increase frequency estimate stability
- Computational cost: Scales linearly with `n_boot`

**`seed`** (default: `NULL`)
- Type: Integer or NULL
- Description: Random seed for reproducible sampling
- If NULL: Uses current RNG state

### Family Selection Arguments

**`families`** (default: `NULL`)
- Type: Character vector or NULL
- Description: Candidate families to evaluate
- If NULL: Uses default panel (18 families across all supports)
- Default panel includes:
  - Count: `PO`, `NBI`, `ZIP`, `ZINBI`, `ZIP2`, `BI`, `BB`
  - Unit: `BE`, `BEINF`, `BEO`, `BEZI`, `BEo`, `BEINF0`
  - Positive: `GA`, `GG`, `IG`, `LOGNO`
  - Real: `NO`, `TF`, `GU`

**`group_by_support`** (default: `TRUE`)
- Type: Logical
- Description: Whether to filter families by empirical support
- TRUE: Families restricted to support-compatible subset
- FALSE: All families in `families` are tested
- Interaction with `transform_mode`: Affects default transformation mode

**`transform_mode`** (default: `NULL`)
- Type: Character or NULL (`"strict"`, `"safe"`, or `NULL`)
- Description: Transformation mode for family comparison
- NULL behavior:
  - `group_by_support = TRUE` → `"strict"`
  - `group_by_support = FALSE` → `"safe"`
- `"strict"`: Domain-preserving, observation exclusion
- `"safe"`: Global affine transformations, all observations retained

**`binom_bd`** (default: `NULL`)
- Type: Numeric scalar, vector, or NULL
- Description: Binomial denominator for BI/BB families
- NULL: Inferred per feature as `max(y)` when all values are integers
- Scalar: Same denominator for all samples
- Vector: Length must equal `ncol(counts_matrix)`

**`filter_beta_inflated`** (default: `TRUE`)
- Type: Logical
- Description: Whether to exclude inflated Beta families without empirical evidence
- Thresholds controlled by `thr_zero` and `thr_one`

**`thr_zero`** (default: `0.005`)
- Type: Numeric in [0, 1]
- Description: Minimum proportion of zeros required for zero-inflated families
- Used when `filter_beta_inflated = TRUE`

**`thr_one`** (default: `0.005`)
- Type: Numeric in [0, 1]
- Description: Minimum proportion of ones required for one-inflated families
- Used when `filter_beta_inflated = TRUE`

### Information Criterion Arguments

**`criterion`** (default: `"GAIC"`)
- Type: Character (`"GAIC"`, `"BIC"`, or `"AIC"`)
- Description: Information criterion for model selection
- Formulas:
  - AIC: `-2 log L + 2·df`
  - BIC: `-2 log L + log(n)·df`
  - GAIC: `-2 log L + k·df` where `k` specified by `gaic_k`

**`gaic_k`** (default: `NULL`)
- Type: Numeric or NULL
- Description: Penalty multiplier for GAIC
- NULL behavior: Uses `log(n_valid)` (equivalent to BIC)
- Ignored if `criterion != "GAIC"`

**`min_n`** (default: `5`)
- Type: Integer ≥ 1
- Description: Minimum valid observations required after common masking
- Features with `n_valid < min_n` are skipped
- Effect: Higher values exclude features with insufficient data

### Output Arguments

**`top_n`** (default: `4`)
- Type: Integer > 0
- Description: Number of top families to return
- Families ranked by frequency across bootstrap samples

**`verbose`** (default: `TRUE`)
- Type: Logical
- Description: Whether to print progress messages and summary
- TRUE: Displays pull-by-pull progress and final report
- FALSE: Silent execution

---

## Return Value

List with eight elements:

### 1. `top_families_overall`
- Type: Character vector (length = `top_n`)
- Description: Most frequent families across all bootstrap samples
- Ordering: Descending by frequency

### 2. `top_families_by_support`
- Type: Named list with elements `count`, `unit`, `positive`, `real`
- Description: Top families within each empirical support category
- Length: Up to `top_n` per category (may be fewer if insufficient data)

### 3. `freq_table_overall`
- Type: Named integer vector (table object)
- Description: Absolute frequencies of each selected family
- Names: Family identifiers
- Values: Selection counts across all bootstrap samples

### 4. `prop_table_overall`
- Type: Named numeric vector
- Description: Proportions of each selected family
- Values: Frequencies normalized by total successful fits
- Range: [0, 1], sum ≤ 1

### 5. `freq_by_support`
- Type: Named list with elements `count`, `unit`, `positive`, `real`
- Description: Frequency tables stratified by support category
- Each element: Named integer vector (table object)

### 6. `prop_by_support`
- Type: Named list with elements `count`, `unit`, `positive`, `real`
- Description: Proportion tables stratified by support category
- Each element: Named numeric vector

### 7. `sampled_results`
- Type: Tibble with columns:
  - `bootstrap` (integer): Pull index (1 to `n_boot`)
  - `feature` (character): Feature identifier
  - `family` (character): Best-fitting family (NA if skipped)
  - `skipped` (logical): Whether feature was excluded
  - `n_valid` (integer): Valid observations after common masking
  - `support` (character): Empirical support classification
- Rows: `n_genes × n_boot` (one per feature per pull)

### 8. `transform_mode`
- Type: Character scalar (`"strict"` or `"safe"`)
- Description: Transformation mode used in analysis
- Value: Resolved from `transform_mode` argument or auto-default

---

## Examples

### Basic Usage

```r
library(perseo)

# Suppose counts_matrix is a simulated gene expression dataset
ff <- find_families(
  counts_matrix = counts_matrix,
  n_genes = 100,
  n_boot = 5,
  top_n = 4,
  criterion = "BIC",
  seed = 123
)

# Examine results
ff$top_families_overall
#> [1] "NBI"   "LOGNO" "GG"    "TF"

ff$freq_table_overall
#> NBI LOGNO    GG    TF    GA    NO 
#> 125    87    76    43    22    12

ff$prop_table_overall
#>      NBI    LOGNO       GG       TF       GA       NO 
#> 0.342466 0.238356 0.208219 0.117808 0.060274 0.032877
```

### Custom Family Panel

```r
# Test only negative binomial variants and normal
custom_families <- c("PO", "NBI", "ZIP", "ZINBI", "NO", "TF")

ff_custom <- find_families(
  counts_matrix = counts_matrix,
  n_genes = 150,
  n_boot = 10,
  families = custom_families,
  criterion = "BIC",
  min_n = 20,
  seed = 456
)
```

### Support-Agnostic Selection

```r
# Test all families regardless of empirical support
# Uses "safe" transformation mode by default
ff_agnostic <- find_families(
  counts_matrix = counts_matrix,
  n_genes = 100,
  n_boot = 5,
  group_by_support = FALSE,
  transform_mode = "safe",  # Explicit (optional)
  criterion = "BIC",
  seed = 789
)

# Compare with support-aware selection
ff_aware <- find_families(
  counts_matrix = counts_matrix,
  n_genes = 100,
  n_boot = 5,
  group_by_support = TRUE,
  transform_mode = "strict",  # Explicit (optional)
  criterion = "BIC",
  seed = 789
)

# Different families may be selected
ff_agnostic$top_families_overall
ff_aware$top_families_overall
```

### Stratified Results by Support

```r
ff <- find_families(
  counts_matrix = counts_matrix,
  n_genes = 200,
  n_boot = 10,
  criterion = "BIC",
  seed = 101
)

# Top families for count data
ff$top_families_by_support$count
#> [1] "NBI" "PO"  "ZIP"

# Top families for positive continuous
ff$top_families_by_support$positive
#> [1] "LOGNO" "GG"    "GA"

# Frequency distribution within count support
ff$freq_by_support$count
#> NBI  PO ZIP 
#> 142  78  23
```

---

## Statistical Details

### Information Criterion Selection

**AIC (Akaike Information Criterion)**
```
AIC = -2 log L(θ̂; y) + 2·df
```
- Penalty: Linear in parameters
- Asymptotic behavior: Minimizes Kullback-Leibler divergence
- Tendency: Favors more complex models

**BIC (Bayesian Information Criterion)**
```
BIC = -2 log L(θ̂; y) + log(n)·df
```
- Penalty: Logarithmic in sample size
- Asymptotic behavior: Consistent model selector
- Tendency: Stronger parsimony than AIC for `n > 7`

**GAIC (Generalized AIC)**
```
GAIC = -2 log L(θ̂; y) + k·df
```
- Penalty: User-specified `k`
- Special cases: `k = 2` (AIC), `k = log(n)` (BIC)
- Flexibility: Allows intermediate penalty strengths

### Effective Sample Size

After applying common mask, effective sample size is:

```
n_valid = Σⱼ 𝟙(mask_common[j])
```

where `mask_common = ⋂_{f ∈ ℱ} mask_f`.

Features with `n_valid < min_n` are excluded to ensure:
- Sufficient degrees of freedom for estimation
- Reliable IC comparisons
- Stable variance estimates

### Empirical Support Classification

Based on `infer_support()` function:

1. **Count**: All values are integers ≥ 0 (tolerance: `|y - round(y)| < 1e-8`)
2. **Unit**: All values in [0, 1]
3. **Positive**: All values > 0 (not count or unit)
4. **Real**: Contains negative values or other

### Bootstrap Aggregation

Family frequencies are computed as:

```
freq(f) = Σ_{b=1}^B Σ_{i=1}^m 𝟙(best_family(yᵢ, b) = f)
```

Proportions:

```
prop(f) = freq(f) / Σ_{f'} freq(f')
```

---

## Computational Considerations

### Time Complexity

For `F` families, `m` features per pull, `B` pulls, `n` samples:
- Per-family model fit: O(n·p²) where p = number of parameters
- Per-feature comparison: O(F·n·p²)
- Total: O(B·m·F·n·p²)

### Memory Requirements

- Primary storage: `counts_matrix` (features × samples)
- Intermediate: Transformed data per family (m × n × F per pull)
- Results: `sampled_results` tibble (m·B rows)

### Parallelization

Not currently parallelized within `find_families()`. For large datasets:
- Can be run on subsets and results combined
- Future implementation may support parallel bootstrap pulls

---

## Relationship to Other Functions

### Upstream

- **Data preparation**: User responsibility (QC, normalization, filtering)
- **Design matrix**: Not used (intercept-only models)

### Downstream

1. **`fit_gamlss_models()`**: Uses selected families for full regression
2. **`run_perseo()`**: Orchestrates `find_families()` → `fit_gamlss_models()` workflow

### Integration Example

```r
# Step 1: Family selection (intercept-only models)
ff <- find_families(
  counts_matrix = counts,
  n_genes = 200,
  n_boot = 10,
  criterion = "BIC",
  seed = 123
)

# Examine selected families
ff$top_families_overall
#> [1] "GG"    "LOGNO" "NBI"   "GA"

# Step 2: Differential expression with selected families and covariates
# WORKFLOW A: Formula + automatic contrasts
de_results <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = "~ condition + age + batch",  # formula with covariates
  metadata = metadata,                           # required for formula
  candidate_families = ff$top_families_overall,  # use selected families
  contrast_variable = "condition",               # auto-generate contrasts
  criterion = "BIC",
  transform_mode = ff$transform_mode,            # use same transform mode
  parallel = TRUE,
  workers = 4
)

# OR WORKFLOW B: Design matrix + manual contrasts
design <- model.matrix(~ condition + age + batch, data = metadata)
C <- matrix(...)  # custom contrast matrix
colnames(C) <- colnames(design)

de_results <- fit_gamlss_models(
  counts_matrix = counts,
  design_matrix = design,                        # pre-built matrix
  candidate_families = ff$top_families_overall,  # use selected families
  contrast_matrix = C,                           # explicit contrasts
  criterion = "BIC",
  transform_mode = ff$transform_mode,
  parallel = TRUE,
  workers = 4
)

# Access results
head(de_results$results)     # coefficient statistics
head(de_results$contrasts)   # contrast results
head(de_results$selection)   # best family per feature
```

**Key Points**:
- `find_families()` uses **intercept-only models** to identify robust families
- `fit_gamlss_models()` uses **full regression models** with covariates
- Same families are tested, but with added complexity (covariates)
- Use same `transform_mode` in both steps for consistency
- See [README](../README.md) for complete workflow documentation

---

## Interpretation Guidelines

### Family Frequency Interpretation

- **High frequency** (>50%): Family fits well across diverse features
- **Moderate frequency** (20-50%): Family suitable for specific feature subsets
- **Low frequency** (<20%): Family rarely optimal, may be redundant

Frequencies do not indicate:
- Absolute model quality (only relative within candidates)
- Feature-level goodness-of-fit (use residual diagnostics)
- Statistical significance (selection is descriptive)

### Support-Stratified Results

When `group_by_support = TRUE`:
- Within-support rankings show family preferences given domain constraints
- Cross-support comparisons are not meaningful (different feature sets)
- Useful for identifying dominant family types per data modality

### Transformation Mode Effects

**Strict mode** (`group_by_support = TRUE` default):
- Conservative: Only domain-compatible observations used
- May exclude features with boundary values
- IC comparisons on consistent validity masks

**Safe mode** (`group_by_support = FALSE` default):
- Inclusive: Global transformations retain all observations
- Allows cross-support family exploration
- IC comparisons on transformed scales (Jacobian-corrected)

---

## Limitations and Caveats

1. **Intercept-only models**: Selection based on marginal distributions
   - Does not account for covariate effects
   - Best family may differ when covariates are included
   - This is intentional: identifies robust families across diverse features
   - Full regression happens in `fit_gamlss_models()` with the selected families
   - **Rationale**: Marginal family selection avoids overfitting to specific covariate structures

2. **Bootstrap sampling**: 
   - Assumes features are exchangeable within support
   - May undersample rare feature types
   - Frequency estimates have sampling variability
   - Use `bootstrap = FALSE` in `run_perseo()` for comprehensive evaluation

3. **Common mask**:
   - Can reduce effective sample size substantially
   - May exclude features with many boundary values
   - Bias toward families with less restrictive domains (in strict mode)
   - Ensures valid IC comparisons (all families on same observations)

4. **IC limitations**:
   - Asymptotic approximations (may be inaccurate for small `n_valid`)
   - Assumes correct model specification within family
   - Does not account for model misspecification uncertainty
   - Relative comparison (not absolute quality measure)

5. **Computational cost**:
   - Scales as O(`n_boot` × `n_genes` × `n_families`)
   - Large family panels may be slow without parallelization
   - Use `parallel = TRUE` in `run_perseo()` for automatic parallelization

---

## References

1. Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for location, scale and shape. *Journal of the Royal Statistical Society: Series C*, 54(3), 507-554.

2. Akaike, H. (1974). A new look at the statistical model identification. *IEEE Transactions on Automatic Control*, 19(6), 716-723.

3. Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461-464.

4. Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference: Understanding AIC and BIC in model selection. *Sociological Methods & Research*, 33(2), 261-304.

---