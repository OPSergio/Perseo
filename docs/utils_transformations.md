# Transformations Utilities (`utils_transformations.R`)

## Overview

This module implements transformation functions for adapting response data to GAMLSS family domains. The transformation system supports two modes:

- **Strict mode**: Domain-preserving transformations with observation exclusion via masking
- **Safe mode**: Global affine transformations allowing reversible domain adaptation

All transformations include Jacobian correction to ensure information criterion comparability across families.

---

## Core Function: `transform_response()`

Unified transformation dispatcher supporting both strict and safe modes.

### Arguments

- `y`: Numeric vector of response values
- `fam`: Character scalar specifying GAMLSS family name (e.g., `"PO"`, `"GA"`, `"BE"`, `"NO"`)
- `mode`: Character scalar, either `"strict"` (default) or `"safe"`
- `eps`: Numeric scalar for epsilon handling in boundary domains (default: `1e-6`)
- `allow_eps`: Logical; whether to apply epsilon adjustments at boundaries (strict mode only)

### Return Value

List with five elements:
- `y`: Numeric vector of transformed response values
- `mask`: Logical vector indicating valid observations
- `logJ_per_obs`: Numeric vector of per-observation log-Jacobian values
- `meta`: List containing transformation metadata (`kind`, `params`)
- `mode_used`: Character scalar indicating applied transformation mode

### Transformation Modes

#### Strict Mode (`mode = "strict"`)

Enforces theoretical domain constraints without data modification:

- **Count families** (`PO`, `NBI`, `ZIP`, `ZINBI`, `ZIP2`, `BI`, `BB`):
  - Transform: Identity `z = y`
  - Validity: `y ≥ 0` and `y ∈ ℤ`
  - Jacobian: `log|∂z/∂y| = 0`

- **Positive continuous families** (`GA`, `GG`, `LOGNO`, `IG`):
  - Transform: Identity `z = y`
  - Validity: `y > 0`
  - Jacobian: `log|∂z/∂y| = 0`

- **Unit interval families** (`BE`, `BEINF`, `BEO`, `BEZI`, `BEo`, `BEINF0`):
  - Transform: Min-max scaling `z = (y - a)/(b - a)` where `a = min(y)`, `b = max(y)`
  - Validity: Depends on inflation variant (`0 < z < 1` for `BE`, `0 ≤ z < 1` for `BEINF`, etc.)
  - Jacobian: `log|∂z/∂y| = -log(b - a)` (constant across observations)

- **Real-valued families** (`NO`, `TF`, `GU`):
  - Transform: Z-score standardization `z = (y - μ)/σ` where `μ = mean(y)`, `σ = sd(y)`
  - Validity: All finite `y`
  - Jacobian: `log|∂z/∂y| = -log(σ)` (constant across observations)

#### Safe Mode (`mode = "safe"`)

Applies global affine transformations `z = ay + b` with `a > 0`:

- **Positive continuous families** (`GA`, `GG`, `LOGNO`, `IG`):
  - If `min(y) ≤ 0`: Apply shift `b = -min(y) + ε`, `a = 1`
  - Otherwise: Identity `a = 1`, `b = 0`
  - Result: `z ∈ (ε, +∞)`
  - Jacobian: `log|∂z/∂y| = log(a) = 0`

- **Unit interval families** (`BE`, `BEINF`, `BEO`, `BEZI`, `BEo`, `BEINF0`):
  - Without epsilon: `a = 1/(max(y) - min(y))`, `b = -min(y) · a`
  - With epsilon (when 0 or 1 excluded): `a = (1 - 2ε)/(max(y) - min(y))`, `b = ε - min(y) · a`
  - Result: `z ∈ [0, 1]` or `z ∈ [ε, 1-ε]` depending on family
  - Jacobian: `log|∂z/∂y| = log(a)` (constant across observations)

- **Real-valued families** (`NO`, `TF`, `GU`):
  - Same as strict mode: Z-score standardization
  - Jacobian: `log|∂z/∂y| = -log(σ)`

- **Count families** (`PO`, `NBI`, `ZIP`, `ZINBI`, `ZIP2`, `BI`, `BB`):
  - Transform: Identity `z = y` (no rounding applied)
  - Validity: All finite `y`
  - Jacobian: `log|∂z/∂y| = 0`

### Mathematical Properties

1. **Jacobian correction**: For any monotonic transformation `z = g(y)`, the log-likelihood on the original scale is:
   ```
   log L(θ; y) = log L(θ; z) + Σᵢ log|∂zᵢ/∂yᵢ|
   ```

2. **Affine transforms**: For `z = ay + b` with `a > 0`:
   - Jacobian: `∂z/∂y = a` (constant)
   - Reversibility: `y = (z - b)/a`

3. **Information criterion adjustment**: IC values become comparable across families by adding `-2 · Σᵢ log|∂zᵢ/∂yᵢ|`

### Examples

```r
# Strict mode: excludes invalid observations
y <- c(-1, 0, 1, 2, 3)
res_strict <- transform_response(y, "GA", mode = "strict")
# res_strict$mask: c(FALSE, FALSE, TRUE, TRUE, TRUE)
# res_strict$y[1:2]: excluded (epsilon or NA)

# Safe mode: global shift to positive domain
res_safe <- transform_response(y, "GA", mode = "safe")
# res_safe$mask: c(TRUE, TRUE, TRUE, TRUE, TRUE)
# res_safe$y: c(1e-6, 1+1e-6, 2+1e-6, 3+1e-6, 4+1e-6)
# res_safe$meta$kind: "affine"
# res_safe$meta$params: list(a = 1, b = 1+1e-6)

# Reversibility in safe mode
y_back <- inverse_transform(res_safe$y, res_safe$meta)
# y_back ≈ y (within numerical tolerance)
```

---

## Legacy Function: `transform_for_family()`

**Status**: Legacy. Use `transform_response()` for new code.

Simplified transformation function without Jacobian correction. Applies observation-wise clipping and rounding.

### Arguments
- `y`: Numeric vector
- `fam`: GAMLSS family name
- `strategy`: `"safe"` (default) or `"strict"`
- `eps`: Epsilon value

### Behavior
- `strategy = "safe"`: Clips/rounds values to fit domain (e.g., `y[y ≤ 0] <- eps` for positive families)
- `strategy = "strict"`: Replaces invalid values with `NA`

### Limitations
- No Jacobian correction
- Observation-wise operations (non-global)
- Not suitable for family comparison

---

## Supporting Function: `transform_for_family_strict()`

**Status**: Internal. Called by `transform_response()` when `mode = "strict"`.

Direct implementation of strict mode transformations. Returns same structure as `transform_response()` but without `mode_used` field.

---

## Function: `inverse_transform()`

Reverses transformations applied by `transform_response()` or `transform_for_family_strict()`.

### Mathematical Formulation

Given transformed data `z` and metadata `meta`, recovers original scale `y`:

- **Identity** (`kind = "identity"`): `y = z`
- **Z-score** (`kind = "zscore"`): `y = μ + σz` where `μ = meta$params$center`, `σ = meta$params$scale`
- **Min-max** (`kind = "minmax"`): `y = a + z(b - a)` where `a = meta$params$min`, `b = meta$params$max`
- **Affine** (`kind = "affine"`): `y = (z - b)/a` where `a = meta$params$a`, `b = meta$params$b`

### Arguments
- `z`: Numeric vector on transformed scale
- `meta`: List with elements `kind` (character) and `params` (list)

### Return Value

Numeric vector on original scale

### Examples

```r
# Z-score inversion
y <- c(-2, -1, 0, 1, 2)
res <- transform_response(y, "NO", mode = "strict")
y_back <- inverse_transform(res$y, res$meta)
all.equal(y, y_back)  # TRUE (within tolerance)

# Affine inversion (SAFE mode)
y <- c(-1, 0, 1, 2, 3)
res <- transform_response(y, "GA", mode = "safe")
y_back <- inverse_transform(res$y, res$meta)
all.equal(y, y_back)  # TRUE

# Min-max inversion
y <- c(0.1, 0.3, 0.5, 0.7, 0.9)
res <- transform_response(y, "BE", mode = "strict")
y_back <- inverse_transform(res$y, res$meta)
all.equal(y, y_back)  # TRUE
```

---

## Function: `family_groups()`

Returns GAMLSS families partitioned by theoretical support.

### Return Value

List with four named elements:
- `count`: Character vector of count families (`"PO"`, `"NBI"`, `"ZIP"`, `"ZINBI"`, `"ZIP2"`, `"BI"`, `"BB"`)
- `unit`: Character vector of unit interval families (`"BE"`, `"BEINF"`, `"BEO"`, `"BEZI"`, `"BEo"`, `"BEINF0"`)
- `positive`: Character vector of positive continuous families (`"GA"`, `"GG"`, `"LOGNO"`, `"IG"`)
- `real`: Character vector of real-valued families (`"NO"`, `"TF"`, `"GU"`)

### Mathematical Domains

- **Count**: `𝕊 = ℤ₊ ∪ {0} = {0, 1, 2, ...}`
- **Unit**: `𝕊 = (0, 1)` or `𝕊 = [0, 1]` depending on inflation variant
- **Positive**: `𝕊 = ℝ₊ = (0, +∞)`
- **Real**: `𝕊 = ℝ = (-∞, +∞)`

### Example

```r
groups <- family_groups()
groups$count     # "PO" "NBI" "ZIP" "ZINBI" "ZIP2" "BI" "BB"
groups$positive  # "GA" "GG" "LOGNO" "IG"
```

---

## Function: `infer_support()`

Determines empirical support of a numeric vector.

### Algorithm

1. Remove non-finite values
2. Check if all values are integers ≥ 0 → `"count"`
3. Else if all values in `[0, 1]` → `"unit"`
4. Else if all values > 0 → `"positive"`
5. Otherwise → `"real"`

### Arguments
- `y`: Numeric vector

### Return Value

Character scalar: `"count"`, `"unit"`, `"positive"`, `"real"`, or `"none"` (if no valid data)

### Examples

```r
infer_support(c(0, 1, 2, 5, 10))          # "count"
infer_support(c(0.1, 0.3, 0.5, 0.9))      # "unit"
infer_support(c(1.2, 3.5, 7.8))           # "positive"
infer_support(c(-2, -1, 0, 1, 2))         # "real"
infer_support(c(NA, NaN, Inf))            # "none"
```

### Limitations

- Threshold for integer detection: `|y - round(y)| < 1e-8`
- Empirical support may differ from theoretical support (e.g., data with small positive values might be count data rounded to avoid zeros)

---

## Function: `jacobian_sum()`

Computes total log-Jacobian contribution over valid observations.

### Mathematical Definition

For transformation `z = g(y)` with Jacobian `J(y) = ∂z/∂y`:

```
jacobian_sum = Σᵢ log|J(yᵢ)| = Σᵢ log|∂zᵢ/∂yᵢ|
```

Used to adjust log-likelihood:

```
log L(θ; y) = log L(θ; z) + jacobian_sum
```

### Arguments
- `logJ_per_obs`: Numeric vector of per-observation log-Jacobians (from `transform_response()`)
- `mask`: Optional logical vector indicating which observations to include (default: all finite values)

### Return Value

Numeric scalar

### Examples

```r
y <- c(1, 2, 3, 4, 5)
res <- transform_response(y, "GA", mode = "strict")

# Sum over all observations
jacobian_sum(res$logJ_per_obs)

# Sum over masked observations only
mask <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
jacobian_sum(res$logJ_per_obs, mask)
```

---

## Workflow Recommendations

### Family Selection Workflow

1. **Data preparation**: Ensure `y` contains numeric finite values
2. **Transform**: Use `transform_response(y, fam, mode = "strict")` for each candidate family
3. **Common mask**: Compute `mask_common = Reduce('&', list_of_masks)`
4. **Subset data**: `z_masked = z[mask_common]` for each family
5. **Fit models**: Fit GAMLSS on `z_masked` with design matrix
6. **Compute IC**: Adjust with `-2 * jacobian_sum(logJ_per_obs, mask_common)`
7. **Select best**: Choose family minimizing adjusted IC

### Exploratory Analysis Workflow

1. **Transform**: Use `transform_response(y, fam, mode = "safe")` to adapt data to family domain
2. **Fit model**: Fit GAMLSS on transformed data
3. **Predictions**: Generate predictions on transformed scale
4. **Back-transform**: Use `inverse_transform(predictions, meta)` to obtain original scale

### Visualization Workflow

1. For direct plotting: Use `transform_response()` with `mode = "safe"` or legacy `transform_for_family()`
2. For diagnostic plots: Use transformed scale directly (e.g., Q-Q plots on z-score scale)
3. For presentation: Back-transform to original scale using `inverse_transform()`

---

## Statistical Notes

### Jacobian Correction Necessity

When comparing models fitted to transformed data, information criteria (AIC, BIC, GAIC) must account for the transformation Jacobian to be valid on the original data scale. Without correction, IC values are incomparable across families using different transformations.

**Proof**: For transformation `z = g(y)` with density `f_Z(z; θ)`, the density on original scale is:

```
f_Y(y; θ) = f_Z(g(y); θ) · |J(y)|
```

Therefore:

```
log f_Y(y; θ) = log f_Z(g(y); θ) + log|J(y)|
```

Summing over observations:

```
log L(θ; y) = log L(θ; z) + Σᵢ log|J(yᵢ)|
```

### Affine Transform Properties

For `z = ay + b` with `a > 0`:

1. **Linearity preservation**: Differences are scaled by `a`
2. **Monotonicity**: Order preservation (strictly increasing)
3. **Bijectivity**: One-to-one mapping
4. **Constant Jacobian**: `∂z/∂y = a` for all `y`

### Common Mask Rationale

Using a common mask across all candidate families ensures:

1. **Fair comparison**: All models fitted to same observations
2. **Consistent sample size**: IC penalties comparable
3. **Valid inference**: Standard errors account for actual data used

Alternative approach (family-specific masks) would confound model fit quality with data availability, biasing selection toward families with less restrictive domains.

---

## Appendix: Family Specifications

### Count Families

| Family | Full Name | Domain | Parameters | Inflation |
|--------|-----------|--------|------------|-----------|
| PO | Poisson | ℤ₊ ∪ {0} | μ | No |
| NBI | Negative Binomial (Type I) | ℤ₊ ∪ {0} | μ, σ | No |
| ZIP | Zero-Inflated Poisson | ℤ₊ ∪ {0} | μ, σ | At 0 |
| ZINBI | Zero-Inflated NBI | ℤ₊ ∪ {0} | μ, σ, ν | At 0 |
| ZIP2 | Zero-Inflated Poisson (alt) | ℤ₊ ∪ {0} | μ, σ | At 0 |
| BI | Binomial | {0, 1, ..., n} | μ, bd | No |
| BB | Beta-Binomial | {0, 1, ..., n} | μ, σ, bd | No |

### Unit Interval Families

| Family | Full Name | Domain | Inflation |
|--------|-----------|--------|-----------|
| BE | Beta | (0, 1) | No |
| BEINF | Beta Inflated | [0, 1] | At 0 and 1 |
| BEO | Beta Inflated at One | (0, 1] | At 1 |
| BEZI | Beta Inflated at Zero | [0, 1) | At 0 |
| BEo | One-Inflated Beta (alt) | (0, 1] | At 1 |
| BEINF0 | Zero-One Inflated Beta (alt) | [0, 1] | At 0 and 1 |

### Positive Continuous Families

| Family | Full Name | Domain | Skewness |
|--------|-----------|--------|----------|
| GA | Gamma | (0, +∞) | Right |
| GG | Generalized Gamma | (0, +∞) | Flexible |
| LOGNO | Log-Normal | (0, +∞) | Right |
| IG | Inverse Gaussian | (0, +∞) | Right |

### Real-Valued Families

| Family | Full Name | Domain | Tails |
|--------|-----------|--------|-------|
| NO | Normal (Gaussian) | ℝ | Light |
| TF | t Family | ℝ | Heavy |
| GU | Gumbel | ℝ | Asymmetric |
