# Transformations Utilities (`utils_transformations.R`)

This module implements variable transformation functions required to make expression data compatible with GAMLSS families.  

The main goal is to ensure that each response vector `y` is projected into the correct domain of each distribution, enabling:

- **Diagnostic transformations (SAFE)** for visualization or robust preprocessing.  
- **Strict transformations (STRICT)** with Jacobian correction for fair family comparison in the selection process (`find_families`).  
- **Helper functions** to infer support, group families, and reverse transformations.  

---

## 1. `transform_for_family()`

*SAFE* transformation for a given vector `y` and a family `fam`.

- **Use cases**:  
  - Visualization.  
  - Robust preprocessing in exploratory pipelines.  

- **Do not use**: for family comparison/selection (since it may apply clipping, rounding, or rescaling that alters likelihood).

### Arguments
- `y`: numeric vector of expression values.  
- `fam`: GAMLSS family name (`"PO"`, `"GA"`, `"BE"`, `"NO"`, etc.).  
- `strategy`: `"safe"` (default) → clipping/rescaling; `"strict"` → replaces invalid values with `NA`.  
- `eps`: small numeric constant for smoothing/clipping.  

### Example
```r
y <- c(-1, 0, 0.5, 2)
transform_for_family(y, "GA")   # => positive values (0.000001, 0, 0.5, 2)
transform_for_family(y, "BE")   # => scaled to [0,1] with clipping
```

---

## 2. `transform_for_family_strict()`

*STRICT* transformation for fair family selection.  

Returns:
- `y`: vector transformed into the working scale `z = g(y)`.  
- `mask`: valid observations according to theoretical domain.  
- `logJ_per_obs`: per-observation log Jacobian.  
- `meta`: metadata of the transformation (`kind`, parameters).  

### Cases
- **Counts (`PO`, `NBI`, `ZIP`, `ZINBI`, `BI`, `BB`)**: identity, valid if integer ≥ 0.  
- **Positive continuous (`GA`, `GG`, `LOGNO`, `IG`)**: identity, valid if > 0.  
- **Unit interval (`BE`, `BEINF`, `BEZI`, `BEO`, etc.)**: min–max scaling to (0,1), Jacobian constant `-log(b - a)`.  
- **Real values (`NO`, `TF`, `GU`)**: z-score standardization, Jacobian `-log(sd)`.  

### Example
```r
res <- transform_for_family_strict(c(0.2, 0.5, 0.9), "BE")
res$y             # scaled to (0,1)
res$mask          # TRUE for valid observations
res$logJ_per_obs  # Jacobian values
```

---

## 3. `inverse_transform()`

Inverse transformation to map data back from the working scale (`z`) to the original scale (`y`).

- **Input**:  
  - `z`: transformed vector.  
  - `meta`: metadata from `transform_for_family_strict()`.  
- **Output**: numeric vector in the original scale.  

### Example
```r
res <- transform_for_family_strict(c(0.2, 0.5, 0.9), "BE")
inverse_transform(res$y, res$meta)  # returns ~ (0.2, 0.5, 0.9)
```

---

## 4. `family_groups()`

Returns a list of families grouped by theoretical support:
- `count`: count distributions (e.g. `PO`, `NBI`).  
- `unit`: proportions (e.g. `BE`, `BEINF`).  
- `positive`: positive continuous (e.g. `GA`, `GG`).  
- `real`: real-valued (e.g. `NO`, `TF`).  

---

## 5. `infer_support()`

Infers the empirical support of a vector `y`. Returns one of:  
- `"count"`  
- `"unit"`  
- `"positive"`  
- `"real"`  
- `"none"` (if no valid data)  

---

## 6. `jacobian_sum()`

Computes the sum of the log-Jacobian over a set of observations.

- **Use case**: correction of log-likelihood in `find_families`.  
- **Input**:  
  - `logJ_per_obs`: vector from `transform_for_family_strict()`.  
  - `mask`: optional logical vector.  
- **Output**: numeric scalar with the sum.  

---

## Recommended Workflow

1. For exploratory analysis or visualization → use `transform_for_family()`.  
2. For family selection → use `transform_for_family_strict()` + `jacobian_sum()`.  
3. For reconstructing simulated data or visualization back on the original scale → use `inverse_transform()`.  
