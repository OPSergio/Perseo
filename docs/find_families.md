# `find_families`: Identify Best-Fitting GAMLSS Families Across Features

## Description

`find_families()` is a utility designed to identify the most plausible **distribution families** for a set of features 
(e.g. genes, OTUs, proteins, or clinical markers) in high-dimensional data.  
Instead of fitting all features at once, it performs **bootstrap sampling** of features and evaluates 
intercept-only GAMLSS models across candidate families.  

The function corrects log-likelihoods back to the original data scale using the **Jacobian of the transformation**, 
and ensures fair comparisons by applying a **common mask** (only observations valid across all tested families).  

It is primarily used to pre-select a **reduced set of families** that will later be fitted in full regression 
models (e.g. `fit_gamlss_models()`).  

---

## Arguments

- **counts_matrix**:  
  Numeric matrix (features × samples).  
  Rows represent features, columns represent samples. Row names must be feature identifiers.

- **n_genes** *(default: 200)*:  
  Number of features to sample in each bootstrap pull.

- **n_boot** *(default: 10)*:  
  Number of bootstrap pulls. Increasing this improves stability but increases runtime.

- **top_n** *(default: 4)*:  
  Number of top frequent families to return.

- **families** *(default: NULL)*:  
  Character vector of candidate families.  
  If `NULL`, a default panel is used, including:
  - **Unit interval**: BE, BEO, BEINF, BEZI, BEo, BEINF0  
  - **Counts**: PO, NBI, ZIP, ZINBI, BI, BB  
  - **Positive continuous**: GA, GG, IG, LOGNO  
  - **Real-valued**: NO, TF, GU  

- **verbose** *(default: TRUE)*:  
  Whether to print progress and summary reports.

- **min_n** *(default: 5)*:  
  Minimum number of valid samples required after applying the common mask.

- **seed** *(default: NULL)*:  
  Random seed for reproducibility.

- **group_by_support** *(default: TRUE)*:  
  Whether to restrict candidate families to those compatible with the **empirical support** of the feature:
  - `count` (0,1,2,...) → Poisson-like  
  - `unit` (0–1 bounded) → Beta-like  
  - `positive` (>0) → Gamma-like  
  - `real` (any real) → Normal-like  

- **binom_bd** *(default: NULL)*:  
  Binomial denominator for BI/BB families.  
  - If `NULL`, inferred per feature as `max(y)` when consistent.  
  - Can also be provided as a scalar or vector of length = number of samples.

---

## Returns

A list with:

- **top_families_overall**:  
  Character vector of the most frequent families across all bootstrap pulls.

- **top_families_by_support**:  
  List with top families per empirical support category (`count`, `unit`, `positive`, `real`).

- **freq_table_overall**:  
  Frequency table of families across all pulls.

- **prop_table_overall**:  
  Proportions of families across all pulls.

- **freq_by_support / prop_by_support**:  
  Tables of family frequencies and proportions stratified by support type.

- **sampled_results**:  
  Tibble with row-level results per feature and bootstrap pull.  
  Includes:
  - `feature` : Feature ID  
  - `family` : Best-fitting family  
  - `skipped` : Whether the feature was skipped  
  - `n_valid` : Valid sample size after masking  
  - `support` : Empirical support class  
  - `bootstrap` : Pull index  

---

## Example

```r
set.seed(123)

# Suppose counts_matrix is a simulated gene expression dataset
ff <- find_families(
  counts_matrix   = counts_matrix,
  n_genes         = 50,     # features per pull
  n_boot          = 5,      # number of pulls
  top_n           = 4,
  verbose         = TRUE,
  min_n           = 10,
  group_by_support = TRUE
)

# Inspect top families
ff$top_families_overall
ff$top_families_by_support
```

---

## Notes

- **Why GAIC / BIC?**  
  AIC tends to favor **more flexible families** (e.g. BEINF over BE, NBI over PO).  
  We use GAIC with `k = log(n)` (equivalent to BIC) to penalize unnecessary complexity.

- **Jacobian correction**  
  Because we transform data differently per family (e.g. min-max for Beta, log for Gamma),  
  log-likelihoods are adjusted by adding the Jacobian term to ensure comparability.

- **Common mask**  
  Ensures that families are compared on the **same subset of observations**, avoiding bias 
  when some transformations drop more data than others.

---