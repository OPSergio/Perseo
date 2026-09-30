# PERSEO API reference

All exported functions, with signatures copied from the source. Defaults are shown inline.

## run_perseo(): end-to-end pipeline (find_families, then fit_gamlss_models)

```r
run_perseo(counts_matrix, design_matrix, contrast_matrix = NULL,
  bootstrap = TRUE, n_genes = 200, n_boot = 10, top_n = 4, families = NULL,
  group_by_support = TRUE, criterion = c("GAIC","BIC","AIC"), gaic_k = NULL, min_n = 5,
  binom_bd = NULL, filter_beta_inflated = TRUE, p_adjust_method = "BH",
  transform_mode = NULL, show_progress = TRUE, seed = NULL, metadata = NULL,
  contrast_variable = NULL, omnibus = FALSE, omnibus_threshold = 0.05,
  omnibus_test = c("Wald","LRT"), parallel = FALSE, workers = NULL,
  thr_zero = 0.005, thr_one = 0.005)
```

| arg | meaning |
|---|---|
| `counts_matrix` | numeric `matrix`, features x samples, rownames = feature IDs |
| `design_matrix` | formula string `"~ a + b"`, a `formula`, or a numeric matrix from `model.matrix()` (samples x coefs). A formula needs `metadata` and may use gamlss smoothers/random effects (`pb(age)`, `random(subject)`); `offset()` is dropped |
| `metadata` | data.frame, `nrow = ncol(counts)`, **same row order as the counts columns** |
| `contrast_variable` | factor column in `metadata`; all pairwise contrasts, named `B_vs_A` |
| `contrast_matrix` | numeric matrix, rows = contrasts (rownames = names), colnames = coefficient names. Takes precedence over `contrast_variable` |
| `bootstrap` | `TRUE`: `n_boot` pulls of `n_genes` random features (needs `nrow >= n_genes`). `FALSE`: every feature votes (slow) |
| `top_n` | number of most frequently winning families passed to the DE step |
| `families` | candidate set for selection (`NULL` = 21 defaults, see below) |
| `group_by_support` | restrict families to each feature's empirical support |
| `criterion`, `gaic_k` | IC. GAIC with `gaic_k = NULL` uses `k = log(n)`, which equals BIC |
| `min_n` | minimum valid observations after the common mask, otherwise the feature is skipped |
| `binom_bd` | BI/BB denominator: `NULL` = per-feature `max(y)`, scalar, or a per-sample vector |
| `filter_beta_inflated`, `thr_zero`, `thr_one` | keep only the families whose point masses match the observed exact 0s/1s (share >= threshold). For unit data: BE without 0/1, BEZI/BEINF0 with zeros, BEINF1 with ones, BEINF with both (see methods.md) |
| `p_adjust_method` | any `p.adjust` method |
| `transform_mode` | `NULL` gives `"strict"` if `group_by_support`, else `"safe"` |
| `omnibus*` | gate contrasts on a joint test of the factor's coefficients (needs contrasts). `LRT` refits the reduced model (~2x slower) |
| `parallel`, `workers` | `future::multisession`. `workers = NULL` uses all cores minus 1. Each worker copies the data (RAM) |

Returns class `perseo_results`: `$family_selection` (find_families output),
`$differential_expression` (fit_gamlss_models output), `$summary` (parameters, `status`,
`models_fitted`, `families_selected`), and `$input_data` (`counts_matrix`, `metadata`,
`contrast_variable`). If no family is selected, `$differential_expression` is `NULL` and
`$summary$status == "failed_family_selection"`.

## fit_gamlss_models(): DE with a known family set

```r
fit_gamlss_models(counts_matrix, design_matrix, metadata = NULL, candidate_families,
  criterion = c("GAIC","BIC","AIC"), gaic_k = NULL, min_n = 5,
  contrast_matrix = NULL, contrast_variable = NULL,
  omnibus = FALSE, omnibus_threshold = 0.05, omnibus_test = c("Wald","LRT"),
  p_adjust = "BH", workers = NULL, parallel = FALSE, show_progress = TRUE,
  progress_label = "Fitting features", transform_mode = "strict",
  group_by_support = FALSE, filter_beta_inflated = TRUE, thr_zero = 0.005, thr_one = 0.005)
```

Returns a list: `results`, `selection`, plus `contrasts` if contrasts were requested and
`omnibus` if `omnibus = TRUE`. The schemas are in SKILL.md. For a numeric design matrix, a
`(Intercept)` column is removed and re-added by the model, so coefficient names are the
matrix column names.

## find_families(): choose families only

```r
find_families(counts_matrix, n_genes = 200, n_boot = 10, top_n = 4, families = NULL,
  criterion = c("GAIC","BIC","AIC"), gaic_k = NULL, min_n = 5, seed = NULL,
  group_by_support = TRUE, binom_bd = NULL, filter_beta_inflated = TRUE,
  thr_zero = 0.005, thr_one = 0.005, bootstrap = TRUE, transform_mode = NULL,
  workers = NULL, parallel = FALSE, show_progress = TRUE)
```

It fits intercept-only models, so it needs no design. It returns `top_families_overall` (chr,
length <= top_n), `top_families_by_support` (list), `freq_table_overall`, `prop_table_overall`,
`freq_by_support`, `prop_by_support`, `sampled_results` (tibble:
`bootstrap feature family skipped n_valid support`) and `transform_mode`.

## Plots (ggplot2 objects)

```r
plot_volcano(x, contrast = NULL, fdr_threshold = 0.05, lfc_threshold = 1,
  label_top = 10, point_size = 1.8, alpha = 0.7, title = NULL)
plot_ma(x, contrast = NULL, counts_matrix = NULL, fdr_threshold = 0.05,
  label_top = 10, point_size = 1.8, alpha = 0.7, title = NULL)
```

- `x` can be run_perseo output, fit_gamlss_models output, or a contrasts tibble.
- `contrast` is required if there is more than one.
- The volcano plot's x axis is `estimate` (the family's link scale; see interpretation.md) and
  its y axis is `-log10(p_adj)`.
- `plot_ma` uses `log2(rowMeans(counts)+1)` on the x axis when `counts_matrix` is given
  (pass it explicitly), otherwise the |estimate| rank.
- Labels use `ggrepel` if installed.
- Save with `ggplot2::ggsave(file, p, width = 8, height = 5.5)`.

## report_perseo(): self-contained HTML report

```r
report_perseo(x, output_file = "perseo_report.html", output_dir = ".",
  title = "PERSEO Analysis Report", open = TRUE, quiet = FALSE)
```

- `x` must be `run_perseo()` output.
- Requires `rmarkdown`, `DT`, `htmltools` and `jsonlite`.
- Returns the path invisibly.
- Use `open = FALSE, quiet = TRUE` in scripts.

## Utilities

| function | purpose |
|---|---|
| `infer_support(y)` | `"count"` (non-negative integers), `"unit"` ([0,1]), `"positive"` (>0), `"zi_positive"` (>=0, non-integer, has zeros), `"real"`, `"none"` |
| `family_groups()` | family names per support (below) |
| `has_insufficient_variation(y)` | `TRUE` if fewer than 2 distinct values, or 2 values with one appearing once (such features are skipped) |
| `transform_response(y, fam, mode = c("strict","safe"), eps = 1e-6, allow_eps = TRUE)` | response transform used for a family: `list(y, mask, logJ_per_obs, meta, mode_used)` |
| `transform_for_family_strict(y, fam, eps, allow_eps)`, `transform_for_family(y, fam, strategy = "safe", eps)` | the strict/safe back ends |
| `inverse_transform(z, meta)` | back to the original scale |
| `jacobian_sum(logJ_per_obs, mask = NULL)` | sum of log-Jacobians over the mask |
| `print(<perseo_results>)` | short text summary (safe to print) |

## Families

Default candidates (`families = NULL`):
`BE BEINF BEZI BEINF0 BEINF1 PO NBI ZIP ZINBI BI BB PIG GA GG IG LOGNO ZAGA ZAIG NO TF GU`.
Any name is looked up in `gamlss.dist` and then on the search path. Unknown names are skipped
silently, so check the spelling of custom `families` (e.g. `BEo`, not `BEO`).

`family_groups()` (used when `group_by_support = TRUE`):

| support | families |
|---|---|
| count | PO NBI ZIP ZINBI ZIP2 BI BB PIG |
| unit | BE BEINF BEZI BEo BEINF0 BEINF1 (BEo: same distribution as BE with shape-parameter mu; not a default) |
| positive | GA GG LOGNO IG NO TF GU |
| zi_positive | ZAGA ZAIG (+ the positive families when zeros are below `thr_zero`) |
| real | NO TF GU |

Typical choices, per the package FAQ:

| data | families |
|---|---|
| RNA-seq counts | `NBI PO ZIP ZINBI` |
| proteomics intensities | `GG GA LOGNO IG` |
| beta values or proportions | `BE BEZI BEINF0 BEINF1 BEINF` |
| normalised or log data | `NO TF` |
