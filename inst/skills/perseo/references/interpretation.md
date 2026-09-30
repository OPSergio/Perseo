# Interpreting PERSEO results

## What `effect` / `estimate` mean

Coefficients are on the default GAMLSS link of the mu parameter for the feature's
`best_family`, fitted to the response after the strict transform. The scale is therefore
**per feature**. Check `de$selection$best_family` (or `de$contrasts$family`) before comparing
magnitudes.

| best_family | response the model sees (strict mode) | link | a coefficient of `b` means |
|---|---|---|---|
| PO NBI PIG ZIP ZINBI ZIP2 | raw counts | log | mean multiplied by `exp(b)` (ln fold change) |
| BI BB | raw counts, denominator = `binom_bd` or `max(y)` | logit | log odds ratio |
| GA GG IG | raw values (exact 0 set to 1e-6) | log | mean multiplied by `exp(b)` |
| LOGNO | raw values | identity on log(y) | geometric mean multiplied by `exp(b)` |
| ZAGA ZAIG | raw values incl. zeros | log | mean of the positive part multiplied by `exp(b)` |
| BE BEZI BEINF0 BEINF1 BEINF | raw proportions (unsupported exact 0/1 nudged by 1e-6) | logit | log odds ratio of the mean of the beta component |
| NO TF GU | per-feature z-score | identity | shift in that feature's SDs |

Notes:

- The log is natural: `log2FC = b / log(2)`.
- `plot_volcano(lfc_threshold = 1)` uses these raw units. For log-link families that is an
  e-fold (about 2.7x) change, not 2x.
- `safe` mode applies a global affine shift or scale first (see methods.md), so effects are
  on that transformed scale.
- To report fold changes, restrict to log-link families or report the family alongside
  each effect.

## Tests and multiple testing

- **Coefficients** (`results`): `stat = effect/se` uses the robust mu covariance
  `(X'WX)^-1`, with a t distribution on the residual df (normal if unavailable). `padj`
  applies `p.adjust` **within each `term`** across features.
- **Contrasts** (`contrasts`): `estimate = C beta` and `se = sqrt(C V C')`, with a two-sided
  z test. `p_adj` applies **within each `contrast`** across features.
- **Omnibus** (`omnibus`): a joint Wald (`beta' V^-1 beta`, chi-sq df = number of
  coefficients) or LRT test over the coefficients that carry non-zero weight in the contrast
  matrix. For `contrast_variable`, those are all of the factor's levels.
  `pass = p_value < omnibus_threshold` uses the **unadjusted** p-value. It only decides which
  features get contrasts; the contrast FDR is computed as usual over the passing features.
- The `(Intercept)` term is always "significant". Drop it before counting hits.

## Missing or NA output: causes

| symptom | cause |
|---|---|
| feature absent from `selection`/`results` | insufficient variation, no eligible family for its support, fewer than `min_n` valid observations, or all fits failed |
| NA `se`/`z`/`p_value` in contrasts | covariance non-finite or fit failed on refit |
| contrasts only for some features | `omnibus = TRUE` and the feature did not pass |
| whole `contrasts` NA for a contrast | contrast matrix colnames do not match coefficient names |
| `differential_expression` is NULL | no family selected in bootstrap (`summary$status`) |
| fewer samples than expected | NA in a formula variable (complete cases only) |

## Useful summaries (small outputs)

```r
de <- res$differential_expression
nrow(de$selection); nrow(counts)                       # fitted vs. total features
sort(table(de$selection$best_family), decreasing = TRUE)
r <- subset(de$results, term != "(Intercept)")
tapply(r$padj < 0.05, r$term, sum, na.rm = TRUE)       # hits per term
c <- de$contrasts
tapply(c$p_adj < 0.05, c$contrast, sum, na.rm = TRUE)  # hits per contrast
head(c[order(c$p_adj), c("feature","family","contrast","estimate","p_adj")], 10)
if (!is.null(de$omnibus)) mean(de$omnibus$pass)        # omnibus pass rate
```

## Reporting checklist

State the following:

- the family set used (`res$summary$families_selected`) and the per-feature winner distribution
- the criterion, transform mode, and FDR method and grouping (per term / per contrast)
- the number of skipped features
- the effect scale caveat
- whether an omnibus gate was used and at which raw threshold
