# PERSEO methods (for methods sections and deep debugging)

## Pipeline

1. **Family selection** (`find_families`): `n_boot` pulls of `n_genes` random features, or all
   features if `bootstrap = FALSE`. Each feature fits intercept-only models for its eligible
   families and votes for the winner. The `top_n` most frequent winners overall become the
   candidates for DE.
2. **DE** (`fit_gamlss_models`): for each feature, fit the design on each candidate family,
   restricted per feature by support when `group_by_support = TRUE`. Pick the winner and
   extract mu coefficients.
3. **Contrasts / omnibus**: refit the winner, then compute `C beta` with the robust covariance.
   Optionally gate the contrasts on an omnibus test.
4. **FDR**: coefficients within each term, contrasts within each contrast.

## Support inference (`infer_support`) and family filtering

- The support is inferred from finite values:
  - count: all non-negative integers
  - unit: all values in [0,1]
  - positive: all values > 0
  - zi_positive: values >= 0, including zeros and non-integers
  - real: anything else
- With `group_by_support = TRUE`, only `family_groups()[[support]]` are eligible.
  zi_positive admits the zero-adjusted families plus the positive families.
- `filter_beta_inflated = TRUE` (applies regardless of `group_by_support`):
  - unit: keep only the beta families whose point masses match the evidence (a share of
    exact 0s >= `thr_zero`, of exact 1s >= `thr_one`): BE/BEo with neither, BEZI/BEINF0 with
    zeros, BEINF1 with ones, BEINF with both. With material 0s/1s this also drops non-unit
    families, as for zi_positive (a family without a point mass would fit them via the nudge)
  - count: drop ZIP/ZINBI/ZIP2 if the share of 0s < `thr_zero`
  - positive: drop ZAIG/ZAGA if the share of 0s < `thr_zero`
  - zi_positive with the share of 0s >= `thr_zero`: keep only the zero-adjusted families
    (a continuous family would otherwise "fit" the zeros via the epsilon nudge)
- Features failing `has_insufficient_variation()` are skipped. `find_families` also skips
  all-zero features.

## Transforms, common mask and Jacobian

**Strict mode** (the default with support grouping) never repairs data. Invalid observations
are masked, per family:

- counts: identity; the mask is non-negative integers
- positive: identity, with exact 0 nudged to 1e-6
- zero-adjusted: identity, y >= 0
- unit: identity. Values outside [0,1] are masked; exact 0/1 are kept for families with a
  point mass there, otherwise nudged by 1e-6 (`allow_eps = TRUE`, the pipeline default) or
  masked (`allow_eps = FALSE`)
- real: z-score; log-Jacobian `-log(sd)`

**Safe mode** applies global affine maps instead:

- positive with min <= 0: shift
- unit: min-max
- real: z-score
- counts: identity

All families of a feature are compared on the **intersection of their masks** (the common
mask). The IC is then

`IC = -2 logLik + penalty * df - 2 * sum(logJ over mask)`

so every family is scored on the original data scale. The penalty is 2 for AIC, `log(n)` for
BIC, and `gaic_k` (default `log(n)`) for GAIC.

## Winner selection

Take the family with the minimum IC. Families within 2 IC units of the minimum are treated as
tied, and the tie goes to the best goodness of fit: the Filliben correlation between sorted
normalised quantile residuals and normal quantiles (`select_by_ic_gof`, `delta = 2`).

## Coefficient covariance

The primary covariance is `(X'WX)^-1` from the fit's mu design and IRLS weights
(`mu_vcov_robust`), with a Moore-Penrose pseudo-inverse if the matrix is singular. It
reproduces the GAMLSS Wald covariance but works inside parallel workers. The fallbacks, in
order, are `vcov(fit, "mu")`, the diagonal of the `summary()` SEs, and `fit$vcov.mu`.
Coefficient tests use a t distribution on the residual df; contrasts use a z test.

## Parallelism and reproducibility

- `parallel = TRUE` sets `future::plan(multisession, workers)` and restores the previous plan
  on exit.
- `future.seed = TRUE` together with `seed` makes the bootstrap reproducible.
- Model fits run with `gamlss.control(c.crit = 0.001)`, and all console output from gamlss
  is captured.
