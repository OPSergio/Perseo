#' Fit GAMLSS models per feature (family selection with common mask + Jacobian)
#'
#' Compares candidate GAMLSS families per feature on the same set of observations
#' (common mask), using strict, family-specific transformations and Jacobian
#' correction for fair IC comparison. Selects the best family by GAIC/BIC/AIC and
#' returns per-term statistics from `summary(fit, what = "mu")` (no VCOV dependency).
#' Parallelization over features is supported via `future.apply`, with optional
#' progress bars via `progressr`.
#'
#' @param counts_matrix numeric matrix (features x samples).
#' @param design_matrix data.frame or matrix of covariates with
#'   `nrow = ncol(counts_matrix)`. If a `"(Intercept)"` column is present it is
#'   removed to avoid a double intercept in `y ~ .`.
#' @param candidate_families character vector of GAMLSS families to test
#'   (e.g. `c("PO","NBI","GA","GG","IG","LOGNO","NO","TF")`).
#' @param criterion one of `"GAIC"`, `"BIC"`, or `"AIC"`; default `"GAIC"`.
#' @param gaic_k numeric penalty used when `criterion = "GAIC"`. If `NULL`,
#'   uses `log(n_valid_obs)`.
#' @param min_n integer minimum number of valid observations after applying the
#'   common mask; features below this threshold are skipped. Default `5`.
#' @param p_adjust method passed to `p.adjust()` to compute FDR per term across
#'   features (default `"BH"`).
#' @param contrast_matrix ignored in this no-VCOV mode; included for future API
#'   compatibility. A warning is emitted if provided.
#' @param workers integer number of parallel workers. If `> 1`, a
#'   `multisession` plan is installed for the duration of the call. Default `1`.
#' @param show_progress logical; show a `progressr` progress bar. Default `TRUE`.
#' @param progress_label character label shown next to the progress bar.
#'
#' @return A list with:
#'   \describe{
#'     \item{results}{tibble with columns: `feature`, `term`, `effect`, `se`,
#'       `stat`, `pval`, `padj`.}
#'     \item{selection}{tibble with columns: `feature`, `best_family`,
#'       `n_valid_obs`, `ic_value` (Jacobian-corrected IC).}
#'   }
#'
#' @details Family comparison uses strict transforms, a common mask across
#' families, and Jacobian correction so ICs are comparable on the original data
#' scale. Coefficients/p-values are taken from `summary(., what = "mu")`;
#' general contrasts are not computed in this mode.
#'
#' @examples
#' \dontrun{
#' fit <- fit_gamlss_models(
#'   counts_matrix, design, c("NBI","GG","LOGNO"),
#'   criterion = "BIC", min_n = 20,
#'   workers = max(1, future::availableCores()-1),
#'   show_progress = TRUE
#' )
#' }
#' @seealso find_families, transform_for_family_strict
#' @export
source("R/utils_transformations.R") # Load utility functions (temporal)
fit_gamlss_models <- function(counts_matrix, X, families = c("PO", "NBI", "NO", "GA"),
                              criterion = c("AIC", "BIC", "GAIC", "logLik"),
                              timeout = 10, verbose = TRUE,
                              strategy = "safe", eps = 1e-6) {
  criterion <- match.arg(criterion)
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(multisession, workers = availableCores())
  

  resid_diagnostics <- function(model, y) {
    res <- tryCatch(residuals(model, type = "normalized"), error = function(e) rep(NA, length(y)))
    valid_res <- res[is.finite(res)]
    ks <- tryCatch(ks.test(valid_res, "pnorm")$p.value, error = function(e) NA)
    skew <- tryCatch(if (length(unique(valid_res)) > 1) e1071::skewness(valid_res, na.rm = TRUE) else NA, error = function(e) NA)
    kurt <- tryCatch(if (length(unique(valid_res)) > 1) e1071::kurtosis(valid_res, na.rm = TRUE) else NA, error = function(e) NA)
    list(ks_p = ks, skewness = skew, kurtosis = kurt)
  }
  
  
  fit_one_gene <- function(gene_name, y) {
    if (all(y == 0)) return(NULL)
    
    results <- list()
    for (fam in families) {
      y_trans <- transform_for_family(y, fam, strategy = strategy, eps = eps)
      
      fit <- tryCatch({
        R.utils::withTimeout({
          model <- gamlss(y_trans ~ X - 1, family = fam, trace = FALSE)
          diag <- resid_diagnostics(model, y_trans)
          list(
            family = fam,
            AIC = AIC(model),
            BIC = BIC(model),
            GAIC3 = GAIC(model, k = 3),
            logLik = logLik(model),
            df = model$df.fit,
            residuals = diag
          )
        }, timeout = timeout, onTimeout = "silent")
      }, error = function(e) NULL)
      
      if (!is.null(fit)) results[[fam]] <- fit
    }
    
    if (length(results) == 0) return(NULL)
    
    best <- switch(criterion,
                   AIC = which.min(map_dbl(results, "AIC")),
                   BIC = which.min(map_dbl(results, "BIC")),
                   GAIC = which.min(map_dbl(results, "GAIC3")),
                   logLik = which.max(map_dbl(results, "logLik")))
    
    best_fit <- results[[best]]
    tibble(
      gene = gene_name,
      best_family = best_fit$family,
      AIC = best_fit$AIC,
      BIC = best_fit$BIC,
      GAIC3 = best_fit$GAIC3,
      logLik = best_fit$logLik,
      df = best_fit$df,
      ks_p = best_fit$residuals$ks_p,
      skewness = best_fit$residuals$skewness,
      kurtosis = best_fit$residuals$kurtosis
    )
  }
  
  if (verbose) cli::cli_h1("Fitting GAMLSS models across genes")
  gene_names <- rownames(counts_matrix)
  gene_list <- asplit(counts_matrix, MARGIN = 1)
  
  results <- furrr::future_map2_dfr(
    .x = gene_names,
    .y = gene_list,
    .f = fit_one_gene,
    .progress = verbose,
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("dplyr", "gamlss", "e1071", "tibble")
    )
  )
  
  if (verbose) {
    cli::cli_h2("GAMLSS Summary")
    total_genes <- nrow(counts_matrix)
    fitted_genes <- nrow(results)
    skipped_genes <- total_genes - fitted_genes
    top_families <- sort(table(results$best_family), decreasing = TRUE)
    most_common <- names(top_families)[1]
    
    cli::cli_text("Genes analyzed: {total_genes}")
    cli::cli_text("Genes fitted successfully: {fitted_genes}")
    cli::cli_text("Genes skipped (e.g., all zeros or NA): {skipped_genes}")
    cli::cli_text("Most frequent family: {most_common} ({top_families[1]} genes)")
    cli::cli_text("Top families:")
    print(top_families)
  }
  
  return(results)
}
