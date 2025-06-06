#' Fit GAMLSS models across multiple genes (enhanced version)
#'
#' Applies GAMLSS models to each gene (row-wise in expression matrix) and selects the best model per gene
#' according to a selected criterion (AIC, BIC, GAIC(k=3), or log-likelihood).
#' Only model diagnostics and scores are stored, not the full model, to keep memory usage low.
#'
#' @param counts_matrix Numeric matrix: expression values (genes x samples)
#' @param X Model matrix
#' @param families Character vector of GAMLSS family names
#' @param criterion Character: selection criterion, one of "AIC", "BIC", "GAIC", "logLik"
#' @param timeout Max time (sec) per model fit
#' @param verbose Logical, whether to print progress messages
#'
#' @return A tibble with gene name, selected family, scores, and diagnostics
#' @export
fit_gamlss_models <- function(counts_matrix, X, families = c("PO", "NBI", "NO", "GA"),
                              criterion = c("AIC", "BIC", "GAIC", "logLik"),
                              timeout = 10, verbose = TRUE) {
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
  

  # Family-aware transformation
  transform_for_family <- function(y, fam) {
    if (fam %in% c("GA", "GG", "LOGNO", "IG")) {
      y[y <= 0] <- NA
      return(y)
    } else if (fam %in% c("PO", "NBI", "ZINBI", "ZIP", "ZIP2")) {
      y[y < 0] <- NA
      return(round(y))
    } else if (fam %in% c("BE", "BEINF", "BEO", "BEZI", "BEo", "BEINF0")) {
      y <- (y - min(y)) / (max(y) - min(y) + 1e-8)
      y[y <= 0] <- NA
      y[y >= 1] <- NA
      return(y)
    } else if (fam %in% c("NO", "TF", "GU")) {
      return(scale(y))
    } else {
      return(y)
    }
  }
  
  fit_one_gene <- function(gene_name, y) {
    if (all(y == 0)) return(NULL)
    
    results <- list()
    for (fam in families) {
      y_trans <- transform_for_family(y, fam)
      
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
