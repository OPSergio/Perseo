#' Run Complete PERSEO Workflow: Family Selection + Differential Expression
#'
#' High-level orchestration function that executes the full PERSEO pipeline:
#' (1) identifies best-fitting GAMLSS families via bootstrap sampling,
#' (2) fits differential expression models with the selected families, and
#' (3) adjusts p-values for multiple testing.
#'
#' This is the recommended entry point for end-to-end analyses.
#'
#' @param counts_matrix Numeric matrix (features × samples) with rownames as feature IDs.
#' @param design_matrix Data frame or matrix with covariates (samples × predictors).
#'   Will be cleaned automatically (intercept column removed, dimension validated).
#' @param contrast_matrix Optional numeric matrix where each row defines a linear
#'   combination of mu coefficients. Column names must match coefficient names from
#'   the design matrix. If provided, contrasts are computed and adjusted for multiple
#'   testing. Use limma::makeContrasts() to build contrast matrices easily.
#' @param bootstrap Logical, whether to use bootstrap sampling for family selection
#'   (default: TRUE). If TRUE, randomly samples n_genes features in n_boot pulls.
#'   If FALSE, evaluates all families on ALL features (full evaluation).
#' @param n_genes Integer, number of features to sample per bootstrap pull for family
#'   selection (default: 200). Ignored if bootstrap = FALSE.
#' @param n_boot Integer, number of bootstrap pulls for family selection (default: 10).
#'   Ignored if bootstrap = FALSE.
#' @param top_n Integer, number of top families to use in differential expression
#'   (default: 4).
#' @param families Character vector of candidate families; if NULL, uses default set.
#' @param group_by_support Logical; if TRUE, restricts families to those matching
#'   empirical support (default: TRUE).
#' @param criterion Character, information criterion: "GAIC", "BIC", or "AIC" (default: "GAIC").
#' @param gaic_k Numeric, GAIC penalty; if NULL, defaults to log(n) (default: NULL).
#' @param min_n Integer, minimum valid observations after common mask (default: 5).
#' @param binom_bd Binomial denominator for BI/BB families; NULL (infer), scalar, or
#'   vector of length = ncol(counts_matrix) (default: NULL).
#' @param filter_beta_inflated Logical, filter inflated Beta families without evidence
#'   (default: TRUE).
#' @param p_adjust_method Character, p-value adjustment method passed to `p.adjust()`
#'   (default: "BH"). Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param transform_mode Character: "strict" or "safe", or NULL for automatic selection.
#'   If NULL (default), automatically sets to "strict" when group_by_support = TRUE,
#'   and "safe" when group_by_support = FALSE. Controls transformation behavior for
#'   both family selection and differential expression.
#' @param verbose Logical, print progress messages (default: TRUE).
#' @param seed Integer or NULL, random seed for reproducibility (default: NULL).
#'
#' @return List with three components:
#' \describe{
#'   \item{family_selection}{Output from `find_families()`: top families overall and
#'     by support, frequency tables, sampled results tibble.}
#'   \item{differential_expression}{Output from `fit_gamlss_models()` with FDR-adjusted
#'     p-values: results tibble with per-term statistics, selection tibble with best
#'     family per feature, and contrasts tibble (if contrast_matrix provided).}
#'   \item{summary}{List with execution summary: number of features tested, families selected,
#'     models fitted, adjustment method used, and transform_mode.}
#' }
#'
#' @details
#' **Workflow:**
#' 1. **Family Selection**: Calls `find_families()` to identify the most plausible
#'    GAMLSS families via bootstrap sampling with Jacobian-corrected IC comparison.
#' 2. **Differential Expression**: Calls `fit_gamlss_models()` with the top families
#'    identified in step 1, fitting models with design matrix covariates.
#' 3. **Multiple Testing Correction**: Adjusts p-values globally across all features
#'    and terms using the specified method (default: Benjamini-Hochberg FDR).
#'
#' **Verbosity Control:**
#' - `verbose = TRUE`: Prints progress for both family selection and model fitting.
#' - `verbose = FALSE`: Suppresses all informational messages (errors still shown).
#'
#' **P-value Adjustment:**
#' Applied globally across all features and model terms. The `p_adjust_method`
#' parameter is passed directly to `stats::p.adjust()`. Common choices:
#' - `"BH"` (default): Benjamini-Hochberg FDR control (recommended for omics).
#' - `"bonferroni"`: Conservative family-wise error rate control.
#' - `"none"`: No adjustment (not recommended for high-dimensional data).
#'
#' @examples
#' \dontrun{
#' # Basic example: Differential expression
#' metadata$tissue_type <- relevel(factor(metadata$tissue_type), ref = "Normal")
#' design <- model.matrix(~ tissue_type + age + gender, data = metadata)
#' counts <- as.matrix(expression_data[, rownames(design)])
#'
#' results <- run_perseo(
#'   counts_matrix = counts,
#'   design_matrix = design,
#'   n_genes = 200,
#'   n_boot = 10,
#'   top_n = 4,
#'   criterion = "BIC",
#'   p_adjust_method = "BH",
#'   verbose = TRUE,
#'   seed = 12345
#' )
#'
#' # Extract significant features (FDR < 0.05)
#' sig_features <- results$differential_expression$results %>%
#'   filter(term != "(Intercept)", p_adj < 0.05)
#'
#' # Advanced: Custom contrasts for multi-level factors
#' # Build contrast matrix manually
#' coef_names <- colnames(design)
#' contrast_matrix <- matrix(0, nrow = 3, ncol = length(coef_names))
#' colnames(contrast_matrix) <- coef_names
#' rownames(contrast_matrix) <- c("Tumor_vs_Normal", "Metastasis_vs_Normal", "Metastasis_vs_Tumor")
#' 
#' contrast_matrix["Tumor_vs_Normal", "tissue_typeTumor"] <- 1
#' contrast_matrix["Metastasis_vs_Normal", "tissue_typeMetastasis"] <- 1
#' contrast_matrix["Metastasis_vs_Tumor", "tissue_typeMetastasis"] <- 1
#' contrast_matrix["Metastasis_vs_Tumor", "tissue_typeTumor"] <- -1
#'
#' results_contrasts <- run_perseo(
#'   counts_matrix = counts,
#'   design_matrix = design,
#'   contrast_matrix = contrast_matrix,
#'   criterion = "BIC",
#'   p_adjust_method = "BH",
#'   seed = 12345
#' )
#'
#' # View contrast results
#' sig_contrasts <- results_contrasts$differential_expression$contrasts %>%
#'   filter(p_adj < 0.05)
#' }
#'
#' @seealso
#' \code{\link{find_families}} for family selection details.
#' \code{\link{fit_gamlss_models}} for differential expression fitting.
#'
#' @export
run_perseo <- function(counts_matrix,
                       design_matrix,
                       contrast_matrix = NULL,
                       bootstrap = TRUE,
                       n_genes = 200,
                       n_boot = 10,
                       top_n = 4,
                       families = NULL,
                       group_by_support = TRUE,
                       criterion = c("GAIC", "BIC", "AIC"),
                       gaic_k = NULL,
                       min_n = 5,
                       binom_bd = NULL,
                       filter_beta_inflated = TRUE,
                       p_adjust_method = "BH",
                       transform_mode = NULL,
                       verbose = TRUE,
                       seed = NULL) {
  
  # ---- Argument validation ----
  criterion <- match.arg(criterion)
  
  if (!p_adjust_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
    stop(
      "p_adjust_method must be one of: 'holm', 'hochberg', 'hommel', ",
      "'bonferroni', 'BH', 'BY', 'fdr', 'none'"
    )
  }
  
  # ---- Step 1: Family Selection ----
  if (verbose) {
    message("\n", strrep("=", 70))
    message(" PERSEO Pipeline - Step 1: Family Selection")
    message(strrep("=", 70))
  }
  
  family_results <- find_families(
    counts_matrix = counts_matrix,
    bootstrap = bootstrap,
    n_genes = n_genes,
    n_boot = n_boot,
    top_n = top_n,
    families = families,
    verbose = verbose,
    min_n = min_n,
    seed = seed,
    group_by_support = group_by_support,
    binom_bd = binom_bd,
    criterion = criterion,
    gaic_k = gaic_k,
    filter_beta_inflated = filter_beta_inflated,
    transform_mode = transform_mode
  )
  
  # Extract selected families
  selected_families <- family_results$top_families_overall
  
  if (length(selected_families) == 0) {
    warning(
      "No families selected in bootstrap analysis. ",
      "Check data quality and consider adjusting parameters."
    )
    return(list(
      family_selection = family_results,
      differential_expression = NULL,
      summary = list(
        n_features_total = nrow(counts_matrix),
        n_samples = ncol(counts_matrix),
        families_selected = character(0),
        models_fitted = 0,
        p_adjust_method = p_adjust_method,
        status = "failed_family_selection"
      )
    ))
  }
  
  if (verbose) {
    message("\nSelected families (top ", top_n, "): ", paste(selected_families, collapse = ", "))
  }
  
  # ---- Step 2: Differential Expression ----
  if (verbose) {
    message("\n", strrep("=", 70))
    message(" PERSEO Pipeline - Step 2: Differential Expression")
    message(strrep("=", 70))
  }
  
  de_results <- fit_gamlss_models(
    counts_matrix = counts_matrix,
    design_matrix = design_matrix,
    candidate_families = selected_families,
    criterion = criterion,
    gaic_k = gaic_k,
    min_n = min_n,
    contrast_matrix = contrast_matrix,
    p_adjust = "none",  # We'll do global adjustment below
    transform_mode = family_results$transform_mode
  )
  
  # ---- Step 3: Multiple Testing Correction ----
  if (verbose) {
    message("\n", strrep("=", 70))
    message(" PERSEO Pipeline - Step 3: Multiple Testing Correction")
    message(strrep("=", 70))
    message("Method: ", p_adjust_method)
  }
  
  # Apply p-value adjustment globally across all features and terms
  if (!is.null(de_results$results) && nrow(de_results$results) > 0) {
    # Extract p-values
    p_values <- de_results$results$pval
    
    # Adjust
    p_adjusted <- stats::p.adjust(p_values, method = p_adjust_method)
    
    # Add to results
    de_results$results$p_adj <- p_adjusted
    
    if (verbose) {
      n_sig_raw <- sum(p_values < 0.05, na.rm = TRUE)
      n_sig_adj <- sum(p_adjusted < 0.05, na.rm = TRUE)
      message(sprintf(
        "Significant tests (p < 0.05): %d raw, %d adjusted",
        n_sig_raw, n_sig_adj
      ))
    }
  } else {
    if (verbose) {
      message("No models fitted successfully.")
    }
  }
  
  # Apply p-value adjustment to contrasts if present
  if (!is.null(de_results$contrasts) && nrow(de_results$contrasts) > 0) {
    p_values_contrasts <- de_results$contrasts$p_value
    p_adjusted_contrasts <- stats::p.adjust(p_values_contrasts, method = p_adjust_method)
    de_results$contrasts$p_adj <- p_adjusted_contrasts
    
    if (verbose) {
      n_sig_raw_c <- sum(p_values_contrasts < 0.05, na.rm = TRUE)
      n_sig_adj_c <- sum(p_adjusted_contrasts < 0.05, na.rm = TRUE)
      message(sprintf(
        "Significant contrasts (p < 0.05): %d raw, %d adjusted",
        n_sig_raw_c, n_sig_adj_c
      ))
    }
  }
  
  # ---- Summary ----
  n_models_fitted <- if (!is.null(de_results$results)) {
    length(unique(de_results$results$feature))
  } else {
    0
  }
  
  summary_info <- list(
    n_features_total = nrow(counts_matrix),
    n_samples = ncol(counts_matrix),
    n_features_sampled = n_genes,
    n_bootstrap_pulls = n_boot,
    families_selected = selected_families,
    n_families_selected = length(selected_families),
    criterion = criterion,
    transform_mode = family_results$transform_mode,
    models_fitted = n_models_fitted,
    p_adjust_method = p_adjust_method,
    status = "completed"
  )
  
  if (verbose) {
    message("\n", strrep("=", 70))
    message(" PERSEO Pipeline - Completed Successfully")
    message(strrep("=", 70))
    message("Features analyzed: ", summary_info$n_features_total)
    message("Models fitted: ", summary_info$models_fitted)
    message("Families used: ", paste(selected_families, collapse = ", "))
    message(strrep("=", 70), "\n")
  }
  
  # ---- Return structured results ----
  structure(
    list(
      family_selection = family_results,
      differential_expression = de_results,
      summary = summary_info
    ),
    class = c("perseo_results", "list")
  )
}


#' Print method for PERSEO results
#'
#' @param x A `perseo_results` object returned by `run_perseo()`.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.perseo_results <- function(x, ...) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(" PERSEO Analysis Results\n")
  cat(strrep("=", 70), "\n\n")
  
  cat("Status:", x$summary$status, "\n")
  cat("Features analyzed:", x$summary$n_features_total, "\n")
  cat("Samples:", x$summary$n_samples, "\n")
  cat("Models fitted:", x$summary$models_fitted, "\n")
  cat("Families selected:", paste(x$summary$families_selected, collapse = ", "), "\n")
  cat("Criterion:", x$summary$criterion, "\n")
  cat("P-value adjustment:", x$summary$p_adjust_method, "\n")
  
  if (!is.null(x$differential_expression$results)) {
    cat("\nDifferential expression results:\n")
    cat("  Total tests:", nrow(x$differential_expression$results), "\n")
    
    if ("p_adj" %in% names(x$differential_expression$results)) {
      n_sig <- sum(x$differential_expression$results$p_adj < 0.05, na.rm = TRUE)
      cat("  Significant (FDR < 0.05):", n_sig, "\n")
    }
  }
  
  cat("\n")
  cat("Access components with:\n")
  cat("  $family_selection     - Family selection results\n")
  cat("  $differential_expression - DE models with adjusted p-values\n")
  cat("  $summary              - Execution summary\n")
  cat(strrep("=", 70), "\n\n")
  
  invisible(x)
}
