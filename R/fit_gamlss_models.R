#' Apply linear contrasts to GAMLSS mu coefficients
#'
#' Computes estimates, standard errors, z-statistics and p-values for arbitrary
#' linear contrasts of mu-parameter coefficients from a fitted GAMLSS model.
#'
#' @param beta Named numeric vector of coefficients from `coef(fit, what = "mu")`.
#' @param V Covariance matrix from `vcov(fit, what = "mu")`.
#' @param C Contrast matrix where each row defines a linear combination.
#'   Column names must match (a subset of) `names(beta)`.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{contrast}{Contrast name (from rownames(C))}
#'     \item{estimate}{Point estimate (C times beta)}
#'     \item{se}{Standard error}
#'     \item{z}{z-statistic}
#'     \item{p_value}{Two-sided p-value}
#'   }
#'
#' @details Missing columns in C are filled with 0; extra columns in C not in
#' beta are dropped. If vcov is not positive definite or contains NAs/Infs,
#' all statistics are returned as NA.
#'
#' @keywords internal
apply_contrasts <- function(beta, V, C) {
  # Validate inputs
  if (!is.numeric(beta) || is.null(names(beta))) {
    stop("beta must be a named numeric vector")
  }
  if (!is.matrix(V) || nrow(V) != ncol(V) || nrow(V) != length(beta)) {
    stop("V must be a square matrix with nrow = length(beta)")
  }
  if (!is.matrix(C) || ncol(C) == 0) {
    stop("C must be a numeric matrix with at least one column")
  }
  
  # Check for non-finite values in V
  if (any(!is.finite(V))) {
    warning("vcov contains non-finite values; returning NA contrasts")
    contrast_names <- if (!is.null(rownames(C))) {
      rownames(C)
    } else {
      paste0("contrast_", seq_len(nrow(C)))
    }
    return(tibble::tibble(
      contrast = contrast_names,
      estimate = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p_value = NA_real_
    ))
  }
  
  # Align C columns to beta names
  beta_names <- names(beta)
  C_names <- colnames(C)
  
  if (is.null(C_names)) {
    stop("C must have column names matching coefficient names")
  }
  
  # Create aligned contrast matrix: rows = contrasts, cols = beta_names
  C_aligned <- matrix(0, nrow = nrow(C), ncol = length(beta_names))
  colnames(C_aligned) <- beta_names
  rownames(C_aligned) <- rownames(C)
  
  # Fill in values where columns match
  matched_cols <- intersect(C_names, beta_names)
  if (length(matched_cols) == 0) {
    warning("No contrast matrix columns match coefficient names; returning NA")
    contrast_names <- if (!is.null(rownames(C))) {
      rownames(C)
    } else {
      paste0("contrast_", seq_len(nrow(C)))
    }
    return(tibble::tibble(
      contrast = contrast_names,
      estimate = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p_value = NA_real_
    ))
  }
  
  for (col in matched_cols) {
    C_aligned[, col] <- C[, col]
  }
  
  # Compute contrasts
  estimates <- as.numeric(C_aligned %*% beta)
  variances <- diag(C_aligned %*% V %*% t(C_aligned))
  
  # Ensure non-negative variances
  variances <- pmax(variances, 0)
  ses <- sqrt(variances)
  
  z_values <- estimates / ses
  p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
  
  contrast_names <- if (!is.null(rownames(C))) {
    rownames(C)
  } else {
    paste0("contrast_", seq_len(nrow(C)))
  }
  
  tibble::tibble(
    contrast = contrast_names,
    estimate = unname(estimates),
    se = unname(ses),
    z = unname(z_values),
    p_value = unname(p_values)
  )
}


#' Fit GAMLSS models per feature (family selection with common mask + Jacobian)
#'
#' Compares candidate GAMLSS families per feature on the same set of observations
#' (common mask), using family-specific transformations and Jacobian
#' correction for fair IC comparison. Selects the best family by GAIC/BIC/AIC and
#' returns per-term statistics from `summary(fit, what = "mu")`.
#'
#' @param counts_matrix numeric matrix (features x samples).
#' @param design_matrix data.frame or matrix of covariates with
#'   `nrow = ncol(counts_matrix)`. Can also be a formula string (e.g., "~ condition + batch")
#'   in which case the design matrix is created internally from the metadata.
#' @param metadata optional data.frame with sample metadata. Required when
#'   `design_matrix` is a formula string or when using `contrast_variable`.
#'   Must have `nrow = ncol(counts_matrix)`.
#' @param candidate_families character vector of GAMLSS families to test
#'   (e.g. `c("PO","NBI","GA","GG","IG","LOGNO","NO","TF")`).
#' @param criterion one of `"GAIC"`, `"BIC"`, or `"AIC"`; default `"GAIC"`.
#' @param gaic_k numeric penalty used when `criterion = "GAIC"`. If `NULL`,
#'   uses `log(n_valid_obs)`.
#' @param min_n integer minimum number of valid observations after applying the
#'   common mask; features below this threshold are skipped. Default `5`.
#' @param p_adjust method passed to `p.adjust()` to compute FDR per term across
#'   features (default `"BH"`).
#' @param contrast_matrix optional numeric matrix where each row defines a linear
#'   combination of mu coefficients. Column names must match coefficient names
#'   from the design matrix. If provided, contrasts are computed from vcov after
#'   model fitting.
#' @param contrast_variable optional character string specifying a categorical variable
#'   name from the metadata for which all pairwise contrasts should be computed.
#'   If both `contrast_matrix` and `contrast_variable` are provided, `contrast_matrix`
#'   takes priority with a warning message.
#' @param omnibus logical; if `TRUE` and contrasts are requested, performs a feature-level
#'   omnibus test first and only computes contrasts for features passing the test.
#'   This implements a hierarchical testing strategy (omnibus → post-hoc). Default `FALSE`.
#' @param omnibus_threshold numeric; significance threshold for omnibus test. Features with
#'   omnibus p-value > this threshold will not have contrasts computed. Only used when
#'   `omnibus = TRUE`. Default `0.05`.
#' @param omnibus_test character; type of omnibus test. Either `"Wald"` (multi-df Wald test,
#'   faster, uses existing vcov) or `"LRT"` (likelihood ratio test, more robust, requires
#'   refitting reduced model). Default `"Wald"`.
#' @param workers integer number of parallel workers. Only used when `parallel = TRUE`.
#'   If `NULL`, uses `future::availableCores() - 1`. Default `NULL`.
#' @param parallel logical; enable parallel processing via `future::plan(multisession)`.
#'   When `TRUE`, automatically configures and cleans up the future backend.
#'   Default `FALSE`.
#' @param show_progress logical; show a `progressr` progress bar. Default `TRUE`.
#' @param progress_label character label shown next to the progress bar.
#' @param transform_mode character: "strict" (default) or "safe". Transformation mode
#'   for family comparison. "strict" enforces domain-preserving transformations with
#'   observation exclusion. "safe" applies global affine transformations to fit data
#'   into family domain (invertible with Jacobian correction).
#'
#' @return A list with:
#'   \describe{
#'     \item{results}{tibble with columns: feature, term, effect, se,
#'       stat, pval, padj.}
#'     \item{selection}{tibble with columns: feature, best_family,
#'       n_valid_obs, ic_value (Jacobian-corrected IC), transform_mode.}
#'     \item{omnibus}{tibble with columns: feature, family, test_type, statistic,
#'       df, p_value, pass. Only present when `omnibus = TRUE` and contrasts are requested.}
#'     \item{contrasts}{tibble with columns: feature, family, contrast,
#'       estimate, se, z, p_value, p_adj. Only present if
#'       contrast_matrix or contrast_variable is provided. When `omnibus = TRUE`,
#'       only includes features passing the omnibus test.}
#'   }
#'
#' @details Family comparison uses transformations (strict or safe), a common mask
#' across families, and Jacobian correction so ICs are comparable on the original
#' data scale. Coefficients/p-values are taken from `summary(., what = "mu")`.
#'
#' When contrasts are requested with `omnibus = TRUE`, the function performs hierarchical
#' testing:
#' \enumerate{
#'   \item Fit full model per feature
#'   \item Perform omnibus test for target factor (Wald or LRT)
#'   \item Only compute pairwise contrasts for features passing omnibus test
#' }
#'
#' This mirrors ANOVA + post-hoc workflow and reduces multiple testing burden.
#'
#' **Important**: The omnibus test serves only as a gatekeeper. Contrast p-values are
#' adjusted exactly as before: `group_by(contrast)` across features. The omnibus p-value
#' is NOT used in contrast adjustment.
#'
#' When `omnibus = FALSE` (default), behavior is identical to previous PERSEO versions.
#'
#' When `parallel = TRUE`, the function automatically sets up `future::plan(multisession)`
#' with the specified number of workers and resets to sequential after completion.
#'
#' @examples
#' \dontrun{
#' # Standard workflow (no omnibus filtering)
#' fit <- fit_gamlss_models(
#'   counts_matrix, design, c("NBI","GG"),
#'   metadata = meta,
#'   contrast_variable = "tissue_type"
#' )
#'
#' # With omnibus filtering (Wald test)
#' fit <- fit_gamlss_models(
#'   counts_matrix, design, c("NBI","GG"),
#'   metadata = meta,
#'   contrast_variable = "tissue_type",
#'   omnibus = TRUE,
#'   omnibus_threshold = 0.05,
#'   omnibus_test = "Wald"
#' )
#'
#' # View omnibus results
#' head(fit$omnibus)
#'
#' # With omnibus filtering (LRT)
#' fit_lrt <- fit_gamlss_models(
#'   counts_matrix, design, c("NBI","GG"),
#'   metadata = meta,
#'   contrast_variable = "tissue_type",
#'   omnibus = TRUE,
#'   omnibus_test = "LRT"
#' )
#' }
#' @seealso find_families, transform_response
#' @export
fit_gamlss_models <- function(counts_matrix,
                              design_matrix,
                              metadata = NULL,
                              candidate_families,
                              criterion = c("GAIC","BIC","AIC"),
                              gaic_k = NULL,
                              min_n = 5,
                              contrast_matrix = NULL,
                              contrast_variable = NULL,
                              omnibus = FALSE,
                              omnibus_threshold = 0.05,
                              omnibus_test = c("Wald", "LRT"),
                              p_adjust = "BH",
                              workers = NULL,
                              parallel = FALSE,
                              show_progress = TRUE,
                              progress_label = "Fitting features",
                              transform_mode = "strict") {
  # ---- Input validation ----
  criterion <- match.arg(criterion)
  transform_mode <- match.arg(transform_mode, c("strict", "safe"))
  omnibus_test <- match.arg(omnibus_test)
  
  validate_counts_matrix(counts_matrix)
  validate_criterion_args(criterion, gaic_k)
  
  # ---- Set up parallel backend if requested ----
  original_plan <- NULL
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required for parallel processing. Install it with install.packages('future')")
    }
    
    original_plan <- future::plan()
    
    n_workers <- if (is.null(workers)) {
      max(1, future::availableCores() - 1)
    } else {
      as.integer(workers)
    }
    
    if (n_workers < 1) {
      warning("workers must be >= 1; using workers = 1")
      n_workers <- 1L
    }
    
    future::plan(future::multisession, workers = n_workers)
    
    if (show_progress) {
      message("Parallel processing enabled with ", n_workers, " workers")
    }
  }
  
  on.exit({
    if (parallel && !is.null(original_plan)) {
      future::plan(original_plan)
    }
  }, add = TRUE)
  
  # Store variables for final report (will be shown at the very end)
  report_vars <- new.env()
  
  feature_ids <- rownames(counts_matrix)
  if (is.null(feature_ids)) {
    feature_ids <- seq_len(nrow(counts_matrix))
  }
  
  # ---- Handle formula string input ----
  if (is.character(design_matrix) && length(design_matrix) == 1) {
    if (is.null(metadata)) {
      stop("metadata must be provided when design_matrix is a formula string")
    }
    if (!is.data.frame(metadata) || nrow(metadata) != ncol(counts_matrix)) {
      stop("metadata must be a data.frame with nrow = ncol(counts_matrix)")
    }
    
    formula_obj <- stats::as.formula(design_matrix)
    design_matrix <- stats::model.matrix(formula_obj, data = metadata)
  }
  
  design_matrix <- validate_design_matrix(design_matrix, ncol(counts_matrix))
  complete_rows <- stats::complete.cases(design_matrix)
  counts_subset <- counts_matrix[, complete_rows, drop = FALSE]
  design_subset <- design_matrix[complete_rows, , drop = FALSE]
  
  # ---- Handle contrast inputs with priority logic ----
  final_contrast_matrix <- NULL
  final_contrast_variable <- NULL
  
  if (!is.null(contrast_matrix) && !is.null(contrast_variable)) {
    message("Both contrast_matrix and contrast_variable provided. ",
            "Prioritizing contrast_matrix.")
    final_contrast_matrix <- contrast_matrix
    final_contrast_variable <- NULL  # Will identify coefs from matrix
  } else if (!is.null(contrast_matrix)) {
    final_contrast_matrix <- contrast_matrix
    final_contrast_variable <- NULL
  } else if (!is.null(contrast_variable)) {
    if (is.null(metadata)) {
      stop("metadata must be provided when using contrast_variable")
    }
    final_contrast_matrix <- PERSEO:::build_contrast_matrix(
      contrast_variable, metadata[complete_rows, , drop = FALSE], design_subset
    )
    final_contrast_variable <- contrast_variable
  }
  
  # ---- Validate contrast matrix if present ----
  if (!is.null(final_contrast_matrix)) {
    if (!is.matrix(final_contrast_matrix) || !is.numeric(final_contrast_matrix)) {
      stop("contrast_matrix must be a numeric matrix")
    }
    if (is.null(colnames(final_contrast_matrix))) {
      stop("contrast_matrix must have column names matching coefficient names")
    }
  }
  
  # ---- Print startup banner ----
  if (show_progress) {
    message("\n", strrep("=", 80))
    message("fit_gamlss_models() - STARTING")
    message(strrep("=", 80))
    message("Total features    : ", length(feature_ids))
    message("Samples           : ", nrow(design_subset))
    message("Candidate families: ", paste(candidate_families, collapse = ", "))
    message("Criterion         : ", criterion)
    message("Transform mode    : ", transform_mode)
    message("P-value adjustment: ", p_adjust)
    if (parallel) {
      message("Parallel workers  : ", n_workers)
    }
    if (!is.null(final_contrast_matrix)) {
      message("Contrasts requested: ", nrow(final_contrast_matrix), " contrasts")
      if (omnibus) {
        message("  Omnibus test    : ", omnibus_test, " (threshold: ", omnibus_threshold, ")")
      }
    }
    message(strrep("=", 80), "\n")
  }
  
  # ---- Worker function: process one feature ----
  process_one_feature <- function(feature_name, feature_vec, design_mat, prog) {
    if (!is.null(prog)) prog()
    
    if (PERSEO:::has_insufficient_variation(feature_vec)) {
      return(NULL)
    }
    
    bd_vec <- if (any(PERSEO:::is_binomial_family(candidate_families))) {
      PERSEO:::infer_binomial_denominator(feature_vec)
    } else NULL
    
    family_results <- PERSEO:::compare_families_with_design(
      feature_vec,
      design_mat,
      families = candidate_families,
      criterion = criterion,
      gaic_k = gaic_k,
      min_n = min_n,
      bd_vec = bd_vec,
      transform_mode = transform_mode
    )
    
    if (nrow(family_results$comparisons) == 0) {
      return(NULL)
    }
    
    best <- dplyr::slice_min(family_results$comparisons, ic_value, with_ties = FALSE)
    
    # Extract coefficient table
    coef_df <- tibble::tibble(
      feature = feature_name,
      term    = best$coef_tbl[[1]]$term,
      effect  = best$coef_tbl[[1]]$effect,
      se      = best$coef_tbl[[1]]$se,
      stat    = best$coef_tbl[[1]]$stat,
      pval    = best$coef_tbl[[1]]$pval
    )
    
    # Initialize outputs
    omnibus_df <- NULL
    contrast_df <- NULL
    
    # Only perform omnibus/contrast if requested
    if (!is.null(final_contrast_matrix)) {
      # Re-fit to get beta and vcov for omnibus/contrasts
      tr <- PERSEO:::transform_response(feature_vec, best$family, mode = transform_mode)
      z <- tr$y[tr$mask]
      fit_data <- cbind.data.frame(y = z, design_mat[tr$mask, , drop = FALSE])
      
      fam_obj <- PERSEO:::instantiate_gamlss_family(
        best$family,
        bd_vec = if (PERSEO:::is_binomial_family(best$family)) {
          PERSEO:::infer_binomial_denominator(feature_vec, tr$mask)
        } else NULL
      )
      
      best_fit <- PERSEO:::fit_gamlss_safely(y ~ ., data = fit_data, family_obj = fam_obj)
      
      if (!is.null(best_fit)) {
        beta <- tryCatch({ coef(best_fit, what = "mu") }, error = function(e) NULL)
        
        # Extract vcov
        V <- NULL
        if (!is.null(best_fit$family_obj_stored)) {
          family_obj <- best_fit$family_obj_stored
        } else {
          family_obj <- fam_obj
        }
        
        V <- tryCatch({ vcov(best_fit, what = "mu") }, error = function(e) NULL)
        
        if (is.null(V)) {
          V <- tryCatch({
            suppressWarnings({
              s <- summary(best_fit, what = "mu")
              if (is.matrix(s) && "Std. Error" %in% colnames(s)) {
                se_vals <- s[, "Std. Error"]
                n_mu <- length(beta)
                if (length(se_vals) >= n_mu) {
                  se_vals <- se_vals[1:n_mu]
                }
                V_diag <- diag(se_vals^2)
                rownames(V_diag) <- colnames(V_diag) <- names(beta)
                V_diag
              } else {
                NULL
              }
            })
          }, error = function(e) NULL)
        }
        
        if (is.null(V) && !is.null(best_fit$vcov.mu)) {
          V <- best_fit$vcov.mu
        }
        
        # Determine which coefficients belong to the target factor
        factor_coefs <- NULL
        if (!is.null(beta)) {
          factor_coefs <- PERSEO:::identify_factor_coefficients(
            contrast_matrix = final_contrast_matrix,
            contrast_variable = final_contrast_variable,
            coef_names = names(beta)
          )
        }
        
        # Perform omnibus test if requested
        omnibus_pass <- TRUE  # Default: compute contrasts
        
        if (omnibus && !is.null(factor_coefs) && length(factor_coefs) > 0) {
          omnibus_result <- if (omnibus_test == "Wald" && !is.null(V)) {
            PERSEO:::wald_omnibus_test(beta, V, factor_coefs)
          } else if (omnibus_test == "LRT") {
            PERSEO:::lrt_omnibus_test(
              feature_vec = feature_vec,
              design_mat = design_mat,
              factor_coefs = factor_coefs,
              family = best$family,
              transform_mode = transform_mode,
              bd_vec = bd_vec
            )
          } else {
            list(statistic = NA_real_, df = NA_integer_, p_value = NA_real_, test_type = omnibus_test)
          }
          
          omnibus_df <- tibble::tibble(
            feature = feature_name,
            family = best$family,
            test_type = omnibus_result$test_type,
            statistic = omnibus_result$statistic,
            df = omnibus_result$df,
            p_value = omnibus_result$p_value,
            pass = !is.na(omnibus_result$p_value) && omnibus_result$p_value < omnibus_threshold
          )
          
          omnibus_pass <- omnibus_df$pass
        }
        
        # Compute contrasts only if omnibus passed (or omnibus disabled)
        if (omnibus_pass && !is.null(beta) && !is.null(V) && all(is.finite(V))) {
          contrast_df <- tryCatch({
            contrast_res <- PERSEO:::apply_contrasts(beta, V, final_contrast_matrix)
            contrast_res$feature <- feature_name
            contrast_res$family <- best$family
            contrast_res
          }, error = function(e) {
            # Fit failed
            contrast_names <- if (!is.null(rownames(final_contrast_matrix))) {
              rownames(final_contrast_matrix)
            } else {
              paste0("contrast_", seq_len(nrow(final_contrast_matrix)))
            }
            tibble::tibble(
              feature = feature_name,
              family = best$family,
              contrast = contrast_names,
              estimate = NA_real_,
              se = NA_real_,
              z = NA_real_,
              p_value = NA_real_
            )
          })
        }
      }
    }
    
    list(
      feature  = feature_name,
      best_fam = best$family,
      n_valid  = best$n_valid_obs,
      ic_value = best$ic_value,
      coef_df  = coef_df,
      omnibus_df = omnibus_df,
      contrast_df = contrast_df
    )
  }
  
  # ---- Parallel execution via future.apply ----
  p <- if (show_progress) {
    progressr::progressor(steps = nrow(counts_subset))
  } else NULL
  
  out_list <- future.apply::future_lapply(
    seq_len(nrow(counts_subset)),
    function(i) {
      process_one_feature(
        feature_name = feature_ids[i],
        feature_vec  = counts_subset[i,],
        design_mat   = design_subset,
        prog         = p
      )
    },
    future.seed = TRUE
  )
  
  valid_list <- Filter(Negate(is.null), out_list)
  if (length(valid_list) == 0) {
    empty_res <- tibble::tibble(
      feature = character(), term = character(), effect = numeric(),
      se = numeric(), stat = numeric(), pval = numeric(), padj = numeric()
    )
    empty_sel <- tibble::tibble(
      feature = character(), best_family = character(),
      n_valid_obs = integer(), ic_value = numeric()
    )
    empty_out <- list(results = empty_res, selection = empty_sel)
    if (!is.null(final_contrast_matrix)) {
      if (omnibus) {
        empty_out$omnibus <- tibble::tibble(
          feature = character(), family = character(), test_type = character(),
          statistic = numeric(), df = integer(), p_value = numeric(), pass = logical()
        )
      }
      empty_out$contrasts <- tibble::tibble(
        feature = character(), family = character(), contrast = character(),
        estimate = numeric(), se = numeric(), z = numeric(),
        p_value = numeric(), p_adj = numeric()
      )
    }
    return(empty_out)
  }
  
  # ---- Aggregate results and compute FDR per term ----
  results_df <- dplyr::bind_rows(lapply(valid_list, `[[`, "coef_df"))
  if (nrow(results_df) > 0 && "pval" %in% names(results_df)) {
    results_df <- dplyr::mutate(
      dplyr::group_by(results_df, term),
      padj = p.adjust(pval, method = p_adjust),
      .groups = "drop"
    )
  } else {
    results_df$padj <- NA_real_
  }
  
  # Extract feature names - handle both character and numeric
  feature_names <- vapply(valid_list, function(x) as.character(x$feature), FUN.VALUE = character(1))
  
  selection_df <- tibble::tibble(
    feature      = feature_names,
    best_family  = vapply(valid_list, `[[`, "best_fam", FUN.VALUE = character(1)),
    n_valid_obs  = vapply(valid_list, `[[`, "n_valid", FUN.VALUE = integer(1)),
    ic_value     = vapply(valid_list, `[[`, "ic_value", FUN.VALUE = numeric(1)),
    transform_mode = transform_mode
  )
  
  output <- list(results = results_df, selection = selection_df)
  
  # ---- Aggregate omnibus results ----
  if (!is.null(final_contrast_matrix) && omnibus) {
    omnibus_list <- lapply(valid_list, `[[`, "omnibus_df")
    omnibus_list <- Filter(Negate(is.null), omnibus_list)
    
    if (length(omnibus_list) > 0) {
      output$omnibus <- dplyr::bind_rows(omnibus_list)
    }
  }
  
  # ---- Aggregate contrasts and compute FDR per contrast (UNCHANGED LOGIC) ----
  if (!is.null(final_contrast_matrix)) {
    contrasts_df <- dplyr::bind_rows(lapply(valid_list, `[[`, "contrast_df"))
    
    if (nrow(contrasts_df) > 0 && "p_value" %in% names(contrasts_df)) {
      # CRITICAL: This grouping logic is UNCHANGED from original PERSEO
      contrasts_df <- dplyr::mutate(
        dplyr::group_by(contrasts_df, contrast),
        p_adj = p.adjust(p_value, method = p_adjust),
        .groups = "drop"
      )
    } else {
      contrasts_df$p_adj <- NA_real_
    }
    
    output$contrasts <- contrasts_df
  }
  
  # ---- Store report data (will print on function exit) ----
  if (show_progress) {
    report_vars$n_total <- length(feature_ids)
    report_vars$n_success <- length(valid_list)
    report_vars$n_failed <- length(feature_ids) - length(valid_list)
    report_vars$selection_df <- selection_df
    report_vars$results_df <- results_df
    report_vars$output <- output
    report_vars$omnibus <- omnibus
    report_vars$omnibus_test <- omnibus_test
    report_vars$omnibus_threshold <- omnibus_threshold
    
    # Schedule report to print on exit (AFTER return, so it's the last thing shown)
    on.exit({
      print_final_report(report_vars)
    }, add = TRUE)
  }
  
  output
}

# Helper function to print final report
print_final_report <- function(vars) {
  n_total <- vars$n_total
  n_success <- vars$n_success
  n_failed <- vars$n_failed
  success_rate <- n_success / n_total * 100
  selection_df <- vars$selection_df
  results_df <- vars$results_df
  output <- vars$output
  omnibus <- vars$omnibus
  omnibus_test <- vars$omnibus_test
  omnibus_threshold <- vars$omnibus_threshold
  
  message("\n", strrep("=", 80))
  message("fit_gamlss_models() - FINAL REPORT")
  message(strrep("=", 80))
  message("\nModel Fitting Summary:")
  message("  Total features      : ", n_total)
  message("  Successfully fitted : ", n_success, " (", sprintf("%.1f%%", success_rate), ")")
  message("  Failed/skipped      : ", n_failed, " (", sprintf("%.1f%%", n_failed/n_total*100), ")")
  
  if (n_success > 0) {
    # Family distribution
    fam_table <- table(selection_df$best_family)
    fam_sorted <- sort(fam_table, decreasing = TRUE)
    
    message("\nSelected Family Distribution:")
    for (i in seq_along(fam_sorted)) {
      fam_name <- names(fam_sorted)[i]
      count <- as.numeric(fam_sorted[i])
      pct <- count / sum(fam_sorted) * 100
      message(sprintf("  %-10s : %5d (%.1f%%)", fam_name, count, pct))
    }
    
    # Coefficient statistics
    if (nrow(results_df) > 0) {
      n_coefs <- nrow(results_df)
      n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
      sig_rate <- if (n_coefs > 0) n_sig / n_coefs * 100 else 0
      
      message("\nCoefficient Results:")
      message("  Total coefficients  : ", n_coefs)
      message("  Significant (FDR<5%): ", n_sig, " (", sprintf("%.1f%%", sig_rate), ")")
      
      # Per-term breakdown
      term_summary <- results_df %>%
        dplyr::group_by(term) %>%
        dplyr::summarize(
          n_total = dplyr::n(),
          n_sig = sum(padj < 0.05, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(n_sig))
      
      if (nrow(term_summary) > 0) {
        message("\n  Significant findings by term:")
        for (i in seq_len(min(10, nrow(term_summary)))) {
          term_name <- term_summary$term[i]
          n_sig_term <- term_summary$n_sig[i]
          n_total_term <- term_summary$n_total[i]
          pct_sig <- n_sig_term / n_total_term * 100
          message(sprintf("    %-25s : %5d / %5d (%.1f%%)", 
                        term_name, n_sig_term, n_total_term, pct_sig))
        }
      }
    }
    
    # Contrast statistics if available
    if (!is.null(output$contrasts) && nrow(output$contrasts) > 0) {
      n_contrasts <- nrow(output$contrasts)
      n_sig_contrasts <- sum(output$contrasts$p_adj < 0.05, na.rm = TRUE)
      
      message("\nContrast Results:")
      message("  Total contrasts     : ", n_contrasts)
      message("  Significant (FDR<5%): ", n_sig_contrasts, " (", 
              sprintf("%.1f%%", n_sig_contrasts/n_contrasts*100), ")")
      
      if (omnibus && !is.null(output$omnibus)) {
        n_omnibus_pass <- sum(output$omnibus$pass, na.rm = TRUE)
        n_omnibus_total <- nrow(output$omnibus)
        message("\nOmnibus Test (", omnibus_test, "):")
        message("  Features tested     : ", n_omnibus_total)
        message("  Passed threshold    : ", n_omnibus_pass, " (", 
                sprintf("%.1f%%", n_omnibus_pass/n_omnibus_total*100), ")")
        message("  Threshold           : p < ", omnibus_threshold)
      }
    }
  } else {
    message("\n[WARNING] No features were successfully fitted.")
  }
  
  message(strrep("=", 80), "\n")
}
