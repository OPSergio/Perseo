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
    return(tibble::tibble(
      contrast = rownames(C) %||% paste0("contrast_", seq_len(nrow(C))),
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
    return(tibble::tibble(
      contrast = rownames(C) %||% paste0("contrast_", seq_len(nrow(C))),
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
  
  tibble::tibble(
    contrast = rownames(C) %||% paste0("contrast_", seq_len(nrow(C))),
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
#'   `nrow = ncol(counts_matrix)`.
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
#'   model fitting. Use limma-style contrast matrices.
#' @param workers integer number of parallel workers. If `> 1`, caller should
#'   configure `future::plan()` before calling. Default `1`.
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
#'     \item{contrasts}{tibble with columns: feature, family, contrast,
#'       estimate, se, z, p_value, p_adj. Only present if
#'       contrast_matrix is not NULL.}
#'   }
#'
#' @details Family comparison uses transformations (strict or safe), a common mask
#' across families, and Jacobian correction so ICs are comparable on the original
#' data scale. Coefficients/p-values are taken from `summary(., what = "mu")`.
#'
#' When `contrast_matrix` is provided, contrasts are computed using the variance-
#' covariance matrix from the best-fitting model. If vcov extraction fails, contrasts
#' will be NA for that feature.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' fit <- fit_gamlss_models(
#'   counts_matrix, design, c("NBI","GG","LOGNO"),
#'   criterion = "BIC", min_n = 20,
#'   workers = 4,
#'   show_progress = TRUE,
#'   transform_mode = "strict"
#' )
#'
#' # With custom contrast matrix (limma-style)
#' # For 3-level factor with columns: (Intercept), StatusB, StatusC
#' C <- matrix(c(0, 1, -1), nrow = 1)  # B vs C
#' colnames(C) <- c("(Intercept)", "StatusB", "StatusC")
#' rownames(C) <- "B_vs_C"
#' 
#' fit <- fit_gamlss_models(
#'   counts_matrix, design, c("NBI","GG"),
#'   contrast_matrix = C
#' )
#' head(fit$contrasts)
#' }
#' @seealso find_families, transform_response
#' @export
fit_gamlss_models <- function(counts_matrix,
                              design_matrix,
                              candidate_families,
                              criterion = c("GAIC","BIC","AIC"),
                              gaic_k = NULL,
                              min_n = 5,
                              contrast_matrix = NULL,
                              p_adjust = "BH",
                              workers = 1,
                              show_progress = TRUE,
                              progress_label = "Fitting features",
                              transform_mode = "strict") {
  # ---- Input validation ----
  criterion <- match.arg(criterion)
  transform_mode <- match.arg(transform_mode, c("strict", "safe"))
  
  validate_counts_matrix(counts_matrix)
  validate_criterion_args(criterion, gaic_k)
  
  feature_ids <- rownames(counts_matrix)
  if (is.null(feature_ids)) {
    feature_ids <- seq_len(nrow(counts_matrix))
  }
  
  design_matrix <- validate_design_matrix(design_matrix, ncol(counts_matrix))
  complete_rows <- stats::complete.cases(design_matrix)
  counts_subset <- counts_matrix[, complete_rows, drop = FALSE]
  design_subset <- design_matrix[complete_rows, , drop = FALSE]
  
  # ---- Validate contrast matrix if provided ----
  if (!is.null(contrast_matrix)) {
    if (!is.matrix(contrast_matrix) || !is.numeric(contrast_matrix)) {
      stop("contrast_matrix must be a numeric matrix")
    }
    if (is.null(colnames(contrast_matrix))) {
      stop("contrast_matrix must have column names matching coefficient names")
    }
  }
  
  # ---- Worker function: process one feature (family selection + coef table + contrasts) ----
  process_one_feature <- function(feature_name, feature_vec, design_mat, prog) {
    if (!is.null(prog)) prog()
    
    if (PERSEO:::has_insufficient_variation(feature_vec)) {
      return(NULL)
    }
    
    # Infer binomial denominator if BI/BB families are candidates
    bd_vec <- if (any(PERSEO:::is_binomial_family(candidate_families))) {
      PERSEO:::infer_binomial_denominator(feature_vec)
    } else NULL
    
    # group_by_support hard-coded to FALSE: test all families
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
    
    # Compute contrasts if requested (re-fit best model to get vcov)
    contrast_df <- NULL
    if (!is.null(contrast_matrix)) {
      contrast_df <- tryCatch({
        # Re-fit best model to extract vcov
        tr <- PERSEO:::transform_response(feature_vec, best$family, mode = transform_mode)
        z <- tr$y[tr$mask]
        
        # Prepare data
        fit_data <- cbind.data.frame(y = z, design_mat[tr$mask, , drop = FALSE])
        
        # Instantiate family
        fam_obj <- PERSEO:::instantiate_gamlss_family(
          best$family,
          bd_vec = if (PERSEO:::is_binomial_family(best$family)) {
            PERSEO:::infer_binomial_denominator(feature_vec, tr$mask)
          } else NULL
        )
        
        # Re-fit with same formula as compare_families_with_design (y ~ .)
        best_fit <- PERSEO:::fit_gamlss_safely(y ~ ., data = fit_data, family_obj = fam_obj)
        
        if (!is.null(best_fit)) {
          beta <- tryCatch({
            coef(best_fit, what = "mu")
          }, error = function(e) NULL)
          
          # Try to extract vcov - use fallback if vcov() fails
          V <- NULL
          
          # Approach 1: Try vcov() with family_obj in scope
          if (!is.null(best_fit$family_obj_stored)) {
            family_obj <- best_fit$family_obj_stored
          } else {
            family_obj <- fam_obj
          }
          
          V <- tryCatch({
            vcov(best_fit, what = "mu")
          }, error = function(e) NULL)
          
          # Approach 2: Extract diagonal vcov from summary if vcov() fails
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
          
          # Approach 3: Try accessing vcov.mu directly
          if (is.null(V) && !is.null(best_fit$vcov.mu)) {
            V <- best_fit$vcov.mu
          }
          
          if (!is.null(beta) && !is.null(V) && all(is.finite(V))) {
            contrast_res <- PERSEO:::apply_contrasts(beta, V, contrast_matrix)
            contrast_res$feature <- feature_name
            contrast_res$family <- best$family
            contrast_res  # Return this
          } else {
            # vcov failed or contains non-finite values
            tibble::tibble(
              feature = feature_name,
              family = best$family,
              contrast = rownames(contrast_matrix) %||% paste0("contrast_", seq_len(nrow(contrast_matrix))),
              estimate = NA_real_,
              se = NA_real_,
              z = NA_real_,
              p_value = NA_real_
            )
          }
        } else {
          # Fit failed
          tibble::tibble(
            feature = feature_name,
            family = best$family,
            contrast = rownames(contrast_matrix) %||% paste0("contrast_", seq_len(nrow(contrast_matrix))),
            estimate = NA_real_,
            se = NA_real_,
            z = NA_real_,
            p_value = NA_real_
          )
        }
      }, error = function(e) {
        # Return NA-filled tibble on error
        tibble::tibble(
          feature = feature_name,
          family = best$family,
          contrast = rownames(contrast_matrix) %||% paste0("contrast_", seq_len(nrow(contrast_matrix))),
          estimate = NA_real_,
          se = NA_real_,
          z = NA_real_,
          p_value = NA_real_
        )
      })
    }
    
    list(
      feature  = feature_name,
      best_fam = best$family,
      n_valid  = best$n_valid_obs,
      ic_value = best$ic_value,
      coef_df  = coef_df,
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
    future.seed = TRUE,
    future.packages = "PERSEO"
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
    if (!is.null(contrast_matrix)) {
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
  
  selection_df <- tibble::tibble(
    feature      = vapply(valid_list, `[[`, "feature", FUN.VALUE = character(1)),
    best_family  = vapply(valid_list, `[[`, "best_fam", FUN.VALUE = character(1)),
    n_valid_obs  = vapply(valid_list, `[[`, "n_valid", FUN.VALUE = integer(1)),
    ic_value     = vapply(valid_list, `[[`, "ic_value", FUN.VALUE = numeric(1)),
    transform_mode = transform_mode
  )
  
  # ---- Aggregate contrasts and compute FDR per contrast ----
  output <- list(results = results_df, selection = selection_df)
  
  if (!is.null(contrast_matrix)) {
    contrasts_df <- dplyr::bind_rows(lapply(valid_list, `[[`, "contrast_df"))
    
    if (nrow(contrasts_df) > 0 && "p_value" %in% names(contrasts_df)) {
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
  
  output
}
