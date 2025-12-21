#' Fit GAMLSS models per feature (family selection with common mask + Jacobian)
#'
#' Compares candidate GAMLSS families per feature on the same set of observations
#' (common mask), using strict, family-specific transformations and Jacobian
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
#' @param contrast_matrix ignored; included for future API compatibility.
#' @param workers integer number of parallel workers. If `> 1`, caller should
#'   configure `future::plan()` before calling. Default `1`.
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
#' scale. Coefficients/p-values are taken from `summary(., what = "mu")`.
#'
#' @examples
#' \dontrun{
#' fit <- fit_gamlss_models(
#'   counts_matrix, design, c("NBI","GG","LOGNO"),
#'   criterion = "BIC", min_n = 20,
#'   workers = 4,
#'   show_progress = TRUE
#' )
#' }
#' @seealso find_families, transform_for_family_strict
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
                              progress_label = "Fitting features") {
  # ---- Input validation ----
  criterion <- match.arg(criterion)
  if (!is.null(contrast_matrix)) {
    warning("contrast_matrix is ignored; returning per-coefficient tests only.")
  }
  
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
  
  # ---- Worker function: process one feature (family selection + coef table) ----
  process_one_feature <- function(feature_name, feature_vec, prog) {
    if (!is.null(prog)) prog()
    
    if (has_insufficient_variation(feature_vec)) {
      return(NULL)
    }
    
    # Infer binomial denominator if BI/BB families are candidates
    bd_vec <- if (any(is_binomial_family(candidate_families))) {
      infer_binomial_denominator(feature_vec)
    } else NULL
    
    # group_by_support hard-coded to FALSE: test all families
    family_results <- compare_families_with_design(
      feature_vec,
      design_subset,
      families = candidate_families,
      criterion = criterion,
      gaic_k = gaic_k,
      min_n = min_n,
      bd_vec = bd_vec
    )
    
    if (nrow(family_results$comparisons) == 0) {
      return(NULL)
    }
    
    best <- dplyr::slice_min(family_results$comparisons, ic_value, with_ties = FALSE)
    
    list(
      feature  = feature_name,
      best_fam = best$family,
      n_valid  = best$n_valid_obs,
      ic_value = best$ic_value,
      coef_df  = tibble::tibble(
        feature = feature_name,
        term    = best$coef_tbl[[1]]$term,
        effect  = best$coef_tbl[[1]]$effect,
        se      = best$coef_tbl[[1]]$se,
        stat    = best$coef_tbl[[1]]$stat,
        pval    = best$coef_tbl[[1]]$pval
      )
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
    return(list(results = empty_res, selection = empty_sel))
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
    ic_value     = vapply(valid_list, `[[`, "ic_value", FUN.VALUE = numeric(1))
  )
  
  list(results = results_df, selection = selection_df)
}
