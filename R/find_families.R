#' Find Best-Fitting GAMLSS Families Across Features (Bootstrap + Jacobian-corrected)
#'
#' Randomly samples features (genes, OTUs, proteins, etc.) in bootstrap pulls,
#' fits intercept-only GAMLSS models with candidate families, and selects the best
#' family per feature using a Jacobian-corrected information criterion (AIC/BIC/GAIC).
#' Comparisons are made on a common mask across candidate families.
#'
#' The function aggregates results across bootstrap pulls and returns the top
#' families (overall and by empirical support) with frequencies and proportions.
#'
#' @param counts_matrix numeric matrix (features x samples) with rownames as feature IDs
#' @param bootstrap logical, whether to use bootstrap sampling (default: TRUE).
#'   If TRUE, randomly samples n_genes features in n_boot pulls.
#'   If FALSE, evaluates all families on ALL features (full evaluation, no sampling).
#' @param n_genes integer, number of features per bootstrap pull (default: 200).
#'   Ignored if bootstrap = FALSE.
#' @param n_boot integer, number of bootstrap pulls (default: 10).
#'   Ignored if bootstrap = FALSE.
#' @param top_n integer, number of top families to return (default: 4)
#' @param families character vector of families to consider; if NULL, a default set is used
#' @param verbose logical, print progress (default: TRUE)
#' @param min_n integer, min valid samples per feature after common mask (default: 5)
#' @param seed integer or NULL, seed for reproducibility (default: NULL)
#' @param group_by_support logical; if TRUE (default), restricts candidate families
#'   to those matching the empirical support of each feature. If FALSE, all families
#'   in `families` are tested regardless of inferred support, allowing exploration
#'   of model fit across support boundaries at the cost of increased computation.
#' @param binom_bd NULL | numeric scalar | numeric vector len = ncol(counts_matrix).
#'   Binomial denominator for BI/BB. If NULL, infer per feature as max(y) when consistent.
#' @param criterion "GAIC", "BIC", or "AIC" (default: "GAIC")
#' @param gaic_k numeric, penalty for GAIC. If NULL, defaults to log(n_valid)
#' @param filter_beta_inflated logical, filter inflated Beta families if no 0/1 evidence (default: TRUE)
#' @param thr_zero numeric threshold for zero inflation evidence (default: 0.005)
#' @param thr_one numeric threshold for one inflation evidence (default: 0.005)
#' @param transform_mode character: "strict" or "safe". Transformation mode for family comparison.
#'   - "strict" (default when group_by_support = TRUE): Conservative, domain-preserving.
#'     Invalid observations excluded via mask. No data repair.
#'   - "safe" (default when group_by_support = FALSE): Global affine transformations
#'     to fit data into family domain. Invertible with Jacobian correction.
#'   If NULL, defaults based on group_by_support setting.
#'
#' @return list(
#'   top_families_overall   = character[],
#'   top_families_by_support= list(count=..., unit=..., positive=..., real=...),
#'   freq_table_overall     = table,
#'   prop_table_overall     = named numeric,
#'   freq_by_support        = named list of tables,
#'   prop_by_support        = named list of named numerics,
#'   sampled_results        = tibble with per-feature decisions across pulls,
#'   transform_mode         = character indicating mode used
#' )
#' @export
find_families <- function(counts_matrix,
                          bootstrap = TRUE,
                          n_genes = 200,
                          n_boot  = 10,
                          top_n   = 4,
                          families = NULL,
                          verbose = TRUE,
                          min_n = 5,
                          seed = NULL,
                          group_by_support = TRUE,
                          binom_bd = NULL,
                          criterion = c("GAIC","BIC","AIC"),
                          gaic_k = NULL,
                          filter_beta_inflated = TRUE,
                          thr_zero = 0.005,
                          thr_one  = 0.005,
                          transform_mode = NULL) {
  
  # ---- Input validation ----
  criterion <- match.arg(criterion)
  if (bootstrap) {
    validate_counts_matrix(counts_matrix, min_features = n_genes)
  } else {
    validate_counts_matrix(counts_matrix, min_features = 1)
  }
  validate_criterion_args(criterion, gaic_k)
  
  # Set default transform_mode based on group_by_support
  if (is.null(transform_mode)) {
    transform_mode <- if (group_by_support) "strict" else "safe"
  } else {
    transform_mode <- match.arg(transform_mode, c("strict", "safe"))
  }
  
  if (is.null(families)) {
    families <- default_candidate_families()
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  feature_ids_all <- rownames(counts_matrix)
  if (is.null(feature_ids_all)) {
    feature_ids_all <- seq_len(nrow(counts_matrix))
  }
  
  n_samples <- ncol(counts_matrix)
  
  # ---- Per-feature evaluation worker ----
  evaluate_feature <- function(feature_id) {
    y <- as.numeric(counts_matrix[feature_id, ])
    
    # Early exit: insufficient variation
    if (has_insufficient_variation(y)) {
      return(list(
        feature = feature_id,
        family = NA_character_,
        skipped = TRUE,
        n_valid = 0L,
        support = infer_support(y)
      ))
    }
    
    # Early exit: all zeros
    if (all(y == 0, na.rm = TRUE)) {
      return(list(
        feature = feature_id,
        family = NA_character_,
        skipped = TRUE,
        n_valid = 0L,
        support = "count"
      ))
    }
    
    # Filter candidate families
    filtered <- filter_candidate_families(
      feature_vec = y,
      candidate_families = families,
      group_by_support = group_by_support,
      filter_beta_inflated = filter_beta_inflated,
      thr_zero = thr_zero,
      thr_one = thr_one
    )
    
    if (length(filtered$families_to_test) == 0) {
      return(list(
        feature = feature_id,
        family = NA_character_,
        skipped = TRUE,
        n_valid = 0L,
        support = filtered$support
      ))
    }
    
    # Compare families and select best
    selection_result <- compare_families_on_feature(
      feature_vec = y,
      families = filtered$families_to_test,
      min_n = min_n,
      criterion = criterion,
      gaic_k = gaic_k,
      bd_vec = filtered$bd_vec,
      transform_mode = transform_mode
    )
    
    list(
      feature = feature_id,
      family = selection_result$best_family,
      skipped = is.na(selection_result$best_family),
      n_valid = selection_result$n_valid,
      support = filtered$support
    )
  }
  
  # ---- Evaluation execution ----
  if (bootstrap) {
    # Bootstrap sampling mode
    if (isTRUE(verbose)) {
      message("Running find_families with ", n_boot, " bootstrap pulls of ", n_genes, " features each...")
    }
    
    results_list <- lapply(seq_len(n_boot), function(b) {
      if (isTRUE(verbose)) message("  Pull ", b, " of ", n_boot)
      pull_ids <- sample(feature_ids_all, n_genes, replace = FALSE)
      future.apply::future_lapply(
        pull_ids, 
        evaluate_feature, 
        future.seed = TRUE,
        future.packages = "PERSEO"
      )
    })
  } else {
    # Full evaluation mode (all features, no bootstrap)
    if (isTRUE(verbose)) {
      message("Running find_families on ALL ", length(feature_ids_all), " features (no bootstrap)...")
    }
    
    results_list <- list(
      future.apply::future_lapply(
        feature_ids_all,
        evaluate_feature,
        future.seed = TRUE,
        future.packages = "PERSEO"
      )
    )
  }
  
  # ---- Aggregate results ----
  sampled_results <- bind_bootstrap_results(results_list)
  summary_info <- summarize_family_frequencies(sampled_results, top_n, verbose)
  
  # Append sampled results and metadata to summary
  summary_info$sampled_results <- sampled_results
  summary_info$transform_mode <- transform_mode
  summary_info
}