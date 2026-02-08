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
#' @param show_progress logical, show progress messages and reports (default: TRUE)
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
#' @param bootstrap logical, whether to use bootstrap sampling (default: TRUE).
#'   If TRUE, randomly samples n_genes features in n_boot pulls.
#'   If FALSE, evaluates all families on ALL features (full evaluation, no sampling).
#' @param transform_mode character: "strict" or "safe". Transformation mode for family comparison.
#'   - "strict" (default when group_by_support = TRUE): Conservative, domain-preserving.
#'     Invalid observations excluded via mask. No data repair.
#'   - "safe" (default when group_by_support = FALSE): Global affine transformations
#'     to fit data into family domain. Invertible with Jacobian correction.
#'   If NULL, defaults based on group_by_support setting.
#' @param seed integer random seed for reproducibility (default NULL).
#' @param workers integer number of parallel workers. Only used when `parallel = TRUE`.
#'   If `NULL`, uses `future::availableCores() - 1`. Default `NULL`.
#' @param parallel logical; enable parallel processing via `future::plan(multisession)`.
#'   When `TRUE`, automatically configures and cleans up the future backend.
#'   Default `FALSE`.
#' @param show_progress logical; show progress messages and reports (default TRUE).
#'
#' @return A list with:
#'   \describe{
#'     \item{top_families_overall}{character vector (length = top_n) of most
#'       frequently selected families across all bootstrap samples.}
#'     \item{top_families_by_support}{named list with top families per support
#'       class (count/unit/positive/real); metadata only if group_by_support=FALSE.}
#'     \item{freq_table_overall}{named integer vector with selection counts.}
#'     \item{prop_table_overall}{named numeric vector with selection proportions.}
#'     \item{freq_by_support}{list of frequency tables by support.}
#'     \item{prop_by_support}{list of proportion tables by support.}
#'     \item{sampled_results}{tibble with one row per feature×bootstrap pull,
#'       including columns: bootstrap, feature, family, skipped (logical),
#'       n_valid (integer), support (character).}
#'     \item{transform_mode}{character indicating which transformation mode was used.}
#'   }
#'
#' @details
#' Performs bootstrap family selection:
#' \enumerate{
#'   \item Sample n_genes features from counts_matrix.
#'   \item For each feature, fit intercept-only GAMLSS models across candidate families.
#'   \item Use strict transformations and a common mask to compare models on the same
#'         observations.
#'   \item Apply Jacobian correction so ICs are comparable on the original scale.
#'   \item Record the best family per feature.
#'   \item Repeat n_boot times.
#'   \item Return top_n most frequently selected families overall and by support.
#' }
#'
#' When `parallel = TRUE`, the function automatically sets up `future::plan(multisession)`
#' with the specified number of workers and resets to sequential after completion.
#'
#' @examples
#' \dontrun{
#' # Sequential execution (default)
#' ff <- find_families(
#'   counts_matrix = counts,
#'   n_genes = 200,
#'   n_boot = 10,
#'   top_n = 4,
#'   criterion = "BIC"
#' )
#'
#' # Parallel execution (automatic setup)
#' ff <- find_families(
#'   counts_matrix = counts,
#'   n_genes = 200,
#'   n_boot = 10,
#'   top_n = 4,
#'   criterion = "BIC",
#'   parallel = TRUE,
#'   workers = 8
#' )
#' }
#' @export
find_families <- function(counts_matrix,
                         n_genes = 200,
                         n_boot = 10,
                         top_n = 4,
                         families = NULL,
                         criterion = c("GAIC", "BIC", "AIC"),
                         gaic_k = NULL,
                         min_n = 5,
                         seed = NULL,
                         group_by_support = FALSE,
                         binom_bd = NULL,
                         filter_beta_inflated = TRUE,
                         thr_zero = 0.005,
                         thr_one = 0.005,
                         bootstrap = TRUE,
                         transform_mode = "strict",
                         workers = NULL,
                         parallel = FALSE,
                         show_progress = TRUE) {
  criterion <- match.arg(criterion)
  
  # ---- Set up parallel backend if requested ----
  original_plan <- NULL
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required for parallel processing. Install it with install.packages('future')")
    }
    
    # Store original plan
    original_plan <- future::plan()
    
    # Determine number of workers
    n_workers <- if (is.null(workers)) {
      max(1, future::availableCores() - 1)
    } else {
      as.integer(workers)
    }
    
    if (n_workers < 1) {
      warning("workers must be >= 1; using workers = 1")
      n_workers <- 1L
    }
    
    # Set up multisession plan
    future::plan(future::multisession, workers = n_workers)
    
    if (show_progress) {
      message("Parallel processing enabled with ", n_workers, " workers")
    }
  }
  
  # Ensure we reset the plan on exit
  on.exit({
    if (parallel && !is.null(original_plan)) {
      future::plan(original_plan)
    }
  }, add = TRUE)
  
  # ---- Input validation ----
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
  
  # ---- Capture internal functions for parallel workers ----
  .extract_mu_coefficients <- extract_mu_coefficients
  
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
      transform_mode = transform_mode,
      extract_coef_fn = .extract_mu_coefficients
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
  has_cli <- requireNamespace("cli", quietly = TRUE)
  
  if (bootstrap) {
    # Bootstrap sampling mode
    if (isTRUE(show_progress)) {
      if (has_cli) {
        cli::cli_rule("find_families() - STARTING")
        cli::cli_alert_info("Mode              : {.val Bootstrap sampling}")
        cli::cli_alert_info("Bootstrap pulls   : {.val {n_boot}}")
        cli::cli_alert_info("Features/pull     : {.val {n_genes}}")
        cli::cli_alert_info("Total evaluations : {.val {n_boot * n_genes}}")
        cli::cli_alert_info("Candidate families: {.val {paste(families, collapse = ', ')}}")
        cli::cli_alert_info("Transform mode    : {.val {transform_mode}}")
        if (parallel) {
          cli::cli_alert_info("Parallel workers  : {.val {n_workers}}")
        }
        cli::cli_rule()
      } else {
        message("\n", strrep("=", 80))
        message("find_families() - STARTING")
        message(strrep("=", 80))
        message("Mode              : Bootstrap sampling")
        message("Bootstrap pulls   : ", n_boot)
        message("Features/pull     : ", n_genes)
        message("Total evaluations : ", n_boot * n_genes)
        message("Candidate families: ", paste(families, collapse = ", "))
        message("Transform mode    : ", transform_mode)
        if (parallel) {
          message("Parallel workers  : ", n_workers)
        }
        message(strrep("=", 80), "\n")
      }
    }
    
    results_list <- lapply(seq_len(n_boot), function(b) {
      if (isTRUE(show_progress)) {
        if (has_cli) {
          cli::cli_progress_step(
            "Bootstrap pull {b}/{n_boot}",
            msg_done = "Bootstrap pull {b}/{n_boot} completed"
          )
        } else {
          message("[", format(Sys.time(), "%H:%M:%S"), "] Bootstrap pull ", b, "/", n_boot, " - sampling ", n_genes, " features...")
        }
      }
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
    if (isTRUE(show_progress)) {
      if (has_cli) {
        cli::cli_rule("find_families() - STARTING")
        cli::cli_alert_info("Mode              : {.val Full evaluation (no bootstrap)}")
        cli::cli_alert_info("Total features    : {.val {length(feature_ids_all)}}")
        cli::cli_alert_info("Candidate families: {.val {paste(families, collapse = ', ')}}")
        cli::cli_alert_info("Transform mode    : {.val {transform_mode}}")
        if (parallel) {
          cli::cli_alert_info("Parallel workers  : {.val {n_workers}}")
        }
        cli::cli_rule()
      } else {
        message("\n", strrep("=", 80))
        message("find_families() - STARTING")
        message(strrep("=", 80))
        message("Mode              : Full evaluation (no bootstrap)")
        message("Total features    : ", length(feature_ids_all))
        message("Candidate families: ", paste(families, collapse = ", "))
        message("Transform mode    : ", transform_mode)
        if (parallel) {
          message("Parallel workers  : ", n_workers)
        }
        message(strrep("=", 80), "\n")
      }
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
  summary_info <- summarize_family_frequencies(sampled_results, top_n, show_progress)
  
  # Append sampled results and metadata to summary
  summary_info$sampled_results <- sampled_results
  summary_info$transform_mode <- transform_mode
  summary_info
}