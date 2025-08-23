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
  criterion <- match.arg(criterion)
  if (!is.null(contrast_matrix)) {
    warning("contrast_matrix is ignored in no-VCOV mode; returning per-coefficient tests only.")
  }
  
  # --- helpers ---
  penalty_value <- function(n_valid_obs) switch(
    criterion,
    "AIC"  = 2,
    "BIC"  = log(n_valid_obs),
    "GAIC" = if (is.null(gaic_k)) log(n_valid_obs) else gaic_k
  )
  get_family_object <- function(family_name) {
    if (exists(family_name, where = asNamespace("gamlss.dist"), inherits = FALSE)) {
      get(family_name, envir = asNamespace("gamlss.dist"))()
    } else if (exists(family_name, mode = "function")) {
      get(family_name)()
    } else NULL
  }
  sanitize_design <- function(dm) {
    dm <- as.data.frame(dm)
    if ("(Intercept)" %in% colnames(dm)) {
      dm <- dm[, setdiff(colnames(dm), "(Intercept)"), drop = FALSE]
    }
    dm
  }
  # Robust extractor of per-coefficient table for mu (Estimate, SE, t, p)
  extract_mu_coef_table <- function(fit) {
    s <- suppressWarnings({
      s_obj <- NULL
      utils::capture.output({ s_obj <- summary(fit, what = "mu") }, file = NULL)
      s_obj
    })
    if (is.data.frame(s) || is.matrix(s)) {
      tab <- as.data.frame(s)
      names(tab) <- sub("Std\\.Error", "Std. Error", names(tab))
      if (!"Pr(>|t|)" %in% names(tab)) {
        idx <- which(tolower(names(tab)) %in% c("p-value","p.value","pr(>|t|)"))
        if (length(idx) == 1) names(tab)[idx] <- "Pr(>|t|)"
      }
      if (!"Estimate"   %in% names(tab)) tab[["Estimate"]]   <- NA_real_
      if (!"Std. Error" %in% names(tab)) tab[["Std. Error"]] <- NA_real_
      if (!"t value"    %in% names(tab)) tab[["t value"]]    <- tab[["Estimate"]]/tab[["Std. Error"]]
      if (!"Pr(>|t|)"   %in% names(tab)) tab[["Pr(>|t|)"]]   <- NA_real_
      term_names <- rownames(tab)
      if (is.null(term_names)) {
        if ("term" %in% names(tab)) term_names <- as.character(tab$term) else term_names <- names(coef(fit, what = "mu"))
      }
      return(tibble::tibble(
        term   = term_names,
        effect = as.numeric(tab[["Estimate"]]),
        se     = as.numeric(tab[["Std. Error"]]),
        stat   = as.numeric(tab[["t value"]]),
        pval   = as.numeric(tab[["Pr(>|t|)"]])
      ))
    }
    if (is.list(s)) {
      candidates <- list(s$coef.table, s$coefficients, s$coef, s$coeftable)
      for (tab in candidates) {
        if (!is.null(tab)) {
          tab <- as.data.frame(tab)
          names(tab) <- sub("Std\\.Error", "Std. Error", names(tab))
          if (!"Pr(>|t|)" %in% names(tab)) {
            idx <- which(tolower(names(tab)) %in% c("p-value","p.value","pr(>|t|)"))
            if (length(idx) == 1) names(tab)[idx] <- "Pr(>|t|)"
          }
          if (!"Estimate"   %in% names(tab)) tab[["Estimate"]]   <- NA_real_
          if (!"Std. Error" %in% names(tab)) tab[["Std. Error"]] <- NA_real_
          if (!"t value"    %in% names(tab)) tab[["t value"]]    <- tab[["Estimate"]]/tab[["Std. Error"]]
          if (!"Pr(>|t|)"   %in% names(tab)) tab[["Pr(>|t|)"]]   <- NA_real_
          term_names <- rownames(tab)
          if (is.null(term_names)) {
            if ("term" %in% names(tab)) term_names <- as.character(tab$term) else term_names <- names(coef(fit, what = "mu"))
          }
          return(tibble::tibble(
            term   = term_names,
            effect = as.numeric(tab[["Estimate"]]),
            se     = as.numeric(tab[["Std. Error"]]),
            stat   = as.numeric(tab[["t value"]]),
            pval   = as.numeric(tab[["Pr(>|t|)"]])
          ))
        }
      }
    }
    est <- coef(fit, what = "mu")
    tibble::tibble(
      term   = names(est),
      effect = as.numeric(est),
      se     = NA_real_,
      stat   = NA_real_,
      pval   = NA_real_
    )
  }
  
  # --- checks & prep (shared across workers) ---
  feature_ids <- rownames(counts_matrix)
  stopifnot(is.matrix(counts_matrix))
  stopifnot(ncol(counts_matrix) == nrow(design_matrix))
  
  design_matrix <- sanitize_design(design_matrix)
  complete_rows <- stats::complete.cases(design_matrix)
  if (!all(complete_rows)) {
    message(sum(!complete_rows), " rows with NA in design; they will be dropped via the common mask.")
  }
  
  # --- per-feature worker function ---
  process_one_feature <- function(feature_id, progressor = NULL) {
    y_raw <- as.numeric(counts_matrix[feature_id, ])
    
    # 1) strict transforms per candidate
    transformed_by_family <- lapply(candidate_families, function(fam) {
      transform_for_family_strict(y_raw, fam)
    })
    names(transformed_by_family) <- candidate_families
    
    # Common mask: domain-valid AND design complete-cases
    family_masks <- lapply(transformed_by_family, `[[`, "mask")
    common_mask  <- Reduce(`&`, family_masks) & complete_rows
    n_valid_obs  <- sum(common_mask, na.rm = TRUE)
    
    if (n_valid_obs < min_n) {
      if (!is.null(progressor)) progressor(message = feature_id)
      return(list(
        selection = tibble::tibble(
          feature     = feature_id,
          best_family = NA_character_,
          n_valid_obs = n_valid_obs,
          ic_value    = NA_real_
        ),
        results = tibble::tibble(
          feature  = feature_id,
          term     = NA_character_,
          effect   = NA_real_,
          se       = NA_real_,
          stat     = NA_real_,
          pval     = NA_real_
        )
      ))
    }
    
    penalty <- penalty_value(n_valid_obs)
    ic_by_family   <- rep(Inf, length(candidate_families)); names(ic_by_family) <- candidate_families
    fits_by_family <- vector("list", length(candidate_families)); names(fits_by_family) <- candidate_families
    
    for (family_name in candidate_families) {
      tr <- transformed_by_family[[family_name]]
      y_trans  <- tr$y[common_mask]
      logJ_sum <- sum(tr$logJ_per_obs[common_mask])
      
      X_masked <- design_matrix[common_mask, , drop = FALSE]
      data_fit <- cbind.data.frame(y = y_trans, X_masked)
      
      fam_obj <- get_family_object(family_name)
      if (is.null(fam_obj)) next
      
      fit_try <- try({
        fit_tmp <- NULL
        utils::capture.output({
          fit_tmp <- gamlss::gamlss(
            y ~ .,
            data   = data_fit,
            family = fam_obj,
            trace  = FALSE
          )
        }, file = NULL)
        fit_tmp
      }, silent = TRUE)
      if (inherits(fit_try, "try-error") || is.null(fit_try)) next
      fit <- fit_try
      fits_by_family[[family_name]] <- fit
      
      log_lik <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
      if (!is.finite(log_lik)) next
      df_fit  <- fit$df.fit
      base_ic <- (-2 * log_lik) + penalty * df_fit
      ic_by_family[family_name] <- base_ic - 2 * logJ_sum
    }
    
    if (all(!is.finite(ic_by_family))) {
      if (!is.null(progressor)) progressor(message = feature_id)
      return(list(
        selection = tibble::tibble(
          feature     = feature_id,
          best_family = NA_character_,
          n_valid_obs = n_valid_obs,
          ic_value    = NA_real_
        ),
        results = tibble::tibble(
          feature  = feature_id,
          term     = NA_character_,
          effect   = NA_real_,
          se       = NA_real_,
          stat     = NA_real_,
          pval     = NA_real_
        )
      ))
    }
    
    best_family <- names(which.min(ic_by_family))
    best_fit    <- fits_by_family[[best_family]]
    coef_tbl    <- extract_mu_coef_table(best_fit)
    
    if (!is.null(progressor)) progressor(message = feature_id)
    
    list(
      selection = tibble::tibble(
        feature     = feature_id,
        best_family = best_family,
        n_valid_obs = n_valid_obs,
        ic_value    = ic_by_family[best_family]
      ),
      results = coef_tbl |> dplyr::mutate(feature = feature_id, .before = 1)
    )
  }
  
  # --- parallel plan (if requested) ---
  if (workers > 1) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
  }
  
  # --- run with/without progress bar ---
  run_apply <- function(fun) {
    if (workers > 1) {
      future.apply::future_lapply(feature_ids, fun, future.seed = TRUE)
    } else {
      lapply(feature_ids, fun)
    }
  }
  
  if (isTRUE(show_progress)) {
    out_list <- progressr::with_progress({
      p <- progressr::progressor(steps = length(feature_ids), message = progress_label)
      run_apply(function(fid) process_one_feature(fid, progressor = p))
    })
  } else {
    out_list <- run_apply(function(fid) process_one_feature(fid, progressor = NULL))
  }
  
  # --- bind & adjust p-values by term ---
  selection <- dplyr::bind_rows(lapply(out_list, `[[`, "selection"))
  results   <- dplyr::bind_rows(lapply(out_list, `[[`, "results"))
  
  if (nrow(results) > 0 && "term" %in% names(results)) {
    results <- results |>
      dplyr::group_by(term) |>
      dplyr::mutate(padj = p.adjust(pval, method = p_adjust)) |>
      dplyr::ungroup()
  } else {
    results$padj <- NA_real_
  }
  
  list(results = results, selection = selection)
}
