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
#' @param n_genes integer, number of features per bootstrap pull (default: 200)
#' @param n_boot integer, number of bootstrap pulls (default: 10)
#' @param top_n integer, number of top families to return (default: 4)
#' @param families character vector of families to consider; if NULL, a default set is used
#' @param verbose logical, print progress (default: TRUE)
#' @param min_n integer, min valid samples per feature after common mask (default: 5)
#' @param seed integer or NULL, seed for reproducibility (default: NULL)
#' @param group_by_support logical, restrict families per feature by inferred support (default: TRUE)
#' @param binom_bd NULL | numeric scalar | numeric vector len = ncol(counts_matrix).
#'   Binomial denominator for BI/BB. If NULL, infer per feature as max(y) when consistent.
#' @param criterion "GAIC", "BIC", or "AIC" (default: "GAIC")
#' @param gaic_k numeric, penalty for GAIC. If NULL, defaults to log(n_valid)
#' @param filter_beta_inflated logical, filter inflated Beta families if no 0/1 evidence (default: TRUE)
#' @param thr_zero numeric threshold for zero inflation evidence (default: 0.005)
#' @param thr_one numeric threshold for one inflation evidence (default: 0.005)
#'
#' @return list(
#'   top_families_overall   = character[],
#'   top_families_by_support= list(count=..., unit=..., positive=..., real=...),
#'   freq_table_overall     = table,
#'   prop_table_overall     = named numeric,
#'   freq_by_support        = named list of tables,
#'   prop_by_support        = named list of named numerics,
#'   sampled_results        = tibble with per-feature decisions across pulls
#' )
#' @export
find_families <- function(counts_matrix,
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
                          thr_one  = 0.005) {
  
  criterion <- match.arg(criterion)
  
  # ---- defaults ----
  if (is.null(families)) {
    families <- c(
      # unit interval (plain + inflated)
      "BE", "BEO", "BEINF", "BEZI", "BEo", "BEINF0",
      # counts incl. binomial-like
      "PO", "NBI", "ZIP", "ZINBI", "BI", "BB",
      # positive continuous
      "GA", "GG", "IG", "LOGNO",
      # real-valued
      "NO", "TF", "GU"
    )
  }
  
  # ---- checks ----
  if (n_genes > nrow(counts_matrix)) {
    stop("n_genes exceeds the number of rows in counts_matrix.")
  }
  if (!is.null(seed)) set.seed(seed)
  feature_ids_all <- rownames(counts_matrix)
  
  # ---- parallel plan ----
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = future::availableCores())
  
  # ---- helpers ----
  is_binom_fam <- function(f) f %in% c("BI","BB")
  
  get_bd_vec <- function(y, common_mask) {
    # Prefer user-supplied binom_bd if present
    if (!is.null(binom_bd)) {
      if (length(binom_bd) == 1L) {
        return(rep(as.numeric(binom_bd), sum(common_mask)))
      } else if (length(binom_bd) == ncol(counts_matrix)) {
        return(as.numeric(binom_bd[common_mask]))
      } else {
        warning("binom_bd has incompatible length; will attempt per-feature inference.")
      }
    }
    # Per-feature inference: bd := max(y[mask]) if consistent
    yy <- y[common_mask]
    bd_guess <- suppressWarnings(max(yy[is.finite(yy)], na.rm = TRUE))
    if (!is.finite(bd_guess) || bd_guess <= 0) return(NULL)
    bd_guess <- round(bd_guess)
    ok <- all(abs(yy - round(yy)) < 1e-8 & yy >= 0 & yy <= bd_guess, na.rm = TRUE)
    if (!ok) return(NULL)
    rep(bd_guess, length(yy))
  }
  
  # core per-feature evaluation
  eval_one_feature <- function(feature_id) {
    y <- as.numeric(counts_matrix[feature_id, ])
    if (!any(is.finite(y)) || length(unique(y[is.finite(y)])) <= 1) {
      return(list(feature = feature_id, family = NA_character_, skipped = TRUE,
                  n_valid = 0L, support = NA_character_))
    }
    if (all(y == 0, na.rm = TRUE)) {
      return(list(feature = feature_id, family = NA_character_, skipped = TRUE,
                  n_valid = 0L, support = "count"))
    }
    
    fams_try <- families
    grp <- infer_support(y)   # from utils_transformations.R
    
    if (isTRUE(group_by_support)) {
      fg <- family_groups()
      fams_try <- intersect(families, fg[[grp]])
      if (length(fams_try) == 0) {
        return(list(feature = feature_id, family = NA_character_, skipped = TRUE,
                    n_valid = 0L, support = grp))
      }
    }
    
    # Prefilter inflated Beta families if no evidence of 0/1
    if (filter_beta_inflated && grp == "unit") {
      y_fin <- y[is.finite(y)]
      p0 <- mean(y_fin == 0)
      p1 <- mean(y_fin == 1)
      if (!is.na(p0) && !is.na(p1) && p0 < thr_zero && p1 < thr_one) {
        fams_try <- setdiff(fams_try, c("BEINF","BEZI","BEINF0","BEO","BEo"))
        if (length(fams_try) == 0) {
          return(list(feature = feature_id, family = NA_character_, skipped = TRUE,
                      n_valid = 0L, support = grp))
        }
      }
    }
    
    # strict transforms for all candidate families
    tr_list <- lapply(fams_try, function(f) transform_for_family_strict(y, f))
    names(tr_list) <- fams_try
    
    # keep families with enough valid obs
    enough <- vapply(tr_list, function(tr) sum(tr$mask & is.finite(tr$y)) >= min_n, logical(1))
    if (!any(enough)) {
      return(list(feature = feature_id, family = NA_character_, skipped = TRUE,
                  n_valid = 0L, support = grp))
    }
    fams_eff <- fams_try[enough]; tr_list <- tr_list[enough]
    
    # common mask across remaining families
    masks <- lapply(tr_list, `[[`, "mask")
    common_mask <- Reduce(`&`, masks)
    n_valid <- sum(common_mask)
    if (n_valid < min_n) {
      return(list(feature = feature_id, family = NA_character_, skipped = TRUE,
                  n_valid = as.integer(n_valid), support = grp))
    }
    
    # choose k for IC
    k_for <- function(nv) switch(criterion,
                                 "AIC"  = 2,
                                 "BIC"  = log(nv),
                                 "GAIC" = if (is.null(gaic_k)) log(nv) else gaic_k
    )
    k <- k_for(n_valid)
    
    # fit intercept-only and compute Jacobian-corrected IC
    ic_corr <- rep(Inf, length(fams_eff)); names(ic_corr) <- fams_eff
    
    need_bd <- any(is_binom_fam(fams_eff))
    bd_vec <- if (need_bd) get_bd_vec(y, common_mask) else NULL
    
    for (i in seq_along(fams_eff)) {
      fam <- fams_eff[i]
      tr  <- tr_list[[fam]]
      z   <- tr$y[common_mask]
      logJ_sum <- sum(tr$logJ_per_obs[common_mask])
      
      fam_fun <- if (exists(fam, where = asNamespace("gamlss.dist"), inherits = FALSE)) {
        get(fam, envir = asNamespace("gamlss.dist"))
      } else if (exists(fam, mode = "function")) {
        get(fam)
      } else NULL
      if (is.null(fam_fun)) next
      
      fam_obj <- tryCatch({
        if (is_binom_fam(fam)) {
          if (is.null(bd_vec)) return(NULL) # skip BI/BB if bd unavailable
          fam_fun(bd = bd_vec)
        } else {
          fam_fun()
        }
      }, error = function(e) NULL)
      if (is.null(fam_obj)) next
      
      fit <- tryCatch({
        suppressWarnings(gamlss::gamlss(y ~ 1,
                                        data = data.frame(y = z),
                                        family = fam_obj,
                                        trace = FALSE))
      }, error = function(e) NULL)
      if (is.null(fit)) next
      
      ll_t <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
      if (!is.finite(ll_t)) next
      df_fit <- fit$df.fit
      
      base_ic <- (-2 * ll_t) + k * df_fit
      ic_corr[i] <- base_ic - 2 * logJ_sum
    }
    
    if (all(!is.finite(ic_corr))) {
      return(list(feature = feature_id, family = NA_character_, skipped = TRUE,
                  n_valid = as.integer(n_valid), support = grp))
    }
    
    best_fam <- names(which.min(ic_corr))
    list(feature = feature_id, family = best_fam, skipped = FALSE,
         n_valid = as.integer(n_valid), support = grp)
  }
  
  # ---- bootstrap pulls ----
  if (isTRUE(verbose)) {
    message("Running find_families with ", n_boot, " bootstrap pulls of ", n_genes, " features each...")
  }
  all_rows <- vector("list", n_boot)
  
  for (b in seq_len(n_boot)) {
    if (isTRUE(verbose)) message("  Pull ", b, " of ", n_boot)
    pull_ids <- sample(feature_ids_all, n_genes, replace = FALSE)
    
    per_feature <- future.apply::future_lapply(pull_ids, eval_one_feature, future.seed = TRUE)
    
    df_b <- tibble::tibble(
      feature = vapply(per_feature, `[[`, "", "feature"),
      family  = vapply(per_feature, function(x) if (is.null(x$family)) NA_character_ else x$family, ""),
      skipped = vapply(per_feature, `[[`, FALSE, "skipped"),
      n_valid = vapply(per_feature,
                       function(x) { nv <- x$n_valid; if (is.null(nv) || is.na(nv)) NA_integer_ else as.integer(nv) },
                       integer(1L)),
      support = vapply(per_feature, function(x) if (is.null(x$support)) NA_character_ else x$support, "")
    ) %>%
      dplyr::mutate(bootstrap = b, .before = 1L)
    
    all_rows[[b]] <- df_b
  }
  
  sampled_results <- dplyr::bind_rows(all_rows)
  
  # ---- aggregate ----
  fitted <- sampled_results %>% dplyr::filter(!is.na(.data$family), !.data$skipped)
  
  fam_freq <- sort(table(fitted$family), decreasing = TRUE)
  fam_prop <- as.numeric(fam_freq) / sum(fam_freq); names(fam_prop) <- names(fam_freq)
  top_families_overall <- names(head(fam_freq, top_n))
  
  supports <- c("count", "unit", "positive", "real")
  freq_by_support <- lapply(supports, function(sup) {
    tab <- table(fitted$family[fitted$support == sup])
    sort(tab, decreasing = TRUE)
  })
  names(freq_by_support) <- supports
  
  prop_by_support <- lapply(freq_by_support, function(tab) {
    if (length(tab) == 0) return(numeric(0))
    pr <- as.numeric(tab) / sum(tab); names(pr) <- names(tab); pr
  })
  
  top_by_support <- lapply(freq_by_support, function(tab) names(head(tab, top_n)))
  
  if (isTRUE(verbose)) {
    message("\n===== Summary Report =====")
    message("Bootstrap pulls       : ", n_boot)
    message("Features per pull     : ", n_genes)
    message("Successful fit rows   : ", nrow(fitted), " (across pulls)")
    if (length(fam_freq) > 0) {
      message("Most frequent family  : ", names(fam_freq)[1], " (", fam_freq[1], ")")
      message("Top families overall  : ", paste(head(names(fam_freq), top_n), collapse = ", "))
    } else {
      message("No successful fits.")
    }
  }
  
  list(
    top_families_overall    = top_families_overall,
    top_families_by_support = top_by_support,
    freq_table_overall      = fam_freq,
    prop_table_overall      = fam_prop,
    freq_by_support         = freq_by_support,
    prop_by_support         = prop_by_support,
    sampled_results         = sampled_results
  )
}
