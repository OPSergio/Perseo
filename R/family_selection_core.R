#' Compare GAMLSS families for a single feature (intercept-only models)
#'
#' Core selection logic: transforms response with all candidate families, builds
#' a common mask, fits intercept-only GAMLSS models, and selects the best family
#' using Jacobian-corrected IC.
#'
#' @param feature_vec Numeric vector with feature values across samples.
#' @param families Character vector of eligible family names (already
#'   filtered by support/inflation heuristics).
#' @param min_n Integer minimum number of valid observations after common mask.
#' @param criterion Character IC criterion ("AIC", "BIC", "GAIC").
#' @param gaic_k Optional numeric GAIC penalty override.
#' @param bd_vec Optional binomial denominator vector.
#' @param transform_mode Character: "strict" or "safe" transformation mode.
#'
#' @return List with comparisons tibble containing family results.
#'
#' @keywords internal
compare_families_on_feature <- function(feature_vec,
                                        families,
                                        criterion = "GAIC",
                                        gaic_k = NULL,
                                        min_n = 5,
                                        bd_vec = NULL,
                                        transform_mode = "strict",
                                        extract_coef_fn = NULL) {
  y <- feature_vec
  candidate_families <- families
  
  # Use default extractor if not provided
  if (is.null(extract_coef_fn)) {
    extract_coef_fn <- extract_mu_coefficients
  }
  
  # Transform with all candidate families
  transforms <- stats::setNames(
    lapply(candidate_families, function(fam) transform_response(y, fam, mode = transform_mode)),
    candidate_families
  )

  # Build common mask across families
  family_masks <- lapply(transforms, `[[`, "mask")
  common_mask <- Reduce(`&`, family_masks)
  n_valid <- sum(common_mask, na.rm = TRUE)

  # Defensive check: n_valid must be finite and >= min_n
  if (!is.finite(n_valid) || n_valid < min_n) {
    return(list(
      comparisons = tibble::tibble(
        family = character(),
        ic_value = numeric(),
        n_valid_obs = integer(),
        coef_tbl = list()
      ),
      best_family = NA_character_,
      n_valid = 0L
    ))
  }

  # Filter families with enough valid observations
  enough_obs <- vapply(
    transforms,
    function(tr) sum(tr$mask & is.finite(tr$y)) >= min_n,
    logical(1)
  )
  
  if (!any(enough_obs)) {
    return(list(
      comparisons = tibble::tibble(
        family = character(),
        ic_value = numeric(),
        n_valid_obs = integer(),
        coef_tbl = list()
      ),
      best_family = NA_character_,
      n_valid = 0L
    ))
  }

  effective_families <- candidate_families[enough_obs]
  transforms <- transforms[enough_obs]

  # Compute IC penalty
  penalty <- compute_ic_penalty(criterion, gaic_k, n_valid)

  # Prepare binomial denominators if needed
  needs_bd <- any(vapply(effective_families, is_binomial_family, logical(1)))
  bd_vec_masked <- if (needs_bd) {
    if (!is.null(bd_vec)) {
      bd_vec[common_mask]
    } else {
      infer_binomial_denominator(y, common_mask)
    }
  } else {
    NULL
  }

  # Fit models and compute IC values
  results_list <- lapply(effective_families, function(fam) {
    tr <- transforms[[fam]]
    z <- tr$y[common_mask]
    logJ_sum <- sum(tr$logJ_per_obs[common_mask])
    
    # Skip binomial families if bd unavailable
    if (is_binomial_family(fam) && is.null(bd_vec_masked)) {
      return(list(ic = Inf, fit = NULL))
    }
    
    # Instantiate family
    fam_obj <- instantiate_gamlss_family(
      fam,
      bd_vec = if (is_binomial_family(fam)) bd_vec_masked else NULL
    )
    
    if (is.null(fam_obj)) {
      return(list(ic = Inf, fit = NULL))
    }
    
    # Fit intercept-only model
    fit <- fit_gamlss_safely(
      y ~ 1,
      data = data.frame(y = z),
      family_obj = fam_obj
    )
    
    if (is.null(fit)) {
      return(list(ic = Inf, fit = NULL))
    }
    
    # Compute Jacobian-corrected IC
    ic <- compute_jacobian_corrected_ic(fit, penalty, logJ_sum)
    list(ic = ic, fit = fit)
  })

  # Build results tibble
  ic_values <- vapply(results_list, `[[`, numeric(1), "ic")
  fits <- lapply(results_list, `[[`, "fit")
  
  valid_results <- is.finite(ic_values)
  if (!any(valid_results)) {
    return(list(
      comparisons = tibble::tibble(
        family = character(),
        ic_value = numeric(),
        n_valid_obs = integer(),
        coef_tbl = list()
      ),
      best_family = NA_character_,
      n_valid = as.integer(n_valid)
    ))
  }
  
  coef_tbls <- lapply(seq_along(effective_families), function(i) {
    if (valid_results[i] && !is.null(fits[[i]])) {
      extract_coef_fn(fits[[i]])
    } else {
      tibble::tibble(term = character(), effect = numeric(), se = numeric(), stat = numeric(), pval = numeric())
    }
  })
  
  # Select best family (minimum IC, near-ties broken by goodness of fit)
  gof_values <- vapply(fits, function(f) if (is.null(f)) -Inf else gof_score(f), numeric(1))
  best_idx <- select_by_ic_gof(ic_values, gof_values)
  best_family <- if (length(best_idx) > 0 && !is.na(best_idx) && is.finite(ic_values[best_idx])) {
    effective_families[best_idx]
  } else {
    NA_character_
  }
  
  list(
    comparisons = tibble::tibble(
      family = effective_families,
      ic_value = ic_values,
      n_valid_obs = rep(as.integer(n_valid), length(effective_families)),
      coef_tbl = coef_tbls
    ),
    best_family = best_family,
    n_valid = as.integer(n_valid)
  )
}


#' Compare GAMLSS families for a single feature (with design matrix)
#'
#' Similar to `compare_families_on_feature()` but fits models with covariates
#' using a design matrix. Returns comparison results with coefficient tables.
#'
#' @param feature_vec Numeric vector with feature values across samples.
#' @param design_df Data frame with covariates (cleaned, no intercept column).
#' @param formula_mu Optional model formula for direct GAMLSS fitting
#'   (e.g., `y ~ condition + random(subject)`). When provided, covariates are
#'   taken from `metadata_df` and `design_df` is ignored for model fitting.
#' @param metadata_df Optional metadata data.frame used with `formula_mu`.
#' @param families Character vector of eligible family names.
#' @param criterion Character IC criterion.
#' @param gaic_k Optional GAIC penalty.
#' @param min_n Integer minimum valid observations.
#' @param bd_vec Optional binomial denominator vector.
#' @param transform_mode Character: "strict" or "safe" transformation mode.
#'
#' @return List with comparisons tibble.
#'
#' @keywords internal
compare_families_with_design <- function(feature_vec,
                                         design_df,
                                         families,
                                         criterion = "GAIC",
                                         gaic_k = NULL,
                                         min_n = 5,
                                         bd_vec = NULL,
                                         transform_mode = "strict",
                                         extract_coef_fn = NULL,
                                         formula_mu = NULL,
                                         metadata_df = NULL) {
  y <- feature_vec
  candidate_families <- families
  design_matrix <- design_df
  use_formula <- !is.null(formula_mu)
  complete_rows <- if (use_formula) {
    if (is.null(metadata_df) || !is.data.frame(metadata_df)) {
      stop("metadata_df must be a data.frame when formula_mu is provided")
    }
    if (nrow(metadata_df) != length(y)) {
      stop("metadata_df must have one row per sample in feature_vec")
    }
    stats::complete.cases(metadata_df)
  } else {
    stats::complete.cases(design_matrix)
  }
  
  # Use default extractor if not provided
  if (is.null(extract_coef_fn)) {
    extract_coef_fn <- extract_mu_coefficients
  }
  
  # Transform with all candidate families
  transforms <- stats::setNames(
    lapply(candidate_families, function(fam) transform_response(y, fam, mode = transform_mode)),
    candidate_families
  )

  # Build common mask: family validity AND design completeness
  family_masks <- lapply(transforms, `[[`, "mask")
  mask_combined <- c(family_masks, list(complete_rows))
  common_mask <- Reduce(`&`, mask_combined)
  n_valid <- sum(common_mask, na.rm = TRUE)

  # Defensive check: n_valid must be finite and >= min_n
  if (!is.finite(n_valid) || n_valid < min_n) {
    return(list(
      comparisons = tibble::tibble(
        family = character(),
        ic_value = numeric(),
        n_valid_obs = integer(),
        coef_tbl = list()
      )
    ))
  }
  # Compute IC penalty
  penalty <- compute_ic_penalty(criterion, gaic_k, n_valid)

  # Prepare masked predictors
  design_masked <- if (use_formula) {
    metadata_df[common_mask, , drop = FALSE]
  } else {
    design_matrix[common_mask, , drop = FALSE]
  }
  
  # Prepare binomial denominators if needed
  needs_bd <- any(vapply(candidate_families, is_binomial_family, logical(1)))
  bd_vec_masked <- if (needs_bd) {
    if (!is.null(bd_vec)) {
      bd_vec[common_mask]
    } else {
      infer_binomial_denominator(y, common_mask)
    }
  } else {
    NULL
  }

  # Fit models and track IC + fitted objects
  results_list <- lapply(candidate_families, function(fam) {
    tr <- transforms[[fam]]
    z <- tr$y[common_mask]
    logJ_sum <- sum(tr$logJ_per_obs[common_mask])
    
    # Combine response and covariates
    fit_data <- cbind.data.frame(y = z, design_masked)
    
    # Instantiate family
    fam_obj <- instantiate_gamlss_family(
      fam,
      bd_vec = if (is_binomial_family(fam)) bd_vec_masked else NULL
    )
    
    if (is.null(fam_obj)) {
      return(list(ic = Inf, fit = NULL))
    }
    
    # Fit with covariates
    fit <- fit_gamlss_safely(
      if (use_formula) formula_mu else y ~ .,
      data = fit_data,
      family_obj = fam_obj
    )
    
    if (is.null(fit)) {
      return(list(ic = Inf, fit = NULL))
    }
    
    ic <- compute_jacobian_corrected_ic(fit, penalty, logJ_sum)
    list(ic = ic, fit = fit)
  })

  # Build results tibble
  ic_values <- vapply(results_list, `[[`, numeric(1), "ic")
  fits <- lapply(results_list, `[[`, "fit")
  
  valid_results <- is.finite(ic_values)
  if (!any(valid_results)) {
    return(list(
      comparisons = tibble::tibble(
        family = character(),
        ic_value = numeric(),
        n_valid_obs = integer(),
        coef_tbl = list()
      )
    ))
  }
  
  coef_tbls <- lapply(seq_along(candidate_families), function(i) {
    if (valid_results[i] && !is.null(fits[[i]])) {
      extract_coef_fn(fits[[i]])
    } else {
      tibble::tibble(term = character(), effect = numeric(), se = numeric(), stat = numeric(), pval = numeric())
    }
  })

  gof_values <- vapply(fits, function(f) if (is.null(f)) -Inf else gof_score(f), numeric(1))

  list(
    comparisons = tibble::tibble(
      family = candidate_families,
      ic_value = ic_values,
      gof = gof_values,
      n_valid_obs = rep(as.integer(n_valid), length(candidate_families)),
      coef_tbl = coef_tbls
    )
  )
}


#' Bind bootstrap results into a single tibble
#'
#' @param results_list List of lists, each containing per-feature results from
#'   one bootstrap pull.
#'
#' @return Tibble with columns: bootstrap, feature, family, skipped, n_valid, support.
#' @keywords internal
bind_bootstrap_results <- function(results_list) {
  all_rows <- lapply(seq_along(results_list), function(b) {
    pull_results <- results_list[[b]]
    
    # Filter out NULL or empty elements
    pull_results <- Filter(function(x) !is.null(x) && length(x) > 0, pull_results)
    
    if (length(pull_results) == 0) {
      return(tibble::tibble(
        bootstrap = integer(0),
        feature = character(0),
        family = character(0),
        skipped = logical(0),
        n_valid = integer(0),
        support = character(0)
      ))
    }
    
    tibble::tibble(
      bootstrap = b,
      feature = vapply(pull_results, function(x) {
        val <- x$feature
        if (is.null(val) || length(val) == 0) return(NA_character_)
        as.character(val)[1]
      }, character(1)),
      family = vapply(pull_results, function(x) {
        val <- x$family
        if (is.null(val) || length(val) == 0) return(NA_character_)
        as.character(val)[1]
      }, character(1)),
      skipped = vapply(pull_results, function(x) {
        val <- x$skipped
        if (is.null(val) || length(val) == 0) return(NA)
        as.logical(val)[1]
      }, logical(1)),
      n_valid = vapply(pull_results, function(x) {
        val <- x$n_valid
        if (is.null(val) || length(val) == 0) return(NA_integer_)
        as.integer(val)[1]
      }, integer(1)),
      support = vapply(pull_results, function(x) {
        val <- x$support
        if (is.null(val) || length(val) == 0) return(NA_character_)
        as.character(val)[1]
      }, character(1))
    )
  })
  
  dplyr::bind_rows(all_rows)
}


#' Summarize family selection frequencies across bootstrap pulls
#'
#' @param sampled_results Tibble with bootstrap results (from `bind_bootstrap_results()`).
#' @param top_n Integer number of top families to return.
#' @param show_progress Logical; print summary report if TRUE.
#'
#' @return List with frequency tables, proportions, and top families.
#' @keywords internal
summarize_family_frequencies <- function(sampled_results, top_n, show_progress = FALSE) {
  fitted <- sampled_results %>%
    dplyr::filter(!is.na(.data$family), !.data$skipped)

  # Overall frequencies
  fam_freq <- sort(table(fitted$family), decreasing = TRUE)
  fam_prop <- if (length(fam_freq) > 0) {
    props <- as.numeric(fam_freq) / sum(fam_freq)
    stats::setNames(props, names(fam_freq))
  } else {
    numeric(0)
  }
  top_overall <- names(head(fam_freq, top_n))

  # By-support frequencies
  supports <- c("count", "unit", "positive", "real")
  freq_by_support <- lapply(supports, function(sup) {
    tab <- table(fitted$family[fitted$support == sup])
    sort(tab, decreasing = TRUE)
  })
  names(freq_by_support) <- supports

  prop_by_support <- lapply(freq_by_support, function(tab) {
    if (length(tab) == 0) return(numeric(0))
    vals <- as.numeric(tab) / sum(tab)
    stats::setNames(vals, names(tab))
  })

  top_by_support <- lapply(freq_by_support, function(tab) {
    names(head(tab, top_n))
  })

  # Enhanced progress reporting
  if (isTRUE(show_progress)) {
    n_boot <- max(sampled_results$bootstrap, na.rm = TRUE)
    n_genes <- nrow(sampled_results) / n_boot
    n_total_eval <- nrow(sampled_results)
    n_skipped <- sum(sampled_results$skipped, na.rm = TRUE)
    n_success <- nrow(fitted)
    success_rate <- n_success / n_total_eval * 100
    
    message("\n", strrep("=", 80))
    message("find_families() - FINAL REPORT")
    message(strrep("=", 80))
    message("\nEvaluation Summary:")
    message("  Bootstrap pulls        : ", n_boot)
    message("  Features per pull      : ", n_genes)
    message("  Total evaluations      : ", n_total_eval)
    message("  Successful fits        : ", n_success, " (", sprintf("%.1f%%", success_rate), ")")
    message("  Skipped (low variation): ", n_skipped, " (", sprintf("%.1f%%", n_skipped/n_total_eval*100), ")")
    
    if (length(fam_freq) > 0) {
      message("\nTop Selected Families (overall):")
      top_display <- head(fam_freq, min(10, length(fam_freq)))
      for (i in seq_along(top_display)) {
        fam_name <- names(top_display)[i]
        count <- as.numeric(top_display[i])
        pct <- count / sum(fam_freq) * 100
        marker <- if (i <= top_n) " *" else ""
        message(sprintf("  %2d. %-10s : %5d (%.1f%%)%s", i, fam_name, count, pct, marker))
      }
      if (top_n < length(top_display)) {
        message("  (* = selected for downstream analysis)")
      }
      
      # Report by support if available
      has_support_data <- any(!is.na(fitted$support))
      if (has_support_data) {
        message("\nTop Families by Data Support:")
        for (sup in supports) {
          tab <- freq_by_support[[sup]]
          if (length(tab) > 0) {
            message(sprintf("  %-8s : %s", sup, paste(head(names(tab), 3), collapse = ", ")))
          }
        }
      }
    } else {
      message("\n[WARNING] No successful fits across any features.")
    }
    
    message(strrep("=", 80), "\n")
  }

  list(
    top_families_overall = top_overall,
    top_families_by_support = top_by_support,
    freq_table_overall = fam_freq,
    prop_table_overall = fam_prop,
    freq_by_support = freq_by_support,
    prop_by_support = prop_by_support
  )
}
