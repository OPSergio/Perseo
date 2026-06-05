#' Filter candidate families based on empirical support and heuristics
#'
#' Pure function that determines which GAMLSS families are eligible for a given
#' feature vector. Handles support-based filtering (when enabled) and beta
#' inflation heuristics.
#'
#' @param y Numeric vector with feature values across samples.
#' @param families Character vector with candidate family names.
#' @param group_by_support Logical; if TRUE, restricts families to those
#'   matching the empirical support inferred from `y`. If FALSE, all families
#'   in `families` are retained (subject only to beta inflation filtering).
#' @param filter_beta_inflated Logical; if TRUE and support is "unit", removes
#'   inflated beta families (BEINF, BEZI, etc.) when there is insufficient
#'   evidence of zero/one inflation.
#' @param thr_zero Numeric threshold for declaring zero inflation (proportion
#'   of exact zeros required).
#' @param thr_one Numeric threshold for declaring one inflation (proportion
#'   of exact ones required).
#'
#' @return List with two elements:
#'   \describe{
#'     \item{families}{Character vector of eligible families (possibly empty).}
#'     \item{support}{Character scalar with inferred support category.}
#'   }
#'
#' @details
#' When `group_by_support = FALSE`, all families are tested regardless of
#' inferred support, allowing exploration of model fit across support boundaries
#' (e.g., fitting count models to continuous data). This increases computation
#' and may cause convergence issues but can be useful for exploratory analysis.
#'
#' @keywords internal
filter_candidate_families <- function(feature_vec,
                                      candidate_families,
                                      group_by_support = TRUE,
                                      filter_beta_inflated = TRUE,
                                      thr_zero = 0.005,
                                      thr_one = 0.005) {
  support <- infer_support(feature_vec)
  eligible_families <- candidate_families
  y <- feature_vec

  # Support-based filtering (preserves group_by_support = FALSE)
  if (isTRUE(group_by_support)) {
    groups <- family_groups()
    if (support == "zi_positive") {
      # Include both ZI families and regular positive families: let IC decide
      # whether modeling zero-inflation is worth the extra parameters
      eligible_families <- intersect(
        eligible_families,
        c(groups[["zi_positive"]], groups[["positive"]])
      )
    } else {
      eligible_families <- intersect(eligible_families, groups[[support]])
    }
  }

  # ZI/ZA inflation filtering (applies regardless of group_by_support)
  if (filter_beta_inflated && length(eligible_families) > 0) {
    y_finite <- y[is.finite(y)]
    if (length(y_finite) > 0) {
      p_zero <- mean(y_finite == 0)
      p_one  <- mean(y_finite == 1)

      # Unit support: drop inflated beta when no 0/1 evidence
      if (identical(support, "unit") && !is.na(p_zero) && !is.na(p_one) &&
          p_zero < thr_zero && p_one < thr_one) {
        eligible_families <- setdiff(eligible_families,
                                     c("BEINF", "BEZI", "BEINF0", "BEO", "BEo"))
      }

      # Count support: drop ZIP/ZINBI when no zero-inflation evidence
      if (identical(support, "count") && !is.na(p_zero) && p_zero < thr_zero) {
        eligible_families <- setdiff(eligible_families, c("ZIP", "ZINBI", "ZIP2"))
      }

      # Positive support: drop ZAIG/ZAGA when no zero-augmentation evidence
      if (identical(support, "positive") && !is.na(p_zero) && p_zero < thr_zero) {
        eligible_families <- setdiff(eligible_families, c("ZAIG", "ZAGA"))
      }
    }
  }

  bd_vec <- if (any(is_binomial_family(eligible_families))) {
    infer_binomial_denominator(feature_vec)
  } else NULL
  
  list(families_to_test = eligible_families, support = support, bd_vec = bd_vec)
}


#' Check if a family requires binomial denominators
#'
#' @param family_name Character scalar with GAMLSS family name.
#'
#' @return Logical; TRUE for BI/BB families, FALSE otherwise.
#' @keywords internal
is_binomial_family <- function(family_name) {
  family_name %in% c("BI", "BB")
}


#' Infer or validate binomial denominators for BI/BB families
#'
#' Handles user-supplied denominators (scalar or vector) with validation and
#' fallback to per-feature inference when denominators are not provided.
#'
#' @param y Numeric vector with original feature counts.
#' @param common_mask Optional logical mask indicating valid observations for fitting.
#' @param binom_bd_user User-supplied binomial denominator: NULL (infer),
#'   numeric scalar (recycle), or numeric vector (length = total samples).
#' @param n_samples Integer total number of samples in the original matrix.
#'
#' @return Numeric vector of denominators with length = sum(common_mask), or
#'   NULL if inference fails or denominators are invalid.
#'
#' @details
#' Inference logic: bd := max(y[mask]) if all masked values are integers in
#' [0, max(y)] with no violations. This is conservative and may fail for
#' irregular data.
#'
#' @keywords internal
infer_binomial_denominator <- function(y,
                                       common_mask = NULL,
                                       binom_bd_user = NULL,
                                       n_samples = NULL) {
  # Default mask to all valid
  if (is.null(common_mask)) {
    common_mask <- is.finite(y)
  }
  
  n_valid <- sum(common_mask)
  if (n_valid == 0) {
    return(NULL)
  }

  # User-supplied denominator handling
  if (!is.null(binom_bd_user)) {
    if (length(binom_bd_user) == 1L) {
      return(rep(as.numeric(binom_bd_user), n_valid))
    }
    if (!is.null(n_samples) && length(binom_bd_user) == n_samples) {
      return(as.numeric(binom_bd_user[common_mask]))
    }
    warning(
      "binom_bd has length ", length(binom_bd_user),
      " but expected 1 or ", n_samples, "; falling back to inference."
    )
  }

  # Per-feature inference
  y_masked <- y[common_mask]
  bd_guess <- suppressWarnings(max(y_masked[is.finite(y_masked)], na.rm = TRUE))
  
  if (!is.finite(bd_guess) || bd_guess <= 0) {
    return(NULL)
  }

  bd_guess <- round(bd_guess)
  is_valid <- all(
    abs(y_masked - round(y_masked)) < 1e-8 &
    y_masked >= 0 &
    y_masked <= bd_guess,
    na.rm = TRUE
  )

  if (!is_valid) {
    return(NULL)
  }

  rep(bd_guess, n_valid)
}


#' Check if feature has insufficient variation for modeling
#'
#' @param y Numeric vector with feature values.
#'
#' @return Logical; TRUE if feature should be skipped.
#' @export
has_insufficient_variation <- function(y) {
  y_finite <- y[is.finite(y)]
  n_unique <- length(unique(y_finite))
  # Constant or nearly constant (only 2 values with one rare)
  if (n_unique < 2) return(TRUE)
  if (n_unique == 2) {
    counts <- table(y_finite)
    return(min(counts) == 1)  # One value appears only once
  }
  FALSE
}
