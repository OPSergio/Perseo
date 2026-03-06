#' Identify coefficients corresponding to a factor variable
#'
#' Given a contrast matrix or variable name, identifies which model coefficients
#' belong to the target factor for omnibus testing.
#'
#' @param contrast_matrix Numeric matrix with contrast definitions (rows = contrasts,
#'   cols = coefficients). Can be NULL if using contrast_variable.
#' @param contrast_variable Character string naming the factor variable. Can be NULL
#'   if using contrast_matrix.
#' @param coef_names Character vector of coefficient names from the fitted model.
#'
#' @return Character vector of coefficient names belonging to the target factor,
#'   excluding the intercept. Returns NULL if no matching coefficients found.
#'
#' @keywords internal
identify_factor_coefficients <- function(contrast_matrix = NULL,
                                        contrast_variable = NULL,
                                        coef_names) {
  # Strategy 1: Extract from contrast_matrix
  if (!is.null(contrast_matrix)) {
    # Find columns with non-zero entries (excluding intercept)
    involved_coefs <- colnames(contrast_matrix)[
      colSums(abs(contrast_matrix)) > 0
    ]
    involved_coefs <- setdiff(involved_coefs, "(Intercept)")
    
    # Validate they exist in model
    matched <- intersect(involved_coefs, coef_names)
    
    if (length(matched) > 0) {
      return(matched)
    }
  }
  
  # Strategy 2: Extract from contrast_variable name
  if (!is.null(contrast_variable)) {
    # Find coefficients starting with variable name
    pattern <- paste0("^", contrast_variable)
    matched <- grep(pattern, coef_names, value = TRUE)
    
    if (length(matched) > 0) {
      return(matched)
    }
  }
  
  # No matches found
  return(NULL)
}


#' Perform Wald omnibus test for factor coefficients
#'
#' Tests the null hypothesis that all coefficients for a factor are jointly zero
#' using a multi-degree-of-freedom Wald test.
#'
#' @param beta Named numeric vector of coefficients.
#' @param V Variance-covariance matrix (must match beta).
#' @param factor_coefs Character vector naming which coefficients to test.
#'
#' @return List with:
#'   \describe{
#'     \item{statistic}{Wald chi-square statistic}
#'     \item{df}{Degrees of freedom (length of factor_coefs)}
#'     \item{p_value}{P-value from chi-square distribution}
#'     \item{test_type}{"Wald"}
#'   }
#'
#' @details Test statistic: W = beta' V^-1 beta ~ chi^2(df)
#'   If V is singular or non-invertible, returns NA values.
#'
#' @keywords internal
wald_omnibus_test <- function(beta, V, factor_coefs) {
  # Validate inputs
  if (!all(factor_coefs %in% names(beta))) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "Wald"
    ))
  }
  
  # Extract relevant subset
  beta_subset <- beta[factor_coefs]
  V_subset <- V[factor_coefs, factor_coefs, drop = FALSE]
  
  # Check for non-finite values
  if (any(!is.finite(beta_subset)) || any(!is.finite(V_subset))) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "Wald"
    ))
  }
  
  # Try to invert V_subset
  V_inv <- tryCatch({
    solve(V_subset)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(V_inv)) {
    # Singular matrix, try pseudo-inverse
    V_inv <- tryCatch({
      MASS::ginv(V_subset)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  if (is.null(V_inv)) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "Wald"
    ))
  }
  
  # Compute Wald statistic: beta' V^-1 beta
  wald_stat <- tryCatch({
    as.numeric(t(beta_subset) %*% V_inv %*% beta_subset)
  }, error = function(e) {
    return(NA_real_)
  })
  
  if (is.na(wald_stat) || wald_stat < 0) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "Wald"
    ))
  }
  
  # Degrees of freedom
  df <- length(factor_coefs)
  
  # P-value from chi-square distribution
  p_value <- pchisq(wald_stat, df = df, lower.tail = FALSE)
  
  list(
    statistic = wald_stat,
    df = df,
    p_value = p_value,
    test_type = "Wald"
  )
}


#' Perform LRT omnibus test for factor coefficients
#'
#' Tests the null hypothesis that a factor has no effect by comparing log-likelihoods
#' of full vs reduced models.
#'
#' @param feature_vec Numeric vector of observations for one feature.
#' @param design_mat Full design matrix.
#' @param factor_coefs Character vector of coefficients to exclude in reduced model.
#' @param family GAMLSS family string.
#' @param transform_mode "strict" or "safe".
#' @param bd_vec Optional binomial denominator vector.
#' @param formula_mu Optional full model formula (response must be `y`) for
#'   direct GAMLSS fitting.
#' @param metadata_df Optional metadata data.frame used with `formula_mu`.
#' @param contrast_variable Optional variable name used to define reduced model
#'   when `formula_mu` is provided.
#'
#' @return List with:
#'   \describe{
#'     \item{statistic}{LRT statistic: 2*(logLik_full - logLik_reduced)}
#'     \item{df}{Degrees of freedom (number of excluded coefficients)}
#'     \item{p_value}{P-value from chi-square distribution}
#'     \item{test_type}{"LRT"}
#'   }
#'
#' @details Fits both models and computes likelihood ratio.
#'   Returns NA if either model fails to converge.
#'
#' @keywords internal
lrt_omnibus_test <- function(feature_vec,
                            design_mat,
                            factor_coefs,
                            family,
                            transform_mode = "strict",
                            bd_vec = NULL,
                            formula_mu = NULL,
                            metadata_df = NULL,
                            contrast_variable = NULL) {
  # Transform response
  tr <- tryCatch({
    transform_response(feature_vec, family, mode = transform_mode)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(tr)) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "LRT"
    ))
  }
  
  y <- tr$y[tr$mask]

  use_formula <- !is.null(formula_mu)
  if (use_formula) {
    if (is.null(metadata_df) || !is.data.frame(metadata_df)) {
      return(list(
        statistic = NA_real_,
        df = NA_integer_,
        p_value = NA_real_,
        test_type = "LRT"
      ))
    }
    if (is.null(contrast_variable) || !is.character(contrast_variable) || length(contrast_variable) != 1) {
      return(list(
        statistic = NA_real_,
        df = NA_integer_,
        p_value = NA_real_,
        test_type = "LRT"
      ))
    }

    metadata_subset <- metadata_df[tr$mask, , drop = FALSE]
    data_full <- cbind.data.frame(y = y, metadata_subset)

    terms_obj <- stats::terms(formula_mu)
    term_labels <- attr(terms_obj, "term.labels")
    has_intercept <- attr(terms_obj, "intercept")
    var_pattern <- paste0("(^|[^[:alnum:]_])", contrast_variable, "([^[:alnum:]_]|$)")
    reduced_terms <- term_labels[!grepl(var_pattern, term_labels)]

    reduced_formula <- stats::reformulate(
      termlabels = reduced_terms,
      response = "y",
      intercept = has_intercept
    )
    environment(reduced_formula) <- environment(formula_mu)
  } else {
    design_subset <- design_mat[tr$mask, , drop = FALSE]

    all_cols <- colnames(design_subset)
    reduced_cols <- setdiff(all_cols, factor_coefs)

    if (length(reduced_cols) == 0) {
      # Can't fit reduced model (would be empty or intercept-only without data)
      return(list(
        statistic = NA_real_,
        df = NA_integer_,
        p_value = NA_real_,
        test_type = "LRT"
      ))
    }

    design_full <- design_subset
    design_reduced <- design_subset[, reduced_cols, drop = FALSE]

    data_full <- cbind.data.frame(y = y, design_full)
    data_reduced <- cbind.data.frame(y = y, design_reduced)
    reduced_formula <- y ~ .
  }
  
  # Instantiate family
  fam_obj <- tryCatch({
    instantiate_gamlss_family(
      family,
      bd_vec = if (is_binomial_family(family)) {
        infer_binomial_denominator(feature_vec, tr$mask)
      } else NULL
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fam_obj)) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "LRT"
    ))
  }
  
  # Fit full model
  fit_full <- fit_gamlss_safely(if (use_formula) formula_mu else y ~ ., data = data_full, family_obj = fam_obj)
  
  # Fit reduced model
  fit_reduced <- fit_gamlss_safely(reduced_formula, data = if (use_formula) data_full else data_reduced, family_obj = fam_obj)
  
  if (is.null(fit_full) || is.null(fit_reduced)) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "LRT"
    ))
  }
  
  # Extract log-likelihoods
  loglik_full <- tryCatch({
    as.numeric(logLik(fit_full))
  }, error = function(e) NA_real_)
  
  loglik_reduced <- tryCatch({
    as.numeric(logLik(fit_reduced))
  }, error = function(e) NA_real_)
  
  if (is.na(loglik_full) || is.na(loglik_reduced)) {
    return(list(
      statistic = NA_real_,
      df = NA_integer_,
      p_value = NA_real_,
      test_type = "LRT"
    ))
  }
  
  # Compute LRT statistic
  lrt_stat <- 2 * (loglik_full - loglik_reduced)
  
  # Handle negative LRT (can happen due to numerical issues)
  if (lrt_stat < 0) {
    lrt_stat <- 0
  }
  
  df <- length(factor_coefs)
  p_value <- pchisq(lrt_stat, df = df, lower.tail = FALSE)
  
  list(
    statistic = lrt_stat,
    df = df,
    p_value = p_value,
    test_type = "LRT"
  )
}
