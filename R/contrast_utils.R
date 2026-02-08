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


#' Build pairwise contrast matrix from categorical variable
#'
#' Creates a contrast matrix for all pairwise comparisons of a
#' categorical variable's levels.
#'
#' @param var_name Character string naming the variable in metadata.
#' @param metadata Data frame containing sample metadata.
#' @param design_matrix Design matrix with coefficient names.
#'
#' @return Numeric matrix with one row per pairwise contrast. Column names
#'   match design matrix coefficients, row names describe comparisons.
#'
#' @keywords internal
build_contrast_matrix <- function(var_name, metadata, design_matrix) {
  if (!var_name %in% names(metadata)) {
    stop("contrast_variable '", var_name, "' not found in metadata")
  }
  
  var_values <- metadata[[var_name]]
  
  if (!is.factor(var_values)) {
    var_values <- as.factor(var_values)
  }
  
  levels_vec <- levels(var_values)
  
  if (length(levels_vec) < 2) {
    stop("contrast_variable '", var_name, "' must have at least 2 levels")
  }
  
  # Get coefficient names from design matrix
  coef_names <- colnames(design_matrix)
  
  # Find which coefficients correspond to this variable
  # Assuming standard R factor coding: first level absorbed in intercept
  var_coef_pattern <- paste0("^", var_name)
  var_coefs <- grep(var_coef_pattern, coef_names, value = TRUE)
  
  if (length(var_coefs) == 0) {
    stop("No coefficients found matching variable '", var_name, "'")
  }
  
  # Generate all pairwise contrasts
  n_levels <- length(levels_vec)
  pairs <- combn(n_levels, 2, simplify = FALSE)
  
  contrast_list <- lapply(pairs, function(pair) {
    level1 <- levels_vec[pair[1]]
    level2 <- levels_vec[pair[2]]
    
    # Create contrast vector
    contrast_vec <- rep(0, length(coef_names))
    names(contrast_vec) <- coef_names
    
    # Determine coefficient names for each level
    # First level (reference) has no coefficient (absorbed in intercept)
    coef1 <- if (pair[1] == 1) {
      "(Intercept)"
    } else {
      paste0(var_name, level1)
    }
    
    coef2 <- if (pair[2] == 1) {
      "(Intercept)"
    } else {
      paste0(var_name, level2)
    }
    
    # Set contrast: level2 - level1
    if (coef1 %in% coef_names) {
      contrast_vec[coef1] <- -1
    }
    if (coef2 %in% coef_names) {
      contrast_vec[coef2] <- 1
    }
    
    list(
      name = paste0(level2, "_vs_", level1),
      vector = contrast_vec
    )
  })
  
  # Build matrix
  contrast_matrix <- do.call(rbind, lapply(contrast_list, `[[`, "vector"))
  rownames(contrast_matrix) <- vapply(contrast_list, `[[`, "name", FUN.VALUE = character(1))
  
  contrast_matrix
}
