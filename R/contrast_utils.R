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
