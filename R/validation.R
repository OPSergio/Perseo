#' Validate counts matrix structure and dimensions
#'
#' @param counts_matrix Numeric matrix with features in rows and samples in columns.
#' @param min_features Optional integer minimum row count to enforce.
#'
#' @return Invisible NULL; throws informative errors if validation fails.
#' @keywords internal
validate_counts_matrix <- function(counts_matrix, min_features = NULL) {
  if (!is.matrix(counts_matrix)) {
    stop("counts_matrix must be a matrix; got ", class(counts_matrix)[1])
  }
  
  if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
    stop("counts_matrix must have at least one row and one column.")
  }
  
  if (!is.null(min_features) && nrow(counts_matrix) < min_features) {
    stop(
      "counts_matrix has ", nrow(counts_matrix), " features but ",
      min_features, " are required."
    )
  }
  
  invisible(NULL)
}


#' Validate and clean design matrix
#'
#' Ensures design matrix matches sample count and removes duplicate intercept
#' columns to avoid collinearity with the formula intercept.
#'
#' @param design_matrix Data frame or matrix with covariates.
#' @param n_samples Expected number of rows (must match sample count).
#'
#' @return Data frame with cleaned design matrix.
#' @keywords internal
validate_design_matrix <- function(design_matrix, n_samples) {
  design_df <- as.data.frame(design_matrix)
  
  if (nrow(design_df) != n_samples) {
    stop(
      "design_matrix has ", nrow(design_df), " rows but counts_matrix has ",
      n_samples, " samples (columns)."
    )
  }
  
  # Remove explicit intercept column to avoid duplication in formula
  if ("(Intercept)" %in% colnames(design_df)) {
    design_df <- design_df[, setdiff(colnames(design_df), "(Intercept)"), drop = FALSE]
  }
  
  design_df
}


#' Validate information criterion arguments
#'
#' @param criterion Character scalar (already match.arg'd).
#' @param gaic_k Optional numeric penalty for GAIC.
#'
#' @return Invisible NULL; throws error if gaic_k is invalid.
#' @keywords internal
validate_criterion_args <- function(criterion, gaic_k = NULL) {
  if (!is.null(gaic_k)) {
    if (!is.numeric(gaic_k) || length(gaic_k) != 1 || gaic_k <= 0) {
      stop("gaic_k must be a positive numeric scalar or NULL.")
    }
  }
  
  invisible(NULL)
}


#' Default candidate GAMLSS families
#'
#' Returns the standard set of families used when user does not provide a
#' custom list. Includes count, unit, positive, and real-valued distributions.
#'
#' @return Character vector with 18 family names.
#' @keywords internal
default_candidate_families <- function() {
  c(
    # Unit interval (plain + inflated)
    "BE", "BEO", "BEINF", "BEZI", "BEo", "BEINF0",
    # Counts (including zero-inflated and binomial)
    "PO", "NBI", "ZIP", "ZINBI", "BI", "BB",
    # Positive continuous
    "GA", "GG", "IG", "LOGNO",
    # Real-valued
    "NO", "TF", "GU"
  )
}
