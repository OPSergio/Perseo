#' Transform expression values for compatibility with GAMLSS family
#'
#' Applies family-specific preprocessing and transformations to the expression vector `y`,
#' ensuring values are in the domain required by each distribution.
#'
#' @param y Numeric vector of expression values.
#' @param fam Character: name of the GAMLSS family.
#' @param strategy Character: either "safe" (default) or "strict". "Safe" clips or rescales values to avoid errors; "strict" replaces invalid values with NA.
#' @param eps Small numeric value used for smoothing and clipping.
#'
#' @return A numeric vector with transformed values (same length as `y`)
#' @export