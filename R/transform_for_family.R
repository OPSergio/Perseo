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
 

transform_for_family <- function(y, fam, strategy = "safe", eps = 1e-6) {
  # Batch A: (0, ∞) – strictly positive continuous (e.g. GA, GG, LOGNO, IG)
  if (fam %in% c("GA", "GG", "LOGNO", "IG")) {
    if (strategy == "safe") {
      y[y <= 0] <- eps
    } else {
      y[y <= 0] <- NA
    }
    return(y)
  }

  # Batch B: [0, ∞) – counts (e.g. PO, NBI, ZIP, ZINBI, ZIP2)
  else if (fam %in% c("PO", "NBI", "ZINBI", "ZIP", "ZIP2")) {
    y <- round(y)
    if (strategy == "safe") {
      y[y < 0] <- 0
    } else {
      y[y < 0] <- NA
    }
    return(y)
  }

  # Batch C: (0,1) or inflated variants – proportions (e.g. BE, BEINF, BEZI)
  else if (fam %in% c("BE", "BEINF", "BEO", "BEZI", "BEo", "BEINF0")) {
    y <- (y - min(y)) / (max(y) - min(y) + eps)

    if (strategy == "safe") {
      # For inflated families, allow true 0 or 1 if model permits it
      allow_zero <- fam %in% c("BEINF", "BEZI", "BEINF0")
      allow_one  <- fam %in% c("BEINF", "BEo", "BEINF0")

      if (!allow_zero) y[y <= 0] <- eps
      if (!allow_one)  y[y >= 1] <- 1 - eps

    } else {
      y[y <= 0] <- NA
      y[y >= 1] <- NA
    }
    return(y)
  }

  # Batch D: ℝ – real-valued (e.g. NO, TF, GU)
  else if (fam %in% c("NO", "TF", "GU")) {
    return(scale(y))
  }

  # If no match, return as-is
  return(y)
}
