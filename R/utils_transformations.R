#' Transform expression values for compatibility with a GAMLSS family (SAFE)
#'
#' Applies family-specific preprocessing to `y` so that it meets the domain
#' required by each distribution. The "safe" version may **clip** or
#' **rescale** (e.g., replace out-of-domain values with valid limits).
#' Use this for visualization/diagnostics. For fair model-family comparison,
#' use \code{transform_for_family_strict()} with Jacobian correction.
#'
#' @param y Numeric vector of expression values.
#' @param fam Character: GAMLSS family name.
#' @param strategy "safe" (default) or "strict". "strict" only replaces
#'        invalid values with NA (for compatibility); comparative selection
#'        should use \code{transform_for_family_strict()}.
#' @param eps Small numeric value for smoothing/clipping.
#'
#' @return Numeric vector transformed (same length as `y`).
#' @export
transform_for_family <- function(y, fam, strategy = "safe", eps = 1e-6) {
  # A: (0, ∞) – strictly positive continuous (GA, GG, LOGNO, IG)
  if (fam %in% c("GA", "GG", "LOGNO", "IG")) {
    if (strategy == "safe") {
      y[y <= 0] <- eps
    } else {
      y[y <= 0] <- NA
    }
    return(y)
  }

  # B: [0, ∞) – counts (PO, NBI, ZIP, ZINBI, ZIP2)
  if (fam %in% c("PO", "NBI", "ZINBI", "ZIP", "ZIP2")) {
    y <- round(y)
    if (strategy == "safe") {
      y[y < 0] <- 0
    } else {
      y[y < 0] <- NA
    }
    return(y)
  }

  # C: (0,1) or inflated – proportions (BE, BEINF, BEO, BEZI, BEo, BEINF0)
  if (fam %in% c("BE", "BEINF", "BEO", "BEZI", "BEo", "BEINF0")) {
    rng <- max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
    y <- (y - min(y, na.rm = TRUE)) / (rng + eps)

    if (strategy == "safe") {
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

  # D: ℝ – real-valued (NO, TF, GU)
  if (fam %in% c("NO", "TF", "GU")) {
    return(as.numeric(scale(y))) # scale() returns a matrix; coerce to numeric.
  }

  # Default: return unchanged
  return(y)
}


#' Strict transform for model selection with Jacobian correction
#'
#' Recommended for *fair family comparison*: only enforces the theoretical
#' domain (no clipping/rounding), returns the working response `z = g(y)`,
#' a validity mask, and the per-observation log-Jacobian \eqn{\log|g'(y)|}.
#' Also includes metadata to invert the transformation.
#'
#' Families:
#' - Counts: \code{PO, NBI, ZIP, ZINBI, ZIP2, BI, BB} → identity (valid if integer ≥ 0)
#' - Positive: \code{GA, GG, LOGNO, IG} → identity (valid if y > 0)
#' - Unit interval and inflated: \code{BE, BEINF, BEO, BEZI, BEo, BEINF0}
#'     → min–max to (0,1), Jacobian \eqn{-\log(b-a)}, validity depends on inflation
#' - Real: \code{NO, TF, GU} → z-score standardization with sd > 0,
#'     Jacobian \eqn{-\log(sd)}.
#'
#' @param y Numeric vector (original response).
#' @param fam Character GAMLSS family name.
#' @param eps Small numeric for epsilon handling in domains that exclude 0/1.
#' @param allow_eps Logical; if TRUE, nudges boundary values slightly inside domain.
#'
#' @return list with:
#'   - y: numeric, response on the working scale z = g(y)
#'   - mask: logical, valid observations (no clipping/rounding)
#'   - logJ_per_obs: numeric, \eqn{\log|g'(y_i)|} per observation
#'   - meta: list(kind = "identity"/"zscore"/"minmax", params = list(...))
#' @export
transform_for_family_strict <- function(y, fam, eps = 1e-6, allow_eps = TRUE) {
  n <- length(y)
  finite_y <- is.finite(y)

  fam_count    <- c("PO","NBI","ZIP","ZINBI","ZIP2","BI","BB")
  fam_unit     <- c("BE","BEINF","BEO","BEZI","BEo","BEINF0")
  fam_positive <- c("GA","GG","LOGNO","IG")
  fam_real     <- c("NO","TF","GU")

  # Counts: identity; valid if integer >= 0
  if (fam %in% fam_count) {
    y_adj <- y
    if (allow_eps) y_adj[y_adj < 0] <- NA
    mask <- finite_y & (y_adj >= 0) & (abs(y_adj - round(y_adj)) < 1e-8)
    z <- y_adj
    logJ <- rep(0, n)
    meta <- list(kind = "identity", params = list())
    return(list(y = z, mask = mask, logJ_per_obs = logJ, meta = meta))
  }

  # Positive continuous: identity; valid if y > 0
  if (fam %in% fam_positive) {
    y_adj <- y
    if (allow_eps) y_adj[y_adj <= 0] <- eps
    mask <- finite_y & (y_adj > 0)
    z <- y_adj
    logJ <- rep(0, n)
    meta <- list(kind = "identity", params = list())
    return(list(y = z, mask = mask, logJ_per_obs = logJ, meta = meta))
  }

  # Unit interval (and inflated variants): min–max if range > 0
  if (fam %in% fam_unit) {
    yy <- y[finite_y]
    a <- suppressWarnings(min(yy, na.rm = TRUE))
    b <- suppressWarnings(max(yy, na.rm = TRUE))
    if (!is.finite(a) || !is.finite(b) || b <= a) {
      return(list(
        y = rep(NA_real_, n),
        mask = rep(FALSE, n),
        logJ_per_obs = rep(-Inf, n),
        meta = list(kind = "minmax", params = list(min = NA_real_, max = NA_real_))
      ))
    }

    z_all <- (y - a) / (b - a)

    allow_zero <- fam %in% c("BEINF", "BEZI", "BEINF0")
    allow_one  <- fam %in% c("BEINF", "BEo", "BEINF0")

    if (allow_eps) {
      if (!allow_zero) z_all[z_all <= 0] <- eps
      if (!allow_one)  z_all[z_all >= 1] <- 1 - eps
    }

    mask <- is.finite(z_all)
    mask <- mask & if (allow_zero) z_all >= 0 else z_all > 0
    mask <- mask & if (allow_one)  z_all <= 1 else z_all < 1

    logJ <- rep(-log(b - a), n)
    meta <- list(kind = "minmax", params = list(min = a, max = b))
    return(list(y = z_all, mask = mask, logJ_per_obs = logJ, meta = meta))
  }

  # Real-valued: z-score standardization
  if (fam %in% fam_real) {
    yy <- y[finite_y]
    s <- sd(yy, na.rm = TRUE)
    m <- mean(yy, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) {
      return(list(
        y = rep(NA_real_, n),
        mask = rep(FALSE, n),
        logJ_per_obs = rep(-Inf, n),
        meta = list(kind = "zscore", params = list(center = NA_real_, scale = NA_real_))
      ))
    }
    z <- (y - m) / s
    mask <- is.finite(z)
    logJ <- rep(-log(s), n)
    meta <- list(kind = "zscore", params = list(center = m, scale = s))
    return(list(y = as.numeric(z), mask = mask, logJ_per_obs = logJ, meta = meta))
  }

  # Default: identity
  mask <- finite_y
  meta <- list(kind = "identity", params = list())
  list(y = y, mask = mask, logJ_per_obs = rep(0, n), meta = meta)
}


#' Inverse of the strict transform (back to original Y-scale)
#'
#' Inverse mapping from the working scale (`z`) back to the original `y`.
#'
#' @param z Numeric vector on the working scale (output of transform_for_family_strict$y).
#' @param meta List(kind, params) as returned by transform_for_family_strict().
#'
#' @return Numeric vector on the original Y-scale.
#' @export
inverse_transform <- function(z, meta) {
  kind <- tryCatch(meta$kind, error = function(e) "identity")
  if (is.null(kind)) kind <- "identity"

  if (kind == "zscore") {
    m <- meta$params$center
    s <- meta$params$scale
    return(m + s * z)
  } else if (kind == "minmax") {
    a <- meta$params$min
    b <- meta$params$max
    return(a + z * (b - a))
  } else {
    # identity or unknown kind
    return(z)
  }
}


#' Family groups by theoretical support
#'
#' @return List with character vectors of families by support: count, unit, positive, real.
#' @export
family_groups <- function() {
  list(
    count    = c("PO","NBI","ZIP","ZINBI","ZIP2","BI","BB"),
    unit     = c("BE","BEINF","BEO","BEZI","BEo","BEINF0"),
    positive = c("GA","GG","LOGNO","IG"),
    real     = c("NO","TF","GU")
  )
}


#' Infer empirical support of a response vector
#'
#' @param y Numeric vector.
#' @return One of "count", "unit", "positive", "real", or "none".
#' @export
infer_support <- function(y) {
  finite_y <- is.finite(y)
  if (!any(finite_y)) return("none")
  yy <- y[finite_y]

  is_count <- all(abs(yy - round(yy)) < 1e-8) && min(yy, na.rm = TRUE) >= 0
  if (is_count) return("count")

  if (all(yy >= 0, na.rm = TRUE) && max(yy, na.rm = TRUE) <= 1) return("unit")
  if (all(yy > 0, na.rm = TRUE)) return("positive")
  "real"
}


#' Sum of log-Jacobian over a mask
#'
#' @param logJ_per_obs Numeric vector from transform_for_family_strict().
#' @param mask Logical vector; if NULL, sums over finite entries.
#' @return Numeric scalar with the sum.
#' @export
jacobian_sum <- function(logJ_per_obs, mask = NULL) {
  if (is.null(mask)) {
    return(sum(logJ_per_obs[is.finite(logJ_per_obs)]))
  }
  sum(logJ_per_obs[mask & is.finite(logJ_per_obs)])
}
