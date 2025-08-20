#' Transform expression values for compatibility with GAMLSS family (SAFE)
#'
#' Aplica preprocesados específicos por familia a `y` para que cumpla el dominio
#' requerido por cada distribución. La versión "safe" puede **recortar** o
#' **re-escalar** (p. ej., sustituir valores fuera de dominio por límites válidos).
#' Úsala para visualización/diagnóstico. Para la selección justa de familias
#' emplea \code{transform_for_family_strict()} con corrección jacobiana.
#'
#' @param y Numeric vector de valores de expresión.
#' @param fam Character: nombre de la familia GAMLSS.
#' @param strategy "safe" (default) o "strict". "strict" aquí sólo reemplaza
#'        valores inválidos por NA (para compatibilidad); la selección real usa
#'        \code{transform_for_family_strict()}.
#' @param eps Numeric pequeño para suavizado/clipping.
#'
#' @return Numeric vector transformado (misma longitud que `y`).
#' @export
transform_for_family <- function(y, fam, strategy = "safe", eps = 1e-6) {
  # A: (0, ∞) – continuas estrictamente positivas (GA, GG, LOGNO, IG)
  if (fam %in% c("GA", "GG", "LOGNO", "IG")) {
    if (strategy == "safe") {
      y[y <= 0] <- eps
    } else {
      y[y <= 0] <- NA
    }
    return(y)
  }
  
  # B: [0, ∞) – conteos (PO, NBI, ZIP, ZINBI, ZIP2)
  if (fam %in% c("PO", "NBI", "ZINBI", "ZIP", "ZIP2")) {
    y <- round(y)
    if (strategy == "safe") {
      y[y < 0] <- 0
    } else {
      y[y < 0] <- NA
    }
    return(y)
  }
  
  # C: (0,1) o infladas – proporciones (BE, BEINF, BEO, BEZI, BEo, BEINF0)
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
  
  # D: ℝ – reales (NO, TF, GU)
  if (fam %in% c("NO", "TF", "GU")) {
    return(as.numeric(scale(y))) # scale() devuelve matriz; lo forzamos a num.
  }
  
  # Si no coincide, devolver tal cual
  return(y)
}


#' Strict transform for model selection with Jacobian correction
#'
#' Uso recomendado para *comparar familias de forma justa*: sólo impone
#' el dominio teórico (sin recortes/rounding), devuelve la respuesta de trabajo
#' `z = g(y)`, una máscara de validez, y el log-Jacobiano por observación
#' \eqn{\log|g'(y)|}. Además, incluye metadatos para invertir la transformación.
#'
#' Familias:
#' - Conteos: \code{PO, NBI, ZIP, ZINBI, ZIP2, BI, BB} → identidad (válido si entero ≥ 0)
#' - Positivas: \code{GA, GG, LOGNO, IG} → identidad (válido si y > 0)
#' - (0,1) y variantes infladas: \code{BE, BEINF, BEO, BEZI, BEo, BEINF0}
#'     → min–max a (0,1), Jacobiano \eqn{-\log(b-a)}, validez según inflado.
#' - Reales: \code{NO, TF, GU} → tipificación (z-score) con sd>0,
#'     Jacobiano \eqn{-\log(sd)}.
#'
#' @param y Numeric vector (respuesta original).
#' @param fam Character nombre de familia GAMLSS.
#'
#' @return list con:
#'   - y: numeric, respuesta en escala de trabajo z = g(y)
#'   - mask: logical, observaciones válidas (sin clipping/rounding)
#'   - logJ_per_obs: numeric, \eqn{\log|g'(y_i)|} por observación
#'   - meta: list(kind = "identity"/"zscore"/"minmax", params = list(...))
#' @export
transform_for_family_strict <- function(y, fam, eps = 1e-6, allow_eps = TRUE) {
  n <- length(y)
  finite_y <- is.finite(y)
  
  fam_count    <- c("PO","NBI","ZIP","ZINBI","ZIP2","BI","BB")
  fam_unit     <- c("BE","BEINF","BEO","BEZI","BEo","BEINF0")
  fam_positive <- c("GA","GG","LOGNO","IG")
  fam_real     <- c("NO","TF","GU")
  
  # Conteos: identidad; válido si entero >= 0
  if (fam %in% fam_count) {
    y_adj <- y
    if (allow_eps) y_adj[y_adj < 0] <- NA
    mask <- finite_y & (y_adj >= 0) & (abs(y_adj - round(y_adj)) < 1e-8)
    z <- y_adj
    logJ <- rep(0, n)
    meta <- list(kind = "identity", params = list())
    return(list(y = z, mask = mask, logJ_per_obs = logJ, meta = meta))
  }
  
  # Positivas: identidad; válido si y > 0
  if (fam %in% fam_positive) {
    y_adj <- y
    if (allow_eps) y_adj[y_adj <= 0] <- eps
    mask <- finite_y & (y_adj > 0)
    z <- y_adj
    logJ <- rep(0, n)
    meta <- list(kind = "identity", params = list())
    return(list(y = z, mask = mask, logJ_per_obs = logJ, meta = meta))
  }
  
  # (0,1) o infladas: min–max si rango > 0
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
  
  # Reales: z-score
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
  
  # Por defecto: identidad
  mask <- finite_y
  meta <- list(kind = "identity", params = list())
  list(y = y, mask = mask, logJ_per_obs = rep(0, n), meta = meta)
}


#' Inverse of the strict transform (back to Y-scale)
#'
#' Inversa de la transformación estricta, útil para llevar simulaciones en
#' escala de trabajo (z) de vuelta a la escala original Y.
#'
#' @param z numeric vector en escala de trabajo (salida de transform_for_family_strict$y)
#' @param meta list(kind, params) tal como devuelve transform_for_family_strict()
#'
#' @return numeric vector en escala original Y
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
    # identity o tipo desconocido
    return(z)
  }
}


#' Family groups by theoretical support
#'
#' @return list con vectores de familias por soporte: count, unit, positive, real
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
#' @param y numeric vector
#' @return one of "count", "unit", "positive", "real", or "none"
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
#' @param logJ_per_obs numeric vector desde transform_for_family_strict()
#' @param mask logical vector; si NULL, suma sobre entradas finitas
#' @return escalar numérico
#' @export
jacobian_sum <- function(logJ_per_obs, mask = NULL) {
  if (is.null(mask)) {
    return(sum(logJ_per_obs[is.finite(logJ_per_obs)]))
  }
  sum(logJ_per_obs[mask & is.finite(logJ_per_obs)])
}
