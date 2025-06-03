#' Extract contrast statistics from a fitted GAMLSS model
#'
#' Computes the estimate, standard error, z-score and p-value for a specific coefficient
#'
#' @param model A fitted GAMLSS model
#' @param coef_name Character: name of the coefficient to test (e.g. "groupTumor")
#' @param verbose Logical: print messages if TRUE
#'
#' @return A named list with: estimate, std_error, z_score, pval
#' @export
contrast_test <- function(model, coef_name) {

  est <- model$mu.coefficients[coef_name]
  
  if (is.null(est)) {
    message("Coeficiente no encontrado: ", coef_name)
    return(NULL)
  }

  vcov_mat <- tryCatch(model$mu.qr$qr, error = function(e) NULL)
  
  if (is.null(vcov_mat)) {
    message("No se pudo extraer la matriz QR")
    return(NULL)
  }

  XtX_inv <- tryCatch({
    chol2inv(qr.R(model$mu.qr))
  }, error = function(e) NULL)
  
  if (is.null(XtX_inv)) {
    message("Fallo al invertir la matriz para el error estándar")
    return(NULL)
  }

  coef_names <- names(model$mu.coefficients)
  idx <- match(coef_name, coef_names)
  
  if (is.na(idx)) {
    message("Nombre de coeficiente no encontrado en QR")
    return(NULL)
  }
  
  se <- sqrt(diag(XtX_inv))[idx]
  z <- est / se
  p <- 2 * pnorm(-abs(z))
  
  return(data.frame(
    term = coef_name,
    estimate = est,
    std_error = se,
    z_value = z,
    p_value = p
  ))
}

