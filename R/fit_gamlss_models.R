#' Fit GAMLSS models using multiple families for a single feature
#'
#' Fits GAMLSS models for a given feature across multiple families.
#' Provides feedback using cli and handles timeouts.
#'
#' @param y Numeric vector: expression values for a single feature (e.g., gene counts)
#' @param X Model matrix
#' @param families Character vector of GAMLSS family names
#' @param timeout Maximum time (seconds) to fit each model (default = 10)
#' @param verbose Logical, whether to print feedback messages (default = TRUE)
#'
#' @return A list of successfully fitted models with family name, AIC, and BIC
#' @export
fit_gamlss_models <- function(y, X, families = c("PO", "NBI", "NO", "GA"), timeout = 10, verbose = TRUE) {
  
  results <- list()
  
  if (verbose) cli::cli_h1("Fitting GAMLSS models for one feature")
  
  for (fam in families) {
    if (verbose) cli::cli_alert_info("Trying family {.strong {fam}}")
    
    start_time <- Sys.time()
    
    try_result <- tryCatch({
      R.utils::withTimeout({
        model <- gamlss::gamlss(y ~ X - 1, family = fam, trace = FALSE)
        list(
          family = fam,
          model = model,
          AIC = AIC(model),
          BIC = BIC(model),
          fit_time = round(as.numeric(Sys.time() - start_time), 2)
        )
      }, timeout = timeout, onTimeout = "silent")
    }, error = function(e) {
      if (verbose) cli::cli_alert_danger("Failed for {.strong {fam}}: {e$message}")
      NULL
    })
    
    if (!is.null(try_result)) {
      if (verbose) cli::cli_alert_success("Model with {.strong {fam}} fit in {try_result$fit_time}s [AIC={round(try_result$AIC, 2)}]")
      results[[fam]] <- try_result
    }
  }
  
  if (verbose && length(results) == 0) {
    cli::cli_alert_warning("No families converged for this feature.")
  }
  
  return(results)
}
