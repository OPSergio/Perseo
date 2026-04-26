#' Compute information criterion penalty term
#'
#' @param n_obs Integer number of valid observations used in the model fit.
#' @param criterion Character scalar: "AIC", "BIC", or "GAIC".
#' @param gaic_k Optional numeric penalty override for GAIC. If NULL and
#'   criterion is GAIC, defaults to log(n_obs).
#'
#' @return Numeric scalar penalty multiplier for model degrees of freedom.
#' @keywords internal
compute_ic_penalty <- function(criterion, gaic_k = NULL, n_obs) {
  stopifnot(length(n_obs) == 1L, n_obs > 0)
  stopifnot(length(criterion) == 1L)
  
  switch(
    criterion,
    "AIC" = 2,
    "BIC" = log(n_obs),
    "GAIC" = if (is.null(gaic_k)) log(n_obs) else gaic_k,
    stop("Unsupported criterion: ", criterion)
  )
}


#' Instantiate a GAMLSS family object by name
#'
#' Safely constructs a GAMLSS family object from its name, searching first in
#' the `gamlss.dist` namespace and then in the global environment.
#'
#' @param family_name Character scalar with the GAMLSS family name (e.g., "NBI").
#' @param bd_vec Optional numeric vector of binomial denominators (for BI/BB).
#'
#' @return A GAMLSS family object or NULL if the constructor cannot be found or
#'   fails to initialize.
#' @keywords internal
instantiate_gamlss_family <- function(family_name, bd_vec = NULL) {
  # Locate family constructor
  fam_constructor <- NULL
  if (exists(family_name, where = asNamespace("gamlss.dist"), inherits = FALSE)) {
    fam_constructor <- get(family_name, envir = asNamespace("gamlss.dist"))
  } else if (exists(family_name, mode = "function")) {
    fam_constructor <- get(family_name)
  }
  
  if (is.null(fam_constructor)) {
    return(NULL)
  }

  # Call constructor with optional binomial denominator
  tryCatch({
    if (!is.null(bd_vec)) {
      fam_constructor(bd = bd_vec)
    } else {
      fam_constructor()
    }
  }, error = function(e) NULL)
}


#' Fit a GAMLSS model while suppressing ALL console output
#'
#' Wraps `gamlss::gamlss()` to completely suppress all printed output (messages,
#' warnings, and cat() calls). Uses capture.output() to catch any remaining
#' output that suppressWarnings/Messages don't catch. Returns NULL on error for
#' safe iteration over multiple families/features.
#'
#' @param formula Model formula passed to `gamlss()`.
#' @param data Data frame containing the response and covariates.
#' @param family_obj Instantiated GAMLSS family object.
#'
#' @return A fitted `gamlss` object or NULL if fitting fails.
#' @keywords internal
fit_gamlss_safely <- function(formula, data, family_obj) {
  tryCatch({
    # Store family_obj BEFORE calling gamlss so it's in scope
    family_obj_stored <- family_obj
    
    # Suppress all GAMLSS output using capture.output
    fit_result <- NULL
    dummy <- capture.output({
      fit_result <- suppressWarnings(suppressMessages(
        gamlss::gamlss(
          formula = formula,
          data = data,
          family = family_obj,
          trace = FALSE,
          control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.001),
          i.control = gamlss::glim.control(trace = FALSE)
        )
      ))
    })
    
    # Store family_obj in fit for vcov() to find it
    if (!is.null(fit_result)) {
      fit_result$family_obj_stored <- family_obj_stored
    }
    
    fit_result
  }, error = function(e) {
    NULL
  })
}


#' Extract log-likelihood from a fitted GAMLSS model
#'
#' @param fit A fitted `gamlss` object.
#'
#' @return Numeric scalar or NA_real_ if extraction fails.
#' @keywords internal
extract_gamlss_loglik <- function(fit) {
  tryCatch(
    as.numeric(stats::logLik(fit)),
    error = function(e) NA_real_
  )
}


#' Compute Jacobian-corrected information criterion
#'
#' Calculates the information criterion with adjustment for the log-Jacobian
#' of the transformation applied to the response. This ensures IC values are
#' comparable across families on the original data scale.
#'
#' @param fit Fitted `gamlss` object.
#' @param penalty Numeric penalty term from `compute_ic_penalty()`.
#' @param logJ_sum Numeric scalar with the sum of log-Jacobian terms over
#'   valid observations.
#'
#' @return Numeric scalar IC value or Inf if logLik extraction fails.
#' @keywords internal
compute_jacobian_corrected_ic <- function(fit, penalty, logJ_sum) {
  log_lik <- extract_gamlss_loglik(fit)
  
  if (!is.finite(log_lik)) {
    return(Inf)
  }

  df_fit <- fit$df.fit
  if (is.null(df_fit) || length(df_fit) != 1L || !is.finite(df_fit)) {
    return(Inf)
  }
  base_ic <- (-2 * log_lik) + (penalty * df_fit)
  base_ic - (2 * logJ_sum)
}


#' Extract coefficient table from GAMLSS fit for mu parameter
#'
#' Robustly extracts coefficient estimates, standard errors, test statistics,
#' and p-values from `summary(fit, what = "mu")` output. Handles various
#' return structures from different GAMLSS versions.
#'
#' @param fit A fitted `gamlss` object.
#'
#' @return Tibble with columns: term, effect, se, stat, pval.
#' @keywords internal
extract_mu_coefficients <- function(fit) {
  # Get summary output while suppressing output
  summary_obj <- NULL
  
  dummy <- capture.output({
    summary_obj <- suppressWarnings(suppressMessages({
      summary(fit, what = "mu")
    }))
  })

  # Helper to coerce various table structures
  coerce_summary_table <- function(tab) {
    tab <- as.data.frame(tab)
    
    # Normalize column names
    names(tab) <- sub("Std\\.Error", "Std. Error", names(tab))
    
    # Ensure p-value column has standard name
    if (!"Pr(>|t|)" %in% names(tab)) {
      p_col_idx <- which(tolower(names(tab)) %in% c("p-value", "p.value", "pr(>|t|)"))
      if (length(p_col_idx) == 1) {
        names(tab)[p_col_idx] <- "Pr(>|t|)"
      }
    }
    
    # Fill missing columns with NA
    if (!"Estimate" %in% names(tab)) tab[["Estimate"]] <- NA_real_
    if (!"Std. Error" %in% names(tab)) tab[["Std. Error"]] <- NA_real_
    if (!"t value" %in% names(tab)) {
      tab[["t value"]] <- tab[["Estimate"]] / tab[["Std. Error"]]
    }
    if (!"Pr(>|t|)" %in% names(tab)) tab[["Pr(>|t|)"]] <- NA_real_
    
    # Extract term names
    term_names <- rownames(tab)
    if (is.null(term_names)) {
      if ("term" %in% names(tab)) {
        term_names <- as.character(tab$term)
      } else {
        term_names <- names(coef(fit, what = "mu"))
      }
    }
    
    # Clean up term names: remove leading X. and fix mangled intercepts
    # Pattern examples: "X.Intercept.", "X.Intercept..1", "X.Intercept..2"
    # Should all become "(Intercept)"
    term_names <- gsub("^X\\.", "", term_names)  # Remove "X." prefix
    term_names <- gsub("\\.+$", "", term_names)  # Remove trailing dots
    term_names <- gsub("\\.\\d+$", "", term_names)  # Remove .1, .2, etc at end
    term_names <- gsub("\\.+", ".", term_names)  # Collapse multiple dots to single dot
    term_names <- gsub("^\\.|\\.$", "", term_names)  # Remove leading/trailing dots again
    term_names <- ifelse(grepl("^Intercept", term_names, ignore.case = TRUE), "(Intercept)", term_names)
    
    tibble::tibble(
      term = term_names,
      effect = as.numeric(tab[["Estimate"]]),
      se = as.numeric(tab[["Std. Error"]]),
      stat = as.numeric(tab[["t value"]]),
      pval = as.numeric(tab[["Pr(>|t|)"]])
    )
  }

  # Try data.frame/matrix format
  if (is.data.frame(summary_obj) || is.matrix(summary_obj)) {
    return(coerce_summary_table(summary_obj))
  }

  # Try list format with nested coefficient table
  if (is.list(summary_obj)) {
    table_candidates <- list(
      summary_obj$coef.table,
      summary_obj$coefficients,
      summary_obj$coef,
      summary_obj$coeftable
    )
    
    for (candidate_table in table_candidates) {
      if (!is.null(candidate_table)) {
        return(coerce_summary_table(candidate_table))
      }
    }
  }

  # Fallback: extract coefficients only (no SE/p-values)
  est <- coef(fit, what = "mu")
  tibble::tibble(
    term = names(est),
    effect = as.numeric(est),
    se = NA_real_,
    stat = NA_real_,
    pval = NA_real_
  )
}
