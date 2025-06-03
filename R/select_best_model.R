#' Select best model from fitted GAMLSS models
#'
#' From a list of fitted models (e.g., output of fit_gamlss_models), selects the best one by AIC or BIC.
#'
#' @param models List of fitted models (as returned by fit_gamlss_models)
#' @param criterion Character: criterion to use for selection ("AIC" or "BIC")
#'
#' @return A list with elements: $best_model, $best_family, $criterion_value
#' @export
select_best_model <- function(models, criterion = c("AIC", "BIC")) {
  criterion <- match.arg(criterion)
  
  if (length(models) == 0) {
    warning("No models to select from.")
    return(NULL)
  }
  
  scores <- purrr::map_dfr(models, function(m) {
    tibble(family = m$family, score = m[[criterion]])
  })
  
  best <- scores %>% arrange(score) %>% slice(1)
  
  return(list(
    best_model = models[[best$family]],
    best_family = best$family,
    criterion_value = best$score
  ))
}
