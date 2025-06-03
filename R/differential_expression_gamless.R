differential_expression_gamlss <- function(counts, design,
                                           families = c("NBI", "GA", "NO", "PO"),
                                           coef_name,
                                           criterion = "AIC",
                                           p_adjust_method = "BH",
                                           timeout = 10,
                                           verbose = TRUE,
                                           workers = 4) {

  future::plan(future::multisession, workers = workers)
  
  genes <- rownames(counts)
  n_genes <- length(genes)

  process_gene <- function(gene, i) {
    y <- as.numeric(counts[gene, ])
    
    fit <- fit_gamlss_models(y = y, X = design, families = families, timeout = timeout, verbose = FALSE)
    best <- select_best_model(fit, criterion = criterion)
    if (is.null(best)) return(NULL)
    
    model <- best$best_model$model
    fam <- best$best_family
    
    res <- contrast_test(model, coef_name)
    if (is.null(res)) return(NULL)
    
    data.frame(
      gene = gene,
      family = fam,
      logFC = res$estimate,
      std_error = res$std_error,
      z_value = res$z_value,
      p_value = res$p_value,
      stringsAsFactors = FALSE
    )
  }

  handlers(global = TRUE)
  progressr::with_progress({
    p <- progressor(steps = n_genes)
    
    results <- future_lapply(seq_along(genes), function(i) {
      p(message = sprintf("Procesando gen %d/%d: %s", i, n_genes, genes[i]))
      process_gene(genes[i], i)
    }, future.seed = TRUE)
  })

  results_df <- do.call(rbind, results)
  if (!is.null(results_df) && nrow(results_df) > 0) {
    results_df$p_adj <- p.adjust(results_df$p_value, method = p_adjust_method)
  }
  
  return(results_df)
}
