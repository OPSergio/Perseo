differential_expression_gamlss_optimized <- function(counts, design,
                                                     families = c("NBI", "GA", "NO", "PO"),
                                                     coef_name,
                                                     criterion = "AIC",
                                                     p_adjust_method = "BH",
                                                     timeout = 10,
                                                     verbose = TRUE,
                                                     workers = 4) {
  load_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  load_package("future.apply")
  load_package("progressr")
  load_package("gamlss")
  
  future::plan(future::multisession, workers = workers)

  genes <- rownames(counts)
  n_genes <- length(genes)

  process_gene <- function(gene, i) {
    y <- as.numeric(counts[gene, ])
    best_model <- NULL
    best_aic <- Inf
    best_family <- NULL
    
    for (fam in families) {
      fit <- tryCatch(
        suppressWarnings(
          gamlss::gamlss(y ~ design - 1, family = fam, trace = FALSE, control = gamlss.control(n.cyc = 20))
        ),
        error = function(e) NULL
      )
      if (!is.null(fit) && isTRUE(fit$converged)) {
        aic <- AIC(fit)
        if (aic < best_aic) {
          best_model <- fit
          best_aic <- aic
          best_family <- fam
        }
      }
    }
    
    if (is.null(best_model)) return(NULL)

    coefs <- coef(best_model)
    vcov_mat <- tryCatch(vcov(best_model, what = "mu"), error = function(e) NULL)
    
    if (is.null(vcov_mat) || !(coef_name %in% names(coefs))) return(NULL)
    
    est <- coefs[coef_name]
    se <- sqrt(vcov_mat[coef_name, coef_name])
    z <- est / se
    p <- 2 * pnorm(-abs(z))
    
    return(data.frame(
      gene = gene,
      family = best_family,
      logFC = est,
      std_error = se,
      z_value = z,
      p_value = p,
      stringsAsFactors = FALSE
    ))
  }

  progressr::handlers(global = TRUE)
  progressr::with_progress({
    p <- progressr::progressor(steps = n_genes)
    
    results <- future.apply::future_lapply(seq_along(genes), function(i) {
      gene <- genes[i]
      if (verbose && i %% 50 == 0) message(sprintf("Procesando gen %d/%d: %s", i, n_genes, gene))
      p()
      process_gene(gene, i)
    }, future.seed = TRUE)
  })

  results_df <- do.call(rbind, results)
  if (!is.null(results_df) && nrow(results_df) > 0) {
    results_df$p_adj <- p.adjust(results_df$p_value, method = p_adjust_method)
  }
  
  return(results_df)
}
