
transform_for_family <- function(y, fam) {
  if (fam %in% c("GA", "GG", "LOGNO", "IG")) {
    y[y <= 0] <- NA
    return(y)
  } else if (fam %in% c("PO", "NBI", "ZINBI", "ZIP", "ZIP2")) {
    y[y < 0] <- NA
    return(round(y))
  } else if (fam %in% c("BE", "BEINF", "BEO", "BEZI", "BEo", "BEINF0")) {
    y <- (y - min(y)) / (max(y) - min(y) + 1e-8)
    y[y <= 0] <- NA
    y[y >= 1] <- NA
    return(y)
  } else if (fam %in% c("NO", "TF", "GU")) {
    return(scale(y))
  } else {
    return(y)
  }
}

find_families <- function(counts_matrix, n_genes = 200, top_n = 4, families = NULL, verbose = TRUE) {
  if (is.null(families)) {
    families <- c("BE", "BEO", "BEINF", "BEZI",
                  "BI", "BB", "NBI",
                  "PO", "ZIP", "ZINBI",
                  "GA", "GG", "IG",
                  "NO", "TF", "GU")
  }
  
  if (n_genes > nrow(counts_matrix)) {
    stop("n_genes exceeds the number of rows in counts_matrix.")
  }
  
  selected_genes <- sample(rownames(counts_matrix), n_genes)
  total_genes <- length(selected_genes)
  
  if (verbose) message(" Running model fitting on ", total_genes, " genes using parallel execution...")
  
  plan(multisession)
  
  results <- future_lapply(selected_genes, function(gene_id) {
    y <- as.numeric(counts_matrix[gene_id, ])
    if (all(y == 0)) return(list(gene = gene_id, family = NA, skipped = TRUE))
    
    best_fit <- NA
    best_aic <- Inf
    
    for (fam in families) {
      y_t <- transform_for_family(y, fam)
      if (any(is.na(y_t))) next
      
      aic_val <- tryCatch({
        suppressWarnings({
          fit <- gamlss(y_t ~ 1, family = fam, trace = FALSE)
          AIC(fit)
        })
      }, error = function(e) NA)
      
      if (!is.na(aic_val) && aic_val < best_aic) {
        best_fit <- fam
        best_aic <- aic_val
      }
    }
    
    list(gene = gene_id, family = best_fit, skipped = FALSE)
  })
  
  plan(sequential)
  
  results_df <- bind_rows(lapply(results, as.data.frame))
  
  fitted <- na.omit(results_df$family)
  skipped <- sum(results_df$skipped)
  
  fam_freq <- sort(table(fitted), decreasing = TRUE)
  top_families <- names(head(fam_freq, top_n))
  
  if (verbose) {
    message("\n ===== Summary Report =====")
    message("Genes analyzed: ", total_genes)
    message("Genes skipped (all 0s): ", skipped)
    message("Genes successfully fitted: ", length(fitted))
    if (length(fam_freq) > 0) {
      message("Most frequent family: ", names(fam_freq)[1],
              " (", fam_freq[1], " genes)")
    } else {
      message("️ No successful fits.")
    }
    message(" Top ", top_n, " families:")
    print(fam_freq[1:min(top_n, length(fam_freq))])
  }
  
  return(top_families)
}
