#' Find Best-Fitting GAMLSS Families Across Genes
#'
#' This function randomly samples a number of genes and fits GAMLSS models with different families
#' to identify which distribution best fits each gene, based on AIC. It reports the top families across all tested genes.
#'
#' @param counts_matrix A numeric matrix of expression/abundance values (genes x samples)
#' @param n_genes Integer: number of genes to sample for fitting (default: 200)
#' @param top_n Integer: number of top frequent families to return (default: 4)
#' @param families Character vector of GAMLSS family names to test. If NULL, a default set is used.
#' @param verbose Logical: whether to print progress messages (default: TRUE)
#' @param strategy Character: transformation strategy ("safe" recommended). Not used yet in base logic.
#'
#' @return A character vector of the top `top_n` families, ranked by how often they best fit the sampled genes.
#'
#' @examples
#' set.seed(42)
#' top_families <- find_families(counts_matrix = my_counts, n_genes = 100)
#'
#' @export

source("R/utils_transformations.R") # Load utility functions (temporal)
find_families <- function(counts_matrix, n_genes = 200, top_n = 4, families = NULL, verbose = TRUE, strategy = "safe", eps = 1e-6) {
  if (is.null(families)) {
    families <- c(
      "BE", "BEO", "BEINF", "BEZI",
      "BI", "BB", "NBI",
      "PO", "ZIP", "ZINBI",
      "GA", "GG", "IG",
      "NO", "TF", "GU"
    )
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
    if (all(y == 0)) {
      return(list(gene = gene_id, family = NA, skipped = TRUE))
    }

    best_fit <- NA
    best_aic <- Inf

    for (fam in families) {
      y_t <- transform_for_family(y, fam, strategy = strategy, eps = eps)
      if (any(is.na(y_t))) next

      aic_val <- tryCatch(
        {
          suppressWarnings({
            fit <- gamlss(y_t ~ 1, family = fam, trace = FALSE)
            AIC(fit)
          })
        },
        error = function(e) NA
      )

      if (!is.na(aic_val) && aic_val < best_aic) {
        best_fit <- fam
        best_aic <- aic_val
      }
    }

    list(gene = gene_id, family = best_fit, skipped = FALSE, y_t = y_t)
  })

  plan(sequential)

  results_df <- tibble::tibble(
    gene = map_chr(results, "gene"),
    family = map_chr(results, "family"),
    skipped = map_lgl(results, "skipped"),
    y_t = map(results, "y_t")  
  )

  fitted <- results_df[!is.na(results_df$family), ]
  skipped <- sum(results_df$skipped)
  
  fam_freq <- sort(table(fitted$family), decreasing = TRUE)
  top_families <- names(head(fam_freq, top_n))

  if (verbose) {
    message("\n ===== Summary Report =====")
    message("Genes analyzed: ", total_genes)
    message("Genes skipped (all 0s): ", skipped)
    message("Genes successfully fitted: ", length(fitted))
    if (length(fam_freq) > 0) {
      message(
        "Most frequent family: ", names(fam_freq)[1],
        " (", fam_freq[1], " genes)"
      )
    } else {
      message("️ No successful fits.")
    }
    message(" Top ", top_n, " families:")
    print(fam_freq[1:min(top_n, length(fam_freq))])
  }

  return(list(
    top_families = top_families,
    results_df = results_df
  ))
}
