test_that("fit_gamlss_models accepts formula string input", {
  skip_on_cran()
  
  set.seed(123)
  counts <- matrix(rpois(300, lambda = 10), nrow = 10, ncol = 30)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    condition = factor(rep(c("A", "B", "C"), each = 10)),
    batch = factor(rep(1:3, 10))
  )
  
  result <- fit_gamlss_models(
    counts,
    design_matrix = "~ condition + batch",
    metadata = metadata,
    candidate_families = c("PO", "NBI"),
    workers = 1,
    show_progress = FALSE
  )
  
  expect_s3_class(result$results, "tbl_df")
  expect_s3_class(result$selection, "tbl_df")
  expect_true(nrow(result$results) > 0)
  expect_true(all(c("feature", "term", "effect", "pval") %in% names(result$results)))
})

test_that("fit_gamlss_models auto-generates contrasts from contrast_variable", {
  skip_on_cran()
  
  set.seed(456)
  counts <- matrix(rpois(300, lambda = 10), nrow = 10, ncol = 30)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    tissue_type = factor(rep(c("liver", "brain", "heart"), each = 10))
  )
  
  design <- model.matrix(~ tissue_type, data = metadata)
  
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO", "NBI"),
    contrast_variable = "tissue_type",
    workers = 1,
    show_progress = FALSE
  )
  
  expect_true("contrasts" %in% names(result))
  expect_s3_class(result$contrasts, "tbl_df")
  expect_true(all(c("feature", "contrast", "estimate", "p_value", "p_adj") %in% names(result$contrasts)))
  
  # Should have 3 pairwise contrasts for 3 levels
  unique_contrasts <- unique(result$contrasts$contrast)
  expect_equal(length(unique_contrasts), 3)
})

test_that("fit_gamlss_models prioritizes contrast_matrix over contrast_variable", {
  skip_on_cran()
  
  set.seed(789)
  counts <- matrix(rpois(200, lambda = 10), nrow = 10, ncol = 20)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    status = factor(rep(c("A", "B"), each = 10))
  )
  
  design <- model.matrix(~ status, data = metadata)
  
  # Custom contrast matrix
  C <- matrix(c(0, 1), nrow = 1)
  colnames(C) <- c("(Intercept)", "statusB")
  rownames(C) <- "B_vs_A"
  
  expect_message(
    result <- fit_gamlss_models(
      counts,
      design_matrix = design,
      metadata = metadata,
      candidate_families = c("PO"),
      contrast_matrix = C,
      contrast_variable = "status",
      workers = 1,
      show_progress = FALSE
    ),
    "Prioritizing contrast_matrix"
  )
  
  expect_true("contrasts" %in% names(result))
  expect_equal(unique(result$contrasts$contrast), "B_vs_A")
})

test_that("fit_gamlss_models requires metadata for formula string", {
  counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  
  expect_error(
    fit_gamlss_models(
      counts,
      design_matrix = "~ condition",
      candidate_families = c("PO"),
      workers = 1,
      show_progress = FALSE
    ),
    "metadata must be provided"
  )
})

test_that("fit_gamlss_models parallelizes correctly", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if(!requireNamespace("PERSEO", quietly = TRUE), "PERSEO must be installed for parallel tests")
  
  set.seed(999)
  counts <- matrix(rpois(500, lambda = 10), nrow = 50, ncol = 10)
  rownames(counts) <- paste0("gene", 1:50)
  
  design <- data.frame(
    condition = factor(rep(c("A", "B"), each = 5))
  )
  
  # Test with workers = 1 (sequential)
  result_seq <- fit_gamlss_models(
    counts,
    design_matrix = model.matrix(~ condition, design),
    candidate_families = c("PO", "NBI"),
    workers = 1,
    show_progress = FALSE
  )
  
  # Test with workers > 1 (parallel)
  future::plan(future::multisession, workers = 2)
  result_par <- fit_gamlss_models(
    counts,
    design_matrix = model.matrix(~ condition, design),
    candidate_families = c("PO", "NBI"),
    workers = 2,
    show_progress = FALSE
  )
  future::plan(future::sequential)
  
  expect_equal(nrow(result_seq$results), nrow(result_par$results))
  expect_equal(nrow(result_seq$selection), nrow(result_par$selection))
})

test_that("fit_gamlss_models works with parallel = TRUE", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if(!requireNamespace("PERSEO", quietly = TRUE), "PERSEO must be installed for parallel tests")
  
  set.seed(999)
  counts <- matrix(rpois(500, lambda = 10), nrow = 50, ncol = 10)
  rownames(counts) <- paste0("gene", 1:50)
  
  design <- data.frame(
    condition = factor(rep(c("A", "B"), each = 5))
  )
  
  # Test with parallel = TRUE (automatic setup)
  result_par <- fit_gamlss_models(
    counts,
    design_matrix = model.matrix(~ condition, design),
    candidate_families = c("PO", "NBI"),
    parallel = TRUE,
    workers = 2,
    show_progress = FALSE
  )
  
  expect_s3_class(result_par$results, "tbl_df")
  expect_s3_class(result_par$selection, "tbl_df")
  expect_true(nrow(result_par$results) > 0)
})

test_that("fit_gamlss_models parallel and sequential give same results", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if(!requireNamespace("PERSEO", quietly = TRUE), "PERSEO must be installed for parallel tests")
  
  set.seed(123)
  counts <- matrix(rpois(200, lambda = 10), nrow = 20, ncol = 10)
  rownames(counts) <- paste0("gene", 1:20)
  
  design <- data.frame(
    condition = factor(rep(c("A", "B"), each = 5))
  )
  
  design_mat <- model.matrix(~ condition, design)
  
  # Sequential
  result_seq <- fit_gamlss_models(
    counts,
    design_matrix = design_mat,
    candidate_families = c("PO", "NBI"),
    parallel = FALSE,
    show_progress = FALSE
  )
  
  # Parallel
  result_par <- fit_gamlss_models(
    counts,
    design_matrix = design_mat,
    candidate_families = c("PO", "NBI"),
    parallel = TRUE,
    workers = 2,
    show_progress = FALSE
  )
  
  expect_equal(nrow(result_seq$results), nrow(result_par$results))
  expect_equal(nrow(result_seq$selection), nrow(result_par$selection))
  expect_equal(sort(result_seq$selection$feature), sort(result_par$selection$feature))
})

test_that("fit_gamlss_models cleans up future plan after execution", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if(!requireNamespace("PERSEO", quietly = TRUE), "PERSEO must be installed for parallel tests")
  
  # Ensure we start with sequential
  future::plan(future::sequential)
  initial_plan <- class(future::plan())[1]
  
  counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  design <- model.matrix(~ factor(rep(c("A", "B"), each = 5)))
  
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    candidate_families = c("PO"),
    parallel = TRUE,
    workers = 2,
    show_progress = FALSE
  )
  
  # Plan should be reset to original (sequential)
  final_plan <- class(future::plan())[1]
  expect_equal(initial_plan, final_plan)
})

test_that("fit_gamlss_models validates workers parameter", {
  counts <- matrix(rpois(50, lambda = 10), nrow = 5, ncol = 10)
  design <- model.matrix(~ factor(rep(c("A", "B"), each = 5)))
  
  expect_warning(
    fit_gamlss_models(
      counts,
      design_matrix = design,
      candidate_families = c("PO"),
      parallel = TRUE,
      workers = 0,  # Invalid
      show_progress = FALSE
    ),
    "workers must be >= 1"
  )
})

test_that("fit_gamlss_models with omnibus = FALSE behaves as before", {
  skip_on_cran()
  
  set.seed(100)
  counts <- matrix(rpois(300, lambda = 10), nrow = 10, ncol = 30)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    tissue = factor(rep(c("A", "B", "C"), each = 10))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  # Without omnibus
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO", "NBI"),
    contrast_variable = "tissue",
    omnibus = FALSE,
    show_progress = FALSE
  )
  
  expect_true("contrasts" %in% names(result))
  expect_false("omnibus" %in% names(result))
  
  # All features should have contrasts
  n_features_with_contrasts <- length(unique(result$contrasts$feature))
  expect_true(n_features_with_contrasts > 0)
})

test_that("fit_gamlss_models with omnibus = TRUE adds omnibus output", {
  skip_on_cran()
  
  set.seed(200)
  counts <- matrix(rpois(300, lambda = 10), nrow = 10, ncol = 30)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    tissue = factor(rep(c("A", "B", "C"), each = 10))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO", "NBI"),
    contrast_variable = "tissue",
    omnibus = TRUE,
    omnibus_threshold = 0.05,
    omnibus_test = "Wald",
    show_progress = FALSE
  )
  
  expect_true("omnibus" %in% names(result))
  expect_s3_class(result$omnibus, "tbl_df")
  
  expect_true(all(c("feature", "family", "test_type", "statistic", 
                    "df", "p_value", "pass") %in% names(result$omnibus)))
  
  expect_equal(unique(result$omnibus$test_type), "Wald")
  expect_true(all(result$omnibus$df == 2))  # 3 levels - 1
})

test_that("fit_gamlss_models omnibus filtering reduces contrast computations", {
  skip_on_cran()
  
  set.seed(300)
  counts <- matrix(rpois(500, lambda = 10), nrow = 20, ncol = 25)
  rownames(counts) <- paste0("gene", 1:20)
  
  metadata <- data.frame(
    tissue = factor(rep(c("A", "B", "C"), length.out = 25))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  # Without omnibus - all features get contrasts
  result_no_omnibus <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO"),
    contrast_variable = "tissue",
    omnibus = FALSE,
    show_progress = FALSE
  )
  
  # With omnibus - only significant features get contrasts
  result_with_omnibus <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO"),
    contrast_variable = "tissue",
    omnibus = TRUE,
    omnibus_threshold = 0.05,
    show_progress = FALSE
  )
  
  n_contrasts_no_omnibus <- length(unique(result_no_omnibus$contrasts$feature))
  n_contrasts_with_omnibus <- length(unique(result_with_omnibus$contrasts$feature))
  
  # With omnibus should have same or fewer features with contrasts
  expect_true(n_contrasts_with_omnibus <= n_contrasts_no_omnibus)
  
  # Only features that passed omnibus should have contrasts
  if (nrow(result_with_omnibus$contrasts) > 0) {
    features_with_contrasts <- unique(result_with_omnibus$contrasts$feature)
    features_passed_omnibus <- result_with_omnibus$omnibus$feature[
      result_with_omnibus$omnibus$pass
    ]
    
    expect_true(all(features_with_contrasts %in% features_passed_omnibus))
  }
})

test_that("fit_gamlss_models omnibus with LRT works", {
  skip_on_cran()
  
  set.seed(400)
  counts <- matrix(rpois(200, lambda = 10), nrow = 10, ncol = 20)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    tissue = factor(rep(c("A", "B"), each = 10))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO"),
    contrast_variable = "tissue",
    omnibus = TRUE,
    omnibus_test = "LRT",
    show_progress = FALSE
  )
  
  expect_true("omnibus" %in% names(result))
  expect_equal(unique(result$omnibus$test_type), "LRT")
  
  # Some features might not have omnibus results if they fail
  if (nrow(result$omnibus) > 0 && any(!is.na(result$omnibus$df))) {
    valid_df <- result$omnibus$df[!is.na(result$omnibus$df)]
    expect_true(all(valid_df == 1))  # 2 levels - 1
  }
})

test_that("fit_gamlss_models omnibus preserves contrast p-value adjustment", {
  skip_on_cran()
  
  set.seed(500)
  counts <- matrix(rpois(300, lambda = 15), nrow = 15, ncol = 20)
  rownames(counts) <- paste0("gene", 1:15)
  
  metadata <- data.frame(
    tissue = factor(rep(c("A", "B", "C"), length.out = 20))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO", "NBI"),
    contrast_variable = "tissue",
    omnibus = TRUE,
    omnibus_threshold = 1.0,  # Allow all to pass
    show_progress = FALSE
  )
  
  if (nrow(result$contrasts) > 0) {
    # Check that p_adj exists
    expect_true("p_adj" %in% names(result$contrasts))
    
    # Check that adjustment is by contrast (not by feature)
    # For each contrast, p_adj should be computed across features
    for (contrast_name in unique(result$contrasts$contrast)) {
      contrast_subset <- result$contrasts[result$contrasts$contrast == contrast_name, ]
      
      # p_adj should differ from p_value (unless all p-values are identical)
      if (length(unique(contrast_subset$p_value[!is.na(contrast_subset$p_value)])) > 1) {
        expect_false(all(contrast_subset$p_adj == contrast_subset$p_value, na.rm = TRUE))
      }
    }
  }
})

test_that("fit_gamlss_models omnibus handles edge cases", {
  skip_on_cran()
  
  set.seed(600)
  counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    tissue = factor(rep(c("A", "B"), each = 5))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  # Very strict threshold - no features should pass
  result_strict <- fit_gamlss_models(
    counts,
    design_matrix = design,
    metadata = metadata,
    candidate_families = c("PO"),
    contrast_variable = "tissue",
    omnibus = TRUE,
    omnibus_threshold = 0.0001,
    show_progress = FALSE
  )
  
  expect_true("omnibus" %in% names(result_strict))
  
  # Contrasts might be empty or only have NA values
  if (nrow(result_strict$contrasts) > 0) {
    features_with_contrasts <- unique(result_strict$contrasts$feature)
    expect_true(length(features_with_contrasts) <= nrow(counts))
  }
})

test_that("fit_gamlss_models omnibus with custom contrast_matrix", {
  skip_on_cran()
  
  set.seed(700)
  counts <- matrix(rpois(200, lambda = 10), nrow = 10, ncol = 20)
  rownames(counts) <- paste0("gene", 1:10)
  
  metadata <- data.frame(
    tissue = factor(rep(c("A", "B", "C"), length.out = 20))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  # Custom contrast: B vs C only
  C <- matrix(c(0, 1, -1), nrow = 1)
  colnames(C) <- c("(Intercept)", "tissueB", "tissueC")
  rownames(C) <- "B_vs_C"
  
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    candidate_families = c("PO"),
    contrast_matrix = C,
    omnibus = TRUE,
    omnibus_test = "Wald",
    show_progress = FALSE
  )
  
  expect_true("omnibus" %in% names(result))
  
  # Should identify tissueB and tissueC as factor coefficients
  expect_true(all(result$omnibus$df == 2))
})

test_that("fit_gamlss_models without contrasts does not create omnibus", {
  skip_on_cran()
  
  set.seed(800)
  counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("gene", 1:10)
  
  design <- model.matrix(~ factor(rep(c("A", "B"), each = 5)))
  
  # No contrasts requested, omnibus should be ignored
  result <- fit_gamlss_models(
    counts,
    design_matrix = design,
    candidate_families = c("PO"),
    omnibus = TRUE,  # This should be ignored
    show_progress = FALSE
  )
  
  expect_false("omnibus" %in% names(result))
  expect_false("contrasts" %in% names(result))
})
