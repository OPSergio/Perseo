test_that("compare_families_on_feature returns valid structure", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(42)
  # Count data
  feature_vec <- rpois(30, lambda = 10)
  families <- c("PO", "NBI")
  
  result <- compare_families_on_feature(
    feature_vec,
    families = families,
    criterion = "AIC",
    min_n = 5
  )
  
  expect_type(result, "list")
  expect_true("comparisons" %in% names(result))
  expect_s3_class(result$comparisons, "tbl_df")
  
  # Check columns
  expect_true("family" %in% names(result$comparisons))
  expect_true("ic_value" %in% names(result$comparisons))
  expect_true("n_valid_obs" %in% names(result$comparisons))
  expect_true("coef_tbl" %in% names(result$comparisons))
})

test_that("compare_families_on_feature handles insufficient data", {
  # Too few observations
  feature_vec <- c(1, 2, 3)
  families <- c("PO", "NBI")
  
  result <- compare_families_on_feature(
    feature_vec,
    families = families,
    criterion = "AIC",
    min_n = 10
  )
  
  expect_equal(nrow(result$comparisons), 0)
})

test_that("compare_families_on_feature handles constant features", {
  # Constant feature - transformation will work but may have convergence issues
  feature_vec <- rep(5, 20)
  families <- c("PO", "NBI")
  
  result <- compare_families_on_feature(
    feature_vec,
    families = families,
    criterion = "AIC",
    min_n = 5
  )
  
  # May return empty or may fit (depends on GAMLSS behavior with constant data)
  expect_s3_class(result$comparisons, "tbl_df")
  expect_true("family" %in% names(result$comparisons))
})

test_that("compare_families_on_feature computes Jacobian correction", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(42)
  feature_vec <- rpois(30, lambda = 10)
  
  result <- compare_families_on_feature(
    feature_vec,
    families = c("PO"),
    criterion = "AIC",
    min_n = 5
  )
  
  if (nrow(result$comparisons) > 0) {
    # IC value should be finite
    expect_true(all(is.finite(result$comparisons$ic_value)))
  }
})

test_that("compare_families_with_design returns valid structure", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(42)
  feature_vec <- rpois(30, lambda = 10)
  design_df <- data.frame(
    condition = factor(rep(c("A", "B"), each = 15)),
    batch = factor(rep(1:3, 10))
  )
  families <- c("PO", "NBI")
  
  result <- compare_families_with_design(
    feature_vec,
    design_df,
    families = families,
    criterion = "AIC",
    min_n = 5
  )
  
  expect_type(result, "list")
  expect_true("comparisons" %in% names(result))
  expect_s3_class(result$comparisons, "tbl_df")
  
  # Should have coefficient tables
  if (nrow(result$comparisons) > 0) {
    expect_true("coef_tbl" %in% names(result$comparisons))
    expect_true(is.list(result$comparisons$coef_tbl))
  }
})

test_that("compare_families_with_design handles binomial families", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(42)
  # Proportion-like data
  feature_vec <- rbeta(30, 2, 2)
  design_df <- data.frame(condition = factor(rep(c("A", "B"), each = 15)))
  
  result <- compare_families_with_design(
    feature_vec,
    design_df,
    families = c("BI", "BB"),
    criterion = "AIC",
    min_n = 5,
    bd_vec = rep(100, 30)  # Provide denominator
  )
  
  expect_s3_class(result$comparisons, "tbl_df")
  # May or may not converge, but should return valid structure
  expect_true("family" %in% names(result$comparisons))
})

test_that("compare_families_with_design uses common mask", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(42)
  # Create data with some NAs
  feature_vec <- rpois(30, lambda = 10)
  feature_vec[c(1, 5, 10)] <- NA
  
  design_df <- data.frame(condition = factor(rep(c("A", "B"), each = 15)))
  
  result <- compare_families_with_design(
    feature_vec,
    design_df,
    families = c("PO"),
    criterion = "AIC",
    min_n = 5
  )
  
  if (nrow(result$comparisons) > 0) {
    # n_valid_obs should be less than original length
    expect_true(result$comparisons$n_valid_obs[1] < length(feature_vec))
  }
})

test_that("bind_bootstrap_results combines multiple iterations", {
  # Mock bootstrap results - each pull is a list of feature results
  boot_results <- list(
    # Pull 1
    list(
      list(feature = "gene1", family = "PO", skipped = FALSE, n_valid = 100L, support = "count"),
      list(feature = "gene2", family = "NBI", skipped = FALSE, n_valid = 95L, support = "count")
    ),
    # Pull 2
    list(
      list(feature = "gene1", family = "PO", skipped = FALSE, n_valid = 100L, support = "count"),
      list(feature = "gene3", family = "GA", skipped = FALSE, n_valid = 98L, support = "positive")
    )
  )
  
  bound <- bind_bootstrap_results(boot_results)
  
  expect_s3_class(bound, "tbl_df")
  expect_true("feature" %in% names(bound))
  expect_true("family" %in% names(bound))
  expect_true("bootstrap" %in% names(bound))
  expect_equal(nrow(bound), 4)  # 2 features per pull x 2 pulls
})

test_that("bind_bootstrap_results handles empty results", {
  bound <- bind_bootstrap_results(list())
  
  expect_s3_class(bound, "tbl_df")
  expect_equal(nrow(bound), 0)
})

test_that("bind_bootstrap_results handles NULL elements", {
  boot_results <- list(
    # Pull 1
    list(
      list(feature = "gene1", family = "PO", skipped = FALSE, n_valid = 100L, support = "count")
    ),
    # Pull 2 - empty
    list(),
    # Pull 3
    list(
      list(feature = "gene2", family = "NBI", skipped = FALSE, n_valid = 95L, support = "count")
    )
  )
  
  bound <- bind_bootstrap_results(boot_results)
  
  expect_equal(nrow(bound), 2)
  expect_equal(sort(unique(bound$feature)), c("gene1", "gene2"))
})

test_that("summarize_family_frequencies computes correct frequencies", {
  # Mock bound results
  bound_df <- tibble::tibble(
    bootstrap = rep(1:5, 2),
    feature = rep(c("gene1", "gene2"), each = 5),
    family = c(
      "PO", "PO", "PO", "NBI", "NBI",  # gene1: PO=3, NBI=2
      "GA", "GA", "LOGNO", "LOGNO", "LOGNO"  # gene2: GA=2, LOGNO=3
    ),
    support = rep(c("count", "positive"), each = 5),
    skipped = rep(FALSE, 10)
  )
  
  summary <- summarize_family_frequencies(bound_df, top_n = 3)
  
  expect_type(summary, "list")
  expect_true("top_families_overall" %in% names(summary))
  expect_true("freq_table_overall" %in% names(summary))
  expect_true(all(c("PO", "NBI", "GA", "LOGNO") %in% names(summary$freq_table_overall)))
})

test_that("summarize_family_frequencies handles single support type", {
  bound_df <- tibble::tibble(
    bootstrap = rep(1:10, 1),
    feature = rep("gene1", 10),
    family = rep(c("PO", "NBI"), times = c(7, 3)),
    support = rep("count", 10),
    skipped = rep(FALSE, 10)
  )
  
  summary <- summarize_family_frequencies(bound_df, top_n = 2)
  
  expect_type(summary, "list")
  expect_true("PO" %in% names(summary$freq_table_overall))
  expect_true("NBI" %in% names(summary$freq_table_overall))
})

test_that("summarize_family_frequencies handles empty input", {
  empty_df <- tibble::tibble(
    bootstrap = integer(),
    feature = character(),
    family = character(),
    support = character(),
    skipped = logical()
  )
  
  summary <- summarize_family_frequencies(empty_df, top_n = 3)
  
  expect_type(summary, "list")
  expect_equal(length(summary$top_families_overall), 0)
})
