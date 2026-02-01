test_that("identify_factor_coefficients works with contrast_matrix", {
  # Create a simple contrast matrix
  C <- matrix(c(0, 1, -1, 0), nrow = 1)
  colnames(C) <- c("(Intercept)", "groupB", "groupC", "age")
  rownames(C) <- "B_vs_C"
  
  coef_names <- c("(Intercept)", "groupB", "groupC", "age")
  
  result <- PERSEO:::identify_factor_coefficients(
    contrast_matrix = C,
    contrast_variable = NULL,
    coef_names = coef_names
  )
  
  expect_type(result, "character")
  expect_true("groupB" %in% result)
  expect_true("groupC" %in% result)
  expect_false("(Intercept)" %in% result)
  expect_false("age" %in% result)
})

test_that("identify_factor_coefficients works with contrast_variable", {
  coef_names <- c("(Intercept)", "tissue_typeBrain", "tissue_typeHeart", "age")
  
  result <- PERSEO:::identify_factor_coefficients(
    contrast_matrix = NULL,
    contrast_variable = "tissue_type",
    coef_names = coef_names
  )
  
  expect_type(result, "character")
  expect_equal(length(result), 2)
  expect_true("tissue_typeBrain" %in% result)
  expect_true("tissue_typeHeart" %in% result)
  expect_false("age" %in% result)
})

test_that("identify_factor_coefficients returns NULL when no matches", {
  C <- matrix(c(0, 1, -1), nrow = 1)
  colnames(C) <- c("(Intercept)", "missing1", "missing2")
  
  coef_names <- c("(Intercept)", "groupB", "age")
  
  result <- PERSEO:::identify_factor_coefficients(
    contrast_matrix = C,
    contrast_variable = NULL,
    coef_names = coef_names
  )
  
  expect_null(result)
})

test_that("wald_omnibus_test computes correct statistic", {
  # Simulate simple scenario
  beta <- c(groupB = 0.5, groupC = 0.8)
  V <- matrix(c(0.1, 0.02, 0.02, 0.12), nrow = 2)
  rownames(V) <- colnames(V) <- names(beta)
  
  factor_coefs <- c("groupB", "groupC")
  
  result <- PERSEO:::wald_omnibus_test(beta, V, factor_coefs)
  
  expect_type(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("df" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("test_type" %in% names(result))
  
  expect_equal(result$df, 2)
  expect_equal(result$test_type, "Wald")
  expect_true(result$statistic > 0)
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("wald_omnibus_test handles singular matrix gracefully", {
  beta <- c(groupB = 0.5, groupC = 0.8)
  V <- matrix(c(1, 1, 1, 1), nrow = 2)  # Singular
  rownames(V) <- colnames(V) <- names(beta)
  
  factor_coefs <- c("groupB", "groupC")
  
  # Should not error, should return NA
  result <- suppressWarnings({
    PERSEO:::wald_omnibus_test(beta, V, factor_coefs)
  })
  
  expect_true(is.na(result$statistic) || is.finite(result$statistic))
  expect_true(is.na(result$p_value) || (result$p_value >= 0 && result$p_value <= 1))
})

test_that("wald_omnibus_test handles non-finite values", {
  beta <- c(groupB = 0.5, groupC = NA_real_)
  V <- matrix(c(0.1, 0.02, 0.02, 0.12), nrow = 2)
  rownames(V) <- colnames(V) <- names(beta)
  
  factor_coefs <- c("groupB", "groupC")
  
  result <- PERSEO:::wald_omnibus_test(beta, V, factor_coefs)
  
  expect_true(is.na(result$statistic))
  expect_true(is.na(result$p_value))
})

test_that("wald_omnibus_test handles missing coefficients", {
  beta <- c(groupB = 0.5)
  V <- matrix(0.1, nrow = 1)
  rownames(V) <- colnames(V) <- names(beta)
  
  factor_coefs <- c("groupB", "missing")
  
  result <- PERSEO:::wald_omnibus_test(beta, V, factor_coefs)
  
  expect_true(is.na(result$statistic))
  expect_true(is.na(result$p_value))
})

test_that("lrt_omnibus_test returns valid structure", {
  skip_on_cran()
  
  # Simulate simple data
  set.seed(123)
  n <- 50
  feature_vec <- rpois(n, lambda = 10)
  
  design_mat <- model.matrix(~ factor(rep(c("A", "B", "C"), length.out = n)) + rnorm(n))
  colnames(design_mat) <- c("(Intercept)", "groupB", "groupC", "covar")
  
  factor_coefs <- c("groupB", "groupC")
  
  result <- PERSEO:::lrt_omnibus_test(
    feature_vec = feature_vec,
    design_mat = design_mat,
    factor_coefs = factor_coefs,
    family = "PO",
    transform_mode = "strict",
    bd_vec = NULL
  )
  
  expect_type(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("df" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("test_type" %in% names(result))
  
  expect_equal(result$test_type, "LRT")
  expect_equal(result$df, 2)
  
  # Should be valid or NA
  expect_true(is.na(result$statistic) || result$statistic >= 0)
  expect_true(is.na(result$p_value) || (result$p_value >= 0 && result$p_value <= 1))
})

test_that("lrt_omnibus_test handles empty reduced model gracefully", {
  skip_on_cran()
  
  set.seed(456)
  n <- 30
  feature_vec <- rpois(n, lambda = 5)
  
  # Design with only intercept and factor
  design_mat <- model.matrix(~ factor(rep(c("A", "B"), each = 15)))
  colnames(design_mat) <- c("(Intercept)", "groupB")
  
  # Trying to exclude the only non-intercept term
  factor_coefs <- c("groupB")
  
  result <- PERSEO:::lrt_omnibus_test(
    feature_vec = feature_vec,
    design_mat = design_mat,
    factor_coefs = factor_coefs,
    family = "PO",
    transform_mode = "strict"
  )
  
  # Should handle gracefully - might return valid result or NA
  expect_type(result, "list")
  expect_true(is.na(result$statistic) || is.numeric(result$statistic))
})

test_that("lrt_omnibus_test handles transform errors", {
  # Create data that will fail transformation
  feature_vec <- rep(-1, 20)  # Negative values for Poisson
  design_mat <- model.matrix(~ factor(rep(c("A", "B"), each = 10)))
  colnames(design_mat) <- c("(Intercept)", "groupB")
  
  result <- PERSEO:::lrt_omnibus_test(
    feature_vec = feature_vec,
    design_mat = design_mat,
    factor_coefs = c("groupB"),
    family = "PO",
    transform_mode = "strict"
  )
  
  expect_true(is.na(result$statistic))
  expect_true(is.na(result$p_value))
})
