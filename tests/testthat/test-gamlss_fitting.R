test_that("compute_ic_penalty returns correct values", {
  n <- 100
  
  # AIC: penalty = 2
  expect_equal(compute_ic_penalty("AIC", NULL, n), 2)
  
  # BIC: penalty = log(n)
  expect_equal(compute_ic_penalty("BIC", NULL, n), log(n))
  
  # GAIC with NULL k: uses log(n)
  expect_equal(compute_ic_penalty("GAIC", NULL, n), log(n))
  
  # GAIC with specified k
  expect_equal(compute_ic_penalty("GAIC", 3.5, n), 3.5)
})

test_that("compute_ic_penalty handles edge cases", {
  # Small sample size
  expect_equal(compute_ic_penalty("AIC", NULL, 5), 2)
  expect_equal(compute_ic_penalty("BIC", NULL, 5), log(5))
  
  # Large sample size
  n_large <- 10000
  expect_equal(compute_ic_penalty("BIC", NULL, n_large), log(n_large))
})

test_that("instantiate_gamlss_family retrieves known families", {
  skip_if_not_installed("gamlss.dist")
  
  # Common families
  po_fam <- instantiate_gamlss_family("PO")
  expect_s3_class(po_fam, "gamlss.family")
  
  nbi_fam <- instantiate_gamlss_family("NBI")
  expect_s3_class(nbi_fam, "gamlss.family")
  
  ga_fam <- instantiate_gamlss_family("GA")
  expect_s3_class(ga_fam, "gamlss.family")
})

test_that("instantiate_gamlss_family returns NULL for unknown families", {
  unknown_fam <- instantiate_gamlss_family("NONEXISTENT_FAMILY_XYZ")
  expect_null(unknown_fam)
})

test_that("compute_jacobian_corrected_ic computes correctly", {
  # Mock fit object
  fit <- list(df.fit = 3)
  class(fit) <- "gamlss"
  attr(fit, "class") <- "gamlss"
  
  # Mock logLik
  loglik_value <- -50
  
  # Temporarily mock extract_gamlss_loglik
  with_mocked_bindings(
    extract_gamlss_loglik = function(f) loglik_value,
    {
      penalty <- 2
      logJ_sum <- 10
      
      # Formula: -2*loglik + penalty*df_fit - 2*logJ_sum
      expected_ic <- (-2 * loglik_value) + (penalty * 3) - (2 * logJ_sum)
      
      result <- compute_jacobian_corrected_ic(fit, penalty, logJ_sum)
      expect_equal(result, expected_ic)
    }
  )
})

test_that("compute_jacobian_corrected_ic handles zero Jacobian", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(123)
  y <- rpois(50, lambda = 5)
  data_fit <- data.frame(y = y)
  family_obj <- instantiate_gamlss_family("PO")
  
  fit <- fit_gamlss_safely(y ~ 1, data_fit, family_obj)
  
  if (!is.null(fit) && !inherits(fit, "try-error")) {
    penalty <- 2
    logJ_sum <- 0  # No transformation
    
    result <- compute_jacobian_corrected_ic(fit, penalty, logJ_sum)
    expect_true(is.finite(result))
  }
})

test_that("compute_jacobian_corrected_ic handles negative Jacobian", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(123)
  y <- rpois(50, lambda = 5)
  data_fit <- data.frame(y = y)
  family_obj <- instantiate_gamlss_family("PO")
  
  fit <- fit_gamlss_safely(y ~ 1, data_fit, family_obj)
  
  if (!is.null(fit) && !inherits(fit, "try-error")) {
    penalty <- 2
    logJ_sum <- -5  # Negative Jacobian sum is mathematically valid
    
    result <- compute_jacobian_corrected_ic(fit, penalty, logJ_sum)
    expect_true(is.finite(result))
  }
})

test_that("fit_gamlss_safely handles successful fits", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  # Simple Poisson data
  set.seed(123)
  y <- rpois(50, lambda = 5)
  data_fit <- data.frame(y = y)
  family_obj <- instantiate_gamlss_family("PO")
  
  fit <- fit_gamlss_safely(y ~ 1, data_fit, family_obj)
  
  expect_s3_class(fit, "gamlss")
  expect_false(inherits(fit, "try-error"))
})

test_that("fit_gamlss_safely returns try-error for bad data", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  # Negative values for Poisson (invalid)
  y <- c(-1, -2, -3, 1, 2)
  data_fit <- data.frame(y = y)
  family_obj <- instantiate_gamlss_family("PO")
  
  fit <- fit_gamlss_safely(y ~ 1, data_fit, family_obj)
  
  # Should return try-error or NULL
  expect_true(inherits(fit, "try-error") || is.null(fit))
})

test_that("extract_gamlss_loglik extracts log-likelihood", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(123)
  y <- rpois(50, lambda = 5)
  data_fit <- data.frame(y = y)
  family_obj <- instantiate_gamlss_family("PO")
  
  fit <- fit_gamlss_safely(y ~ 1, data_fit, family_obj)
  
  if (!inherits(fit, "try-error")) {
    loglik <- extract_gamlss_loglik(fit)
    expect_type(loglik, "double")
    expect_true(is.finite(loglik))
    expect_true(loglik < 0)  # Log-likelihood typically negative
  }
})

test_that("extract_gamlss_loglik returns NA for failed fits", {
  # Mock a failed fit
  fake_fit <- structure(list(), class = "fake_gamlss")
  
  loglik <- extract_gamlss_loglik(fake_fit)
  expect_true(is.na(loglik) || !is.finite(loglik))
})

test_that("extract_mu_coefficients returns tibble with correct columns", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(123)
  y <- rpois(50, lambda = 5)
  x <- rnorm(50)
  data_fit <- data.frame(y = y, x = x)
  family_obj <- instantiate_gamlss_family("PO")
  
  fit <- fit_gamlss_safely(y ~ x, data_fit, family_obj)
  
  if (!inherits(fit, "try-error")) {
    coef_tbl <- extract_mu_coefficients(fit)
    
    expect_s3_class(coef_tbl, "tbl_df")
    expect_true("term" %in% names(coef_tbl))
    expect_true("effect" %in% names(coef_tbl))
    expect_true("se" %in% names(coef_tbl))
    expect_true("stat" %in% names(coef_tbl))
    expect_true("pval" %in% names(coef_tbl))
    
    # Should have intercept and x coefficient
    expect_equal(nrow(coef_tbl), 2)
    expect_true("(Intercept)" %in% coef_tbl$term)
    expect_true("x" %in% coef_tbl$term)
  }
})

test_that("extract_mu_coefficients handles intercept-only models", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  
  set.seed(123)
  y <- rpois(50, lambda = 5)
  data_fit <- data.frame(y = y)
  family_obj <- instantiate_gamlss_family("PO")
  
  fit <- fit_gamlss_safely(y ~ 1, data_fit, family_obj)
  
  if (!inherits(fit, "try-error")) {
    coef_tbl <- extract_mu_coefficients(fit)
    
    expect_equal(nrow(coef_tbl), 1)
    expect_equal(coef_tbl$term[1], "(Intercept)")
    expect_true(is.finite(coef_tbl$effect[1]))
  }
})
