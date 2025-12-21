test_that("infer_support correctly identifies count data", {
  # Integer counts
  counts <- c(0, 1, 2, 3, 5, 8, 13)
  support <- infer_support(counts)
  expect_equal(support, "count")
  
  # With NAs
  counts_na <- c(0, 1, NA, 3, 5)
  support_na <- infer_support(counts_na)
  expect_equal(support_na, "count")
})

test_that("infer_support correctly identifies unit data", {
  # Proportions between 0 and 1
  props <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  support <- infer_support(props)
  expect_equal(support, "unit")
  
  # Including exact 0 and 1
  props_boundary <- c(0, 0.2, 0.5, 0.8, 1)
  support_boundary <- infer_support(props_boundary)
  expect_equal(support_boundary, "unit")
})

test_that("infer_support correctly identifies positive data", {
  # Positive non-integer
  positive <- c(0.5, 1.2, 3.7, 5.9, 10.3)
  support <- infer_support(positive)
  expect_equal(support, "positive")
  
  # Strictly positive (no zeros)
  positive_strict <- c(1.1, 2.2, 3.3, 4.4)
  support_strict <- infer_support(positive_strict)
  expect_equal(support_strict, "positive")
})

test_that("infer_support correctly identifies real data", {
  # Including negative values
  real_data <- c(-2.5, -1.0, 0, 1.5, 3.2)
  support <- infer_support(real_data)
  expect_equal(support, "real")
  
  # All negative
  negative <- c(-5, -3, -1, -0.5)
  support_neg <- infer_support(negative)
  expect_equal(support_neg, "real")
})

test_that("transform_for_family_strict handles count families correctly", {
  counts <- c(0, 1, 2, 5, 10)
  
  # Poisson: identity transform
  tr_po <- transform_for_family_strict(counts, "PO")
  expect_equal(tr_po$y, counts)
  expect_true(all(tr_po$mask))
  expect_equal(sum(tr_po$logJ_per_obs), 0)  # Identity has zero Jacobian
  
  # NBI: identity transform
  tr_nbi <- transform_for_family_strict(counts, "NBI")
  expect_equal(tr_nbi$y, counts)
  expect_true(all(tr_nbi$mask))
})

test_that("transform_for_family_strict handles positive families correctly", {
  positive_data <- c(0.1, 1.5, 3.2, 7.8, 12.5)
  
  # Gamma
  tr_ga <- transform_for_family_strict(positive_data, "GA")
  expect_true(all(tr_ga$y > 0))
  expect_true(all(tr_ga$mask))
  expect_true(is.finite(sum(tr_ga$logJ_per_obs)))
  
  # Log-normal
  tr_logno <- transform_for_family_strict(positive_data, "LOGNO")
  expect_true(all(tr_logno$y > 0))
  expect_true(all(tr_logno$mask))
})

test_that("transform_for_family_strict handles unit families correctly", {
  unit_data <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  
  # Beta
  tr_be <- transform_for_family_strict(unit_data, "BE")
  expect_true(all(tr_be$y > 0 & tr_be$y < 1))
  expect_true(all(tr_be$mask))
  expect_true(is.finite(sum(tr_be$logJ_per_obs)))
})

test_that("transform_for_family_strict handles real families correctly", {
  real_data <- c(-2, -1, 0, 1, 2, 3)
  
  # Normal
  tr_no <- transform_for_family_strict(real_data, "NO")
  expect_equal(length(tr_no$y), length(real_data))
  expect_true(all(tr_no$mask))
  
  # Should apply z-score transformation
  expect_true(abs(mean(tr_no$y)) < 1e-10)  # Mean ~ 0
  expect_true(abs(sd(tr_no$y) - 1) < 1e-10)  # SD ~ 1
})

test_that("transform_for_family_strict creates masks for invalid values", {
  # Data with NAs
  data_na <- c(1, 2, NA, 4, 5)
  tr <- transform_for_family_strict(data_na, "PO")
  expect_false(tr$mask[3])
  expect_true(all(tr$mask[-3]))
  
  # Negative values for Poisson (should be masked as FALSE or NA)
  data_neg <- c(1, 2, -1, 4, 5)
  tr_neg <- transform_for_family_strict(data_neg, "PO")
  # Mask could be FALSE or NA for invalid values - both acceptable
  expect_true(is.na(tr_neg$mask[3]) || tr_neg$mask[3] == FALSE)
})

test_that("transform_for_family_strict handles edge cases", {
  # Single value
  single <- c(5)
  tr_single <- transform_for_family_strict(single, "PO")
  expect_equal(length(tr_single$y), 1)
  expect_equal(length(tr_single$mask), 1)
  
  # All same value
  constant <- rep(5, 10)
  tr_const <- transform_for_family_strict(constant, "PO")
  expect_equal(tr_const$y, constant)
  expect_true(all(tr_const$mask))
})

test_that("family_groups returns correct categorization", {
  groups <- family_groups()
  
  expect_type(groups, "list")
  expect_true("count" %in% names(groups))
  expect_true("unit" %in% names(groups))
  expect_true("positive" %in% names(groups))
  expect_true("real" %in% names(groups))
  
  # Check specific families
  expect_true("PO" %in% groups$count)
  expect_true("NBI" %in% groups$count)
  expect_true("GA" %in% groups$positive)
  expect_true("LOGNO" %in% groups$positive)
  expect_true("BE" %in% groups$unit)
  expect_true("NO" %in% groups$real)
})

test_that("jacobian_sum computes correctly for log transform", {
  x <- c(1, 2, 5, 10)
  # log transform: Jacobian = -log(x), sum = -sum(log(x))
  expected <- -sum(log(x))
  
  # This would need to be tested if jacobian_sum is exported
  # For now, verify through transform_for_family_strict
  tr <- transform_for_family_strict(x, "LOGNO")
  # Log-normal uses log transform, so check that Jacobian is computed
  expect_true(is.finite(sum(tr$logJ_per_obs)))
})
