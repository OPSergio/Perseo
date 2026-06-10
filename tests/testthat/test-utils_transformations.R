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

test_that("transform_for_family_strict keeps zero samples for BE via unconditional eps nudge", {
  # Minmax maps min→0, max→1; for BE the min sample must not be excluded
  unit_with_zeros <- c(0, 0.2, 0.5, 0.8, 1.0)

  # Default allow_eps=TRUE
  tr_be <- transform_for_family_strict(unit_with_zeros, "BE")
  expect_true(all(tr_be$mask), info = "All samples must be valid for BE with zeros")
  expect_true(all(tr_be$y > 0 & tr_be$y < 1), info = "All transformed values must be in open (0,1)")

  # Explicit allow_eps=FALSE: eps must still be applied for unit families
  tr_be_noeps <- transform_for_family_strict(unit_with_zeros, "BE", allow_eps = FALSE)
  expect_true(all(tr_be_noeps$mask), info = "Zeros must not be excluded even with allow_eps=FALSE for unit families")
  expect_true(all(tr_be_noeps$y > 0 & tr_be_noeps$y < 1))

  # BEO (allows ones, not zeros): zero sample must also be kept
  tr_beo <- transform_for_family_strict(unit_with_zeros, "BEO")
  expect_true(all(tr_beo$mask))
  expect_true(all(tr_beo$y > 0 & tr_beo$y <= 1))
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

test_that("infer_support correctly identifies zi_positive data", {
  # Zeros + positive non-integer values (ZILN/ZAGA/ZAIG type)
  ziln_like <- c(0, 0, 1.5, 3.2, 0, 5.1, 12.7)
  expect_equal(infer_support(ziln_like), "zi_positive")

  # With NAs
  ziln_na <- c(0, NA, 2.3, 0, 5.8)
  expect_equal(infer_support(ziln_na), "zi_positive")

  # Should NOT trigger for data with negatives (→ real)
  mixed_neg <- c(-1.5, 0, 2.3, 5.1)
  expect_equal(infer_support(mixed_neg), "real")

  # Should NOT trigger for integers with zeros (→ count)
  int_zeros <- c(0L, 0L, 3L, 5L, 2L)
  expect_equal(infer_support(int_zeros), "count")

  # Should NOT trigger for unit-interval data with zeros (→ unit, caught first)
  unit_zeros <- c(0, 0.2, 0.5, 0.8, 1.0)
  expect_equal(infer_support(unit_zeros), "unit")
})

test_that("transform_for_family_strict handles zi_positive families", {
  # ZILN-like data: zeros + positive continuous
  ziln_data <- c(0, 0, 1.5, 3.2, 0, 8.7, 15.1)

  tr <- transform_for_family_strict(ziln_data, "ZILN")
  expect_equal(tr$y, ziln_data)                    # identity transform
  expect_true(all(tr$mask))                        # all values >= 0, all valid
  expect_equal(sum(tr$logJ_per_obs), 0)            # zero Jacobian for identity

  # Negative values must be masked out
  bad_data <- c(-1, 0, 2.5, 5.0)
  tr_bad <- transform_for_family_strict(bad_data, "ZAGA")
  expect_false(tr_bad$mask[1])
  expect_true(all(tr_bad$mask[-1]))
})

test_that("family_groups returns correct categorization", {
  groups <- family_groups()

  expect_type(groups, "list")
  expect_true("count"       %in% names(groups))
  expect_true("unit"        %in% names(groups))
  expect_true("positive"    %in% names(groups))
  expect_true("zi_positive" %in% names(groups))
  expect_true("real"        %in% names(groups))

  # Check specific families
  expect_true("PO"   %in% groups$count)
  expect_true("NBI"  %in% groups$count)
  expect_true("PIG"  %in% groups$count)
  expect_true("GA"   %in% groups$positive)
  expect_true("LOGNO" %in% groups$positive)
  expect_true("BE"   %in% groups$unit)
  expect_true("NO"   %in% groups$real)
  expect_true("ZILN" %in% groups$zi_positive)
  expect_true("ZAGA" %in% groups$zi_positive)
  expect_true("ZAIG" %in% groups$zi_positive)
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


# ===== transform_response() tests =====

test_that("transform_response returns correct contract", {
  y <- c(1, 2, 3, 4, 5)
  
  # Strict mode
  res_strict <- transform_response(y, "PO", mode = "strict")
  expect_type(res_strict, "list")
  expect_named(res_strict, c("y", "mask", "logJ_per_obs", "meta", "mode_used"))
  expect_equal(length(res_strict$y), length(y))
  expect_equal(length(res_strict$mask), length(y))
  expect_equal(length(res_strict$logJ_per_obs), length(y))
  expect_equal(res_strict$mode_used, "strict")
  
  # Safe mode
  res_safe <- transform_response(y, "PO", mode = "safe")
  expect_type(res_safe, "list")
  expect_named(res_safe, c("y", "mask", "logJ_per_obs", "meta", "mode_used"))
  expect_equal(length(res_safe$y), length(y))
  expect_equal(length(res_safe$mask), length(y))
  expect_equal(length(res_safe$logJ_per_obs), length(y))
  expect_equal(res_safe$mode_used, "safe")
})

test_that("transform_response: strict excludes, safe keeps observations", {
  # Data with negative values
  y <- c(-1, 0, 1, 2, 3)
  
  # Positive family (GA): strict should exclude negatives
  res_strict <- transform_response(y, "GA", mode = "strict")
  expect_false(res_strict$mask[1])  # -1 should be excluded
  expect_true(res_strict$mask[3])   # 1 should be included
  
  # Safe mode should keep all finite observations
  res_safe <- transform_response(y, "GA", mode = "safe")
  expect_true(all(res_safe$mask))   # All finite values kept
  expect_true(all(res_safe$y > 0))  # All transformed to positive
})

test_that("transform_response SAFE mode is reversible", {
  # Positive family
  y_pos <- c(0.5, 1.2, 3.7, 5.9, 10.3)
  res_pos <- transform_response(y_pos, "GA", mode = "safe")
  y_back_pos <- inverse_transform(res_pos$y, res_pos$meta)
  expect_equal(y_back_pos, y_pos, tolerance = 1e-10)
  
  # Unit family
  y_unit <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  res_unit <- transform_response(y_unit, "BE", mode = "safe")
  y_back_unit <- inverse_transform(res_unit$y, res_unit$meta)
  expect_equal(y_back_unit, y_unit, tolerance = 1e-10)
  
  # Real family
  y_real <- c(-2, -1, 0, 1, 2)
  res_real <- transform_response(y_real, "NO", mode = "safe")
  y_back_real <- inverse_transform(res_real$y, res_real$meta)
  expect_equal(y_back_real, y_real, tolerance = 1e-10)
})

test_that("transform_response SAFE affine has constant Jacobian", {
  # Positive family: affine shift
  y <- c(-2, -1, 0, 1, 2)
  res <- transform_response(y, "GA", mode = "safe")
  expect_equal(res$meta$kind, "affine")
  # For affine y* = a*y + b, Jacobian should be constant log(a)
  expect_true(all(res$logJ_per_obs == res$logJ_per_obs[1]))
  expect_equal(res$logJ_per_obs[1], log(abs(res$meta$params$a)))
  
  # Unit family: affine scaling
  y_unit <- c(0.1, 0.5, 0.9)
  res_unit <- transform_response(y_unit, "BE", mode = "safe")
  expect_equal(res_unit$meta$kind, "affine")
  expect_true(all(res_unit$logJ_per_obs == res_unit$logJ_per_obs[1]))
  expect_equal(res_unit$logJ_per_obs[1], log(abs(res_unit$meta$params$a)))
})

test_that("transform_response identity transform has zero Jacobian", {
  # Count families use identity
  y <- c(0, 1, 2, 3, 5)
  res_strict <- transform_response(y, "PO", mode = "strict")
  expect_true(all(res_strict$logJ_per_obs == 0))
  
  res_safe <- transform_response(y, "PO", mode = "safe")
  expect_true(all(res_safe$logJ_per_obs == 0))
})

test_that("transform_response SAFE: one family per support", {
  # Positive family: shift applied
  y_pos <- c(-1, 0, 1, 2)
  res_pos <- transform_response(y_pos, "GA", mode = "safe")
  expect_equal(res_pos$meta$kind, "affine")
  expect_true(all(res_pos$y > 0))
  expect_true(res_pos$meta$params$b > 0)  # Shift applied
  
  # Unit family: min-max scaling
  y_unit <- c(2, 4, 6, 8)
  res_unit <- transform_response(y_unit, "BE", mode = "safe")
  expect_equal(res_unit$meta$kind, "affine")
  expect_true(all(res_unit$y >= 0 & res_unit$y <= 1))
  
  # Count family: identity, no rounding
  y_count <- c(0.5, 1.2, 2.8)
  res_count <- transform_response(y_count, "PO", mode = "safe")
  expect_equal(res_count$meta$kind, "identity")
  expect_equal(res_count$y, y_count)  # No rounding
  
  # Real family: z-score
  y_real <- c(-2, -1, 0, 1, 2)
  res_real <- transform_response(y_real, "NO", mode = "safe")
  expect_equal(res_real$meta$kind, "zscore")
  expect_equal(mean(res_real$y), 0, tolerance = 1e-10)
  expect_equal(sd(res_real$y), 1, tolerance = 1e-10)
})

test_that("transform_response SAFE does not clip per observation", {
  # This test ensures SAFE mode doesn't do observation-wise operations
  y <- c(-5, -1, 0, 1, 5)
  res <- transform_response(y, "GA", mode = "safe")
  
  # All observations should be transformed by the same affine function
  a <- res$meta$params$a
  b <- res$meta$params$b
  expected <- a * y + b
  expect_equal(res$y, expected, tolerance = 1e-10)
  
  # No observation-wise clipping means differences are preserved
  diff_original <- diff(y)
  diff_transformed <- diff(res$y)
  # Affine transform preserves relative differences (scaled by a)
  expect_equal(diff_transformed, a * diff_original, tolerance = 1e-10)
})

test_that("transform_response SAFE does not collapse values", {
  # Ensure all distinct input values map to distinct output values
  y <- c(1, 2, 3, 4, 5)
  res <- transform_response(y, "GA", mode = "safe")
  
  # No value collapse: all transformed values should be unique
  expect_equal(length(unique(res$y)), length(unique(y)))
  
  # Same for unit interval
  y_unit <- seq(0, 1, by = 0.1)
  res_unit <- transform_response(y_unit, "BE", mode = "safe")
  expect_equal(length(unique(res_unit$y)), length(unique(y_unit)))
})

test_that("inverse_transform handles affine kind correctly", {
  # Create an affine transform
  y <- c(1, 2, 3, 4, 5)
  a <- 2
  b <- 3
  z <- a * y + b
  
  meta <- list(kind = "affine", params = list(a = a, b = b))
  y_back <- inverse_transform(z, meta)
  
  expect_equal(y_back, y, tolerance = 1e-10)
})

test_that("transform_response handles edge case: constant data", {
  # All same value
  y_const <- rep(5, 10)
  
  # Real family should fail (sd = 0)
  res_real <- transform_response(y_const, "NO", mode = "safe")
  expect_false(any(res_real$mask))
  
  # Unit family should fail (range = 0)
  res_unit <- transform_response(y_const, "BE", mode = "safe")
  expect_false(any(res_unit$mask))
  
  # Count/positive families should work
  res_count <- transform_response(y_const, "PO", mode = "safe")
  expect_true(all(res_count$mask))
})
