test_that("apply_contrasts computes correct estimates and SEs", {
  # Synthetic example: 3 coefficients
  beta <- c("(Intercept)" = 1.0, "groupB" = 0.5, "groupC" = -0.3)
  V <- matrix(c(
    0.04, 0.01, 0.01,
    0.01, 0.09, 0.02,
    0.01, 0.02, 0.16
  ), nrow = 3, ncol = 3)
  rownames(V) <- colnames(V) <- names(beta)
  
  # Contrast: B - C = groupB - groupC = 0.5 - (-0.3) = 0.8
  C <- matrix(c(0, 1, -1), nrow = 1)
  colnames(C) <- names(beta)
  rownames(C) <- "B_vs_C"
  
  result <- apply_contrasts(beta, V, C)
  
  expect_equal(nrow(result), 1)
  expect_equal(result$contrast, "B_vs_C")
  expect_equal(result$estimate, 0.8)
  
  # SE = sqrt(Var(B) + Var(C) - 2*Cov(B,C))
  # = sqrt(0.09 + 0.16 - 2*0.02) = sqrt(0.21) ≈ 0.458
  expected_se <- sqrt(0.09 + 0.16 - 2 * 0.02)
  expect_equal(result$se, expected_se, tolerance = 1e-6)
  expect_equal(result$z, 0.8 / expected_se, tolerance = 1e-6)
  expect_true(result$p_value > 0 && result$p_value < 1)
})

test_that("apply_contrasts handles missing columns gracefully", {
  beta <- c("(Intercept)" = 1.0, "groupB" = 0.5)
  V <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2, ncol = 2)
  rownames(V) <- colnames(V) <- names(beta)
  
  # Contrast includes a column not in beta
  C <- matrix(c(0, 1, -1), nrow = 1)
  colnames(C) <- c("(Intercept)", "groupB", "groupC")
  rownames(C) <- "B_vs_C"
  
  result <- apply_contrasts(beta, V, C)
  
  # Should treat groupC as 0
  expect_equal(result$estimate, 0.5)  # Just groupB
})

test_that("apply_contrasts returns NA when vcov is non-finite", {
  beta <- c("(Intercept)" = 1.0, "groupB" = 0.5)
  V <- matrix(c(0.04, Inf, Inf, 0.09), nrow = 2, ncol = 2)
  rownames(V) <- colnames(V) <- names(beta)
  
  C <- matrix(c(0, 1), nrow = 1)
  colnames(C) <- names(beta)
  rownames(C) <- "test"
  
  expect_warning(result <- apply_contrasts(beta, V, C), "non-finite")
  expect_true(is.na(result$estimate))
})

test_that("fit_gamlss_models returns contrasts tibble with custom matrix", {
  skip_if_not_installed("gamlss")
  
  # Simulate small dataset with 3-level factor
  set.seed(123)
  n_samples <- 30
  n_features <- 5
  
  status <- factor(rep(c("A", "B", "C"), each = 10))
  design <- model.matrix(~ status)
  
  # Simulate count data
  counts <- matrix(
    rpois(n_features * n_samples, lambda = 50),
    nrow = n_features, ncol = n_samples
  )
  rownames(counts) <- paste0("gene_", 1:n_features)
  
  # Custom contrast matrix: B vs C
  C <- matrix(c(0, 1, -1), nrow = 1)
  colnames(C) <- c("(Intercept)", "statusB", "statusC")
  rownames(C) <- "B_vs_C"
  
  # Fit with custom contrasts
  result <- fit_gamlss_models(
    counts_matrix = counts,
    design_matrix = design,
    candidate_families = c("NBI", "PO"),
    contrast_matrix = C,
    workers = 1,
    show_progress = FALSE
  )
  
  expect_true("contrasts" %in% names(result))
  expect_s3_class(result$contrasts, "tbl_df")
  expect_true(all(c("feature", "family", "contrast", "estimate", "se", "z", "p_value", "p_adj") %in% names(result$contrasts)))
  
  # Should have 1 contrast per feature
  n_contrasts_per_feature <- result$contrasts %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  expect_true(all(n_contrasts_per_feature$n == 1))
  expect_true(all(result$contrasts$contrast == "B_vs_C"))
  
  # CRITICAL: Verify numerical values are valid
  expect_true(all(is.finite(result$contrasts$estimate) | is.na(result$contrasts$estimate)))
  expect_true(all(result$contrasts$se >= 0 | is.na(result$contrasts$se)))
  expect_true(all(is.finite(result$contrasts$z) | is.na(result$contrasts$z)))
  expect_true(all(result$contrasts$p_value >= 0 & result$contrasts$p_value <= 1 | is.na(result$contrasts$p_value)))
  expect_true(all(result$contrasts$p_adj >= 0 & result$contrasts$p_adj <= 1 | is.na(result$contrasts$p_adj)))
  
  # Verify z = estimate / se relationship (when both are finite)
  finite_idx <- is.finite(result$contrasts$estimate) & is.finite(result$contrasts$se) & result$contrasts$se > 0
  if (any(finite_idx)) {
    expected_z <- result$contrasts$estimate[finite_idx] / result$contrasts$se[finite_idx]
    expect_equal(result$contrasts$z[finite_idx], expected_z, tolerance = 1e-10)
  }
})

test_that("fit_gamlss_models works without contrasts (backward compatibility)", {
  skip_if_not_installed("gamlss")
  
  set.seed(456)
  n_samples <- 20
  n_features <- 3
  
  design <- data.frame(
    intercept = rep(1, n_samples),
    treatment = c(rep(0, 10), rep(1, 10))
  )
  
  counts <- matrix(
    rpois(n_features * n_samples, lambda = 30),
    nrow = n_features, ncol = n_samples
  )
  rownames(counts) <- paste0("gene_", 1:n_features)
  
  result <- fit_gamlss_models(
    counts_matrix = counts,
    design_matrix = design,
    candidate_families = c("PO", "NBI"),
    workers = 1,
    show_progress = FALSE
  )
  
  # Should NOT have contrasts element
  expect_false("contrasts" %in% names(result))
  
  # Should have results and selection
  expect_true("results" %in% names(result))
  expect_true("selection" %in% names(result))
})

test_that("fit_gamlss_models handles multiple custom contrasts", {
  skip_if_not_installed("gamlss")
  
  set.seed(789)
  n_samples <- 30
  n_features <- 3
  
  status <- factor(rep(c("A", "B", "C"), each = 10))
  design <- model.matrix(~ status)
  
  counts <- matrix(
    rpois(n_features * n_samples, lambda = 40),
    nrow = n_features, ncol = n_samples
  )
  rownames(counts) <- paste0("gene_", 1:n_features)
  
  # Multiple contrasts: B vs C and (B+C)/2 vs A
  C_custom <- rbind(
    "B_vs_C" = c(0, 1, -1),
    "BC_avg_vs_A" = c(0, 0.5, 0.5)
  )
  colnames(C_custom) <- c("(Intercept)", "statusB", "statusC")
  
  result <- fit_gamlss_models(
    counts_matrix = counts,
    design_matrix = design,
    candidate_families = c("PO", "NBI"),
    contrast_matrix = C_custom,
    workers = 1,
    show_progress = FALSE
  )
  
  expect_true("contrasts" %in% names(result))
  
  # Should have 2 contrasts per feature
  n_contrasts_per_feature <- result$contrasts %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  expect_true(all(n_contrasts_per_feature$n == 2))
  expect_true(all(c("B_vs_C", "BC_avg_vs_A") %in% result$contrasts$contrast))
  
  # CRITICAL: Verify estimates match manual computation
  # For each feature, verify that contrast = c1*beta1 + c2*beta2
  for (feat in unique(result$contrasts$feature)) {
    feat_results <- result$results[result$results$feature == feat, ]
    feat_contrasts <- result$contrasts[result$contrasts$feature == feat, ]
    
    # Get coefficients from results
    beta_B <- feat_results$effect[feat_results$term == "statusB"]
    beta_C <- feat_results$effect[feat_results$term == "statusC"]
    
    if (length(beta_B) == 1 && length(beta_C) == 1) {
      # Verify B_vs_C = beta_B - beta_C
      b_vs_c_contrast <- feat_contrasts$estimate[feat_contrasts$contrast == "B_vs_C"]
      if (is.finite(b_vs_c_contrast)) {
        expect_equal(b_vs_c_contrast, beta_B - beta_C, tolerance = 1e-10)
      }
      
      # Verify BC_avg_vs_A = 0.5*beta_B + 0.5*beta_C
      bc_avg_contrast <- feat_contrasts$estimate[feat_contrasts$contrast == "BC_avg_vs_A"]
      if (is.finite(bc_avg_contrast)) {
        expect_equal(bc_avg_contrast, 0.5 * beta_B + 0.5 * beta_C, tolerance = 1e-10)
      }
    }
  }
  
  # Verify SE > 0 for all finite contrasts
  finite_contrasts <- result$contrasts[is.finite(result$contrasts$estimate), ]
  if (nrow(finite_contrasts) > 0) {
    expect_true(all(finite_contrasts$se > 0))
  }
  
  # Verify p-values are two-sided from z-scores
  finite_idx <- is.finite(result$contrasts$z)
  if (any(finite_idx)) {
    expected_p <- 2 * pnorm(abs(result$contrasts$z[finite_idx]), lower.tail = FALSE)
    expect_equal(result$contrasts$p_value[finite_idx], expected_p, tolerance = 1e-10)
  }
  
  # Verify p_adj is grouped by contrast
  b_vs_c_pvals <- result$contrasts$p_value[result$contrasts$contrast == "B_vs_C"]
  b_vs_c_padj <- result$contrasts$p_adj[result$contrasts$contrast == "B_vs_C"]
  if (all(is.finite(b_vs_c_pvals))) {
    expect_equal(b_vs_c_padj, p.adjust(b_vs_c_pvals, method = "BH"), tolerance = 1e-10)
  }
})

test_that("contrast estimates are consistent with coefficient estimates", {
  skip_if_not_installed("gamlss")
  
  # Use a realistic dataset where models converge reliably
  set.seed(100)
  n_samples <- 90
  status <- factor(rep(c("A", "B", "C"), each = 30))
  design <- model.matrix(~ status)
  
  # Multiple features with varying expression
  n_features <- 10
  counts <- matrix(
    rpois(n_features * n_samples, lambda = rep(c(100, 150, 200), each = 30)),
    nrow = n_features, ncol = n_samples, byrow = TRUE
  )
  rownames(counts) <- paste0("gene_", 1:n_features)
  
  # Define multiple contrasts
  C <- rbind(
    "B_vs_A" = c(0, 1, 0),      # Same as coefficient statusB
    "C_vs_A" = c(0, 0, 1),      # Same as coefficient statusC
    "B_vs_C" = c(0, 1, -1)      # statusB - statusC
  )
  colnames(C) <- c("(Intercept)", "statusB", "statusC")
  
  # Fit model with contrasts (use stable families)
  result <- fit_gamlss_models(
    counts_matrix = counts,
    design_matrix = design,
    candidate_families = c("PO", "NBI"),
    contrast_matrix = C,
    workers = 1,
    show_progress = FALSE
  )
  
  # For each feature, verify contrast consistency
  for (feat in unique(result$contrasts$feature)) {
    feat_results <- result$results[result$results$feature == feat, ]
    feat_contrasts <- result$contrasts[result$contrasts$feature == feat, ]
    
    # Get coefficients
    beta_B <- feat_results$effect[feat_results$term == "statusB"]
    beta_C <- feat_results$effect[feat_results$term == "statusC"]
    
    if (length(beta_B) == 1 && length(beta_C) == 1) {
      # Verify B_vs_A matches statusB coefficient
      b_vs_a <- feat_contrasts$estimate[feat_contrasts$contrast == "B_vs_A"]
      if (is.finite(b_vs_a)) {
        expect_equal(b_vs_a, beta_B, tolerance = 1e-10,
                     info = paste("Feature:", feat, "- B_vs_A should equal statusB"))
      }
      
      # Verify C_vs_A matches statusC coefficient
      c_vs_a <- feat_contrasts$estimate[feat_contrasts$contrast == "C_vs_A"]
      if (is.finite(c_vs_a)) {
        expect_equal(c_vs_a, beta_C, tolerance = 1e-10,
                     info = paste("Feature:", feat, "- C_vs_A should equal statusC"))
      }
      
      # Verify B_vs_C = statusB - statusC
      b_vs_c <- feat_contrasts$estimate[feat_contrasts$contrast == "B_vs_C"]
      if (is.finite(b_vs_c)) {
        expect_equal(b_vs_c, beta_B - beta_C, tolerance = 1e-10,
                     info = paste("Feature:", feat, "- B_vs_C should equal statusB - statusC"))
      }
    }
  }
  
  # Verify all contrasts have valid statistics when estimate is finite
  finite_contrasts <- result$contrasts[is.finite(result$contrasts$estimate), ]
  if (nrow(finite_contrasts) > 0) {
    # All finite contrasts should have positive SE
    expect_true(all(finite_contrasts$se > 0))
    
    # z-statistic should match estimate/se
    expected_z <- finite_contrasts$estimate / finite_contrasts$se
    expect_equal(finite_contrasts$z, expected_z, tolerance = 1e-12)
    
    # p-value should match 2*pnorm(-|z|)
    expected_p <- 2 * pnorm(abs(finite_contrasts$z), lower.tail = FALSE)
    expect_equal(finite_contrasts$p_value, expected_p, tolerance = 1e-12)
  }
})


test_that("contrasts work with 5 different families and vcov extraction", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("dplyr")
  
  set.seed(456)
  n_samples <- 30
  n_features <- 10  # 2 per family ideally
  
  # 3-level factor
  status <- factor(rep(c("A", "B", "C"), each = 10))
  design <- model.matrix(~ status)
  
  # Simulate data that favors different families:
  # Features 1-2: count data (NBI/PO)
  # Features 3-4: low count data (PO/ZIP)
  # Features 5-6: binomial-like (0-10 range, BI family)
  # Features 7-8: proportions (0-1 range, BETA)
  # Features 9-10: positive continuous (GA/LOGNO)
  
  counts <- matrix(0, nrow = n_features, ncol = n_samples)
  
  # Features 1-2: overdispersed counts (favor NBI)
  counts[1,] <- rnbinom(n_samples, size = 2, mu = 100)
  counts[2,] <- rnbinom(n_samples, size = 3, mu = 80)
  
  # Features 3-4: low counts (favor PO)
  counts[3,] <- rpois(n_samples, lambda = 5)
  counts[4,] <- rpois(n_samples, lambda = 8)
  
  # Features 5-6: binomial-like counts (0-10)
  counts[5,] <- rbinom(n_samples, size = 10, prob = 0.3)
  counts[6,] <- rbinom(n_samples, size = 10, prob = 0.5)
  
  # Features 7-8: beta-like proportions (0-1)
  counts[7,] <- rbeta(n_samples, shape1 = 2, shape2 = 5)
  counts[8,] <- rbeta(n_samples, shape1 = 3, shape2 = 3)
  
  # Features 9-10: positive continuous (favor GA)
  counts[9,] <- rgamma(n_samples, shape = 2, rate = 0.5)
  counts[10,] <- rgamma(n_samples, shape = 3, rate = 1)
  
  rownames(counts) <- paste0("feat_", 1:n_features)
  
  # Define B vs C contrast
  C <- matrix(c(0, 1, -1), nrow = 1)
  colnames(C) <- c("(Intercept)", "statusB", "statusC")
  rownames(C) <- "B_vs_C"
  
  # Use 5 different families
  result <- fit_gamlss_models(
    counts_matrix = counts,
    design_matrix = design,
    candidate_families = c("NBI", "PO", "BI", "BETA", "GA"),
    contrast_matrix = C,
    workers = 1,
    show_progress = FALSE
  )
  
  # ===== DEBUGGING CHECKS =====
  
  # 1. Contrasts should exist
  expect_true("contrasts" %in% names(result))
  expect_s3_class(result$contrasts, "tbl_df")
  
  # Some features may be filtered or fail to fit - expect at least 50% success
  expect_true(nrow(result$contrasts) >= n_features / 2,
              info = sprintf("Only %d/%d features returned contrasts", nrow(result$contrasts), n_features))
  
  # 2. CRITICAL: Estimates should NOT all be NA
  n_na_estimates <- sum(is.na(result$contrasts$estimate))
  expect_true(n_na_estimates < nrow(result$contrasts),
              info = sprintf("All %d estimates are NA - vcov extraction likely failing", nrow(result$contrasts)))
  
  # 3. CRITICAL: SEs should NOT all be NA
  n_na_ses <- sum(is.na(result$contrasts$se))
  expect_true(n_na_ses < nrow(result$contrasts),
              info = sprintf("All %d SEs are NA - vcov likely not extracted", nrow(result$contrasts)))
  
  # 4. For non-NA rows, verify vcov was used correctly: z = estimate / se
  valid_rows <- !is.na(result$contrasts$estimate) & 
                !is.na(result$contrasts$se) & 
                result$contrasts$se > 0
  
  if (sum(valid_rows) > 0) {
    valid_contrasts <- result$contrasts[valid_rows, ]
    expected_z <- valid_contrasts$estimate / valid_contrasts$se
    expect_equal(valid_contrasts$z, expected_z, tolerance = 1e-10,
                 info = "z-statistic doesn't match estimate/se - vcov may be corrupted")
  } else {
    fail("No valid contrast rows - all vcov extractions failed")
  }
  
  # 5. Verify multiple families were actually selected (at least 2)
  n_families_used <- length(unique(result$selection$best_family))
  expect_true(n_families_used >= 2,
              info = sprintf("Only %d families selected, expected >= 2", n_families_used))
  
  # 6. Verify p-values are valid
  valid_p <- !is.na(result$contrasts$p_value)
  if (sum(valid_p) > 0) {
    expect_true(all(result$contrasts$p_value[valid_p] >= 0 & 
                    result$contrasts$p_value[valid_p] <= 1),
                info = "Some p-values outside [0,1] range")
  }
  
  # 7. Print diagnostic info for debugging
  message("\n=== DIAGNOSTIC INFO ===")
  message(sprintf("Features with valid contrasts: %d/%d", sum(valid_rows), n_features))
  message(sprintf("Families used: %s", paste(unique(result$selection$best_family), collapse = ", ")))
  message(sprintf("NA estimates: %d, NA SEs: %d", n_na_estimates, n_na_ses))
})

