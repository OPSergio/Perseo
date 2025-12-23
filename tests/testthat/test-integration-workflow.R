test_that("PERSEO end-to-end workflow on realistic synthetic omics data", {

  # -----------------------------
  # Dependencies
  # -----------------------------
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")

  # -----------------------------
  # Debug utilities (prints)
  # -----------------------------
  debug_on <- identical(Sys.getenv("PERSEO_TEST_DEBUG"), "1") ||
    isTRUE(getOption("perseo.debug", FALSE))

  dbg <- function(...) if (debug_on) message(sprintf(...))
  dbg_section <- function(title) if (debug_on) {
    message("\n", paste0(rep("=", 70), collapse=""))
    message(title)
    message(paste0(rep("=", 70), collapse=""))
  }

  # -----------------------------
  # Reproducibility
  # -----------------------------
  set.seed(20250101)

  dbg_section("1) Synthetic data generation: realistic multi-support omics matrix")

  # Design: 2 groups, batch effect, continuous covariate
  n_samples <- 48
  condition <- factor(rep(c("A", "B"), each = n_samples / 2))
  batch <- factor(rep(rep(1:3, each = 8), times = 2))  # 3 batches balanced by condition
  age <- scale(rnorm(n_samples, mean = 55, sd = 10))[,1]  # continuous covariate

  design <- data.frame(condition = condition, batch = batch, age = age)

  # Features: create a mixed matrix (features x samples)
  # Blocks:
  # - count NB (RNA-seq like)
  # - zero-inflated count (metagenomics-like)
  # - unit interval beta (fractions/proportions)
  # - positive continuous gamma/lognormal (proteomics intensity-like)
  # - real (normalized logFC-like)
  n_nb  <- 120
  n_zic <- 60
  n_unit <- 40
  n_pos <- 40
  n_real <- 40
  n_features <- n_nb + n_zic + n_unit + n_pos + n_real

  feature_ids <- paste0("feat_", seq_len(n_features))

  # True differential set (inject true signal in multiple supports)
  # Keep it moderate: a subset differential, not all.
  n_de_nb   <- 20
  n_de_zic  <- 10
  n_de_unit <- 8
  n_de_pos  <- 8
  n_de_real <- 8

  de_ids <- c(
    feature_ids[1:n_de_nb],
    feature_ids[(n_nb + 1):(n_nb + n_de_zic)],
    feature_ids[(n_nb + n_zic + 1):(n_nb + n_zic + n_de_unit)],
    feature_ids[(n_nb + n_zic + n_unit + 1):(n_nb + n_zic + n_unit + n_de_pos)],
    feature_ids[(n_nb + n_zic + n_unit + n_pos + 1):(n_nb + n_zic + n_unit + n_pos + n_de_real)]
  )
  de_flag <- feature_ids %in% de_ids

  dbg("Samples: %d | Features: %d", n_samples, n_features)
  dbg("DE features injected: %d", sum(de_flag))

  # Helper indices
  idx_A <- which(condition == "A")
  idx_B <- which(condition == "B")

  # --- NB counts (RNA-seq like) ---
  # Baseline mean varies per feature; dispersion varies.
  nb_base_mu <- rgamma(n_nb, shape = 2, rate = 0.2)  # mean ~10
  nb_size <- runif(n_nb, min = 3, max = 10)          # overdispersion

  nb_mat <- matrix(0, nrow = n_nb, ncol = n_samples)
  for (i in seq_len(n_nb)) {
    mu_A <- nb_base_mu[i] * exp(0.15 * age[idx_A])  # mild covariate effect
    mu_B <- nb_base_mu[i] * exp(0.15 * age[idx_B])

    # Inject DE: fold-change ~ 1.5x in B for DE subset
    if (i <= n_de_nb) mu_B <- mu_B * 1.5

    nb_mat[i, idx_A] <- rnbinom(length(idx_A), size = nb_size[i], mu = mu_A)
    nb_mat[i, idx_B] <- rnbinom(length(idx_B), size = nb_size[i], mu = mu_B)
  }

  # --- Zero-inflated counts (metagenomics-like) ---
  # base NB + dropout zeros
  zic_base_mu <- rgamma(n_zic, shape = 1.5, rate = 0.15)  # mean ~10
  zic_size <- runif(n_zic, 2, 8)
  zic_mat <- matrix(0, nrow = n_zic, ncol = n_samples)

  for (i in seq_len(n_zic)) {
    mu_A <- zic_base_mu[i] * exp(0.10 * age[idx_A])
    mu_B <- zic_base_mu[i] * exp(0.10 * age[idx_B])

    # DE: B lower abundance (0.6x)
    if (i <= n_de_zic) mu_B <- mu_B * 0.6

    xA <- rnbinom(length(idx_A), size = zic_size[i], mu = mu_A)
    xB <- rnbinom(length(idx_B), size = zic_size[i], mu = mu_B)

    # dropout probability: batch-dependent + mild condition shift for DE
    p0_A <- pmin(pmax(0.15 + as.numeric(batch[idx_A]) * 0.03, 0.05), 0.6)
    p0_B <- pmin(pmax(0.15 + as.numeric(batch[idx_B]) * 0.03, 0.05), 0.6)
    if (i <= n_de_zic) p0_B <- pmin(p0_B + 0.10, 0.8)  # more zeros in B for DE

    xA[runif(length(xA)) < p0_A] <- 0
    xB[runif(length(xB)) < p0_B] <- 0

    zic_mat[i, idx_A] <- xA
    zic_mat[i, idx_B] <- xB
  }

  # --- Unit interval (beta-like proportions) ---
  unit_mat <- matrix(NA_real_, nrow = n_unit, ncol = n_samples)
  for (i in seq_len(n_unit)) {
    # baseline proportion depends on age + batch
    eta_A <- -0.2 + 0.3 * age[idx_A] + (as.numeric(batch[idx_A]) - 2) * 0.1
    eta_B <- -0.2 + 0.3 * age[idx_B] + (as.numeric(batch[idx_B]) - 2) * 0.1

    # DE: shift mean upward in B
    if (i <= n_de_unit) eta_B <- eta_B + 0.6

    mu_A <- plogis(eta_A)
    mu_B <- plogis(eta_B)

    # concentration: moderate
    phi <- runif(1, 10, 30)
    aA <- mu_A * phi; bA <- (1 - mu_A) * phi
    aB <- mu_B * phi; bB <- (1 - mu_B) * phi

    unit_mat[i, idx_A] <- rbeta(length(idx_A), aA, bA)
    unit_mat[i, idx_B] <- rbeta(length(idx_B), aB, bB)

    # Add sparse 0/1 inflation for a subset (plausible tech artifacts)
    if (i %% 10 == 0) {
      zA <- runif(length(idx_A)) < 0.03
      zB <- runif(length(idx_B)) < 0.03
      unit_mat[i, idx_A][zA] <- 0
      unit_mat[i, idx_B][zB] <- 1
    }
  }

  # --- Positive continuous (gamma/lognormal-like) ---
  pos_mat <- matrix(NA_real_, nrow = n_pos, ncol = n_samples)
  for (i in seq_len(n_pos)) {
    # baseline intensity
    base <- rlnorm(1, meanlog = 2.0, sdlog = 0.3)
    # batch effect multiplicative
    batch_mult_A <- exp((as.numeric(batch[idx_A]) - 2) * 0.15)
    batch_mult_B <- exp((as.numeric(batch[idx_B]) - 2) * 0.15)

    mu_A <- base * batch_mult_A * exp(0.10 * age[idx_A])
    mu_B <- base * batch_mult_B * exp(0.10 * age[idx_B])

    # DE: B higher by 1.4x
    if (i <= n_de_pos) mu_B <- mu_B * 1.4

    # gamma-ish by lognormal noise
    pos_mat[i, idx_A] <- mu_A * rlnorm(length(idx_A), 0, 0.25)
    pos_mat[i, idx_B] <- mu_B * rlnorm(length(idx_B), 0, 0.25)
  }

  # --- Real-valued (approximately normal) ---
  real_mat <- matrix(NA_real_, nrow = n_real, ncol = n_samples)
  for (i in seq_len(n_real)) {
    base <- rnorm(1, 0, 1)
    batch_shift_A <- (as.numeric(batch[idx_A]) - 2) * 0.2
    batch_shift_B <- (as.numeric(batch[idx_B]) - 2) * 0.2

    yA <- base + batch_shift_A + 0.15 * age[idx_A] + rnorm(length(idx_A), 0, 0.7)
    yB <- base + batch_shift_B + 0.15 * age[idx_B] + rnorm(length(idx_B), 0, 0.7)

    # DE: mean shift in B (Cohen ~0.6–0.8)
    if (i <= n_de_real) yB <- yB + 0.7

    real_mat[i, idx_A] <- yA
    real_mat[i, idx_B] <- yB
  }

  counts_matrix <- rbind(nb_mat, zic_mat, unit_mat, pos_mat, real_mat)
  rownames(counts_matrix) <- feature_ids

  # Inject missingness + a couple outliers (realistic)
  # Missingness: ~2% random NAs
  na_mask <- matrix(runif(n_features * n_samples) < 0.02, nrow = n_features)
  counts_matrix[na_mask] <- NA

  # Outliers: a few extreme positives in pos block
  if (n_pos > 0) {
    out_rows <- (n_nb + n_zic + n_unit + 1):(n_nb + n_zic + n_unit + min(3, n_pos))
    out_cols <- sample(seq_len(n_samples), 3)
    counts_matrix[out_rows, out_cols] <- counts_matrix[out_rows, out_cols] * 8
  }

  dbg("Matrix NA rate: %.2f%%", mean(!is.finite(counts_matrix)) * 100)

  # -----------------------------
  # 2) run_perseo: complete pipeline
  # -----------------------------
  dbg_section("2) run_perseo: complete end-to-end pipeline")

  # Keep runtime reasonable in CI: small-ish bootstrap.
  results <- run_perseo(
    counts_matrix = counts_matrix,
    design_matrix = design,
    n_genes = 160,
    n_boot  = 4,
    top_n   = 6,
    verbose = debug_on,
    min_n   = 10,
    seed    = 123,
    group_by_support = FALSE,  # KEY: the 'breakthrough' mode must stay valid
    p_adjust_method = "BH"
  )

  # Validate structure
  expect_s3_class(results, "perseo_results")
  expect_type(results, "list")
  expect_true(all(c("family_selection", "differential_expression", "summary") %in% names(results)))

  # Family selection results
  fams <- results$family_selection
  expect_type(fams, "list")
  expect_true("top_families_overall" %in% names(fams))
  expect_true(length(fams$top_families_overall) > 0)
  expect_true(all(is.character(fams$top_families_overall)))

  dbg("Top families overall: %s", paste(fams$top_families_overall, collapse = ", "))

  # Strong structural invariants
  expect_true("sampled_results" %in% names(fams))
  expect_s3_class(fams$sampled_results, "tbl_df")
  expect_true(all(c("bootstrap","feature","family","skipped","n_valid","support") %in% names(fams$sampled_results)))

  # -----------------------------
  # 3) Differential expression results
  # -----------------------------
  dbg_section("3) Differential expression validation")

  fit <- results$differential_expression

  expect_type(fit, "list")
  expect_true(all(c("selection", "results") %in% names(fit)))
  expect_s3_class(fit$selection, "tbl_df")
  expect_s3_class(fit$results, "tbl_df")

  dbg("Selection rows: %d", nrow(fit$selection))
  dbg("Models rows:    %d", nrow(fit$results))

  # Strong structure checks
  expect_true(all(c("feature","best_family","n_valid_obs","ic_value") %in% names(fit$selection)))
  expect_true(all(c("feature","term","effect","se","pval","padj") %in% names(fit$results)))

  # Sanity checks: no illegal p-values
  expect_true(all(is.na(fit$results$pval) | (fit$results$pval >= 0 & fit$results$pval <= 1)))
  expect_true(all(is.na(fit$results$padj) | (fit$results$padj >= 0 & fit$results$padj <= 1)))

  # Require some non-empty inference
  non_na_p <- fit$results |> dplyr::filter(!is.na(.data$pval))
  expect_true(nrow(non_na_p) > 0)

  # Validate summary
  expect_equal(results$summary$status, "completed")
  expect_equal(results$summary$p_adjust_method, "BH")
  expect_true(results$summary$models_fitted > 0)

  # -----------------------------
  # 4) Scientific validation (robust, not brittle)
  # -----------------------------
  dbg_section("4) Scientific validation: enrichment of DE among top hits (non-brittle)")

  # Focus on the main term: conditionB (or similar)
  # We don't assume exact naming: detect the condition term robustly.
  condition_terms <- non_na_p |> dplyr::filter(grepl("^condition", .data$term))
  expect_true(nrow(condition_terms) > 0)

  # Aggregate per-feature: keep the smallest p-value among condition terms (defensive)
  per_feature_p <- condition_terms |>
    dplyr::group_by(.data$feature) |>
    dplyr::summarise(p = min(.data$pval, na.rm = TRUE),
                     padj = min(.data$p_adj, na.rm = TRUE),
                     .groups = "drop")

  # Attach truth
  truth_tbl <- tibble::tibble(feature = feature_ids, is_de = de_flag)
  per_feature_p <- dplyr::left_join(per_feature_p, truth_tbl, by = "feature")

  # Remove features not present (skipped etc.)
  per_feature_p <- per_feature_p |> dplyr::filter(!is.na(.data$is_de), is.finite(.data$p))

  dbg("Features with valid condition p-values: %d", nrow(per_feature_p))
  dbg("DE among those: %d", sum(per_feature_p$is_de))

  # Non-brittle requirement:
  # DE features should be enriched among top-k hits.
  # Choose k based on available size, but not too small.
  k <- min(30, nrow(per_feature_p))
  topk <- per_feature_p |> dplyr::arrange(.data$p) |> dplyr::slice_head(n = k)

  frac_de_topk <- mean(topk$is_de)
  frac_de_all  <- mean(per_feature_p$is_de)

  dbg("Enrichment check: frac(DE) top-%d = %.3f vs background = %.3f", k, frac_de_topk, frac_de_all)

  # Expect enrichment: top hits should contain more DE than background by a margin.
  # Margin chosen to be strong but robust.
  expect_true(frac_de_topk >= frac_de_all + 0.10)

  # Secondary check: DE p-values tend to be smaller than non-DE (stochastic but stable)
  de_p   <- per_feature_p$p[per_feature_p$is_de]
  nde_p  <- per_feature_p$p[!per_feature_p$is_de]

  if (length(de_p) > 5 && length(nde_p) > 5) {
    dbg("Median p(DE)=%.4f | Median p(non-DE)=%.4f", median(de_p), median(nde_p))
    expect_true(median(de_p) < median(nde_p))
  }

  # -----------------------------
  # 5) Diagnostics dump (only if debug)
  # -----------------------------
  if (debug_on) {
    dbg_section("5) Debug dump: top families, top hits, and basic summaries")

    top_fam_tab <- fams$freq_table_overall
    message("Top family frequency table (overall):")
    print(head(top_fam_tab, 10))

    message("\nTop 15 hits (by p):")
    print(per_feature_p |> dplyr::arrange(.data$p) |> dplyr::slice_head(n = 15))

    message("\nSelection head:")
    print(dplyr::as_tibble(fit$selection) |> dplyr::slice_head(n = 10))

    message("\nResults head:")
    print(dplyr::as_tibble(fit$results) |> dplyr::slice_head(n = 10))
  }
})
