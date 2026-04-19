test_that("build_contrast_matrix creates pairwise contrasts", {
  metadata <- data.frame(
    tissue = factor(c("A", "A", "B", "B", "C", "C"))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  C <- PERSEO:::build_contrast_matrix("tissue", metadata, design)
  
  expect_true(is.matrix(C))
  expect_equal(nrow(C), 3)  # 3 pairwise comparisons for 3 levels
  expect_equal(ncol(C), ncol(design))
  expect_true(all(rownames(C) %in% c("B_vs_A", "C_vs_A", "C_vs_B")))
})

test_that("build_contrast_matrix validates inputs", {
  metadata <- data.frame(
    tissue = factor(c("A", "A", "B", "B"))
  )
  
  design <- model.matrix(~ tissue, data = metadata)
  
  expect_error(
    PERSEO:::build_contrast_matrix("missing_var", metadata, design),
    "not found in metadata"
  )
  
  metadata_single <- data.frame(
    tissue = factor(c("A", "A", "A", "A"))
  )
  
  expect_error(
    PERSEO:::build_contrast_matrix("tissue", metadata_single, design),
    "must have at least 2 levels"
  )
})

test_that("build_contrast_matrix handles non-factor variables", {
  metadata <- data.frame(
    group = c("X", "X", "Y", "Y", "Z", "Z")
  )
  
  design <- model.matrix(~ group, data = metadata)
  
  C <- PERSEO:::build_contrast_matrix("group", metadata, design)
  
  expect_true(is.matrix(C))
  expect_equal(nrow(C), 3)
})

test_that("build_contrast_matrix creates correct contrast structure", {
  metadata <- data.frame(
    status = factor(c(rep("Control", 5), rep("Treatment", 5)))
  )

  design <- model.matrix(~ status, data = metadata)

  C <- PERSEO:::build_contrast_matrix("status", metadata, design)

  expect_equal(nrow(C), 1)  # Only 1 pairwise comparison
  expect_equal(rownames(C), "Treatment_vs_Control")

  # In treatment coding, Treatment_vs_Control = beta_Treatment only.
  # The intercept (reference mean) cancels: mu_T - mu_C = beta_T.
  expect_equal(C[1, "(Intercept)"], 0)
  expect_equal(C[1, "statusTreatment"], 1)
})

test_that("build_contrast_matrix produces correct vectors for 3 groups", {
  metadata <- data.frame(
    cond = factor(rep(c("A", "B", "C"), each = 4))
  )
  design <- model.matrix(~ cond, data = metadata)
  # colnames(design): "(Intercept)", "condB", "condC"

  C <- PERSEO:::build_contrast_matrix("cond", metadata, design)

  expect_equal(nrow(C), 3)
  expect_setequal(rownames(C), c("B_vs_A", "C_vs_A", "C_vs_B"))

  # B_vs_A: mu_B - mu_A = beta_B  →  [0, 1, 0]
  expect_equal(unname(C["B_vs_A", ]), c(0, 1, 0))
  # C_vs_A: mu_C - mu_A = beta_C  →  [0, 0, 1]
  expect_equal(unname(C["C_vs_A", ]), c(0, 0, 1))
  # C_vs_B: mu_C - mu_B = beta_C - beta_B  →  [0, -1, 1]
  expect_equal(unname(C["C_vs_B", ]), c(0, -1, 1))
})
