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
  
  # Check contrast coefficients
  expect_equal(C[1, "(Intercept)"], -1)
  expect_equal(C[1, "statusTreatment"], 1)
})
