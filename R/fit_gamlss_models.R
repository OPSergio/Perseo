#' Fit GAMLSS models per feature (family selection with common mask + Jacobian)
#'
#' Compares candidate GAMLSS families per feature on the same set of observations
#' (common mask), using family-specific transformations and Jacobian correction 
#' for fair IC comparison. Selects the best family by GAIC/BIC/AIC and returns 
#' per-term statistics from `summary(fit, what = "mu")`.
#'
#' @section Input Workflows:
#' 
#' \strong{Workflow A: Formula + Metadata + Automatic Contrasts}
#' 
#' When you want PERSEO to build the design matrix and automatically generate 
#' pairwise contrasts:
#' 
#' \preformatted{
#' fit_gamlss_models(
#'   counts_matrix = counts,
#'   design_matrix = "~ condition + batch",  # formula string
#'   metadata = metadata,                     # REQUIRED
#'   candidate_families = c("NBI", "GG"),
#'   contrast_variable = "condition"          # auto-generate contrasts
#' )
#' }
#' 
#' \strong{Requirements:}
#' \itemize{
#'   \item \code{design_matrix} MUST be a formula (string or formula object)
#'   \item \code{metadata} MUST be provided (data.frame with nrow = ncol(counts_matrix))
#'   \item \code{contrast_variable} names a column in metadata
#'   \item Row names of metadata should match column names of counts_matrix
#' }
#' 
#' \strong{Workflow B: Design Matrix + Manual Contrasts}
#' 
#' When you have a pre-built design matrix and want explicit control over contrasts:
#' 
#' \preformatted{
#' design <- model.matrix(~ condition + batch, data = metadata)
#' C <- matrix(...)  # Custom contrast matrix
#' colnames(C) <- colnames(design)
#' 
#' fit_gamlss_models(
#'   counts_matrix = counts,
#'   design_matrix = design,        # pre-built matrix
#'   candidate_families = c("NBI", "GG"),
#'   contrast_matrix = C             # REQUIRED for contrasts
#' )
#' }
#' 
#' \strong{Requirements:}
#' \itemize{
#'   \item \code{design_matrix} is a numeric matrix (nrow = ncol(counts_matrix))
#'   \item If contrasts are needed, \code{contrast_matrix} MUST be provided
#'   \item \code{contrast_variable} is NOT allowed (automatic generation not supported)
#'   \item Column names of contrast_matrix must match coefficient names
#' }
#' 
#' \strong{Invalid Combinations:}
#' \itemize{
#'   \item design_matrix = matrix + contrast_variable (without metadata)
#'   \item Automatic contrast generation requires formula + metadata
#' }
#'
#' @param counts_matrix Numeric matrix (features × samples). Feature IDs in row names, 
#'   sample IDs in column names.
#' @param design_matrix Either:
#'   \itemize{
#'     \item A formula string (e.g., \code{"~ condition + batch"}) - requires \code{metadata}
#'     \item A formula object - requires \code{metadata}
#'     \item A pre-built model matrix from \code{model.matrix()} with nrow = ncol(counts_matrix)
#'   }
#' @param metadata Data.frame with sample metadata. Required when:
#'   \itemize{
#'     \item \code{design_matrix} is a formula (string or object), OR
#'     \item Using \code{contrast_variable} for automatic contrast generation
#'   }
#'   Must have \code{nrow = ncol(counts_matrix)}. Row names should match column 
#'   names of counts_matrix for proper alignment.
#' @param candidate_families Character vector of GAMLSS families to test
#'   (e.g., \code{c("PO","NBI","GA","GG","IG","LOGNO","NO","TF")}).
#' @param criterion Model selection criterion: \code{"GAIC"}, \code{"BIC"}, or 
#'   \code{"AIC"}. Default \code{"GAIC"}.
#' @param gaic_k Numeric penalty for GAIC. If \code{NULL}, uses \code{log(n_valid_obs)} 
#'   (equivalent to BIC). Ignored if \code{criterion != "GAIC"}.
#' @param min_n Integer minimum number of valid observations after common mask. 
#'   Features below this threshold are skipped. Default \code{5}.
#' @param p_adjust Method passed to \code{p.adjust()} to compute FDR per term 
#'   across features. Default \code{"BH"}.
#' @param contrast_matrix Numeric matrix where each row defines a linear combination 
#'   of mu coefficients. Column names must match coefficient names from design matrix. 
#'   Use this for custom contrasts with pre-built design matrices (Workflow B). 
#'   Mutually exclusive with \code{contrast_variable}.
#' @param contrast_variable Character string naming a categorical variable in 
#'   \code{metadata} for automatic pairwise contrast generation. Only valid with 
#'   formula-based design (Workflow A). Requires \code{metadata}. Mutually 
#'   exclusive with \code{contrast_matrix}.
#' @param omnibus Logical; if \code{TRUE} and contrasts are requested, performs a 
#'   feature-level omnibus test before computing contrasts. Only features passing 
#'   the omnibus test (p < \code{omnibus_threshold}) will have contrasts computed. 
#'   This implements hierarchical testing (omnibus → post-hoc). Default \code{FALSE}.
#' @param omnibus_threshold Numeric; significance threshold for omnibus test. 
#'   Default \code{0.05}.
#' @param omnibus_test Character; type of omnibus test: \code{"Wald"} (fast, uses 
#'   vcov from fitted model) or \code{"LRT"} (robust, refits reduced model). 
#'   Default \code{"Wald"}.
#' @param workers Integer number of parallel workers. Only used when 
#'   \code{parallel = TRUE}. If \code{NULL}, uses \code{future::availableCores() - 1}. 
#'   Default \code{NULL}.
#' @param parallel Logical; enable parallel processing via \code{future::plan(multisession)}. 
#'   When \code{TRUE}, automatically configures and cleans up the future backend. 
#'   Default \code{FALSE}.
#' @param show_progress Logical; show progress messages and final report. 
#'   Default \code{TRUE}.
#' @param progress_label Character label for progress bar.
#' @param transform_mode Character: \code{"strict"} (default) or \code{"safe"}. 
#'   Transformation mode for family comparison. See Details.
#'
#' @return A list with:
#'   \describe{
#'     \item{results}{Tibble with one row per feature × coefficient. Columns: feature, 
#'       term, effect, se, stat, pval, padj. FDR adjustment (padj) is performed 
#'       \strong{by term} across all features.}
#'     \item{selection}{Tibble with one row per feature. Columns: feature, best_family, 
#'       n_valid_obs, ic_value (Jacobian-corrected), transform_mode.}
#'     \item{omnibus}{(Optional) Tibble with omnibus test results. Only present when 
#'       \code{omnibus = TRUE} and contrasts are requested. Columns: feature, family, 
#'       test_type, statistic, df, p_value, pass.}
#'     \item{contrasts}{(Optional) Tibble with contrast results. Only present if 
#'       \code{contrast_matrix} or \code{contrast_variable} is provided. Columns: 
#'       feature, family, contrast, estimate, se, z, p_value, p_adj. FDR adjustment 
#'       (p_adj) is performed \strong{by contrast} across features. When 
#'       \code{omnibus = TRUE}, only includes features where \code{omnibus$pass == TRUE}.}
#'   }
#'
#' @details
#' 
#' \subsection{Model Selection Strategy}{
#'   For each feature, PERSEO fits all candidate families on the \strong{same 
#'   observations} (common mask) using family-specific transformations and Jacobian 
#'   correction. This ensures information criteria (BIC/AIC/GAIC) are comparable 
#'   across families on the original data scale. The best family is selected per 
#'   feature, and coefficient statistics are extracted from \code{summary(fit, what = "mu")}.
#' }
#' 
#' \subsection{Transformation Modes}{
#'   \describe{
#'     \item{\strong{strict}}{(Default) Domain-preserving transformations. Invalid 
#'       observations are excluded via masking. Conservative approach that respects 
#'       theoretical support constraints.}
#'     \item{\strong{safe}}{Global affine transformations (y* = a·y + b) that are 
#'       invertible with Jacobian correction. No observation exclusion. Useful for 
#'       exploratory analysis across support boundaries.}
#'   }
#' }
#' 
#' \subsection{Hierarchical Testing (Omnibus + Contrasts)}{
#'   When \code{omnibus = TRUE} and contrasts are requested, the function implements 
#'   a two-stage hierarchical testing strategy:
#'   \enumerate{
#'     \item \strong{Omnibus Test}: For each feature, test whether the target factor 
#'       has \emph{any} effect (multi-df test).
#'     \item \strong{Contrasts}: Only for features passing the omnibus test 
#'       (p < \code{omnibus_threshold}), compute pairwise contrasts.
#'   }
#'   
#'   This mirrors the ANOVA + post-hoc workflow and reduces the multiple testing 
#'   burden, increasing power to detect specific pairwise differences.
#'   
#'   \strong{Important}: The omnibus test acts only as a gatekeeper. Contrast p-values 
#'   are still adjusted \strong{by contrast} across features using the standard 
#'   FDR procedure (\code{group_by(contrast) |> mutate(p_adj = p.adjust(p_value))}). 
#'   The omnibus p-value is NOT incorporated into contrast adjustment.
#'   
#'   When \code{omnibus = FALSE} (default), all features receive contrast computations 
#'   without filtering.
#' }
#' 
#' \subsection{Omnibus Test Options}{
#'   \describe{
#'     \item{\strong{Wald}}{(Default) Fast multi-df Wald test using variance-covariance 
#'       matrix: W = β'V⁻¹β ~ χ²(df). No model refitting required. Best for large 
#'       samples (n > 30/group) and standard families.}
#'     \item{\strong{LRT}}{Likelihood ratio test requiring reduced model refit: 
#'       LRT = 2(logLik_full - logLik_reduced) ~ χ²(df). More robust, especially 
#'       for small samples (n < 30/group) or complex families. Approximately 2× 
#'       slower than Wald.}
#'   }
#' }
#' 
#' \subsection{Parallel Processing}{
#'   When \code{parallel = TRUE}, the function automatically configures 
#'   \code{future::plan(multisession)} with the specified number of workers and 
#'   resets to sequential plan upon completion. No manual setup required.
#' }
#'
#' @examples
#' \dontrun{
#' # ============================================================================
#' # Workflow A: Formula + Metadata + Automatic Contrasts
#' # ============================================================================
#' 
#' # Prepare data
#' metadata <- data.frame(
#'   sample_id = colnames(counts_matrix),
#'   tissue_type = factor(c(rep("Normal", 50), rep("Tumor", 50))),
#'   age = rnorm(100, 60, 10),
#'   batch = factor(rep(1:4, each = 25))
#' )
#' 
#' # Simple: formula + auto-contrasts
#' fit_auto <- fit_gamlss_models(
#'   counts_matrix = counts_matrix,
#'   design_matrix = "~ tissue_type + age + batch",  # formula string
#'   metadata = metadata,                             # REQUIRED
#'   candidate_families = c("NBI", "GG", "LOGNO"),
#'   contrast_variable = "tissue_type",               # auto-generate all pairwise
#'   criterion = "BIC",
#'   parallel = TRUE,
#'   workers = 4
#' )
#' 
#' # View results
#' head(fit_auto$results)    # Coefficient statistics
#' head(fit_auto$selection)  # Best family per feature
#' head(fit_auto$contrasts)  # Pairwise contrasts (auto-generated)
#' 
#' # ============================================================================
#' # Workflow B: Design Matrix + Manual Contrasts
#' # ============================================================================
#' 
#' # Build design matrix manually
#' metadata$tissue_type <- relevel(metadata$tissue_type, ref = "Normal")
#' design <- model.matrix(~ tissue_type + age + batch, data = metadata)
#' 
#' # Define custom contrast matrix
#' C <- matrix(0, nrow = 1, ncol = ncol(design))
#' colnames(C) <- colnames(design)
#' rownames(C) <- "Tumor_vs_Normal"
#' C["Tumor_vs_Normal", "tissue_typeTumor"] <- 1
#' 
#' # Fit with explicit contrast matrix
#' fit_manual <- fit_gamlss_models(
#'   counts_matrix = counts_matrix,
#'   design_matrix = design,                          # pre-built matrix
#'   candidate_families = c("NBI", "GG", "LOGNO"),
#'   contrast_matrix = C,                             # explicit contrasts
#'   criterion = "BIC",
#'   parallel = TRUE,
#'   workers = 4
#' )
#' 
#' head(fit_manual$contrasts)
#' 
#' # ============================================================================
#' # Hierarchical Testing: Omnibus + Contrasts
#' # ============================================================================
#' 
#' # For multi-level factors (3+ levels), use omnibus filtering
#' metadata_multi <- metadata
#' metadata_multi$tissue_type <- factor(
#'   rep(c("Normal", "Luminal", "Basal"), length.out = 100)
#' )
#' 
#' fit_hierarchical <- fit_gamlss_models(
#'   counts_matrix = counts_matrix,
#'   design_matrix = "~ tissue_type + age + batch",
#'   metadata = metadata_multi,
#'   candidate_families = c("NBI", "GG"),
#'   contrast_variable = "tissue_type",
#'   omnibus = TRUE,                 # Enable hierarchical testing
#'   omnibus_threshold = 0.05,       # Filter threshold
#'   omnibus_test = "Wald",          # "Wald" or "LRT"
#'   criterion = "BIC",
#'   parallel = TRUE,
#'   workers = 4
#' )
#' 
#' # View omnibus results
#' head(fit_hierarchical$omnibus)
#' 
#' # Contrasts only computed for features where omnibus$pass == TRUE
#' head(fit_hierarchical$contrasts)
#' 
#' # Count features passing omnibus filter
#' sum(fit_hierarchical$omnibus$pass, na.rm = TRUE)
#' }
#' @seealso \code{\link{find_families}}, \code{\link{transform_response}}, 
#'   \code{\link{run_perseo}}
#' @export
fit_gamlss_models <- function(counts_matrix,
                              design_matrix,
                              metadata = NULL,
                              candidate_families,
                              criterion = c("GAIC","BIC","AIC"),
                              gaic_k = NULL,
                              min_n = 5,
                              contrast_matrix = NULL,
                              contrast_variable = NULL,
                              omnibus = FALSE,
                              omnibus_threshold = 0.05,
                              omnibus_test = c("Wald", "LRT"),
                              p_adjust = "BH",
                              workers = NULL,
                              parallel = FALSE,
                              show_progress = TRUE,
                              progress_label = "Fitting features",
                              transform_mode = "strict",
                              group_by_support = FALSE,
                              filter_beta_inflated = TRUE,
                              thr_zero = 0.005,
                              thr_one = 0.005) {
  # ---- Input validation ----
  criterion <- match.arg(criterion)
  transform_mode <- match.arg(transform_mode, c("strict", "safe"))
  omnibus_test <- match.arg(omnibus_test)
  
  validate_counts_matrix(counts_matrix)
  validate_criterion_args(criterion, gaic_k)
  
  # ---- Set up parallel backend if requested ----
  original_plan <- NULL
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required for parallel processing. Install it with install.packages('future')")
    }
    
    original_plan <- future::plan()
    
    n_workers <- if (is.null(workers)) {
      max(1, future::availableCores() - 1)
    } else {
      as.integer(workers)
    }
    
    if (n_workers < 1) {
      warning("workers must be >= 1; using workers = 1")
      n_workers <- 1L
    }
    
    future::plan(future::multisession, workers = n_workers)
    
    if (show_progress) {
      message("Parallel processing enabled with ", n_workers, " workers")
    }
  }
  
  on.exit({
    if (parallel && !is.null(original_plan)) {
      future::plan(original_plan)
    }
  }, add = TRUE)
  
  # Store variables for final report (will be shown at the very end)
  report_vars <- new.env()
  
  feature_ids <- rownames(counts_matrix)
  if (is.null(feature_ids)) {
    feature_ids <- seq_len(nrow(counts_matrix))
  }
  
  # ---- Track original design_matrix input type for validation ----
  is_formula_input <- (is.character(design_matrix) && length(design_matrix) == 1) || inherits(design_matrix, "formula")

  formula_mu <- NULL
  design_subset <- NULL
  metadata_subset <- NULL

  # ---- Handle formula vs matrix input ----
  if (is_formula_input) {
    if (is.null(metadata)) {
      stop("metadata must be provided when design_matrix is a formula")
    }
    if (!is.data.frame(metadata) || nrow(metadata) != ncol(counts_matrix)) {
      stop("metadata must be a data.frame with nrow = ncol(counts_matrix)")
    }

    formula_obj <- if (is.character(design_matrix)) {
      stats::as.formula(design_matrix)
    } else {
      design_matrix
    }

    terms_obj <- stats::terms(formula_obj)
    term_labels <- attr(terms_obj, "term.labels")
    has_intercept <- attr(terms_obj, "intercept")
    formula_mu <- stats::reformulate(
      termlabels = term_labels,
      response = "y",
      intercept = has_intercept
    )
    environment(formula_mu) <- environment(formula_obj)

    formula_vars <- intersect(all.vars(formula_obj), names(metadata))
    metadata_for_cc <- if (length(formula_vars) > 0) {
      metadata[, formula_vars, drop = FALSE]
    } else {
      data.frame(.row_id = seq_len(nrow(metadata)))
    }

    complete_rows <- stats::complete.cases(metadata_for_cc)
    counts_subset <- counts_matrix[, complete_rows, drop = FALSE]
    metadata_subset <- metadata[complete_rows, , drop = FALSE]
  } else {
    design_matrix <- validate_design_matrix(design_matrix, ncol(counts_matrix))
    complete_rows <- stats::complete.cases(design_matrix)
    counts_subset <- counts_matrix[, complete_rows, drop = FALSE]
    design_subset <- design_matrix[complete_rows, , drop = FALSE]
  }
  
  # ---- Validate contrast input combinations (CRITICAL) ----
  # Enforce clear contract for automatic contrast generation
  if (!is.null(contrast_variable) && is.null(contrast_matrix)) {
    # User wants automatic contrast generation via contrast_variable
    # This REQUIRES that we can build the design matrix from metadata
    
    if (!is_formula_input && is.null(metadata)) {
      stop(
        "Invalid input combination:\n",
        "  When contrast_variable is provided and contrast_matrix is NULL,\n",
        "  you must EITHER:\n",
        "    1) Provide design_matrix as a formula (string or formula object) with metadata, OR\n",
        "    2) Provide both a design matrix AND metadata (ensuring row alignment).\n",
        "\n",
        "  Current inputs:\n",
        "    - design_matrix: ", if (is_formula_input) "formula" else "numeric matrix", "\n",
        "    - metadata: NULL\n",
        "    - contrast_variable: '", contrast_variable, "'\n",
        "    - contrast_matrix: NULL\n",
        "\n",
        "  Recommended workflows:\n",
        "    Workflow A (Formula + automatic contrasts):\n",
        "      fit_gamlss_models(\n",
        "        counts_matrix = counts,\n",
        "        design_matrix = \"~ condition + batch\",  # formula string\n",
        "        metadata = metadata,                     # required\n",
        "        contrast_variable = \"condition\",        # auto-generate contrasts\n",
        "        ...\n",
        "      )\n",
        "\n",
        "    Workflow B (Design matrix + manual contrasts):\n",
        "      design <- model.matrix(~ condition + batch, data = metadata)\n",
        "      C <- build_contrast_matrix(...)  # or manually create\n",
        "      fit_gamlss_models(\n",
        "        counts_matrix = counts,\n",
        "        design_matrix = design,      # pre-built matrix\n",
        "        contrast_matrix = C,         # explicit specification\n",
        "        ...\n",
        "      )\n"
      )
    }
    
    # If metadata is provided but design was a matrix, validate row alignment
    if (!is_formula_input && !is.null(metadata)) {
      if (nrow(metadata) != nrow(design_matrix)) {
        stop(
          "Row count mismatch:\n",
          "  design_matrix has ", nrow(design_matrix), " rows\n",
          "  metadata has ", nrow(metadata), " rows\n",
          "  These must match when using contrast_variable with a pre-built design matrix."
        )
      }
    }
  }
  
  # ---- Handle contrast inputs with priority logic ----
  final_contrast_matrix <- NULL
  final_contrast_variable <- NULL
  
  if (!is.null(contrast_matrix) && !is.null(contrast_variable)) {
    message("Both contrast_matrix and contrast_variable provided. ",
            "Prioritizing contrast_matrix.")
    final_contrast_matrix <- contrast_matrix
    final_contrast_variable <- NULL  # Will identify coefs from matrix
  } else if (!is.null(contrast_matrix)) {
    final_contrast_matrix <- contrast_matrix
    final_contrast_variable <- NULL
  } else if (!is.null(contrast_variable)) {
    if (is.null(metadata)) {
      stop("metadata must be provided when using contrast_variable")
    }
    if (is_formula_input) {
      final_contrast_matrix <- NULL
    } else {
      final_contrast_matrix <- build_contrast_matrix(
        contrast_variable, metadata[complete_rows, , drop = FALSE], design_subset
      )
    }
    final_contrast_variable <- contrast_variable
  }

  has_contrasts_requested <- !is.null(final_contrast_matrix) || !is.null(final_contrast_variable)
  
  # ---- Validate contrast matrix if present ----
  if (!is.null(final_contrast_matrix)) {
    if (!is.matrix(final_contrast_matrix) || !is.numeric(final_contrast_matrix)) {
      stop("contrast_matrix must be a numeric matrix")
    }
    if (is.null(colnames(final_contrast_matrix))) {
      stop("contrast_matrix must have column names matching coefficient names")
    }
  }
  
  # ---- Print startup banner with cli ----
  if (show_progress) {
    # Check if cli is available, fallback to base messages
    has_cli <- requireNamespace("cli", quietly = TRUE)
    
    if (has_cli) {
      cli::cli_rule("fit_gamlss_models() - STARTING")
      cli::cli_alert_info("Total features    : {.val {length(feature_ids)}}")
      cli::cli_alert_info("Samples           : {.val {ncol(counts_subset)}}")
      cli::cli_alert_info("Candidate families: {.val {paste(candidate_families, collapse = ', ')}}")
      cli::cli_alert_info("Criterion         : {.val {criterion}}")
      cli::cli_alert_info("Transform mode    : {.val {transform_mode}}")
      cli::cli_alert_info("P-value adjustment: {.val {p_adjust}}")
      if (parallel) {
        cli::cli_alert_info("Parallel workers  : {.val {n_workers}}")
      }
      if (has_contrasts_requested) {
        n_contr <- if (is.null(final_contrast_matrix)) "auto" else nrow(final_contrast_matrix)
        cli::cli_alert_info("Contrasts requested: {.val {n_contr}} contrasts")
        if (omnibus) {
          cli::cli_alert_info("  Omnibus test    : {.val {omnibus_test}} (threshold: {.val {omnibus_threshold}})")
        }
      }
      cli::cli_rule()
    } else {
      # Fallback to base messages
      message("\n", strrep("=", 80))
      message("fit_gamlss_models() - STARTING")
      message(strrep("=", 80))
      message("Total features    : ", length(feature_ids))
      message("Samples           : ", ncol(counts_subset))
      message("Candidate families: ", paste(candidate_families, collapse = ", "))
      message("Criterion         : ", criterion)
      message("Transform mode    : ", transform_mode)
      message("P-value adjustment: ", p_adjust)
      if (parallel) {
        message("Parallel workers  : ", n_workers)
      }
      if (has_contrasts_requested) {
        n_contr <- if (is.null(final_contrast_matrix)) "auto" else as.character(nrow(final_contrast_matrix))
        message("Contrasts requested: ", n_contr, " contrasts")
        if (omnibus) {
          message("  Omnibus test    : ", omnibus_test, " (threshold: ", omnibus_threshold, ")")
        }
      }
      message(strrep("=", 80), "\n")
    }
  }
  
  # ---- Capture internal functions for parallel workers ----
  # Workers don't have access to PERSEO::: unless functions are captured locally
  .has_insufficient_variation <- has_insufficient_variation
  .is_binomial_family <- is_binomial_family
  .infer_binomial_denominator <- infer_binomial_denominator
  .filter_candidate_families <- filter_candidate_families
  .compare_families_with_design <- compare_families_with_design
  .transform_response <- transform_response
  .instantiate_gamlss_family <- instantiate_gamlss_family
  .fit_gamlss_safely <- fit_gamlss_safely
  .extract_mu_coefficients <- extract_mu_coefficients
  .identify_factor_coefficients <- identify_factor_coefficients
  .wald_omnibus_test <- wald_omnibus_test
  .lrt_omnibus_test <- lrt_omnibus_test
  .apply_contrasts <- apply_contrasts

  # ---- Worker function: process one feature ----
  process_one_feature <- function(feature_name, feature_vec, design_mat, metadata_df, prog) {
    if (!is.null(prog)) prog()

    if (.has_insufficient_variation(feature_vec)) {
      return(NULL)
    }

    # Per-feature support filtering: restrict families to those compatible with
    # this feature's empirical support before computing the common mask
    feature_families <- if (group_by_support) {
      filtered <- .filter_candidate_families(
        feature_vec          = feature_vec,
        candidate_families   = candidate_families,
        group_by_support     = TRUE,
        filter_beta_inflated = filter_beta_inflated,
        thr_zero             = thr_zero,
        thr_one              = thr_one
      )
      filtered$families_to_test
    } else {
      candidate_families
    }

    if (length(feature_families) == 0) {
      return(NULL)
    }

    bd_vec <- if (any(.is_binomial_family(feature_families))) {
      .infer_binomial_denominator(feature_vec)
    } else NULL

    family_results <- .compare_families_with_design(
      feature_vec,
      design_mat,
      families = feature_families,
      criterion = criterion,
      gaic_k = gaic_k,
      min_n = min_n,
      bd_vec = bd_vec,
      transform_mode = transform_mode,
      extract_coef_fn = .extract_mu_coefficients,
      formula_mu = formula_mu,
      metadata_df = metadata_df
    )
    
    if (nrow(family_results$comparisons) == 0) {
      return(NULL)
    }
    
    cmp <- family_results$comparisons
    best <- cmp[select_by_ic_gof(cmp$ic_value, cmp$gof), ]
    
    # Extract coefficient table
    coef_df <- tibble::tibble(
      feature = feature_name,
      term    = best$coef_tbl[[1]]$term,
      effect  = best$coef_tbl[[1]]$effect,
      se      = best$coef_tbl[[1]]$se,
      stat    = best$coef_tbl[[1]]$stat,
      pval    = best$coef_tbl[[1]]$pval
    )
    
    # Initialize outputs
    omnibus_df <- NULL
    contrast_df <- NULL
    
    # Only perform omnibus/contrast if requested
    if (has_contrasts_requested) {
      # Re-fit to get beta and vcov for omnibus/contrasts
      tr <- .transform_response(feature_vec, best$family, mode = transform_mode)
      z <- tr$y[tr$mask]
      fit_data <- if (is_formula_input) {
        cbind.data.frame(y = z, metadata_df[tr$mask, , drop = FALSE])
      } else {
        cbind.data.frame(y = z, design_mat[tr$mask, , drop = FALSE])
      }
      
      fam_obj <- .instantiate_gamlss_family(
        best$family,
        bd_vec = if (.is_binomial_family(best$family)) {
          .infer_binomial_denominator(feature_vec, tr$mask)
        } else NULL
      )
      
      best_fit <- .fit_gamlss_safely(
        if (is_formula_input) formula_mu else y ~ .,
        data = fit_data,
        family_obj = fam_obj
      )
      
      if (!is.null(best_fit)) {
        beta <- tryCatch({ coef(best_fit, what = "mu") }, error = function(e) NULL)
        
        # Extract vcov
        V <- NULL
        if (!is.null(best_fit$family_obj_stored)) {
          family_obj <- best_fit$family_obj_stored
        } else {
          family_obj <- fam_obj
        }
        
        V <- tryCatch({ vcov(best_fit, what = "mu") }, error = function(e) NULL)
        
        if (is.null(V)) {
          V <- tryCatch({
            suppressWarnings({
              s <- summary(best_fit, what = "mu")
              if (is.matrix(s) && "Std. Error" %in% colnames(s)) {
                se_vals <- s[, "Std. Error"]
                n_mu <- length(beta)
                if (length(se_vals) >= n_mu) {
                  se_vals <- se_vals[1:n_mu]
                }
                V_diag <- diag(se_vals^2)
                rownames(V_diag) <- colnames(V_diag) <- names(beta)
                V_diag
              } else {
                NULL
              }
            })
          }, error = function(e) NULL)
        }
        
        if (is.null(V) && !is.null(best_fit$vcov.mu)) {
          V <- best_fit$vcov.mu
        }
        
        # Determine which coefficients belong to the target factor
        contrast_matrix_current <- final_contrast_matrix
        if (is.null(contrast_matrix_current) && !is.null(final_contrast_variable) && !is.null(beta)) {
          contrast_matrix_current <- tryCatch({
            coef_name_stub <- matrix(0, nrow = nrow(metadata_df), ncol = length(beta))
            colnames(coef_name_stub) <- names(beta)
            build_contrast_matrix(final_contrast_variable, metadata_df, coef_name_stub)
          }, error = function(e) NULL)
        }

        factor_coefs <- NULL
        if (!is.null(beta)) {
          factor_coefs <- .identify_factor_coefficients(
            contrast_matrix = contrast_matrix_current,
            contrast_variable = final_contrast_variable,
            coef_names = names(beta)
          )
        }
        
        # Perform omnibus test if requested
        omnibus_pass <- TRUE  # Default: compute contrasts
        
        if (omnibus && !is.null(factor_coefs) && length(factor_coefs) > 0) {
          omnibus_result <- if (omnibus_test == "Wald" && !is.null(V)) {
            .wald_omnibus_test(beta, V, factor_coefs)
          } else if (omnibus_test == "LRT") {
            .lrt_omnibus_test(
              feature_vec = feature_vec,
              design_mat = design_mat,
              factor_coefs = factor_coefs,
              family = best$family,
              transform_mode = transform_mode,
              bd_vec = bd_vec,
              formula_mu = if (is_formula_input) formula_mu else NULL,
              metadata_df = if (is_formula_input) metadata_df else NULL,
              contrast_variable = final_contrast_variable
            )
          } else {
            list(statistic = NA_real_, df = NA_integer_, p_value = NA_real_, test_type = omnibus_test)
          }
          
          omnibus_df <- tibble::tibble(
            feature = feature_name,
            family = best$family,
            test_type = omnibus_result$test_type,
            statistic = omnibus_result$statistic,
            df = omnibus_result$df,
            p_value = omnibus_result$p_value,
            pass = !is.na(omnibus_result$p_value) && omnibus_result$p_value < omnibus_threshold
          )
          
          omnibus_pass <- omnibus_df$pass
        }
        
        # Compute contrasts only if omnibus passed (or omnibus disabled)
        if (omnibus_pass && !is.null(beta) && !is.null(V) && all(is.finite(V)) && !is.null(contrast_matrix_current)) {
          contrast_df <- tryCatch({
            contrast_res <- .apply_contrasts(beta, V, contrast_matrix_current)
            contrast_res$feature <- feature_name
            contrast_res$family <- best$family
            contrast_res
          }, error = function(e) {
            # Fit failed
            contrast_names <- if (!is.null(rownames(contrast_matrix_current))) {
              rownames(contrast_matrix_current)
            } else {
              paste0("contrast_", seq_len(nrow(contrast_matrix_current)))
            }
            tibble::tibble(
              feature = feature_name,
              family = best$family,
              contrast = contrast_names,
              estimate = NA_real_,
              se = NA_real_,
              z = NA_real_,
              p_value = NA_real_
            )
          })
        }
      }
    }
    
    list(
      feature  = feature_name,
      best_fam = best$family,
      n_valid  = best$n_valid_obs,
      ic_value = best$ic_value,
      coef_df  = coef_df,
      omnibus_df = omnibus_df,
      contrast_df = contrast_df
    )
  }
  
  # ---- Execute fitting (sequential or parallel) ----
  has_cli <- requireNamespace("cli", quietly = TRUE)
  n_features <- nrow(counts_subset)
  
  if (parallel) {
    # Parallel execution
    out_list <- future.apply::future_lapply(
      seq_len(n_features),
      function(i) {
        process_one_feature(
          feature_name = feature_ids[i],
          feature_vec  = counts_subset[i,],
          design_mat   = design_subset,
          metadata_df  = metadata_subset,
          prog         = NULL
        )
      },
      future.seed = TRUE
    )
  } else {
    # Sequential execution with progress reporting
    start_time <- Sys.time()
    out_list <- vector("list", n_features)
    
    # Report every 10% or every 500 features, whichever is smaller
    report_interval <- min(500, max(1, floor(n_features / 10)))
    
    for (i in seq_len(n_features)) {
      out_list[[i]] <- process_one_feature(
        feature_name = feature_ids[i],
        feature_vec  = counts_subset[i,],
        design_mat   = design_subset,
        metadata_df  = metadata_subset,
        prog         = NULL
      )
      
      # Progress reporting
      if (show_progress && (i %% report_interval == 0 || i == n_features)) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        pct <- round(100 * i / n_features, 1)
        rate <- i / elapsed
        eta_sec <- (n_features - i) / rate
        
        if (has_cli) {
          cli::cli_alert_info(
            "Processed {i}/{n_features} ({pct}%) | Elapsed: {round(elapsed, 1)}s | ETA: {round(eta_sec, 1)}s"
          )
        } else {
          message(sprintf(
            "[%s] Processed %d/%d (%.1f%%) | Elapsed: %.1fs | ETA: %.1fs",
            format(Sys.time(), "%H:%M:%S"), i, n_features, pct, elapsed, eta_sec
          ))
        }
      }
    }
  }
  
  valid_list <- Filter(Negate(is.null), out_list)
  if (length(valid_list) == 0) {
    empty_res <- tibble::tibble(
      feature = character(), term = character(), effect = numeric(),
      se = numeric(), stat = numeric(), pval = numeric(), padj = numeric()
    )
    empty_sel <- tibble::tibble(
      feature = character(), best_family = character(),
      n_valid_obs = integer(), ic_value = numeric()
    )
    empty_out <- list(results = empty_res, selection = empty_sel)
    if (has_contrasts_requested) {
      if (omnibus) {
        empty_out$omnibus <- tibble::tibble(
          feature = character(), family = character(), test_type = character(),
          statistic = numeric(), df = integer(), p_value = numeric(), pass = logical()
        )
      }
      empty_out$contrasts <- tibble::tibble(
        feature = character(), family = character(), contrast = character(),
        estimate = numeric(), se = numeric(), z = numeric(),
        p_value = numeric(), p_adj = numeric()
      )
    }
    return(empty_out)
  }
  
  # ---- Aggregate results and compute FDR per term ----
  results_df <- dplyr::bind_rows(lapply(valid_list, `[[`, "coef_df"))
  if (nrow(results_df) > 0 && "pval" %in% names(results_df)) {
    results_df <- dplyr::mutate(
      dplyr::group_by(results_df, term),
      padj = p.adjust(pval, method = p_adjust),
      .groups = "drop"
    )
  } else {
    results_df$padj <- NA_real_
  }
  
  # Extract feature names - handle both character and numeric
  feature_names <- vapply(valid_list, function(x) as.character(x$feature), FUN.VALUE = character(1))
  
  selection_df <- tibble::tibble(
    feature      = feature_names,
    best_family  = vapply(valid_list, `[[`, "best_fam", FUN.VALUE = character(1)),
    n_valid_obs  = vapply(valid_list, `[[`, "n_valid", FUN.VALUE = integer(1)),
    ic_value     = vapply(valid_list, `[[`, "ic_value", FUN.VALUE = numeric(1)),
    transform_mode = transform_mode
  )
  
  output <- list(results = results_df, selection = selection_df)
  
  # ---- Aggregate omnibus results ----
  if (has_contrasts_requested && omnibus) {
    omnibus_list <- lapply(valid_list, `[[`, "omnibus_df")
    omnibus_list <- Filter(Negate(is.null), omnibus_list)
    
    if (length(omnibus_list) > 0) {
      output$omnibus <- dplyr::bind_rows(omnibus_list)
    }
  }
  
  # ---- Aggregate contrasts and compute FDR per contrast (UNCHANGED LOGIC) ----
  if (has_contrasts_requested) {
    contrasts_df <- dplyr::bind_rows(lapply(valid_list, `[[`, "contrast_df"))
    
    if (nrow(contrasts_df) > 0 && "p_value" %in% names(contrasts_df)) {
      # CRITICAL: This grouping logic is UNCHANGED from original PERSEO
      contrasts_df <- dplyr::mutate(
        dplyr::group_by(contrasts_df, contrast),
        p_adj = p.adjust(p_value, method = p_adjust),
        .groups = "drop"
      )
    } else {
      contrasts_df$p_adj <- NA_real_
    }
    
    output$contrasts <- contrasts_df
  }
  
  # ---- Store report data (will print on function exit) ----
  if (show_progress) {
    report_vars$n_total <- length(feature_ids)
    report_vars$n_success <- length(valid_list)
    report_vars$n_failed <- length(feature_ids) - length(valid_list)
    report_vars$selection_df <- selection_df
    report_vars$results_df <- results_df
    report_vars$output <- output
    report_vars$omnibus <- omnibus
    report_vars$omnibus_test <- omnibus_test
    report_vars$omnibus_threshold <- omnibus_threshold
    
    # Schedule report to print on exit (AFTER return, so it's the last thing shown)
    on.exit({
      print_final_report(report_vars)
    }, add = TRUE)
  }
  
  output
}

# Helper function to print final report
print_final_report <- function(vars) {
  n_total <- vars$n_total
  n_success <- vars$n_success
  n_failed <- vars$n_failed
  success_rate <- n_success / n_total * 100
  selection_df <- vars$selection_df
  results_df <- vars$results_df
  output <- vars$output
  omnibus <- vars$omnibus
  omnibus_test <- vars$omnibus_test
  omnibus_threshold <- vars$omnibus_threshold
  
  has_cli <- requireNamespace("cli", quietly = TRUE)
  
  if (has_cli) {
    cli::cli_rule("fit_gamlss_models() - FINAL REPORT")
  } else {
    message("\n", strrep("=", 80))
    message("fit_gamlss_models() - FINAL REPORT")
    message(strrep("=", 80))
  }
  if (has_cli) {
    cli::cli_h2("Model Fitting Summary")
    cli::cli_ul(c(
      "Total features      : {.val {n_total}}",
      "Successfully fitted : {.val {n_success}} ({.val {sprintf('%.1f%%', success_rate)}})",
      "Failed/skipped      : {.val {n_failed}} ({.val {sprintf('%.1f%%', n_failed/n_total*100)}})"
    ))
  } else {
    message("\nModel Fitting Summary:")
    message("  Total features      : ", n_total)
    message("  Successfully fitted : ", n_success, " (", sprintf("%.1f%%", success_rate), ")")
    message("  Failed/skipped      : ", n_failed, " (", sprintf("%.1f%%", n_failed/n_total*100), ")")
  }
  
  if (n_success > 0) {
    # Family distribution
    fam_table <- table(selection_df$best_family)
    fam_sorted <- sort(fam_table, decreasing = TRUE)
    
    if (has_cli) {
      cli::cli_h2("Selected Family Distribution")
      fam_items <- vapply(seq_along(fam_sorted), function(i) {
        fam_name <- names(fam_sorted)[i]
        count <- as.numeric(fam_sorted[i])
        pct <- count / sum(fam_sorted) * 100
        sprintf("%s: %d (%.1f%%)", fam_name, count, pct)
      }, character(1))
      cli::cli_ul(fam_items)
    } else {
      message("\nSelected Family Distribution:")
      for (i in seq_along(fam_sorted)) {
        fam_name <- names(fam_sorted)[i]
        count <- as.numeric(fam_sorted[i])
        pct <- count / sum(fam_sorted) * 100
        message(sprintf("  %-10s : %5d (%.1f%%)", fam_name, count, pct))
      }
    }
    
    # Coefficient statistics — exclude intercept (always significant by definition)
    if (nrow(results_df) > 0) {
      results_no_intercept <- results_df %>%
        dplyr::filter(term != "(Intercept)")
      n_coefs <- nrow(results_no_intercept)
      n_sig <- sum(results_no_intercept$padj < 0.05, na.rm = TRUE)
      sig_rate <- if (n_coefs > 0) n_sig / n_coefs * 100 else 0

      if (has_cli) {
        cli::cli_h2("Coefficient Results (excluding intercept)")
        cli::cli_ul(c(
          "Total coefficients  : {.val {n_coefs}}",
          "Significant (FDR<5%): {.val {n_sig}} ({.val {sprintf('%.1f%%', sig_rate)}})"
        ))
      } else {
        message("\nCoefficient Results (excluding intercept):")
        message("  Total coefficients  : ", n_coefs)
        message("  Significant (FDR<5%): ", n_sig, " (", sprintf("%.1f%%", sig_rate), ")")
      }

      # Per-term breakdown (intercept excluded)
      term_summary <- results_no_intercept %>%
        dplyr::group_by(term) %>%
        dplyr::summarize(
          n_total = dplyr::n(),
          n_sig = sum(padj < 0.05, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(n_sig))
      
      if (nrow(term_summary) > 0) {
        if (has_cli) {
          cli::cli_h3("Significant findings by term (top 10)")
          term_items <- vapply(seq_len(min(10, nrow(term_summary))), function(i) {
            term_name <- term_summary$term[i]
            n_sig_term <- term_summary$n_sig[i]
            n_total_term <- term_summary$n_total[i]
            pct_sig <- n_sig_term / n_total_term * 100
            sprintf("%s: %d / %d (%.1f%%)", term_name, n_sig_term, n_total_term, pct_sig)
          }, character(1))
          cli::cli_ul(term_items)
        } else {
          message("\n  Significant findings by term:")
          for (i in seq_len(min(10, nrow(term_summary)))) {
            term_name <- term_summary$term[i]
            n_sig_term <- term_summary$n_sig[i]
            n_total_term <- term_summary$n_total[i]
            pct_sig <- n_sig_term / n_total_term * 100
            message(sprintf("    %-25s : %5d / %5d (%.1f%%)", 
                          term_name, n_sig_term, n_total_term, pct_sig))
          }
        }
      }
    }
    
    # Contrast statistics if available
    if (!is.null(output$contrasts) && nrow(output$contrasts) > 0) {
      n_contrasts <- nrow(output$contrasts)
      n_sig_contrasts <- sum(output$contrasts$p_adj < 0.05, na.rm = TRUE)
      
      if (has_cli) {
        cli::cli_h2("Contrast Results")
        cli::cli_ul(c(
          "Total contrasts     : {.val {n_contrasts}}",
          "Significant (FDR<5%): {.val {n_sig_contrasts}} ({.val {sprintf('%.1f%%', n_sig_contrasts/n_contrasts*100)}})"
        ))
        
        if (omnibus && !is.null(output$omnibus)) {
          n_omnibus_pass <- sum(output$omnibus$pass, na.rm = TRUE)
          n_omnibus_total <- nrow(output$omnibus)
          cli::cli_h2("Omnibus Test ({omnibus_test})")
          cli::cli_ul(c(
            "Features tested     : {.val {n_omnibus_total}}",
            "Passed threshold    : {.val {n_omnibus_pass}} ({.val {sprintf('%.1f%%', n_omnibus_pass/n_omnibus_total*100)}})",
            "Threshold           : p < {.val {omnibus_threshold}}"
          ))
        }
      } else {
        message("\nContrast Results:")
        message("  Total contrasts     : ", n_contrasts)
        message("  Significant (FDR<5%): ", n_sig_contrasts, " (", 
                sprintf("%.1f%%", n_sig_contrasts/n_contrasts*100), ")")
        
        if (omnibus && !is.null(output$omnibus)) {
          n_omnibus_pass <- sum(output$omnibus$pass, na.rm = TRUE)
          n_omnibus_total <- nrow(output$omnibus)
          message("\nOmnibus Test (", omnibus_test, "):")
          message("  Features tested     : ", n_omnibus_total)
          message("  Passed threshold    : ", n_omnibus_pass, " (", 
                  sprintf("%.1f%%", n_omnibus_pass/n_omnibus_total*100), ")")
          message("  Threshold           : p < ", omnibus_threshold)
        }
      }
    }
  } else {
    if (has_cli) {
      cli::cli_alert_warning("No features were successfully fitted.")
    } else {
      message("\n[WARNING] No features were successfully fitted.")
    }
  }
  
  if (has_cli) {
    cli::cli_rule()
  } else {
    message(strrep("=", 80), "\n")
  }
}
