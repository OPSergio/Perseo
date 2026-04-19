
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [1.0.0] - 2026-02-08

First major release of PERSEO.

### Added

- **Enhanced user experience (UX)**: Major improvements to console output and progress reporting
- **Hierarchical omnibus testing**: `fit_gamlss_models()` and `run_perseo()` now support omnibus tests
  - New parameters: `omnibus`, `omnibus_test`, `omnibus_threshold`
  - Two test types: Wald (fast, vcov-based) and LRT (robust, model refitting)
  - Omnibus test acts as gatekeeper before computing pairwise contrasts
  - Reduces multiple testing burden for multi-level factors (3+ groups)
  - Returns `omnibus` table with test statistics, p-values, and pass/fail indicators
  - Mirrors ANOVA + post-hoc workflow familiar to statisticians
- **Automatic contrast generation**: `fit_gamlss_models()` accepts `contrast_variable` parameter
  - Auto-generates all pairwise contrasts for a categorical variable
  - Requires `metadata` parameter to identify factor levels
  - Simplifies workflow: no manual contrast matrix construction needed
  - Example: `contrast_variable = "tissue_type"` creates all tissue comparisons automatically
- **Formula-based design specification**: `fit_gamlss_models()` now accepts formula strings
  - Alternative to manual `model.matrix()` construction
  - Example: `design_matrix = "~ tissue_type + age + batch"`
  - Requires `metadata` parameter for variable resolution
  - Automatically handles factor encoding and reference levels
- **Built-in parallelization**: Automatic parallel processing configuration
  - New parameters: `parallel` and `workers` in all main functions
  - No manual `future::plan()` setup required
  - Automatic cleanup: resets to sequential plan after completion
  - Memory-efficient: optimal worker allocation based on system resources
- **Enhanced reporting system**: Improved user feedback and progress visibility
  - New S3 print method for `perseo_fit` objects (from `fit_gamlss_models()`)
  - Comprehensive summary reports with family distributions, significance rates, omnibus results
  - Automatic report generation when `show_progress = TRUE`
  - Report prints AFTER function return to avoid terminal clutter
  - Detailed startup banners in `find_families()` and `fit_gamlss_models()`
- **Optional bootstrap mode**: `find_families()` and `run_perseo()` now accept `bootstrap` parameter
  - `bootstrap = TRUE` (default): Fast bootstrap sampling of features
  - `bootstrap = FALSE`: Full evaluation of all families on ALL features (comprehensive but slower)
  - Allows users to choose between speed and exhaustive family selection
- **Custom contrast matrices**: `fit_gamlss_models()` now accepts `contrast_matrix` parameter
  - Contrast specification for arbitrary linear combinations of coefficients
  - Enables extraction of contrasts between non-reference levels (e.g., B-C when A is reference)
  - Returns contrasts tibble with estimates, standard errors, z-statistics, and p-values
  - Full FDR correction across features for contrast p-values
  - New helper function `apply_contrasts()` for robust contrast computation
- **High-level orchestration function**: `run_perseo()`
  - Complete end-to-end pipeline: family selection → differential expression → p-value adjustment
  - Verbosity control for progress messages
  - Global multiple testing correction with `p.adjust()` (BH, bonferroni, etc.)
  - Clean structured output with S3 class `perseo_results`
  - Custom print method for user-friendly summary display
- **Complete modularization**: New module files created
  - `validation.R`: Input validation and default value provision
  - `family_filtering.R`: Pure functions for family eligibility determination
  - `gamlss_fitting.R`: GAMLSS model fitting plumbing and IC computation
  - `family_selection_core.R`: Core family comparison and bootstrap aggregation
- **Comprehensive test suite**: 122 unit tests across 5 test files
  - `test-validation.R` (26 tests)
  - `test-family_filtering.R` (8 tests)
  - `test-utils_transformations.R` (42 tests)
  - `test-gamlss_fitting.R` (33 tests)
  - `test-family_selection_core.R` (13 tests)
- **Integration tests**: End-to-end workflow validation with realistic synthetic omics data
- **Package infrastructure**: DESCRIPTION and NAMESPACE files for R package compliance
- **Documentation**: Comprehensive testing guide with examples and best practices
- Added representations and interactive report.


### Changed

- **Major refactoring**: Eliminated ~90% code duplication between `find_families()` and `fit_gamlss_models()`
  - `find_families()`: Reduced from 280 to ~100 lines
  - `fit_gamlss_models()`: Reduced from 320 to ~176 lines
- **Defensive programming**: Robust handling of NA, NULL, and edge cases
  - `compare_families_on_feature()`: Defensive n_valid calculation with `is.finite()` checks
  - `compare_families_with_design()`: Same defensive improvements
  - `bind_bootstrap_results()`: Robust field extraction with NULL/empty filtering
- **Return structure improvements**: 
  - `compare_families_on_feature()` now returns `best_family` and `n_valid` fields
  - Consistent empty result handling across all comparison functions

### Fixed

- **Variance-covariance extraction**: Implemented 3-tier fallback strategy for robust vcov computation
  - Primary: Direct `vcov(fit, what = "mu")` call
  - Fallback 1: Diagonal vcov from `summary()` standard errors (conservative approximation)
  - Fallback 2: Extract `fit$vcov.mu` directly from model object
  - Addresses gamlss internal bug where `vcov()` searches for `family_obj` in parent frame
  - Fixed dimension mismatch for multi-parameter families (extract only mu coefficients, not sigma)
- **Output suppression**: Fixed console flooding during tests
  - Corrected `extract_mu_coefficients()` to use `capture.output(..., file = character())`
  - Added explicit sink cleanup in `fit_gamlss_safely()` to prevent dangling connections
  - All GAMLSS fitting now properly suppressed during test execution
- **Critical bug**: Fixed NA handling in common mask computation that caused "missing value where TRUE/FALSE needed" errors
- **Bootstrap aggregation**: Fixed vapply length errors when processing NULL or malformed results
- **Test suite corrections**: 
  - Fixed all syntax errors (missing `})` parentheses)
  - Aligned tests with actual function signatures
  - Removed non-existent parameters from test calls
  - Made tests permissive where multiple valid outcomes exist

### Documentation

- Updated `testing_guide.md` with real test count (122 tests)
- Added detailed notes on function signature requirements
- Documented defensive programming patterns
- Added "Real Implementation Focus" principle to test design

---

## [0.1.0] - 2025-07-31
Initial functional release of the preprocessing and GAMLSS model visualization pipeline.
