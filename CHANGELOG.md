
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

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
