
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Initial design of `fit_gamlss_models()` function to fit multiple GAMLSS families per gene.
- Parallelization support using `furrr` to speed up model fitting.
- Automatic summary report of fitted models and most frequent families.
- `transform_for_family()` function to apply family-specific transformations (BE, TF, GA, etc.).
- Support for "safe" or "strict" strategies in transformations.

### Changed

- `find_families()` was updated to keep transformed `y_t` for further visualization.

### Fixed

### Deprecated

### Removed

### Security

---

## [0.1.0] - 2025-07-31
Initial functional release of the preprocessing and GAMLSS model visualization pipeline.
