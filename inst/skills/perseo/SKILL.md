---
name: perseo
description: Differential expression / differential abundance on omics feature-by-sample matrices (RNA-seq counts, proteomics or metabolomics intensities, methylation beta values, microbiome, any numeric features) with the PERSEO R package, which picks the best GAMLSS distribution family per feature. Use when the user mentions PERSEO or GAMLSS-based DE, asks for per-feature distribution/family selection, or wants to run, debug, plot or interpret run_perseo / fit_gamlss_models / find_families results.
license: GPL-3
metadata:
  package: PERSEO
  package-version: "1.0.0"
  repository: https://github.com/OPSergio/Perseo
---

# PERSEO

R package `PERSEO` (uppercase). Per feature: fit candidate GAMLSS families on a common set of
observations, pick the lowest Jacobian-corrected IC (near-ties within 2 IC units broken by
goodness of fit), then report Wald tests per coefficient, optional pairwise contrasts and an
optional omnibus gate. Install: `remotes::install_github("OPSergio/Perseo")` (R >= 4.2).

## Fastest path: the CLI runner (no R code to write, compact output)

`scripts/perseo_run.R` reads files, aligns metadata to samples, runs the pipeline silently,
writes every table to disk and prints a summary of about 30 lines. Locate it with
`Rscript -e 'cat(system.file("skills/perseo/scripts/perseo_run.R", package="PERSEO"))'`.

```sh
R=scripts/perseo_run.R
# 1. Validate inputs only (seconds): dims, alignment, NA drops, levels, contrast names, support
Rscript $R --counts counts.rds --metadata meta.csv --formula "~ condition + batch" \
  --contrast condition --ref Control --check
# 2. Pilot on 300 random features to catch problems cheaply
Rscript $R ... --subset 300 --out pilot
# 3. Full run (long: run in background, output goes to files)
Rscript $R ... --workers 8 --out perseo_out > perseo_out.log 2>&1
```

Inputs: `.rds`, `.csv`, `.tsv` (optionally `.gz`). The counts file has features in rows and
samples in columns, with feature IDs in the first column or in rownames. The metadata file
has one row per sample; sample IDs come from its rownames/first column or from any column
that matches the counts column names. Run `--help` for every flag (`--families`,
`--skip-selection`, `--omnibus`, `--criterion`, `--report`, ...). Outputs are written to
`--out`: `result.rds`, `results.csv`, `selection.csv`, `contrasts.csv`, `omnibus.csv`,
`family_selection.csv` and `warnings.log`.

## R API (when you must write R)

```r
library(PERSEO)
counts <- as.matrix(counts)                      # base numeric matrix, features x samples
meta   <- meta[colnames(counts), , drop = FALSE] # REQUIRED: rows in the same order as columns
meta$condition <- relevel(factor(meta$condition), ref = "Control")

res <- run_perseo(
  counts_matrix = counts, design_matrix = "~ condition + batch", metadata = meta,
  contrast_variable = "condition",   # all pairwise contrasts of this factor
  criterion = "BIC", seed = 1,
  n_genes = min(200, nrow(counts)),  # bootstrap sample size must be <= nrow(counts)
  parallel = TRUE, workers = 8,
  show_progress = FALSE              # ALWAYS in agent runs: otherwise hundreds of log lines
)
saveRDS(res, "perseo_result.rds")    # fitting is expensive; never refit to re-inspect
de <- res$differential_expression
```

Lower level: `find_families()` (family choice only), then `fit_gamlss_models(candidate_families=)`
(DE with known families, returns the `de` list directly). Full signatures are in `references/api.md`.

## Output schema (`de <- res$differential_expression`)

| table | one row per | columns |
|---|---|---|
| `de$results` | feature x coefficient | `feature term effect se stat pval padj` |
| `de$selection` | fitted feature | `feature best_family n_valid_obs ic_value transform_mode` |
| `de$contrasts` | feature x contrast | `feature family contrast estimate se z p_value p_adj` |
| `de$omnibus` | feature (if `omnibus=TRUE`) | `feature family test_type statistic df p_value pass` |

Also available: `res$family_selection$top_families_overall` (the families used) and
`res$summary` (the run parameters).

## Pitfalls (verified against the source)

1. **Column names differ between tables:** `results` uses `pval`/`padj`, while `contrasts` and
   `omnibus` use `p_value`/`p_adj`. `filter(results, p_adj < .05)` fails.
2. **Metadata is matched by position, not by rownames.** Reorder it with
   `meta[colnames(counts), ]` first.
3. `counts_matrix` must be a base numeric `matrix`. For a data.frame, sparse `Matrix` or
   SummarizedExperiment, call `as.matrix(...)` or `assay(se)` first.
4. `bootstrap = TRUE` (the default) errors unless `nrow(counts) >= n_genes` (default 200).
5. **No library-size normalisation, and `offset()` terms in the formula are silently dropped.**
   Sequencing depth must be handled explicitly by the analyst (e.g. normalised input or a depth
   covariate). Tell the user this; don't choose for them.
6. Samples with `NA` in any formula variable are dropped. Constant or near-constant features, and
   features whose every fit fails, are silently absent from the outputs. Compare
   `nrow(de$selection)` with `nrow(counts)`.
7. Contrast names are `<later>_vs_<earlier>` in `levels()` order, e.g. `Treat_vs_Control`.
   Set the reference with `relevel()` and use syntactic level names (`make.names`).
8. `contrast_variable` needs `metadata`. If `contrast_matrix` is also given, the matrix wins.
   Its column names must equal the coefficient names (`colnames(model.matrix(...))`). With
   interactions, auto-contrasts use only main-effect columns, which means they are evaluated
   at the reference level of the other factors.
9. `omnibus = TRUE` requires contrasts. It gates on the **raw** omnibus p-value
   (`p_value < omnibus_threshold`), and contrasts are then computed only for features where
   `pass == TRUE`.
10. **`effect`/`estimate` scale depends on `best_family`.** It is natural log (not log2) for
    count and GA/GG/IG/LOGNO families, logit for BI/BB and beta families, and SD units for
    NO/TF/GU (z-scored). Features with different families are not on a common scale, and
    `plot_volcano(lfc_threshold = 1)` means 1 unit on that scale. See
    `references/interpretation.md`.
11. FDR for coefficients is BH **within each term** across features, and for contrasts **within
    each contrast**. Exclude `(Intercept)` when counting hits.
12. `plot_ma(run_perseo_output)` does not find the counts by itself; pass `counts_matrix=`.
    `report_perseo()` accepts only `run_perseo()` output, needs `rmarkdown`, `DT`,
    `htmltools` and `jsonlite`, and should be called with `open = FALSE` in scripts.
13. `fit_gamlss_models` defaults to `group_by_support = FALSE` (unlike `run_perseo` and
    `find_families`), so pass `TRUE` to reproduce `run_perseo`. Never combine
    `group_by_support = FALSE` with `transform_mode = "strict"`: the common mask across
    families of other supports becomes empty and every feature is skipped.

## Token discipline

- Always set `show_progress = FALSE`, and redirect long runs to a log (`> run.log 2>&1`). Read
  the log only with `tail`/`grep`.
- Never print whole tibbles. Print counts (`table()`, `sum(padj < .05)`) and `head(arrange(x, p_adj), 10)`
  with a few selected columns, or write CSVs and summarise them.
- Pilot on a subset (`counts[sample(nrow(counts), 300), ]`) before the full run. Runtime scales
  with features x families, so pre-filter uninformative features when the user agrees.
- Reuse `saveRDS` results and `--skip-selection --families` to rerun DE without the bootstrap.

## References (load only when needed)

- `references/api.md`: every exported function, all arguments and defaults, family sets.
- `references/interpretation.md`: effect scales per family, testing and FDR scheme, NA and
  skip causes, how to report results.
- `references/methods.md`: the algorithm (support inference, filtering rules, transforms,
  Jacobian, IC and GoF tie-break, vcov, bootstrap), for methods sections and deep debugging.
