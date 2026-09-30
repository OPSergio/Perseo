#!/usr/bin/env Rscript
# PERSEO runner for agents and scripts: reads files, aligns metadata to samples,
# runs the pipeline silently, writes every table to --out and prints a compact
# summary to stdout. Run with --help for usage.

usage <- '
Usage: Rscript perseo_run.R --counts FILE --formula "~ x + y" [options]

Inputs (.rds, .csv, .tsv, optionally .gz; CSV/TSV use the first column as row IDs)
  --counts FILE         features x samples (matrix or data.frame)
  --metadata FILE       one row per sample; IDs from rownames or any column matching
                        the counts column names. Required with --formula.
  --formula STR         design, e.g. "~ condition + batch"
  --contrast VAR        factor in metadata: all pairwise contrasts
  --ref LEVEL           reference level for --contrast
Model
  --families A,B,..     candidate families (default: PERSEO defaults)
  --skip-selection      skip find_families(); fit --families directly
  --criterion X         BIC | AIC | GAIC (default BIC)
  --n-genes N --n-boot N --top-n N   family selection (default 200 / 10 / 4)
  --omnibus             gate contrasts on an omnibus test
  --omnibus-test X      Wald | LRT (default Wald)
  --padj-method X       p.adjust method (default BH)
  --min-n N             minimum valid observations per feature (default 5)
Run
  --check               validate inputs and print the plan; do not fit
  --subset N            pilot on N random features
  --workers N           parallel workers (default 1 = sequential)
  --seed N              default 1
  --out DIR             output directory (default perseo_out)
  --fdr X               threshold used in the summary (default 0.05)
  --top N               top features printed per contrast/term (default 5)
  --report              also render report_perseo() HTML into --out
'

parse_args <- function(argv) {
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a, call. = FALSE)
    key <- sub("^--", "", a)
    if (grepl("=", key, fixed = TRUE)) {
      val <- sub("^[^=]*=", "", key)
      key <- sub("=.*$", "", key)
    } else if (i < length(argv) && !startsWith(argv[i + 1L], "--")) {
      i <- i + 1L
      val <- argv[i]
    } else {
      val <- "TRUE"
    }
    out[[gsub("-", "_", key)]] <- val
    i <- i + 1L
  }
  out
}

arg_num <- function(args, key, default) {
  if (is.null(args[[key]])) return(default)
  v <- suppressWarnings(as.numeric(args[[key]]))
  if (is.na(v)) stop("--", gsub("_", "-", key), " must be numeric", call. = FALSE)
  v
}

read_input <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
  switch(ext,
    rds = readRDS(path),
    csv = utils::read.csv(path, row.names = 1, check.names = FALSE),
    tsv = , txt = , tab = utils::read.delim(path, row.names = 1, check.names = FALSE),
    stop("Unsupported file type: ", path, " (use .rds, .csv or .tsv)", call. = FALSE)
  )
}

as_counts_matrix <- function(x) {
  if (is.data.frame(x)) {
    num <- vapply(x, is.numeric, logical(1))
    if (!all(num)) {
      stop("counts has non-numeric columns: ",
           paste(utils::head(names(x)[!num], 5), collapse = ", "), call. = FALSE)
    }
  }
  m <- tryCatch(as.matrix(x), error = function(e) {
    stop("Cannot convert counts (class ", class(x)[1], ") to a matrix. ",
         "Extract the assay first (e.g. SummarizedExperiment::assay(se)).", call. = FALSE)
  })
  storage.mode(m) <- "double"
  if (is.null(rownames(m))) rownames(m) <- paste0("feature_", seq_len(nrow(m)))
  if (is.null(colnames(m))) stop("counts has no column (sample) names", call. = FALSE)
  m
}

align_metadata <- function(meta, counts) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  ids <- colnames(counts)
  if (!all(ids %in% rownames(meta))) {
    hit <- names(meta)[vapply(meta, function(v) all(ids %in% as.character(v)), logical(1))]
    if (length(hit) > 0) {
      meta <- meta[!duplicated(as.character(meta[[hit[1]]])), , drop = FALSE]
      rownames(meta) <- as.character(meta[[hit[1]]])
    } else if (all(rownames(counts) %in% rownames(meta))) {
      stop("counts looks like samples x features; transpose it (features must be rows)",
           call. = FALSE)
    } else {
      miss <- setdiff(ids, rownames(meta))
      stop(length(miss), " count columns not found in metadata (e.g. ",
           paste(utils::head(miss, 3), collapse = ", "), ")", call. = FALSE)
    }
  }
  meta[ids, , drop = FALSE]
}

fmt <- function(x) formatC(x, digits = 3, format = "g")

count_line <- function(tab, total) {
  paste(sprintf("%s=%d/%d", names(tab), as.integer(tab), as.integer(total[names(tab)])),
        collapse = " ; ")
}

main <- function(args) {
  if (!is.null(args$help) || length(args) == 0) {
    cat(usage)
    return(invisible(0))
  }
  if (is.null(args$counts)) stop("--counts is required", call. = FALSE)
  if (is.null(args$formula)) stop("--formula is required", call. = FALSE)
  if (is.null(args$metadata)) stop("--metadata is required with --formula", call. = FALSE)

  suppressPackageStartupMessages(library(PERSEO))
  seed <- arg_num(args, "seed", 1)
  fdr <- arg_num(args, "fdr", 0.05)
  top <- arg_num(args, "top", 5)
  out_dir <- if (is.null(args$out)) "perseo_out" else args$out
  contrast <- args$contrast
  families <- if (!is.null(args$families)) trimws(strsplit(args$families, ",")[[1]]) else NULL

  counts <- as_counts_matrix(read_input(args$counts))
  meta <- align_metadata(read_input(args$metadata), counts)

  f <- stats::as.formula(args$formula)
  vars <- all.vars(f)
  missing_vars <- setdiff(vars, names(meta))
  if (length(missing_vars) > 0) {
    stop("formula variables not in metadata: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }
  for (v in vars) if (is.character(meta[[v]])) meta[[v]] <- factor(meta[[v]])

  if (!is.null(contrast)) {
    if (!contrast %in% vars) stop("--contrast '", contrast, "' is not in the formula", call. = FALSE)
    meta[[contrast]] <- factor(meta[[contrast]])
    if (!is.null(args$ref)) {
      if (!args$ref %in% levels(meta[[contrast]])) {
        stop("--ref '", args$ref, "' is not a level of ", contrast, ": ",
             paste(levels(meta[[contrast]]), collapse = ", "), call. = FALSE)
      }
      meta[[contrast]] <- stats::relevel(meta[[contrast]], ref = args$ref)
    }
  }

  set.seed(seed)
  if (!is.null(args$subset)) {
    n_sub <- min(nrow(counts), as.integer(arg_num(args, "subset", nrow(counts))))
    counts <- counts[sort(sample(nrow(counts), n_sub)), , drop = FALSE]
  }

  if (!is.null(args$check)) {
    complete <- if (length(vars) > 0) {
      stats::complete.cases(meta[, vars, drop = FALSE])
    } else {
      rep(TRUE, nrow(meta))
    }
    idx <- if (nrow(counts) > 500) sample(nrow(counts), 500) else seq_len(nrow(counts))
    support <- table(vapply(idx, function(i) infer_support(counts[i, ]), character(1)))
    cat(sprintf("counts: %d features x %d samples | range %s..%s | NA: %d\n",
                nrow(counts), ncol(counts), fmt(min(counts, na.rm = TRUE)),
                fmt(max(counts, na.rm = TRUE)), sum(is.na(counts))))
    cat(sprintf("metadata aligned: %d samples | dropped for NA in formula vars: %d\n",
                nrow(meta), sum(!complete)))
    cat("support (sample of", length(idx), "features):",
        paste(names(support), as.integer(support), sep = "=", collapse = " "), "\n")
    for (v in vars) {
      x <- meta[[v]][complete]
      if (is.factor(x)) {
        cat(sprintf("%s: factor %s\n", v,
                    paste(names(table(x)), as.integer(table(x)), sep = "=", collapse = " ")))
      } else {
        cat(sprintf("%s: %s [%s, %s]\n", v, class(x)[1], fmt(min(x)), fmt(max(x))))
      }
    }
    if (!is.null(contrast)) {
      lv <- levels(droplevels(meta[[contrast]][complete]))
      pr <- utils::combn(length(lv), 2)
      cat("contrasts:", paste0(lv[pr[2, ]], "_vs_", lv[pr[1, ]], collapse = ", "), "\n")
      bad <- lv[make.names(lv) != lv]
      if (length(bad) > 0) cat("WARNING non-syntactic levels (rename with make.names):",
                               paste(bad, collapse = ", "), "\n")
    }
    n_genes <- arg_num(args, "n_genes", 200)
    if (is.null(args$skip_selection) && nrow(counts) < n_genes) {
      cat(sprintf("note: n_genes will be capped to %d (nrow)\n", nrow(counts)))
    }
    cat("check OK\n")
    return(invisible(0))
  }

  common <- list(
    counts_matrix = counts, design_matrix = args$formula, metadata = meta,
    contrast_variable = contrast,
    criterion = if (is.null(args$criterion)) "BIC" else toupper(args$criterion),
    min_n = arg_num(args, "min_n", 5),
    omnibus = !is.null(args$omnibus),
    omnibus_test = if (identical(tolower(args$omnibus_test), "lrt")) "LRT" else "Wald",
    parallel = arg_num(args, "workers", 1) > 1,
    workers = if (arg_num(args, "workers", 1) > 1) as.integer(args$workers) else NULL,
    show_progress = FALSE
  )
  padj_method <- if (is.null(args$padj_method)) "BH" else args$padj_method

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  res <- withCallingHandlers(
    if (!is.null(args$skip_selection)) {
      if (is.null(families)) stop("--skip-selection needs --families", call. = FALSE)
      do.call(fit_gamlss_models, c(common, list(
        candidate_families = families, p_adjust = padj_method,
        group_by_support = TRUE, transform_mode = "strict"
      )))
    } else {
      do.call(run_perseo, c(common, list(
        families = families, p_adjust_method = padj_method, seed = seed,
        n_genes = min(arg_num(args, "n_genes", 200), nrow(counts)),
        n_boot = arg_num(args, "n_boot", 10), top_n = arg_num(args, "top_n", 4)
      )))
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  minutes <- (proc.time()[["elapsed"]] - t0) / 60

  saveRDS(res, file.path(out_dir, "result.rds"))
  writeLines(warns, file.path(out_dir, "warnings.log"))
  de <- if (inherits(res, "perseo_results")) res$differential_expression else res
  if (is.null(de)) {
    cat("FAILED:", res$summary$status, "- no family selected. warnings:", length(warns), "\n")
    return(invisible(1))
  }

  write_tbl <- function(x, name) {
    if (!is.null(x) && nrow(x) > 0) {
      utils::write.csv(x, file.path(out_dir, name), row.names = FALSE)
      name
    }
  }
  fam_sel <- if (inherits(res, "perseo_results")) res$family_selection else NULL
  files <- c(
    "result.rds",
    write_tbl(de$results, "results.csv"),
    write_tbl(de$selection, "selection.csv"),
    write_tbl(de$contrasts, "contrasts.csv"),
    write_tbl(de$omnibus, "omnibus.csv"),
    if (!is.null(fam_sel)) write_tbl(
      data.frame(family = names(fam_sel$freq_table_overall),
                 n = as.integer(fam_sel$freq_table_overall)), "family_selection.csv")
  )
  if (!is.null(args$report) && inherits(res, "perseo_results")) {
    report_perseo(res, output_dir = out_dir, open = FALSE, quiet = TRUE)
    files <- c(files, "perseo_report.html")
  }

  sel <- de$selection
  cat(sprintf("PERSEO %s | %d features x %d samples | fitted %d (skipped %d) | %.1f min\n",
              as.character(utils::packageVersion("PERSEO")), nrow(counts), ncol(counts),
              nrow(sel), nrow(counts) - nrow(sel), minutes))
  if (!is.null(fam_sel)) {
    cat("families used:", paste(fam_sel$top_families_overall, collapse = ","), "\n")
  }
  bf <- sort(table(sel$best_family), decreasing = TRUE)
  cat("best_family:", paste(names(bf), as.integer(bf), sep = "=", collapse = " "), "\n")

  r <- de$results[de$results$term != "(Intercept)", ]
  if (nrow(r) > 0) {
    hits <- tapply(r$padj < fdr, r$term, sum, na.rm = TRUE)
    hits <- sort(hits, decreasing = TRUE)
    cat(sprintf("terms padj<%s: %s\n", fdr, count_line(hits, table(r$term))))
  }

  show_top <- function(x, p_col, est_col) {
    x <- x[order(x[[p_col]]), ]
    x <- utils::head(x[!is.na(x[[p_col]]), ], top)
    for (i in seq_len(nrow(x))) {
      cat(sprintf("  %s %s est=%s %s=%s\n", x$feature[i], x$family[i],
                  fmt(x[[est_col]][i]), p_col, fmt(x[[p_col]][i])))
    }
  }

  cn <- de$contrasts
  if (!is.null(cn) && nrow(cn) > 0) {
    for (k in unique(cn$contrast)) {
      x <- cn[cn$contrast == k, ]
      sig <- !is.na(x$p_adj) & x$p_adj < fdr
      cat(sprintf("contrast %s: %d/%d p_adj<%s (up %d, down %d, NA %d)\n", k, sum(sig),
                  nrow(x), fdr, sum(sig & x$estimate > 0), sum(sig & x$estimate < 0),
                  sum(is.na(x$p_adj))))
      show_top(x, "p_adj", "estimate")
    }
  } else if (nrow(r) > 0) {
    r$family <- sel$best_family[match(r$feature, sel$feature)]
    for (k in utils::head(names(hits), 3)) {
      cat(sprintf("term %s:\n", k))
      show_top(r[r$term == k, ], "padj", "effect")
    }
  }
  if (!is.null(de$omnibus)) {
    cat(sprintf("omnibus (%s, raw p<0.05): %d/%d passed\n", de$omnibus$test_type[1],
                sum(de$omnibus$pass, na.rm = TRUE), nrow(de$omnibus)))
  }
  cat("effect scale: link of best_family (log for counts/GA/GG/IG/LOGNO, logit for beta/BI, SD units for NO/TF/GU)\n")
  cat("files in", out_dir, ":", paste(files, collapse = ", "), "\n")
  cat("warnings:", length(warns),
      if (length(warns) > 0) paste0("(", length(unique(warns)), " unique; see warnings.log)"), "\n")
  invisible(0)
}

status <- tryCatch(
  main(parse_args(commandArgs(trailingOnly = TRUE))),
  error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    1
  }
)
quit(save = "no", status = as.integer(status))
