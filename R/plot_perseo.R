#' Internal ggplot2 theme for PERSEO plots
#'
#' @keywords internal
theme_perseo <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(colour = "#CCCCCC", fill = NA, linewidth = 0.5),
      strip.background   = ggplot2::element_rect(fill = "#F5F5F5", colour = "#CCCCCC"),
      strip.text         = ggplot2::element_text(face = "bold"),
      legend.position    = "right",
      plot.title         = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle      = ggplot2::element_text(colour = "#555555", size = base_size - 1),
      axis.title         = ggplot2::element_text(face = "bold"),
      plot.caption       = ggplot2::element_text(colour = "#888888", size = base_size - 3)
    )
}

#' PERSEO colour palette for significance categories
#' @keywords internal
.perseo_sig_colours <- c(
  "Not sig"    = "#AAAAAA",
  "LFC only"   = "#4E9FD4",
  "FDR only"   = "#F4A460",
  "Sig & LFC"  = "#C0392B"
)

#' Extract contrast tibble from a PERSEO result object
#' @keywords internal
.extract_contrasts <- function(x) {
  if (is.data.frame(x) && "contrast" %in% names(x)) {
    return(x)
  }
  if ("contrasts" %in% names(x)) {
    return(x$contrasts)
  }
  if ("differential_expression" %in% names(x) &&
      "contrasts" %in% names(x$differential_expression)) {
    return(x$differential_expression$contrasts)
  }
  stop("Cannot find contrasts in the supplied object. ",
       "Pass the output of fit_gamlss_models() or run_perseo(), ",
       "or a contrasts tibble directly.")
}

#' Extract counts matrix from a PERSEO result object (for MA plot)
#' @keywords internal
.extract_counts <- function(x) {
  if ("counts_matrix" %in% names(x)) return(x$counts_matrix)
  if ("differential_expression" %in% names(x) &&
      "counts_matrix" %in% names(x$differential_expression)) {
    return(x$differential_expression$counts_matrix)
  }
  NULL
}


#' Volcano plot for PERSEO differential expression results
#'
#' Displays log fold change (estimate) against -log10 adjusted p-value for
#' pairwise contrasts computed by \code{fit_gamlss_models()} or
#' \code{run_perseo()}. Points are coloured by significance category.
#'
#' @param x Output of \code{fit_gamlss_models()}, \code{run_perseo()}, or a
#'   contrasts tibble (must contain columns \code{contrast}, \code{estimate},
#'   \code{p_adj}, and \code{feature}).
#' @param contrast Character string naming which contrast to plot. If the data
#'   contains only one contrast it is selected automatically.
#' @param fdr_threshold Numeric FDR threshold (default: 0.05).
#' @param lfc_threshold Numeric |LFC| threshold for the "LFC only" / "Sig & LFC"
#'   categories (default: 1, i.e. 2-fold on the log scale).
#' @param label_top Integer, number of top significant features to label
#'   (default: 10). Set to 0 to suppress labels. Requires \pkg{ggrepel}.
#' @param point_size Numeric, base point size (default: 1.8).
#' @param alpha Numeric, point transparency (default: 0.7).
#' @param title Character plot title; if NULL a default is constructed.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' result <- fit_gamlss_models(counts, design, candidate_families = c("PO"),
#'                             contrast_matrix = C, show_progress = FALSE)
#' plot_volcano(result, contrast = "B_vs_A")
#' }
#'
#' @export
plot_volcano <- function(x,
                         contrast      = NULL,
                         fdr_threshold = 0.05,
                         lfc_threshold = 1,
                         label_top     = 10,
                         point_size    = 1.8,
                         alpha         = 0.7,
                         title         = NULL) {

  df <- .extract_contrasts(x)

  if (is.null(contrast)) {
    available <- unique(df$contrast)
    if (length(available) == 1L) {
      contrast <- available
    } else {
      stop("Multiple contrasts found: ", paste(available, collapse = ", "),
           ". Please specify one via the 'contrast' argument.")
    }
  }

  df <- df[df$contrast == contrast, ]
  if (nrow(df) == 0L) {
    stop("No rows found for contrast '", contrast, "'.")
  }

  if (!"p_adj" %in% names(df)) {
    stop("Column 'p_adj' not found. Run fit_gamlss_models() with p_adjust != 'none'.")
  }

  df$neg_log10_padj <- -log10(pmax(df$p_adj, 1e-300))

  df$sig_cat <- dplyr::case_when(
    !is.finite(df$estimate) | !is.finite(df$p_adj) ~ "Not sig",
    df$p_adj < fdr_threshold & abs(df$estimate) >= lfc_threshold ~ "Sig & LFC",
    df$p_adj < fdr_threshold                                      ~ "FDR only",
    abs(df$estimate) >= lfc_threshold                             ~ "LFC only",
    TRUE                                                          ~ "Not sig"
  )
  df$sig_cat <- factor(df$sig_cat,
                       levels = c("Not sig", "LFC only", "FDR only", "Sig & LFC"))

  n_sig <- sum(df$sig_cat == "Sig & LFC", na.rm = TRUE)
  plot_title <- title %||%
    paste0("Volcano Plot — ", contrast)
  plot_sub <- paste0("FDR < ", fdr_threshold, "  |  |LFC| \u2265 ", lfc_threshold,
                     "  |  Significant & LFC: ", n_sig, " features")

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = .data$estimate,
                      y = .data$neg_log10_padj,
                      colour = .data$sig_cat)) +
    ggplot2::geom_point(size = point_size, alpha = alpha, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "#555555", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold),
                        linetype = "dashed", colour = "#555555", linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = .perseo_sig_colours, name = NULL,
                                 drop = FALSE) +
    ggplot2::labs(x = "Log Fold Change (estimate)",
                  y = expression(-log[10](p[adj])),
                  title = plot_title,
                  subtitle = plot_sub,
                  caption = "PERSEO \u2022 GAMLSS-based differential expression") +
    theme_perseo()

  if (label_top > 0L) {
    top_df <- df[df$sig_cat == "Sig & LFC", ]
    if (nrow(top_df) == 0L) {
      top_df <- df[df$p_adj < fdr_threshold & is.finite(df$p_adj), ]
    }
    if (nrow(top_df) > 0L) {
      top_df <- top_df[order(top_df$p_adj), ]
      top_df <- utils::head(top_df, label_top)

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data    = top_df,
          mapping = ggplot2::aes(label = .data$feature),
          size    = 3,
          colour  = "#222222",
          box.padding       = 0.4,
          point.padding     = 0.2,
          max.overlaps      = 20,
          show.legend       = FALSE
        )
      } else {
        p <- p + ggplot2::geom_text(
          data    = top_df,
          mapping = ggplot2::aes(label = .data$feature),
          size    = 3,
          colour  = "#222222",
          vjust   = -0.5,
          show.legend = FALSE
        )
      }
    }
  }

  p
}


#' MA plot for PERSEO differential expression results
#'
#' Displays mean expression (log2 scale, derived from the counts matrix when
#' available, otherwise from the contrast estimate rank) against log fold change.
#' Points are coloured by significance.
#'
#' @param x Output of \code{fit_gamlss_models()}, \code{run_perseo()}, or a
#'   contrasts tibble. For the x-axis to show true mean expression, \code{x}
#'   should be the full \code{run_perseo()} result or you can supply
#'   \code{counts_matrix} directly.
#' @param contrast Character string naming which contrast to plot.
#' @param counts_matrix Optional numeric matrix (features \eqn{\times} samples)
#'   used to compute per-feature mean expression for the x-axis. If NULL,
#'   the feature rank is used instead.
#' @param fdr_threshold Numeric FDR threshold for colouring (default: 0.05).
#' @param label_top Integer, number of top features to label (default: 10).
#' @param point_size Numeric, base point size (default: 1.8).
#' @param alpha Numeric, point transparency (default: 0.7).
#' @param title Character plot title; if NULL a default is constructed.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' result <- run_perseo(counts, "~ group", metadata = meta,
#'                      contrast_variable = "group", seed = 42)
#' plot_ma(result, contrast = "B_vs_A", counts_matrix = counts)
#' }
#'
#' @export
plot_ma <- function(x,
                    contrast      = NULL,
                    counts_matrix = NULL,
                    fdr_threshold = 0.05,
                    label_top     = 10,
                    point_size    = 1.8,
                    alpha         = 0.7,
                    title         = NULL) {

  df <- .extract_contrasts(x)

  if (is.null(counts_matrix)) {
    counts_matrix <- .extract_counts(x)
  }

  if (is.null(contrast)) {
    available <- unique(df$contrast)
    if (length(available) == 1L) {
      contrast <- available
    } else {
      stop("Multiple contrasts found: ", paste(available, collapse = ", "),
           ". Please specify one via the 'contrast' argument.")
    }
  }

  df <- df[df$contrast == contrast, ]
  if (nrow(df) == 0L) stop("No rows found for contrast '", contrast, "'.")

  if (!"p_adj" %in% names(df)) {
    stop("Column 'p_adj' not found.")
  }

  # Compute mean expression (A axis)
  if (!is.null(counts_matrix) && is.matrix(counts_matrix)) {
    feat_means <- rowMeans(counts_matrix, na.rm = TRUE)
    mean_df <- data.frame(
      feature  = names(feat_means),
      mean_expr = log2(feat_means + 1),
      stringsAsFactors = FALSE
    )
    df <- merge(df, mean_df, by = "feature", all.x = TRUE)
    x_label <- expression(log[2](Mean~expression + 1))
    x_col   <- "mean_expr"
  } else {
    df$feature_rank <- rank(-abs(df$estimate), ties.method = "first")
    x_label <- "Feature rank (by |LFC|)"
    x_col   <- "feature_rank"
  }

  df$significant <- !is.na(df$p_adj) & df$p_adj < fdr_threshold

  sig_colours <- c("TRUE" = "#C0392B", "FALSE" = "#AAAAAA")

  plot_title <- title %||% paste0("MA Plot — ", contrast)
  n_sig <- sum(df$significant, na.rm = TRUE)
  plot_sub <- paste0("FDR < ", fdr_threshold,
                     "  |  Significant: ", n_sig, " features")

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x     = .data[[x_col]],
                      y     = .data$estimate,
                      colour = .data$significant)) +
    ggplot2::geom_point(size = point_size, alpha = alpha, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                        colour = "#333333", linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = sig_colours,
                                 labels = c("TRUE" = paste0("FDR < ", fdr_threshold),
                                            "FALSE" = "Not sig"),
                                 name = NULL) +
    ggplot2::labs(x       = x_label,
                  y       = "Log Fold Change (estimate)",
                  title   = plot_title,
                  subtitle = plot_sub,
                  caption  = "PERSEO \u2022 GAMLSS-based differential expression") +
    theme_perseo()

  if (label_top > 0L) {
    top_df <- df[df$significant & is.finite(df$p_adj), ]
    if (nrow(top_df) > 0L) {
      top_df <- top_df[order(top_df$p_adj), ]
      top_df <- utils::head(top_df, label_top)

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data    = top_df,
          mapping = ggplot2::aes(label = .data$feature),
          size    = 3,
          colour  = "#222222",
          box.padding   = 0.4,
          point.padding = 0.2,
          max.overlaps  = 20,
          show.legend   = FALSE
        )
      } else {
        p <- p + ggplot2::geom_text(
          data    = top_df,
          mapping = ggplot2::aes(label = .data$feature),
          size    = 3,
          colour  = "#222222",
          vjust   = -0.5,
          show.legend = FALSE
        )
      }
    }
  }

  p
}

# Null-coalescing operator (avoid importing rlang::`%||%`)
`%||%` <- function(a, b) if (!is.null(a)) a else b
