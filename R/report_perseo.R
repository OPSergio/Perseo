#' Generate an interactive HTML report from a PERSEO analysis
#'
#' Renders a self-contained HTML report from the output of \code{run_perseo()}.
#' The report includes KPI summaries, family selection diagnostics, differential
#' expression tables, interactive volcano and MA plots (powered by
#' \pkg{ggiraph}), and a pairwise contrast switcher.
#'
#' @param x Output of \code{run_perseo()}.
#' @param output_file Character, filename for the HTML report
#'   (default: \code{"perseo_report.html"}).
#' @param output_dir Character, directory where the report is written
#'   (default: current working directory \code{"."}).
#' @param title Character, report title shown in the header
#'   (default: \code{"PERSEO Analysis Report"}).
#' @param open Logical, open the report in the default browser after rendering
#'   (default: \code{TRUE}).
#' @param quiet Logical, suppress rmarkdown render messages (default: \code{FALSE}).
#'
#' @return Invisibly returns the path to the generated HTML file.
#'
#' @details
#' Requires the following packages: \pkg{rmarkdown}, \pkg{ggiraph}, \pkg{DT},
#' \pkg{htmltools}. Install them with
#' \code{install.packages(c("rmarkdown", "ggiraph", "DT", "htmltools"))}.
#'
#' @examples
#' \dontrun{
#' results <- run_perseo(counts, design, contrast_variable = "group",
#'                       metadata = meta, seed = 42)
#' report_perseo(results, output_file = "my_analysis.html", open = TRUE)
#' }
#'
#' @export
report_perseo <- function(x,
                          output_file = "perseo_report.html",
                          output_dir  = ".",
                          title       = "PERSEO Analysis Report",
                          open        = TRUE,
                          quiet       = FALSE) {

  for (pkg in c("rmarkdown", "DT", "htmltools", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for report_perseo(). ",
           "Install it with: install.packages('", pkg, "')")
    }
  }

  if (!inherits(x, "perseo_results") &&
      !all(c("family_selection", "differential_expression") %in% names(x))) {
    stop("'x' must be the output of run_perseo().")
  }

  template <- system.file("templates", "perseo_report.Rmd", package = "PERSEO")
  if (!nzchar(template)) {
    stop("Report template not found. Reinstall PERSEO.")
  }

  output_dir  <- normalizePath(output_dir, mustWork = FALSE)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  out_path <- file.path(output_dir, output_file)

  rmarkdown::render(
    input       = template,
    output_file = output_file,
    output_dir  = output_dir,
    params      = list(result = x, report_title = title),
    envir       = new.env(parent = globalenv()),
    quiet       = quiet
  )

  if (open && interactive()) {
    utils::browseURL(out_path)
  }

  message("Report written to: ", out_path)
  invisible(out_path)
}
