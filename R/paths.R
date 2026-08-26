# =============================================================================
# Shared path resolution for susieR2.0-paper
# =============================================================================
# Replaces the hardcoded "/Users/alexmccreight/..." assignment that used to sit
# at the top of every figure script.
#
# Usage (scripts are run from the repository root):
#     source("R/paths.R")
#     fig1_dir <- fig_dir(1)
#
# Figure data that ships with the repo is reached via paper_root(). Raw
# simulation output is far too large to vendor (~22 GB) and stays outside the
# repo; reach it with benchmark_root(), overridable via the environment
# variable SUSIER2_BENCHMARK_ROOT.
# =============================================================================

.susieR2_cache <- new.env(parent = emptyenv())

# Walk upwards from `start` looking for the repo marker file.
.susieR2_find_root <- function(start = NULL) {
  if (is.null(start)) {
    # Under `Rscript path/to/script.R`, prefer the script's own location so the
    # scripts also work when invoked from outside the repo.
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    start <- if (length(file_arg) > 0) {
      dirname(normalizePath(file_arg[1], mustWork = FALSE))
    } else {
      getwd()
    }
  }

  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, ".susieR2-paper-root"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop("Could not locate the susieR2.0-paper repository root.\n",
           "  Searched upwards from: ", normalizePath(start, mustWork = FALSE), "\n",
           "  Expected to find the marker file '.susieR2-paper-root'.\n",
           "  Run scripts from inside the repository, e.g.:\n",
           "    Rscript Main_Figures/figure_1/figure_1_panel_C/plot_panel_C.R",
           call. = FALSE)
    }
    d <- parent
  }
}

#' Absolute path to the repository root.
paper_root <- function() {
  if (is.null(.susieR2_cache$paper_root)) {
    .susieR2_cache$paper_root <- .susieR2_find_root()
  }
  .susieR2_cache$paper_root
}

#' Absolute path to the susieR2.0-benchmark tree holding raw simulation output.
#'
#' Only the handful of `extract_*` scripts need this. Override the location with
#' the SUSIER2_BENCHMARK_ROOT environment variable.
#'
#' @param required If TRUE (default), error when the directory is absent.
benchmark_root <- function(required = TRUE) {
  p <- Sys.getenv("SUSIER2_BENCHMARK_ROOT", unset = "")
  if (!nzchar(p)) {
    p <- file.path(dirname(paper_root()), "susieR2.0-benchmark")
  }
  p <- path.expand(p)
  if (required && !dir.exists(p)) {
    stop("Benchmark data directory not found: ", p, "\n",
         "  This script needs the raw simulation output, which is not vendored\n",
         "  into this repository because of its size (~22 GB).\n",
         "  Point at it with:  export SUSIER2_BENCHMARK_ROOT=/path/to/susieR2.0-benchmark",
         call. = FALSE)
  }
  normalizePath(p, mustWork = FALSE)
}

#' Absolute path to the development susieR package source.
#'
#' Used by scripts that fit models with the in-development package rather than
#' the installed one. Override with the SUSIER_ROOT environment variable.
susier_root <- function(required = TRUE) {
  p <- Sys.getenv("SUSIER_ROOT", unset = "")
  if (!nzchar(p)) p <- file.path(dirname(paper_root()), "susieR")
  p <- path.expand(p)
  if (required && !dir.exists(p)) {
    stop("susieR source directory not found: ", p, "\n",
         "  Point at it with:  export SUSIER_ROOT=/path/to/susieR", call. = FALSE)
  }
  normalizePath(p, mustWork = FALSE)
}

#' Directory for main figure `n`, e.g. fig_dir(1).
fig_dir <- function(n) file.path(paper_root(), "Main_Figures", paste0("figure_", n))

#' Directory for supplementary figure `n`, e.g. supp_dir(10).
supp_dir <- function(n) file.path(paper_root(), "Supplementary_Figures", paste0("S", n))
