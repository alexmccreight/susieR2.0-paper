# Paths and small I/O helpers shared by every script. Scripts are run from the
# repository root and start with source("R/paths.R").
#
# Data that ships with the repository is reached through paper_root(), fig_dir()
# and supp_dir(). Raw simulation output is too large to include and lives outside
# the repository; benchmark_root() finds it (see SUSIER2_BENCHMARK_ROOT).

.susieR2_cache <- new.env(parent = emptyenv())

# Walk upwards from `start` looking for the repo marker file.
.susieR2_find_root <- function(start = NULL) {
  if (is.null(start)) {
    # Prefer the script's own location, so scripts also run from outside the repo.
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
#' Needed only by the `extract_*` scripts. Defaults to ../susieR2.0-benchmark;
#' override with the SUSIER2_BENCHMARK_ROOT environment variable.
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

#' Absolute path to the susieR package source, for scripts that fit models with
#' pkgload. Defaults to ../susieR; override with the SUSIER_ROOT environment variable.
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

#' Directory for supplementary figure `n`, e.g. supp_dir(6).
supp_dir <- function(n) file.path(paper_root(), "Supplementary_Figures", paste0("S", n))

#' The i-th command-line argument, or `default` when it is absent or empty.
cli_arg <- function(i, default) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= i && nzchar(args[i])) args[i] else default
}

#' Save `plot` as <dir>/<name>.pdf and <dir>/<name>.png (300 dpi, white background).
save_figure <- function(plot, dir, name, width, height, ...) {
  for (ext in c("pdf", "png")) {
    out <- file.path(dir, paste0(name, ".", ext))
    ggplot2::ggsave(out, plot, width = width, height = height, units = "in",
                    dpi = 300, bg = "white", ...)
    cat("Saved:", out, "\n")
  }
}
