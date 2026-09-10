#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S10: build the three-method Top-N panels from vendored data
# =============================================================================
# Figure 1 panels C-F draw SuSiE and SuSiE-inf from
# Main_Figures/figure_1/data/complex_S{1,2}_finemapping_summary.rds, which also
# contain SuSiE-ash. This script draws the same four panels with all three
# methods (plot_topn_panel() from R/topn_panels.R) and saves the plot objects
# that plot_S10.R assembles, so S10 is regenerable rather than a snapshot.
#
# Input:  Main_Figures/figure_1/data/complex_S{1,2}_finemapping_summary.rds
# Output: data/panel_{C,D,E,F}_plot.rds
#           C = Complex S1 Power   D = Complex S1 FDR
#           E = Complex S2 Power   F = Complex S2 FDR
# =============================================================================

library(ggplot2)

source("R/paths.R")
source("R/aesthetics.R")
source("R/topn_panels.R")

fig1_data <- file.path(fig_dir(1), "data")
out_dir   <- file.path(supp_dir(10), "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summ <- list(
  S1 = readRDS(file.path(fig1_data, "complex_S1_finemapping_summary.rds")),
  S2 = readRDS(file.path(fig1_data, "complex_S2_finemapping_summary.rds"))
)

panels <- list(
  C = plot_topn_panel(summ$S1, metric = "Power", methods = ALL_METHODS),
  D = plot_topn_panel(summ$S1, metric = "FDR",   methods = ALL_METHODS),
  E = plot_topn_panel(summ$S2, metric = "Power", methods = ALL_METHODS),
  F = plot_topn_panel(summ$S2, metric = "FDR",   methods = ALL_METHODS)
)

for (k in names(panels)) {
  out <- file.path(out_dir, sprintf("panel_%s_plot.rds", k))
  saveRDS(panels[[k]], out)
  cat(sprintf("Saved: %s\n", out))
}

cat("\nDone!\n")
