#!/usr/bin/env Rscript

# Figure 1, panel D: 95% CS FDR of SuSiE and SuSiE-inf against the number N
# of top-|beta| variants treated as causal, for oligogenic effects on a
# moderate infinitesimal background (simulation scenario "complex S1").
# The summary also holds SuSiE-ash, which plot_topn_panel() leaves out by
# default.
#
# Input:  data/complex_S1_finemapping_summary.rds  (extract_complex_ldext.R)
# Output: panel_D.pdf, panel_D.png, panel_D_plot.rds

source("R/paths.R")
source("R/aesthetics.R")
source("R/topn_panels.R")

fig1_dir   <- fig_dir(1)
script_dir <- file.path(fig1_dir, "figure_1_panel_D")

complex_data <- readRDS(file.path(fig1_dir, "data", "complex_S1_finemapping_summary.rds"))

panel_D <- plot_topn_panel(complex_data, metric = "FDR")
save_figure(panel_D, script_dir, "panel_D", width = 6, height = 5)

plot_file <- file.path(script_dir, "panel_D_plot.rds")
saveRDS(panel_D, plot_file)
cat("Saved:", plot_file, "\n")
