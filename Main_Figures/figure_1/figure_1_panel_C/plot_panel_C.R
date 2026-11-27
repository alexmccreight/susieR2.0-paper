#!/usr/bin/env Rscript

# Figure 1, panel C: 95% CS power of SuSiE and SuSiE-inf against the number N
# of top-|beta| variants treated as causal, for oligogenic effects on a
# moderate infinitesimal background (simulation scenario "complex S1").
# The summary also holds SuSiE-ash, which plot_topn_panel() leaves out by
# default.
#
# Input:  data/complex_S1_finemapping_summary.rds  (extract_complex_ldext.R)
# Output: panel_C.pdf, panel_C.png, panel_C_plot.rds

source("R/paths.R")
source("R/aesthetics.R")
source("R/topn_panels.R")

fig1_dir   <- fig_dir(1)
script_dir <- file.path(fig1_dir, "figure_1_panel_C")

complex_data <- readRDS(file.path(fig1_dir, "data", "complex_S1_finemapping_summary.rds"))

panel_C <- plot_topn_panel(complex_data, metric = "Power")
save_figure(panel_C, script_dir, "panel_C", width = 6, height = 5)

plot_file <- file.path(script_dir, "panel_C_plot.rds")
saveRDS(panel_C, plot_file)
cat("Saved:", plot_file, "\n")
