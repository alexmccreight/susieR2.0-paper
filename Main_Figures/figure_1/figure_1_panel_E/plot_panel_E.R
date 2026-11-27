#!/usr/bin/env Rscript

# Figure 1, panel E: 95% CS power of SuSiE and SuSiE-inf against the number N
# of top-|beta| variants treated as causal, for oligogenic effects on an
# extensive infinitesimal background (simulation scenario "complex S2").
# The summary also holds SuSiE-ash, which plot_topn_panel() leaves out by
# default.
#
# Input:  data/complex_S2_finemapping_summary.rds  (extract_complex_ldext.R)
# Output: panel_E.pdf, panel_E.png, panel_E_plot.rds

source("R/paths.R")
source("R/aesthetics.R")
source("R/topn_panels.R")

fig1_dir   <- fig_dir(1)
script_dir <- file.path(fig1_dir, "figure_1_panel_E")

complex_data <- readRDS(file.path(fig1_dir, "data", "complex_S2_finemapping_summary.rds"))

panel_E <- plot_topn_panel(complex_data, metric = "Power")
save_figure(panel_E, script_dir, "panel_E", width = 6, height = 5)

plot_file <- file.path(script_dir, "panel_E_plot.rds")
saveRDS(panel_E, plot_file)
cat("Saved:", plot_file, "\n")
