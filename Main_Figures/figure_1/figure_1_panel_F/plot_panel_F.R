#!/usr/bin/env Rscript

# =============================================================================
# Figure 1, Panel F — Complex S2 FDR vs Top N causal variants
# =============================================================================
# Line plot with points and error bars showing 95% CS FDR against the number
# of top-|beta| variants treated as causal (N = 3, 5, 10, 15, 20, 25).
# The panel is drawn by plot_topn_panel() in R/topn_panels.R, shared with the
# three-method version in Supplementary Figure S10.
#
# Input:  data/complex_S2_finemapping_summary.rds  (from extract_complex_ldext.R; un-extended credible sets)
# Output: panel_F_plot.rds, panel_F.pdf, panel_F.png
# =============================================================================

library(ggplot2)

source("R/paths.R")
source("R/aesthetics.R")
source("R/topn_panels.R")

fig1_dir   <- fig_dir(1)
script_dir <- file.path(fig1_dir, "figure_1_panel_F")
data_path  <- file.path(fig1_dir, "data", "complex_S2_finemapping_summary.rds")

cat("Loading complex fine-mapping summary (S2)...\n")
if (!file.exists(data_path)) stop("Data file not found: ", data_path)
complex_data <- readRDS(data_path)

# Main figure shows SuSiE + SuSiE-inf only. The three-method version that
# retains SuSiE-ash is Supplementary Figure S10 (built by
# Supplementary_Figures/S10/build_S10_panels.R from the same summary file).
cat("Creating Panel F (Complex S2 FDR)...\n")
panel_F <- plot_topn_panel(complex_data, metric = "FDR", methods = MAIN_METHODS)

output_pdf <- file.path(script_dir, "panel_F.pdf")
ggsave(output_pdf, panel_F, width = 6, height = 5, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(script_dir, "panel_F.png")
ggsave(output_png, panel_F, width = 6, height = 5, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

plots_path <- file.path(script_dir, "panel_F_plot.rds")
saveRDS(panel_F, plots_path)
cat(sprintf("Saved: %s\n", plots_path))

cat("\nDone!\n")
