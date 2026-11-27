#!/usr/bin/env Rscript

# Figure 2: fine-mapping and gene expression prediction in ROSMAP.
#   A   number of credible sets for SuSiE and SuSiE-inf, by size
#   B   mean cross-validated R2 of the prediction methods
#   C   signals recovered by both methods or by one only
#   D   PIP-weighted AlphaGenome scores by cross-method group
#   E   the LMOD1 locus: PIPs and annotation tracks
#
# A-D share the top row; E spans the bottom row.
#
# Input:  figure_2_panel_A/panel_A_plots.rds,
#         figure_2_panel_{B-E}/panel_{B-E}_plot.rds  (plot_panel_{A-E}.R)
# Output: susieR2_figure2.pdf, susieR2_figure2.png

library(cowplot)

source("R/paths.R")

fig2_dir <- fig_dir(2)

read_panel <- function(folder, file) {
  readRDS(file.path(fig2_dir, paste0("figure_2_panel_", folder), file))
}

panel_A <- read_panel("A", "panel_A_plots.rds")$A
panel_B <- read_panel("B", "panel_B_plot.rds")
panel_C <- read_panel("C", "panel_C_plot.rds")
panel_D <- read_panel("D", "panel_D_plot.rds")
panel_E <- read_panel("E", "panel_E_plot.rds")

row1 <- plot_grid(
  panel_A, panel_B, panel_C, panel_D,
  nrow = 1,
  rel_widths = c(1.85, 1.45, 1.45, 1.95),
  align = "h", axis = "bt",
  labels = c("A", "B", "C", "D"),
  label_size = 20, label_fontface = "bold"
)

row2 <- plot_grid(panel_E, nrow = 1, labels = "E",
                  label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1.05, 1.65))

save_figure(final_figure, fig2_dir, "susieR2_figure2", width = 18, height = 16.5)
