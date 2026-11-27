#!/usr/bin/env Rscript

# Supplementary Figure S10: the real-data comparison of main Figure 2 (ROSMAP
# eQTL fine-mapping) with all three methods, SuSiE, SuSiE-ash and SuSiE-inf.
#   A   credible set counts by size bin
#   B   TWAS cross-validation R^2
#   C   cross-method signal concordance
#   D   AlphaGenome CS scores by concordance group
#   E   locus example
#
# The panels are drawn by the Main_Figures/figure_2 panel scripts run with all
# three methods, which write into data/ (run_all.R lists the commands). Panel B
# is the panel of Figure 2 itself: SuSiE-ash is not among the prediction
# methods it compares.
#
# Input:  data/panel_A_plots.rds, data/panel_C_plot.rds, data/panel_D_plot.rds,
#         data/panel_E_plot.rds
#         Main_Figures/figure_2/figure_2_panel_B/panel_B_plot.rds
# Output: S10.pdf, S10.png

library(ggplot2)
library(cowplot)

source("R/paths.R")

s10_dir  <- supp_dir(10)
data_dir <- file.path(s10_dir, "data")

plot_A <- readRDS(file.path(data_dir, "panel_A_plots.rds"))$A
plot_B <- readRDS(file.path(fig_dir(2), "figure_2_panel_B", "panel_B_plot.rds"))
plot_C <- readRDS(file.path(data_dir, "panel_C_plot.rds"))
plot_D <- readRDS(file.path(data_dir, "panel_D_plot.rds"))
plot_E <- readRDS(file.path(data_dir, "panel_E_plot.rds"))

row1 <- plot_grid(
  plot_A, plot_B, plot_C,
  nrow = 1,
  rel_widths = c(2, 1.2, 2),
  align = "h", axis = "bt",
  labels = c("A", "B", "C"),
  label_size = 20, label_fontface = "bold"
)

row2 <- plot_grid(
  plot_D, plot_E,
  nrow = 1,
  rel_widths = c(1.2, 2),
  labels = c("D", "E"),
  label_size = 20, label_fontface = "bold"
)

final_figure <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.6))

save_figure(final_figure, s10_dir, "S10", width = 18, height = 16)
