#!/usr/bin/env Rscript

# Supplementary Figure S6: 95% CS power and FDR against the number of top
# variants treated as causal, for SuSiE, SuSiE-ash and SuSiE-inf (the
# three-method version of Figure 1C-F).
#   A, B  power and FDR, moderate infinitesimal background (Complex S1)
#   C, D  power and FDR, extensive infinitesimal background (Complex S2)
# Panels carry no titles: the y-axis label distinguishes power from FDR.
#
# Input:  data/panel_{C,D,E,F}_plot.rds  (build_S6_panels.R)
# Output: S6.pdf, S6.png

library(ggplot2)
library(cowplot)

source("R/paths.R")

s6_dir <- supp_dir(6)

# The files are named after the corresponding Figure 1 panels.
panel_files <- c(A = "panel_C_plot.rds",   # Complex S1 Power
                 B = "panel_D_plot.rds",   # Complex S1 FDR
                 C = "panel_E_plot.rds",   # Complex S2 Power
                 D = "panel_F_plot.rds")   # Complex S2 FDR

# Larger text than the main-figure panel theme, for a standalone figure.
supp_theme <- theme(
  axis.title  = element_text(face = "bold", size = 13),
  axis.text   = element_text(color = "black", size = 11),
  plot.margin = margin(8, 8, 8, 8)
)

panels <- lapply(panel_files, function(f) {
  readRDS(file.path(s6_dir, "data", f)) + supp_theme
})

# Shared legend, taken from the first panel with its legend switched back on
# (the panel theme sets legend.position = "none").
legend <- get_legend(
  panels$A +
    theme(legend.position  = "bottom",
          legend.title     = element_blank(),
          legend.text      = element_text(size = 13),
          legend.key.size  = unit(0.7, "cm"))
)

grid_2x2 <- plot_grid(
  plotlist = panels, nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 18, label_fontface = "bold"
)

final_figure <- plot_grid(grid_2x2, legend, ncol = 1, rel_heights = c(1, 0.07))

save_figure(final_figure, s6_dir, "S6", width = 12, height = 10)
