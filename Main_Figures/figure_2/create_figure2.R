#!/usr/bin/env Rscript

# =============================================================================
# FIGURE 2: REAL DATA APPLICATION — ROSMAP eQTL Fine-Mapping
# =============================================================================
# Combined figure for susieR 2.0 paper (Figure 2)
# Loads pre-built panel plot objects from figure_2/ subdirectories
#
# Layout:
#   Row 1:  [ A1 | A2 ]  [ B (UpSet) ]     [ C (TWAS) ]
#   Row 2:  [ D (Enrichment) ] [ E (ASH gene) ]  [ F (INF gene) ]
#
# Data: ROSMAP eQTL fine-mapping, 3 brain tissues (DLPFC, AC, PCC), 2,654 genes
# =============================================================================

library(ggplot2)
library(cowplot)

cat("Creating Figure 2: Real Data Application...\n\n")

# =============================================================================
# PATHS
# =============================================================================

base_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark"
fig2_dir <- file.path(base_dir, "final_scripts/figure_2")
output_dir <- file.path(base_dir, "final_scripts/figure_2")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD PANEL PLOT OBJECTS
# =============================================================================

cat("Loading panel plot objects...\n")

# Panel A: single grouped bar chart
panel_A_plots <- readRDS(file.path(fig2_dir, "figure_2_panel_A", "panel_A_plots.rds"))
plot_A <- panel_A_plots$A
cat("  Panel A loaded.\n")

# Panel B: UpSet plot
plot_B <- readRDS(file.path(fig2_dir, "figure_2_panel_B", "panel_B_plot.rds"))
cat("  Panel B loaded.\n")

# Panel C: TWAS prediction performance
plot_C <- readRDS(file.path(fig2_dir, "figure_2_panel_C", "panel_C_plot.rds"))
cat("  Panel C loaded.\n")

# Panel D: Functional enrichment
plot_D <- readRDS(file.path(fig2_dir, "figure_2_panel_D", "panel_D_plot.rds"))
cat("  Panel D loaded.\n")

# Panel E: Gene example (stacked PIP + annotation tracks)
plot_E <- readRDS(file.path(fig2_dir, "figure_2_panel_E", "panel_E_plot.rds"))
cat("  Panel E loaded.\n")

# =============================================================================
# ROW 1: [A1 | A2] + [B] + [C]
# =============================================================================
cat("\nAssembling Row 1...\n")

row1 <- plot_grid(
  plot_A, plot_B, plot_C,
  nrow = 1,
  rel_widths = c(2, 2, 1.2),
  align = "h", axis = "bt",
  labels = c("A", "B", "C"),
  label_size = 20, label_fontface = "bold"
)

cat("Row 1 complete.\n")

# =============================================================================
# ROW 2: [D] + [E] + [F]
# =============================================================================
cat("Assembling Row 2...\n")

row2 <- plot_grid(
  plot_D, plot_E,
  nrow = 1,
  rel_widths = c(1.2, 2),
  labels = c("D", "E"),
  label_size = 20, label_fontface = "bold"
)

cat("Row 2 complete.\n")

# =============================================================================
# COMBINE ALL ROWS
# =============================================================================
cat("\nCombining rows into final Figure 2...\n")

final_figure <- plot_grid(
  row1, row2,
  ncol = 1,
  rel_heights = c(1, 1.6)
)

# =============================================================================
# SAVE
# =============================================================================

output_pdf <- file.path(output_dir, "susieR2_figure2.pdf")
ggsave(output_pdf, final_figure,
       width = 18, height = 16, units = "in",
       bg = "white", limitsize = FALSE)
cat(sprintf("\nPDF saved to: %s\n", output_pdf))

output_png <- file.path(output_dir, "susieR2_figure2.png")
ggsave(output_png, final_figure,
       width = 18, height = 16, units = "in", dpi = 300,
       bg = "white", limitsize = FALSE)
cat(sprintf("PNG saved to: %s\n", output_png))

cat("\nDone!\n")
