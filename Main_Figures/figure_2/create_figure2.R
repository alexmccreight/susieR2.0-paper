#!/usr/bin/env Rscript

# =============================================================================
# FIGURE 2: REAL DATA APPLICATION — ROSMAP eQTL Fine-Mapping
# =============================================================================
# Combined figure for susieR 2.0 paper (Figure 2)
# Loads pre-built panel plot objects from figure_2/ subdirectories
#
# Layout:
#   Row 1:  [ A (CS counts) ] [ B (TWAS CV R2) ] [ C (concordance) ] [ D (enrichment) ]
#   Row 2:  [ E (locus example, full width) ]
#
# Data: ROSMAP eQTL fine-mapping, 3 brain tissues (DLPFC, AC, PCC), 2,654 genes
# =============================================================================

library(ggplot2)
library(cowplot)

cat("Creating Figure 2: Real Data Application...\n\n")

# =============================================================================
# PATHS
# =============================================================================

source("R/paths.R")

fig2_dir   <- fig_dir(2)
output_dir <- fig2_dir

# =============================================================================
# LOAD PANEL PLOT OBJECTS
# =============================================================================

cat("Loading panel plot objects...\n")

# Panel A: single grouped bar chart
panel_A_plots <- readRDS(file.path(fig2_dir, "figure_2_panel_A", "panel_A_plots.rds"))
plot_A <- panel_A_plots$A
cat("  Panel A loaded.\n")

# Concordance panel (object plot_B, from figure_2_panel_B/; labelled "C")
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
# ROW 1: [A] + [B] + [C] + [D]  (labels; objects are plot_A, plot_C, plot_B, plot_D)
# =============================================================================
cat("\nAssembling Row 1...\n")

# Panels are laid out and LABELLED in reading order A, B, C, D.
#
# The plot objects keep the variable names of the directories they come from,
# so the mapping is deliberately crossed here:
#     plot_A  (figure_2_panel_A) -> label "A"   CS counts by size bin
#     plot_C  (figure_2_panel_C) -> label "B"   TWAS cross-validation R^2
#     plot_B  (figure_2_panel_B) -> label "C"   cross-method signal concordance
#     plot_D  (figure_2_panel_D) -> label "D"   AlphaGenome CS scores
# Putting concordance (C) next to the AlphaGenome scores (D) keeps the two
# panels that share the Consensus / SuSiE / SuSiE-inf colouring adjacent.
#
# NOTE: this ordering and relabelling are specific to the main figure.
# Supplementary S11 keeps the original three-method layout and its own
# A B C / D E labels — see plot_S11.R.
row1 <- plot_grid(
  plot_A, plot_C, plot_B, plot_D,
  nrow = 1,
  rel_widths = c(1.85, 1.45, 1.45, 1.95),
  align = "h", axis = "bt",
  labels = c("A", "B", "C", "D"),
  label_size = 20, label_fontface = "bold"
)

cat("Row 1 complete.\n")

# =============================================================================
# ROW 2: [E] (full width)
# =============================================================================
cat("Assembling Row 2...\n")

# Panel E gets its own full-width row: it is a 12-track vertical stack and was
# being squeezed beside D.
row2 <- plot_grid(
  plot_E,
  nrow = 1,
  labels = c("E"),
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
  rel_heights = c(1.05, 1.65)
)

# =============================================================================
# SAVE
# =============================================================================

output_pdf <- file.path(output_dir, "susieR2_figure2.pdf")
ggsave(output_pdf, final_figure,
       width = 18, height = 16.5, units = "in",
       bg = "white", limitsize = FALSE)
cat(sprintf("\nPDF saved to: %s\n", output_pdf))

output_png <- file.path(output_dir, "susieR2_figure2.png")
ggsave(output_png, final_figure,
       width = 18, height = 16.5, units = "in", dpi = 300,
       bg = "white", limitsize = FALSE)
cat(sprintf("PNG saved to: %s\n", output_png))

cat("\nDone!\n")
