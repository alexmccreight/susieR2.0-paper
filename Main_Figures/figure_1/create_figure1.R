#!/usr/bin/env Rscript

# =============================================================================
# FIGURE 1: SuSiE 2.0 Overview & Simulation Benchmarks
# =============================================================================
# Combined figure for susieR 2.0 paper (Figure 1)
# Loads pre-built panel plot objects from figure_1/ subdirectories
#
# Layout:
#   Row 1:  [ A: Software Architecture (full width) ]              ~55%
#   Row 2:  [ B: Runtime (placeholder) ] [ C: Complex S1 Power ] [ D: Complex S1 FDR ] [ E: Complex S2 Power ] [ F: Complex S2 FDR ]  ~45%
# =============================================================================

library(ggplot2)
library(cowplot)
library(png)
library(grid)

cat("Creating Figure 1: SuSiE 2.0 Overview & Simulation Benchmarks...\n\n")

# =============================================================================
# PATHS
# =============================================================================

base_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark"
fig1_dir <- file.path(base_dir, "final_scripts/figure_1")
output_dir <- file.path(base_dir, "final_scripts/figure_1")

# =============================================================================
# LOAD PANEL PLOT OBJECTS
# =============================================================================

cat("Loading panel plot objects...\n")

# Panel A: Architecture diagram (loaded as PNG — placeholder, will be expanded)
panel_A_path <- file.path(fig1_dir, "figure_1_panel_A", "panel_A.png")
panel_A_img <- readPNG(panel_A_path)

img_aspect <- ncol(panel_A_img) / nrow(panel_A_img)
cat(sprintf("  Panel A native aspect ratio: %.3f (w/h)\n", img_aspect))

panel_A <- ggdraw() +
  draw_grob(rasterGrob(panel_A_img,
                        width = unit(1, "npc"), height = unit(1, "npc"),
                        interpolate = TRUE))
cat("  Panel A loaded (from PNG).\n")

# Panel B: Runtime Benchmark
panel_B <- readRDS(file.path(fig1_dir, "figure_1_panel_B", "panel_B_plot.rds"))
cat("  Panel B loaded (Runtime Benchmark).\n")

# Panel C: Complex S1 Power
panel_C <- readRDS(file.path(fig1_dir, "figure_1_panel_C", "panel_C_plot.rds"))
cat("  Panel C loaded (Complex S1 Power).\n")

# Panel D: Complex S1 FDR
panel_D <- readRDS(file.path(fig1_dir, "figure_1_panel_D", "panel_D_plot.rds"))
cat("  Panel D loaded (Complex S1 FDR).\n")

# Panel E: Complex S2 Power
panel_E <- readRDS(file.path(fig1_dir, "figure_1_panel_E", "panel_E_plot.rds"))
cat("  Panel E loaded (Complex S2 Power).\n")

# Panel F: Complex S2 FDR
panel_F <- readRDS(file.path(fig1_dir, "figure_1_panel_F", "panel_F_plot.rds"))
cat("  Panel F loaded (Complex S2 FDR).\n")

# ── Increase font sizes for all ggplot panels ──
bigger_theme <- theme(
  axis.title = element_text(face = "bold", size = 18),
  axis.text  = element_text(color = "black", size = 15),
  strip.text = element_text(size = 16, face = "bold")
)
for (p_name in c("panel_B", "panel_C", "panel_D", "panel_E", "panel_F")) {
  p <- get(p_name)
  if (inherits(p, "gg")) {
    p$layers <- lapply(p$layers, function(layer) {
      geom_class <- class(layer$geom)[1]
      if (geom_class == "GeomText" && !is.null(layer$aes_params$size)) {
        layer$aes_params$size <- layer$aes_params$size * 2.2
      } else if (geom_class == "GeomPoint" && !is.null(layer$aes_params$size)) {
        layer$aes_params$size <- layer$aes_params$size * 1.6
      }
      layer
    })
    assign(p_name, p + bigger_theme)
  }
}
cat("  Font sizes increased for panels B-F.\n")

# =============================================================================
# ROW 1: [A: Architecture (full width, centerpiece)]
# =============================================================================
cat("\nAssembling Row 1...\n")

row1 <- panel_A

cat("Row 1 complete.\n")

# =============================================================================
# ROW 2: [B: Runtime (placeholder)] [C: Complex S1 Power] [D: Complex S1 FDR] [E: Complex S2 Power] [F: Complex S2 FDR]
# =============================================================================
cat("Assembling Row 2...\n")

row2 <- plot_grid(
  panel_B, panel_C, panel_D, panel_E, panel_F,
  nrow = 1,
  rel_widths = c(1, 1, 1, 1, 1),
  labels = c("B", "C", "D", "E", "F"),
  label_size = 16, label_fontface = "bold"
)

cat("Row 2 complete.\n")

# =============================================================================
# COMBINE ALL ROWS
# =============================================================================
cat("\nCombining rows into final Figure 1...\n")

# Architecture (55%) + data panels (45%)
final_figure <- plot_grid(
  row1, row2,
  ncol = 1,
  rel_heights = c(75, 25),
  labels = c("A", ""),
  label_size = 16, label_fontface = "bold"
)

# =============================================================================
# SAVE
# =============================================================================

fig_width <- 20
# Panel A gets 55% of total height. For no distortion:
# fig_width / (0.55 * fig_height) = img_aspect
# => fig_height = fig_width / (0.55 * img_aspect)
fig_height <- fig_width / (0.75 * img_aspect)
cat(sprintf("\nComputed figure dimensions: %.1f x %.1f in\n", fig_width, fig_height))

output_pdf <- file.path(output_dir, "susieR2_figure1.pdf")
ggsave(output_pdf, final_figure,
       width = fig_width, height = fig_height, units = "in",
       bg = "white", limitsize = FALSE)
cat(sprintf("PDF saved to: %s\n", output_pdf))

output_png <- file.path(output_dir, "susieR2_figure1.png")
ggsave(output_png, final_figure,
       width = fig_width, height = fig_height, units = "in", dpi = 300,
       bg = "white", limitsize = FALSE)
cat(sprintf("PNG saved to: %s\n", output_png))

cat("\nDone!\n")
