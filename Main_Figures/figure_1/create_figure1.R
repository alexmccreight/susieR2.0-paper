#!/usr/bin/env Rscript

# Figure 1: SuSiE 2.0 architecture and performance.
#   A      architecture diagram (TikZ, rasterised to panel_A.png)
#   B      fold speedup of SuSiE-inf over its original implementation
#   C, D   95% CS power and FDR, moderate infinitesimal background
#   E, F   95% CS power and FDR, extensive infinitesimal background
#
# Panel A spans the full width and takes 75% of the height; B-F share the row
# below. The figure height is derived from the aspect ratio of panel_A.png so
# that the diagram is not distorted.
#
# Input:  figure_1_panel_A/panel_A.png
#         figure_1_panel_{B-F}/panel_{B-F}_plot.rds  (plot_panel_{B-F}.R)
# Output: susieR2_figure1.pdf, susieR2_figure1.png

library(ggplot2)
library(cowplot)
library(png)
library(grid)

source("R/paths.R")

fig1_dir <- fig_dir(1)

# ---- Panel A ----

panel_A_img <- readPNG(file.path(fig1_dir, "figure_1_panel_A", "panel_A.png"))
img_aspect  <- ncol(panel_A_img) / nrow(panel_A_img)

panel_A <- ggdraw() +
  draw_grob(rasterGrob(panel_A_img, width = unit(1, "npc"), height = unit(1, "npc")))

# ---- Panels B-F ----

# The panels are saved at their standalone size; text and points are enlarged
# for the composite. The legend settings only affect panel B, the one panel
# with a legend.
bigger_theme <- theme(
  axis.title  = element_text(face = "bold", size = 18),
  axis.text   = element_text(color = "black", size = 15),
  legend.text = element_text(size = 15),
  legend.key.size = unit(0.6, "cm")
)

panel_ids <- c("B", "C", "D", "E", "F")
panels <- lapply(panel_ids, function(id) {
  p <- readRDS(file.path(fig1_dir, paste0("figure_1_panel_", id),
                         paste0("panel_", id, "_plot.rds")))
  p$layers <- lapply(p$layers, function(layer) {
    if (class(layer$geom)[1] == "GeomPoint") {
      layer$aes_params$size <- layer$aes_params$size * 1.6
    }
    layer
  })
  p + bigger_theme
})

# ---- Assemble ----

a_share <- 0.75  # fraction of the figure height given to panel A

row2 <- plot_grid(plotlist = panels, nrow = 1, labels = panel_ids,
                  label_size = 16, label_fontface = "bold")

final_figure <- plot_grid(panel_A, row2, ncol = 1,
                          rel_heights = c(a_share, 1 - a_share),
                          labels = c("A", ""), label_size = 16, label_fontface = "bold")

# Panel A keeps its aspect ratio when fig_width / (a_share * fig_height) = img_aspect.
fig_width  <- 20
fig_height <- fig_width / (a_share * img_aspect)

save_figure(final_figure, fig1_dir, "susieR2_figure1",
            width = fig_width, height = fig_height, limitsize = FALSE)
