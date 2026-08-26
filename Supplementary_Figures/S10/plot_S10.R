#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S10: Complex S1/S2 Power and FDR -- all three methods
# =============================================================================
# These four panels were Figure 1C-1F before the manuscript revision. Main
# Figure 1 now shows SuSiE and SuSiE-inf only; this supplement preserves the
# original three-method comparison including SuSiE-ash.
#
# Layout (2 x 2):
#   A: Complex S1 Power    B: Complex S1 FDR
#   C: Complex S2 Power    D: Complex S2 FDR
#
# Input:  data/panel_{C,D,E,F}_plot.rds
#           Snapshots of the pre-revision Figure 1 panel objects, taken before
#           the panel scripts were changed to drop SuSiE-ash. Because
#           Main_Figures/figure_1/data/complex_S{1,2}_finemapping_summary.rds
#           still contains all three methods, these are also regenerable: run
#           the Figure 1 panel scripts with MAIN_METHODS replaced by
#           ALL_METHODS (see R/aesthetics.R).
# Output: S10.pdf, S10.png
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(cowplot)
})

source("R/paths.R")
source("R/aesthetics.R")

s10_dir  <- supp_dir(10)
data_dir <- file.path(s10_dir, "data")

# =============================================================================
# Load the snapshotted panel objects
# =============================================================================

panel_files <- c(A = "panel_C_plot.rds",   # Complex S1 Power
                 B = "panel_D_plot.rds",   # Complex S1 FDR
                 C = "panel_E_plot.rds",   # Complex S2 Power
                 D = "panel_F_plot.rds")   # Complex S2 FDR

panel_titles <- c(A = "Complex S1: 95% CS Power",
                  B = "Complex S1: 95% CS FDR",
                  C = "Complex S2: 95% CS Power",
                  D = "Complex S2: 95% CS FDR")

panels <- lapply(names(panel_files), function(k) {
  p <- file.path(data_dir, panel_files[[k]])
  if (!file.exists(p)) stop("Missing snapshot: ", p)
  readRDS(p)
})
names(panels) <- names(panel_files)

# Sanity check: this supplement exists to show all three methods.
for (k in names(panels)) {
  m <- unique(as.character(panels[[k]]$data$Method))
  if (!any(grepl("ash", m, ignore.case = TRUE))) {
    stop("Panel ", k, " has no SuSiE-ash series; the snapshot in data/ looks ",
         "like a post-revision (two-method) object.")
  }
}

# =============================================================================
# Restyle for a standalone supplementary figure
# =============================================================================

supp_theme <- theme(
  axis.title  = element_text(face = "bold", size = 13),
  axis.text   = element_text(color = "black", size = 11),
  plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
  plot.margin = margin(8, 8, 8, 8)
)

panels <- lapply(names(panels), function(k) {
  panels[[k]] + ggtitle(panel_titles[[k]]) + supp_theme
})
names(panels) <- names(panel_files)

# Shared legend, taken from the first panel with its legend switched back on
# (the panel theme sets legend.position = "none" for the main-figure layout).
legend <- get_legend(
  panels$A +
    theme(legend.position  = "bottom",
          legend.title     = element_blank(),
          legend.text      = element_text(size = 13),
          legend.key.size  = unit(0.7, "cm"))
)

# =============================================================================
# Assemble + save
# =============================================================================

grid_2x2 <- plot_grid(
  panels$A, panels$B, panels$C, panels$D,
  nrow = 2, ncol = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 18, label_fontface = "bold"
)

final_figure <- plot_grid(grid_2x2, legend, ncol = 1, rel_heights = c(1, 0.07))

out_pdf <- file.path(s10_dir, "S10.pdf")
ggsave(out_pdf, final_figure, width = 12, height = 10, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s10_dir, "S10.png")
ggsave(out_png, final_figure, width = 12, height = 10, units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
