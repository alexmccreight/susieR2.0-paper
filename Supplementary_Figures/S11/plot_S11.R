#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S11: Real-data comparison including SuSiE-ash
# =============================================================================
# This was Figure 2 before the manuscript revision. Main Figure 2 now shows
# SuSiE and SuSiE-inf only; this supplement preserves the full three-method
# version.
#
# Layout matches the original Figure 2 exactly:
#   Row 1: [ A (CS counts) ] [ B (TWAS CV R2) ] [ C (concordance) ]
#   Row 2: [ D (AlphaGenome enrichment) ] [ E (locus example) ]
#
# Row 1 places TWAS R2 before concordance so the ordering matches main
# Figure 2. The plot objects keep the variable names of the directories they
# came from, so the mapping below is deliberately crossed:
#     plot_A -> label "A"   CS counts by size bin
#     plot_C -> label "B"   TWAS cross-validation R^2
#     plot_B -> label "C"   cross-method signal concordance
#
# Input:  data/panel_{A,B,C,D,E}_plot*.rds
#           Panels A-C retain the pre-revision Figure 2 plot objects.
#
#           Panel D (data/panel_D_plot.rds) is rebuilt from source by
#           run_all.R: prepare_panel_D_data.R and plot_panel_D.R are run on the
#           three-method scores/assignments in Main_Figures/figure_2/data/.
#           Panel B is also regenerable:
#             Rscript .../alphagenome_cs_group_comparison.R standard,ash,inf
#             Rscript .../prepare_panel_B_data.R <3-method assignments> <out>
#           Panel E is rebuilt from the vendored fit objects, gene models and
#           annotation BED files by the existing figure_2_panel_E scripts:
#             prepare_panel_E_data.R --include-ash
#             plot_panel_E.R --include-ash
#           run_all.R executes both before assembling this figure.
# Output: S11.pdf, S11.png
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(cowplot)
})

source("R/paths.R")
source("R/aesthetics.R")

s11_dir  <- supp_dir(11)
data_dir <- file.path(s11_dir, "data")

# =============================================================================
# Load the panel objects
# =============================================================================

need <- function(f) {
  p <- file.path(data_dir, f)
  if (!file.exists(p)) stop("Missing snapshot: ", p)
  readRDS(p)
}

plot_A <- need("panel_A_plots.rds")$A
plot_B <- need("panel_B_plot.rds")
plot_C <- need("panel_C_plot.rds")
plot_D <- need("panel_D_plot.rds")
plot_E <- need("panel_E_plot.rds")
cat("Loaded panels A-E.\n")

# Sanity check: this supplement exists to show all three methods, so panel A
# must still carry a SuSiE-ash series.
a_methods <- unique(as.character(plot_A$data$method))
if (!any(grepl("ash", a_methods, ignore.case = TRUE))) {
  stop("Panel A has no SuSiE-ash series; the snapshot in data/ looks like a ",
       "post-revision (two-method) object.")
}
cat(sprintf("  Panel A methods: %s\n", paste(sort(a_methods), collapse = ", ")))

# =============================================================================
# Assemble (identical geometry to the original create_figure2.R)
# =============================================================================

row1 <- plot_grid(
  plot_A, plot_C, plot_B,
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

# =============================================================================
# Save
# =============================================================================

out_pdf <- file.path(s11_dir, "S11.pdf")
ggsave(out_pdf, final_figure, width = 18, height = 16, units = "in",
       bg = "white", limitsize = FALSE)
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s11_dir, "S11.png")
ggsave(out_png, final_figure, width = 18, height = 16, units = "in",
       dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
