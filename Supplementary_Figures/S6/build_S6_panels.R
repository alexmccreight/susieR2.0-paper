#!/usr/bin/env Rscript

# Supplementary Figure S6: the four panels, drawn with SuSiE, SuSiE-ash and
# SuSiE-inf by plot_topn_panel() from the summary tables behind Figure 1C-F
# (which show SuSiE and SuSiE-inf only). Each output file carries the letter of
# the Figure 1 panel it corresponds to.
#
# Input:  Main_Figures/figure_1/data/complex_S{1,2}_finemapping_summary.rds
#           complex_S1 = moderate, complex_S2 = extensive infinitesimal background
# Output: data/panel_{C,D,E,F}_plot.rds
#           C = Complex S1 Power   D = Complex S1 FDR   (S6 panels A, B)
#           E = Complex S2 Power   F = Complex S2 FDR   (S6 panels C, D)

source("R/paths.R")
source("R/aesthetics.R")
source("R/topn_panels.R")

fig1_data <- file.path(fig_dir(1), "data")
out_dir   <- file.path(supp_dir(6), "data")

summ <- list(
  S1 = readRDS(file.path(fig1_data, "complex_S1_finemapping_summary.rds")),
  S2 = readRDS(file.path(fig1_data, "complex_S2_finemapping_summary.rds"))
)

panels <- list(
  C = plot_topn_panel(summ$S1, metric = "Power", methods = ALL_METHODS),
  D = plot_topn_panel(summ$S1, metric = "FDR",   methods = ALL_METHODS),
  E = plot_topn_panel(summ$S2, metric = "Power", methods = ALL_METHODS),
  F = plot_topn_panel(summ$S2, metric = "FDR",   methods = ALL_METHODS)
)

for (k in names(panels)) {
  out <- file.path(out_dir, sprintf("panel_%s_plot.rds", k))
  saveRDS(panels[[k]], out)
  cat("Saved:", out, "\n")
}
