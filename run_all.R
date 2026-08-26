#!/usr/bin/env Rscript

# =============================================================================
# Rebuild every figure that can be built from vendored data
# =============================================================================
# Run from the repository root:
#     Rscript run_all.R
#
# Skipped deliberately:
#   - extract_* / fit_models scripts, which need raw simulation output outside
#     this repository (see SUSIER2_BENCHMARK_ROOT in README.md)
#   - Figure 2 panel E's prepare_gene_models.R, which queries the Ensembl REST API
#   - Figure 1 panel A, which is a TikZ diagram requiring XeLaTeX
# =============================================================================

source("R/paths.R")

steps <- c(
  # ---- Figure 1 (SuSiE + SuSiE-inf) ----
  "Main_Figures/figure_1/figure_1_panel_B/plot_panel_B.R",
  "Main_Figures/figure_1/figure_1_panel_C/plot_panel_C.R",
  "Main_Figures/figure_1/figure_1_panel_D/plot_panel_D.R",
  "Main_Figures/figure_1/figure_1_panel_E/plot_panel_E.R",
  "Main_Figures/figure_1/figure_1_panel_F/plot_panel_F.R",
  "Main_Figures/figure_1/create_figure1.R",

  # ---- Figure 2 concordance analysis (SuSiE + SuSiE-inf) ----
  # Regenerates alphagenome_cs_group_{assignments,scores}__standard_inf.*
  "Main_Figures/figure_2/data/data_generation/alphagenome_scoring/alphagenome_cs_group_comparison.R",

  # ---- Figure 2 panels ----
  "Main_Figures/figure_2/figure_2_panel_A/prepare_panel_A_data.R",
  "Main_Figures/figure_2/figure_2_panel_A/plot_panel_A.R",
  "Main_Figures/figure_2/figure_2_panel_B/prepare_panel_B_data.R",
  "Main_Figures/figure_2/figure_2_panel_B/plot_panel_B.R",
  "Main_Figures/figure_2/figure_2_panel_C/prepare_panel_C_data.R",
  "Main_Figures/figure_2/figure_2_panel_C/plot_panel_C.R",
  "Main_Figures/figure_2/figure_2_panel_D/prepare_panel_D_data.R",
  "Main_Figures/figure_2/figure_2_panel_D/plot_panel_D.R",
  "Main_Figures/figure_2/figure_2_panel_E/prepare_panel_E_data.R",
  "Main_Figures/figure_2/figure_2_panel_E/plot_panel_E.R",
  "Main_Figures/figure_2/create_figure2.R",

  # ---- Supplementary figures ----
  "Supplementary_Figures/S1/plot_S1.R",
  "Supplementary_Figures/S2/plot_S2.R",
  "Supplementary_Figures/S3/plot_S3.R",
  "Supplementary_Figures/S4/plot_S4.R",
  "Supplementary_Figures/S5/plot_S5.R",
  "Supplementary_Figures/S6/plot_S6.R",
  "Supplementary_Figures/S7/plot_S7.R",
  "Supplementary_Figures/S8/plot_S8.R",
  "Supplementary_Figures/S9/plot_S9.R",
  "Supplementary_Figures/S10/plot_S10.R",
  "Supplementary_Figures/S11/plot_S11.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
failed <- character()

for (s in steps) {
  cat(sprintf("\n=== %s ===\n", s))
  status <- system2(rscript, s, stdout = NULL, stderr = NULL)
  if (!identical(status, 0L)) {
    failed <- c(failed, s)
    cat("  FAILED\n")
  } else {
    cat("  ok\n")
  }
}

cat("\n============================================================\n")
if (length(failed) == 0) {
  cat(sprintf("All %d steps completed.\n", length(steps)))
} else {
  cat(sprintf("%d of %d steps FAILED:\n", length(failed), length(steps)))
  for (f in failed) cat("  ", f, "\n")
  quit(status = 1)
}
