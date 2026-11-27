#!/usr/bin/env Rscript

# Rebuild every figure from the data in this repository.
# Run from the repository root:  Rscript run_all.R
#
# Not run here, because they need inputs from outside the repository (see
# README.md): the extract_*, fit_models and simulate_* scripts, Figure 2E's
# prepare_gene_models.R (Ensembl REST API), and the Figure 1A diagram (XeLaTeX).
#
# A step may carry command-line arguments after the script path; they are
# passed to Rscript unchanged.

source("R/paths.R")

steps <- c(
  # ---- Figure 1 ----
  "Main_Figures/figure_1/figure_1_panel_B/plot_panel_B.R",
  "Main_Figures/figure_1/figure_1_panel_C/plot_panel_C.R",
  "Main_Figures/figure_1/figure_1_panel_D/plot_panel_D.R",
  "Main_Figures/figure_1/figure_1_panel_E/plot_panel_E.R",
  "Main_Figures/figure_1/figure_1_panel_F/plot_panel_F.R",
  "Main_Figures/figure_1/create_figure1.R",

  # ---- Figure 2 ----
  # Cross-method matching for SuSiE and SuSiE-inf: writes
  # alphagenome_cs_group_{assignments,scores}__standard_inf.*
  "Main_Figures/figure_2/data/data_generation/alphagenome_scoring/alphagenome_cs_group_comparison.R",

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
  "Supplementary_Figures/S6/build_S6_panels.R",
  "Supplementary_Figures/S6/plot_S6.R",
  "Supplementary_Figures/S7/plot_S7.R",
  "Supplementary_Figures/S8/plot_S8.R",
  "Supplementary_Figures/S9/plot_S9.R",

  # S10: the Figure 2 scripts run with all three methods, writing to
  # Supplementary_Figures/S10/data. The cross-method matching then writes
  # alphagenome_cs_group_{assignments,scores}.* (no suffix). Panel B is the
  # panel of Figure 2 itself.
  paste("Main_Figures/figure_2/data/data_generation/alphagenome_scoring/alphagenome_cs_group_comparison.R",
        "standard,ash,inf"),
  "Main_Figures/figure_2/figure_2_panel_A/prepare_panel_A_data.R --include-ash",
  "Main_Figures/figure_2/figure_2_panel_A/plot_panel_A.R --include-ash",
  paste("Main_Figures/figure_2/figure_2_panel_C/prepare_panel_C_data.R",
        "Main_Figures/figure_2/data/alphagenome_cs_group_assignments.rds",
        "Supplementary_Figures/S10/data/panel_C_data.rds"),
  paste("Main_Figures/figure_2/figure_2_panel_C/plot_panel_C.R",
        "Supplementary_Figures/S10/data/panel_C_data.rds",
        "Supplementary_Figures/S10/data/panel_C"),
  paste("Main_Figures/figure_2/figure_2_panel_D/prepare_panel_D_data.R",
        "Main_Figures/figure_2/data/alphagenome_cs_group_scores.csv",
        "Main_Figures/figure_2/data/alphagenome_cs_group_assignments.rds",
        "Supplementary_Figures/S10/data/panel_D_data.rds"),
  paste("Main_Figures/figure_2/figure_2_panel_D/plot_panel_D.R",
        "Supplementary_Figures/S10/data/panel_D_data.rds",
        "Supplementary_Figures/S10/data"),
  "Main_Figures/figure_2/figure_2_panel_E/prepare_panel_E_data.R --include-ash",
  "Main_Figures/figure_2/figure_2_panel_E/plot_panel_E.R --include-ash",
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
