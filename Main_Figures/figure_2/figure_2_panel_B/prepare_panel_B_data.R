#!/usr/bin/env Rscript

# Figure 2, panel B (data): cross-validated R2 of six prediction methods in
# ROSMAP.
#
# Keeps SuSiE, SuSiE-inf, Elastic Net, LASSO, BayesR and BayesL (SuSiE-ash and
# mr.ash are left out), renames them for display and orders them by decreasing
# mean R2.
#
# Input:  data/twas_cv_all_metrics_combined.rds
# Output: figure_2_panel_B/panel_B_data.rds

source("R/paths.R")

fig2_dir    <- fig_dir(2)
input_file  <- file.path(fig2_dir, "data", "twas_cv_all_metrics_combined.rds")
output_file <- file.path(fig2_dir, "figure_2_panel_B", "panel_B_data.rds")

method_labels <- c(
  "susie"     = "SuSiE",
  "susie_inf" = "SuSiE-inf",
  "enet"      = "Elastic Net",
  "lasso"     = "LASSO",
  "bayes_r"   = "BayesR",
  "bayes_l"   = "BayesL"
)

d <- readRDS(input_file)
d <- d[d$method %in% names(method_labels), ]
d$method <- method_labels[d$method]

mean_rsq <- tapply(d$rsq, d$method, mean, na.rm = TRUE)
method_order <- names(sort(mean_rsq, decreasing = TRUE))
d$method <- factor(d$method, levels = method_order)

saveRDS(d, output_file)
cat("Saved:", output_file, "\n")
