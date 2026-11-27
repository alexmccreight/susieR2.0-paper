#!/usr/bin/env Rscript

# Figure 1, panel B: fold speedup of SuSiE-inf over its original
# eigendecomposition-based implementation, across variant-to-sample ratios
# (p/n) at p = 5,000.
#
# Speedup is a ratio of median run times over replicates. The reference for
# both bars is the original implementation run on the LD matrix (time_old_LD):
#   LD Matrix        median(time_old_LD) / median(time_new_LD)   susie_rss(R)
#   Genotype Matrix  median(time_old_LD) / median(time_new_X)    susie(X, y)
#
# Input:  ../data/panel_B/susie_inf_comparison_data/benchmark_susieinf_*_nrep*.rds
# Output: panel_B.pdf, panel_B.png, panel_B_plot.rds

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("R/paths.R")
source("R/aesthetics.R")

fig1_dir    <- fig_dir(1)
script_dir  <- file.path(fig1_dir, "figure_1_panel_B")
results_dir <- file.path(fig1_dir, "data", "panel_B", "susie_inf_comparison_data")

# ---- Median run times ----

files <- list.files(results_dir, pattern = "^benchmark_susieinf_.*nrep[0-9]+\\.rds$",
                    full.names = TRUE)

df_all <- bind_rows(lapply(files, function(f) {
  d <- readRDS(f)
  reps <- d$replicates[!sapply(d$replicates, is.null)]
  med <- function(field) median(sapply(reps, `[[`, field), na.rm = TRUE)
  data.frame(
    n = d$parameters$n,
    p = d$parameters$p,
    median_old_LD = med("time_old_LD"),
    median_new_X  = med("time_new_X"),
    median_new_LD = med("time_new_LD")
  )
}))

# ---- Fold speedups ----

df_fold <- df_all %>%
  filter(p == 5000) %>%
  mutate(
    pn_ratio = as.integer(round(p / n)),
    fold_X   = median_old_LD / median_new_X,
    fold_LD  = median_old_LD / median_new_LD
  ) %>%
  select(pn_ratio, fold_LD, fold_X) %>%
  pivot_longer(cols = c(fold_LD, fold_X), names_to = "method", values_to = "fold") %>%
  mutate(
    method   = factor(method, levels = c("fold_LD", "fold_X"),
                      labels = c("LD Matrix", "Genotype Matrix")),
    pn_label = factor(pn_ratio)
  )

# ---- Plot ----

input_colors <- c("LD Matrix" = "#DF8A2C", "Genotype Matrix" = "#2E6B9E")

p_fold <- ggplot(df_fold, aes(x = pn_label, y = fold, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = input_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Feature-Sample Ratio", y = "Speedup (fold)") +
  theme_panel +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    plot.margin = margin(5, 10, 5, 5)
  )

save_figure(p_fold, script_dir, "panel_B", width = 4.5, height = 4)

plot_file <- file.path(script_dir, "panel_B_plot.rds")
saveRDS(p_fold, plot_file)
cat("Saved:", plot_file, "\n")
