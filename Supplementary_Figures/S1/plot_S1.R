#!/usr/bin/env Rscript

# Supplementary Figure S1: fold speedup of susieR 2.0 over susieR 1.0 for the
# lambda-regularized RSS model (susie_rss_lambda), at p = 5,000 variants across
# feature-sample ratios p/n.
#
# Speedup = median run time of susieR 1.0 (LD matrix input) divided by the
# median run time of susieR 2.0 given the same LD matrix or the genotype
# matrix; medians are over the replicates of each condition.
#
# Input:  data/rss_lambda_comparison_data/*.rds  (data_generation/)
# Output: S1.pdf, S1.png

library(ggplot2)

source("R/paths.R")

s1_dir <- supp_dir(1)
files  <- list.files(file.path(s1_dir, "data", "rss_lambda_comparison_data"),
                     pattern = "[.]rds$", full.names = TRUE)

# ---- Fold speedups ----

speedup <- do.call(rbind, lapply(files, function(f) {
  d <- readRDS(f)
  median_time <- function(key) median(sapply(d$replicates, `[[`, key))
  data.frame(
    n     = d$parameters$n,
    p     = d$parameters$p,
    input = c("LD Matrix", "Genotype Matrix"),
    fold  = median_time("time_1.0_R") /
      c(median_time("time_2.0_R"), median_time("time_2.0_X"))
  )
}))

# The folder also holds conditions with p < 5,000, which the figure does not show.
speedup <- speedup[speedup$p == 5000, ]
speedup$pn_ratio <- factor(as.integer(speedup$p / speedup$n))
speedup$input    <- factor(speedup$input, levels = c("LD Matrix", "Genotype Matrix"))

# ---- Plot ----

input_colors <- c(
  "LD Matrix"       = "#DF8A2C",
  "Genotype Matrix" = "#2E6B9E"
)

p_fold <- ggplot(speedup, aes(x = pn_ratio, y = fold, fill = input)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red",
             linewidth = 0.8) +
  scale_fill_manual(values = input_colors, name = NULL) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    breaks = seq(0, 18, by = 2)
  ) +
  labs(
    x = "Feature-Sample Ratio",
    y = "Speedup (fold)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(color = "black", size = 13),
    legend.position  = c(0.5, 0.97),
    legend.justification = c(0.5, 1),
    legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.15, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    plot.margin      = margin(8, 12, 8, 8)
  )

save_figure(p_fold, s1_dir, "S1", width = 6.5, height = 5)
