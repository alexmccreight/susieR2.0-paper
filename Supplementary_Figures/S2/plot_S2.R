#!/usr/bin/env Rscript

# Supplementary Figure S2: fine-mapping in the sparse simulations.
#   A, B   95% CS power and FDR against total PVE (K = 1-5 causal variants)
#   C-G    variant-level ROC curves, one panel per K
#   H-J    PIP calibration for SuSiE, SuSiE-ash and SuSiE-inf, pooled over K
#
# Input:  data/finemapping_data.rds      (extract_finemapping_data.R)
#         data/pip_roc_data.rds,
#         data/pip_calibration_data.rds  (extract_pip_data.R)
# Output: S2.pdf, S2.png

library(ggplot2)
library(cowplot)

source("R/paths.R")
source("R/aesthetics.R")

s2_dir   <- supp_dir(2)
data_dir <- file.path(s2_dir, "data")

bar_data   <- readRDS(file.path(data_dir, "finemapping_data.rds"))
roc_curves <- readRDS(file.path(data_dir, "pip_roc_data.rds"))$curves
cal_data   <- readRDS(file.path(data_dir, "pip_calibration_data.rds"))

bar_data$Method <- factor(bar_data$Method, levels = ALL_METHODS)
# SuSiE is the last level so that its curve is drawn on top.
roc_curves$Method <- factor(roc_curves$Method,
                            levels = c("SuSiE-ash", "SuSiE-inf", "SuSiE"))

theme_bar <- theme_classic(base_size = 14) +
  theme(
    axis.title  = element_text(face = "bold", size = 18),
    axis.text   = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

theme_curve <- theme_classic(base_size = 14) +
  theme(
    axis.title  = element_text(face = "bold", size = 16),
    axis.text   = element_text(color = "black", size = 13),
    plot.title  = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

# ---- A, B: power and FDR ----

create_bar_panel <- function(data, metric = "Power") {
  if (metric == "Power") {
    y_label  <- "95% CS Power"
    y_mean   <- "Power_mean"
    y_se     <- "Power_se"
    y_limits <- c(0.3, 1.0)
  } else {
    y_label  <- "95% CS FDR"
    y_mean   <- "FDR_mean"
    y_se     <- "FDR_se"
    y_limits <- c(0, 0.1)
  }

  p <- ggplot(data, aes(x = factor(Total_PVE), y = .data[[y_mean]], fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
    geom_errorbar(
      aes(ymin = .data[[y_mean]] - .data[[y_se]],
          ymax = .data[[y_mean]] + .data[[y_se]]),
      position = position_dodge(width = 0.8),
      width = 0.15, linewidth = 0.5
    ) +
    scale_fill_manual(values = method_colors) +
    scale_x_discrete(expand = expansion(add = c(0.3, 0.3))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Total PVE", y = y_label) +
    theme_bar +
    coord_cartesian(ylim = y_limits)

  for (K in 1:5) {
    p <- p + annotate("text", x = K + 0.05, y = y_limits[2] * 0.99,
                      label = paste0("K = ", K), size = 4.5, fontface = "bold")
  }

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted",
                        color = "red", linewidth = 0.8)
  }

  p
}

panel_A <- create_bar_panel(bar_data, metric = "Power")
panel_B <- create_bar_panel(bar_data, metric = "FDR")

# ---- C-G: ROC curves ----

roc_panels <- lapply(1:5, function(K) {
  ggplot(roc_curves[roc_curves$K == K, ], aes(x = FPR, y = TPR, color = Method)) +
    geom_vline(xintercept = 0.05, linetype = "dotted", color = "red", linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = method_colors) +
    labs(x = "Variant-level FPR", y = "Variant-level TPR", title = paste0("K = ", K)) +
    theme_curve +
    coord_cartesian(xlim = c(0, 0.10), ylim = c(0.7, 1))
})

# ---- H-J: PIP calibration ----

cal_panels <- lapply(ALL_METHODS, function(method) {
  ggplot(cal_data[cal_data$Method == method, ], aes(x = bin_mid, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey60", linewidth = 0.5) +
    geom_errorbar(aes(ymin = observed - observed_se, ymax = observed + observed_se),
                  width = 0.02, linewidth = 0.5, color = method_colors[method]) +
    geom_point(size = 3.5, color = method_colors[method]) +
    labs(x = "Expected PIP", y = "Observed Frequency") +
    theme_curve +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
})

# ---- Assemble ----

legend <- get_legend(
  ggplot(bar_data, aes(x = factor(Total_PVE), y = Power_mean, fill = Method)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = method_colors) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14))
)

row1 <- plot_grid(panel_A, panel_B, nrow = 1,
                  labels = c("A", "B"), label_size = 20, label_fontface = "bold")
row2 <- plot_grid(plotlist = roc_panels, nrow = 1,
                  labels = c("C", "D", "E", "F", "G"), label_size = 20, label_fontface = "bold")
row3 <- plot_grid(plotlist = cal_panels, nrow = 1,
                  labels = c("H", "I", "J"), label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(row1, row2, row3, legend, ncol = 1,
                          rel_heights = c(1, 0.85, 0.85, 0.06))

save_figure(final_figure, s2_dir, "S2", width = 20, height = 17)
