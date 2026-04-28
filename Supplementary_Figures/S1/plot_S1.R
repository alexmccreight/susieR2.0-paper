#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S1: Sparse Fine-Mapping (8 panels)
# =============================================================================
# Row 1: Panel A (Power bar plot) + Panel B (FDR bar plot)
# Row 2: Panel C (ROC K=1) + Panel D (ROC K=2) + Panel E (ROC K=3) + Panel F (ROC K=4) + Panel G (ROC K=5)
# Row 3: Panel H (PIP cal SuSiE-ash) + Panel I (PIP cal SuSiE-inf) + Panel J (PIP cal SuSiE)
#
# Run extract_finemapping_data.R, extract_pip_roc.R, and
# extract_pip_calibration.R first.
#
# Input:  data/finemapping_data.rds, data/pip_roc_data.rds,
#         data/pip_calibration_data.rds
# Output: S1.pdf / .png
# =============================================================================

library(ggplot2)
library(cowplot)

# =============================================================================
# PATHS
# =============================================================================

s1_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S1"
data_dir <- file.path(s1_dir, "data")

bar_data_file <- file.path(data_dir, "finemapping_data.rds")
roc_data_file <- file.path(data_dir, "pip_roc_data.rds")
cal_data_file <- file.path(data_dir, "pip_calibration_data.rds")

if (!file.exists(bar_data_file)) stop("Run extract_finemapping_data.R first.")
if (!file.exists(roc_data_file)) stop("Run extract_pip_roc.R first.")
if (!file.exists(cal_data_file)) stop("Run extract_pip_calibration.R first.")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

bar_data <- readRDS(bar_data_file)
bar_data <- bar_data[bar_data$K %in% 1:5, ]
bar_data$Method <- factor(bar_data$Method,
                          levels = c("SuSiE", "SuSiE-ash", "SuSiE-inf"))

roc_raw <- readRDS(roc_data_file)
roc_curves <- roc_raw$curves
roc_auroc  <- roc_raw$auroc
roc_curves$Method <- factor(roc_curves$Method,
                            levels = c("SuSiE-ash", "SuSiE-inf", "SuSiE"))
roc_auroc$Method <- factor(roc_auroc$Method,
                           levels = c("SuSiE-ash", "SuSiE-inf", "SuSiE"))

cal_data <- readRDS(cal_data_file)

cat(sprintf("  Bar data: %d rows\n", nrow(bar_data)))
cat(sprintf("  ROC data: %d rows\n", nrow(roc_curves)))
cat(sprintf("  Cal data: %d rows\n", nrow(cal_data)))

# =============================================================================
# SHARED AESTHETICS
# =============================================================================

method_colors <- c(
  "SuSiE"     = "#4A90E2",
  "SuSiE-inf" = "#7CB342",
  "SuSiE-ash" = "#E53935"
)

theme_bar <- theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold", size = 18),
    axis.text    = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin  = margin(10, 10, 10, 10)
  )

theme_roc <- theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(color = "black", size = 13),
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none",
    plot.margin  = margin(10, 10, 10, 10)
  )

# =============================================================================
# PANELS A & B: Bar plots (Power / FDR)
# =============================================================================

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

  k_labels <- data.frame(Total_PVE = c(0.03, 0.06, 0.09, 0.12, 0.15), K = 1:5)

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

  for (i in seq_len(nrow(k_labels))) {
    p <- p + annotate("text", x = i + 0.05, y = y_limits[2] * 0.99,
                       label = paste0("K = ", k_labels$K[i]),
                       size = 4.5, fontface = "bold")
  }

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted",
                         color = "red", linewidth = 0.8)
  }

  return(p)
}

cat("Creating bar panels...\n")
panel_A <- create_bar_panel(bar_data, metric = "Power")
panel_B <- create_bar_panel(bar_data, metric = "FDR")

# =============================================================================
# PANELS C, D, E: PIP ROC curves (one per K)
# =============================================================================

create_roc_panel <- function(curves, auroc, K_val) {
  df <- curves[curves$K == K_val, ]
  auc_df <- auroc[auroc$K == K_val, ]

  # Build AUROC annotation labels
  auc_labels <- sprintf("%s: %.3f", auc_df$Method, auc_df$AUROC)
  auc_colors <- method_colors[as.character(auc_df$Method)]

  p <- ggplot(df, aes(x = FPR, y = TPR, color = Method)) +
    geom_vline(xintercept = 0.05, linetype = "dotted", color = "red", linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = method_colors) +
    labs(x = "Variant-level FPR", y = "Variant-level TPR", title = paste0("K = ", K_val)) +
    theme_roc +
    coord_cartesian(xlim = c(0, 0.10), ylim = c(0.7, 1))

  return(p)
}

cat("Creating ROC panels...\n")
panel_C <- create_roc_panel(roc_curves, roc_auroc, 1)
panel_D <- create_roc_panel(roc_curves, roc_auroc, 2)
panel_E <- create_roc_panel(roc_curves, roc_auroc, 3)
panel_F <- create_roc_panel(roc_curves, roc_auroc, 4)
panel_G <- create_roc_panel(roc_curves, roc_auroc, 5)

# =============================================================================
# PANELS F, G, H: PIP Calibration (one per method, pooled across K)
# =============================================================================

theme_cal <- theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(color = "black", size = 13),
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none",
    plot.margin  = margin(10, 10, 10, 10)
  )

create_cal_panel <- function(data, method_name, color) {
  df <- data[data$Method == method_name, ]

  ggplot(df, aes(x = bin_mid, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey60", linewidth = 0.5) +
    geom_errorbar(aes(ymin = observed - observed_se, ymax = observed + observed_se),
                  width = 0.02, linewidth = 0.5, color = color) +
    geom_point(size = 3.5, color = color) +
    labs(x = "Expected PIP", y = "Observed Frequency") +
    theme_cal +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

cat("Creating calibration panels...\n")
panel_H <- create_cal_panel(cal_data, "SuSiE",     method_colors["SuSiE"])
panel_I <- create_cal_panel(cal_data, "SuSiE-ash", method_colors["SuSiE-ash"])
panel_J <- create_cal_panel(cal_data, "SuSiE-inf", method_colors["SuSiE-inf"])

# =============================================================================
# SHARED LEGEND
# =============================================================================

legend_plot <- ggplot(bar_data, aes(x = factor(Total_PVE), y = Power_mean, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = method_colors) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))

legend <- get_legend(legend_plot)

# =============================================================================
# COMBINE: Row 1 = A+B, Row 2 = C+D+E+F+G, Row 3 = H+I+J, Bottom = legend
# =============================================================================

cat("Assembling figure...\n")

row1 <- plot_grid(
  panel_A, panel_B,
  nrow = 1, rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 20, label_fontface = "bold"
)

row2 <- plot_grid(
  panel_C, panel_D, panel_E, panel_F, panel_G,
  nrow = 1, rel_widths = c(1, 1, 1, 1, 1),
  labels = c("C", "D", "E", "F", "G"),
  label_size = 20, label_fontface = "bold"
)

row3 <- plot_grid(
  panel_H, panel_I, panel_J,
  nrow = 1, rel_widths = c(1, 1, 1),
  labels = c("H", "I", "J"),
  label_size = 20, label_fontface = "bold"
)

final_figure <- plot_grid(
  row1, row2, row3, legend,
  ncol = 1,
  rel_heights = c(1, 0.85, 0.85, 0.06)
)

# =============================================================================
# SAVE
# =============================================================================

fig_width  <- 20
fig_height <- 17

out_pdf <- file.path(s1_dir, "S1.pdf")
ggsave(out_pdf, final_figure, width = fig_width, height = fig_height,
       units = "in", bg = "white")
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s1_dir, "S1.png")
ggsave(out_png, final_figure, width = fig_width, height = fig_height,
       units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
