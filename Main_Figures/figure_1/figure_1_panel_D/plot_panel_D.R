#!/usr/bin/env Rscript

# =============================================================================
# Figure 1, Panel D — Complex Simulation FDR (S1: nOligo=5, nInf=all)
# =============================================================================
# Line plot with points and error bars showing 95% CS FDR vs Top N
# causal variants (3, 5, 7, 9, 11, 23). Includes axis break and FDR=0.05 line.
#
# Input:  complex_finemapping_summary.rds (from complex_results_S1/)
# Output: panel_D_plot.rds, panel_D.pdf, panel_D.png
# =============================================================================

library(ggplot2)
library(dplyr)

# =============================================================================
# Paths
# =============================================================================

fig1_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_1"
script_dir <- file.path(fig1_dir, "figure_1_panel_D")
data_path  <- file.path(fig1_dir, "data", "complex_S1_finemapping_summary.rds")

# =============================================================================
# Shared aesthetics
# =============================================================================

method_colors <- c(
  "SuSiE"     = "#4A90E2",
  "SuSiE-inf" = "#7CB342",
  "SuSiE-ash" = "#E53935"
)

theme_panel <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# =============================================================================
# Load data
# =============================================================================

cat("Loading complex fine-mapping summary (S1)...\n")
if (!file.exists(data_path)) stop("Data file not found: ", data_path)
complex_data <- readRDS(data_path)
complex_data$Method <- factor(complex_data$Method,
                              levels = c("SuSiE", "SuSiE-ash", "SuSiE-inf"))

# =============================================================================
# Create plot
# =============================================================================

cat("Creating Panel D (Complex S1 FDR)...\n")

y_limits <- c(0, 0.35)

complex_data <- complex_data %>%
  mutate(
    N_position = case_when(
      N == 25 ~ 1, N == 20 ~ 2, N == 15 ~ 3,
      N == 10 ~ 4, N == 5 ~ 5, N == 3 ~ 6,
      TRUE ~ NA_real_
    )
  )

panel_D <- ggplot(complex_data, aes(x = N_position, y = FDR_mean,
                                    color = Method, group = Method)) +
  geom_line(linewidth = 1.0) +
  geom_errorbar(
    aes(ymin = FDR_mean - FDR_se, ymax = FDR_mean + FDR_se),
    width = 0.25, linewidth = 0.5
  ) +
  geom_point(size = 3.5) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red", linewidth = 0.8) +
  scale_color_manual(values = method_colors, breaks = c("SuSiE", "SuSiE-ash", "SuSiE-inf")) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c("25", "20", "15", "10", "5", "3")
  ) +
  labs(x = "Top N as Causal Variants", y = "95% CS FDR") +
  theme_panel +
  coord_cartesian(ylim = y_limits, xlim = c(0.5, 6.5))

# =============================================================================
# Save
# =============================================================================

output_pdf <- file.path(script_dir, "panel_D.pdf")
ggsave(output_pdf, panel_D, width = 6, height = 5, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(script_dir, "panel_D.png")
ggsave(output_png, panel_D, width = 6, height = 5, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

plots_path <- file.path(script_dir, "panel_D_plot.rds")
saveRDS(panel_D, plots_path)
cat(sprintf("Saved: %s\n", plots_path))

cat("\nDone!\n")
