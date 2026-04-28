#!/usr/bin/env Rscript

# =============================================================================
# Figure 1, Panel C — Complex Simulation Power (S1: nOligo=5, nInf=all)
# =============================================================================
# Line plot with points and error bars showing 95% CS Power vs Top N
# causal variants (3, 5, 7, 9, 11, 23). Includes axis break.
#
# Input:  complex_finemapping_summary.rds (from complex_results_S1/)
# Output: panel_C_plot.rds, panel_C.pdf, panel_C.png
# =============================================================================

library(ggplot2)
library(dplyr)

# =============================================================================
# Paths
# =============================================================================

fig1_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_1"
script_dir <- file.path(fig1_dir, "figure_1_panel_C")
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

cat("Creating Panel C (Complex S1 Power)...\n")

y_limits <- c(0, 0.9)

complex_data <- complex_data %>%
  mutate(
    N_position = case_when(
      N == 3 ~ 1, N == 5 ~ 2, N == 10 ~ 3,
      N == 15 ~ 4, N == 20 ~ 5, N == 25 ~ 6,
      TRUE ~ NA_real_
    )
  )

panel_C <- ggplot(complex_data, aes(x = N_position, y = Power_mean,
                                    color = Method, group = Method)) +
  geom_line(linewidth = 1.0) +
  geom_errorbar(
    aes(ymin = Power_mean - Power_se, ymax = Power_mean + Power_se),
    width = 0.25, linewidth = 0.5
  ) +
  geom_point(size = 3.5) +
  scale_color_manual(values = method_colors, breaks = c("SuSiE", "SuSiE-ash", "SuSiE-inf")) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c("3", "5", "10", "15", "20", "25")
  ) +
  labs(x = "Top N as Causal Variants", y = "95% CS Power") +
  theme_panel +
  coord_cartesian(ylim = y_limits, xlim = c(0.5, 6.5))

# =============================================================================
# Save
# =============================================================================

output_pdf <- file.path(script_dir, "panel_C.pdf")
ggsave(output_pdf, panel_C, width = 6, height = 5, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(script_dir, "panel_C.png")
ggsave(output_png, panel_C, width = 6, height = 5, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

plots_path <- file.path(script_dir, "panel_C_plot.rds")
saveRDS(panel_C, plots_path)
cat(sprintf("Saved: %s\n", plots_path))

cat("\nDone!\n")
