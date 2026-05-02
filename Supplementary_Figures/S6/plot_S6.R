#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S6: Per-method runtime by scenario
# =============================================================================
# Single-panel bar plot of mean per-fit wall time (seconds) for SuSiE,
# SuSiE-ash, SuSiE-inf, dodged within each scenario on the x-axis. Error bars
# = SE across replicates.
#
# Input:  data/s6_runtime_data.rds  (built by extract_S6_data.R)
# Output: S6.pdf, S6.png
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
})

# =============================================================================
# Paths
# =============================================================================

s6_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S6"
data_dir <- file.path(s6_dir, "data")
data_path <- file.path(data_dir, "s6_runtime_data.rds")

if (!file.exists(data_path)) {
  stop("Run extract_S6_data.R first to produce ", data_path)
}

# =============================================================================
# Load
# =============================================================================

cat("Loading data...\n")
dat <- readRDS(data_path)
runtime_data    <- dat$runtime_data
method_levels   <- dat$meta$method_levels      # SuSiE, SuSiE-ash, SuSiE-inf
scenario_levels <- dat$meta$scenario_levels    # Sparse, Complex, Complex S1, Complex S2

# Defensive: re-apply factor levels in case the saved rds was modified.
runtime_data$method   <- factor(runtime_data$method,   levels = method_levels)
runtime_data$scenario <- factor(runtime_data$scenario, levels = scenario_levels)

cat(sprintf("Total rows: %d\n", nrow(runtime_data)))

# =============================================================================
# Aggregate: mean +/- SE per (scenario, method)
# =============================================================================

agg <- runtime_data %>%
  group_by(scenario, method) %>%
  summarise(
    mean_time = mean(elapsed_time, na.rm = TRUE),
    se_time   = sd(elapsed_time,   na.rm = TRUE) / sqrt(sum(!is.na(elapsed_time))),
    n         = sum(!is.na(elapsed_time)),
    .groups   = "drop"
  )

cat("\n=== Aggregated mean elapsed_time (s) ===\n")
print(as.data.frame(agg))

# =============================================================================
# Aesthetics (match S1 / S5 conventions)
# =============================================================================

method_colors <- c(
  "SuSiE"     = "#4A90E2",
  "SuSiE-ash" = "#E53935",
  "SuSiE-inf" = "#7CB342"
)

# Display labels for x-axis (descriptive simulation names; matches S5)
scenario_labels <- c(
  "Sparse"     = "Sparse",
  "Complex"    = "Oligogenic Effects on a\nPolygenic Background",
  "Complex S1" = "Oligogenic Effects on a Moderate\nInfinitesimal Background",
  "Complex S2" = "Oligogenic Effects on an Extensive\nInfinitesimal Background"
)

theme_s6 <- theme_classic(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold", size = 18),
    axis.text        = element_text(color = "black", size = 15),
    axis.text.x      = element_text(face = "bold", size = 11, lineheight = 0.9),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(face = "bold", size = 14),
    plot.margin      = margin(10, 10, 10, 10)
  )

# =============================================================================
# Bar plot: scenario on x, method as fill, dodged with SE error bars
# =============================================================================

panel <- ggplot(agg, aes(x = scenario, y = mean_time, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_time - se_time),
        ymax = mean_time + se_time),
    position = position_dodge(width = 0.8),
    width    = 0.18,
    linewidth = 0.5
  ) +
  scale_fill_manual(values = method_colors) +
  scale_x_discrete(labels = scenario_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Scenario", y = "Mean Runtime (s)") +
  theme_s6

# =============================================================================
# Save
# =============================================================================

output_pdf <- file.path(s6_dir, "S6.pdf")
ggsave(output_pdf, panel,
       width = 14, height = 6, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(s6_dir, "S6.png")
ggsave(output_png, panel,
       width = 14, height = 6, units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

cat("\nDone!\n")
