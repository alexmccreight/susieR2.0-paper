#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel C — Plotting (Dot Plot)
# =============================================================================
# Creates a dot plot of mean cross-validated R² across 6 TWAS prediction methods,
# aggregated across all three brain tissues (DLPFC, AC, PCC).
#
# Input:  panel_C_data.rds (from prepare_panel_C_data.R)
# Output: panel_C.pdf, panel_C.png, panel_C_plot.rds
# =============================================================================

library(ggplot2)

# =============================================================================
# Paths
# =============================================================================

source("R/paths.R")
source("R/aesthetics.R")

script_dir <- file.path(fig_dir(2), "figure_2_panel_C")
data_path  <- file.path(script_dir, "panel_C_data.rds")

# =============================================================================
# Shared aesthetics
# =============================================================================

# SuSiE and SuSiE-inf take their colours from the shared palette in
# R/aesthetics.R so they cannot drift from panels A, B, D and E. The remaining
# prediction methods are specific to this panel.
method_colors <- c(
  method_colors[MAIN_METHODS],   # SuSiE = #4A90E2, SuSiE-inf = #7CB342
  prediction_colors              # Elastic Net, LASSO, BayesR, BayesL
)

# =============================================================================
# Load data
# =============================================================================

cat("Loading panel C data...\n")
d <- readRDS(data_path)

# =============================================================================
# Build dot plot
# =============================================================================

cat("Creating dot plot...\n")

# Compute mean R² and SE per method
mean_df <- aggregate(rsq ~ method, data = d, FUN = function(x) {
  x <- x[!is.na(x)]
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
})
mean_df <- data.frame(
  method   = mean_df$method,
  mean_rsq = mean_df$rsq[, "mean"],
  se       = mean_df$rsq[, "se"]
)

panel_C <- ggplot(mean_df, aes(x = method, y = mean_rsq, color = method)) +
  geom_point(size = 8) +
  geom_errorbar(aes(ymin = mean_rsq - se, ymax = mean_rsq + se),
                width = 0.2, linewidth = 1.2) +
  scale_color_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = NULL, y = expression(bold("Mean CV" ~ R^2))) +
  theme_classic(base_size = 20) +
  theme(
    axis.title.y  = element_text(face = "bold", size = 20),
    axis.text.x   = element_text(color = "black", size = 17, angle = 45,
                                 hjust = 1, face = "bold"),
    axis.text.y   = element_text(color = "black", size = 18),
    legend.position = "none",
    plot.margin   = margin(10, 10, 5, 5)
  )

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

# Save plot object for use by create_figure2.R
plot_path <- file.path(script_dir, "panel_C_plot.rds")
saveRDS(panel_C, plot_path)
cat(sprintf("Saved: %s\n", plot_path))

cat("\nDone!\n")
