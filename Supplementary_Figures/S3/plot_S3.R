#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S3: TWAS Prediction Figure (2 rows x 4 columns)
# =============================================================================
# Row 1: R² boxplots     (Sparse | Complex 1 | Complex 2 | Complex 3)
# Row 2: Prop replicates  (Sparse | Complex 1 | Complex 2 | Complex 3)
# =============================================================================

library(ggplot2)
library(cowplot)
library(dplyr)

script_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S3"
data_dir   <- file.path(script_dir, "data")

cat("Creating supplementary figure S3: TWAS prediction...\n")

# =============================================================================
# Load extracted data
# =============================================================================
d <- readRDS(file.path(data_dir, "prediction_data.rds"))

method_levels <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")
method_colors <- c(
  "SuSiE"     = "#4A90E2",
  "SuSiE-inf" = "#7CB342",
  "SuSiE-ash" = "#E53935"
)

theme_s3 <- theme_classic(base_size = 14) +
  theme(
    axis.title      = element_text(face = "bold", size = 18),
    axis.text       = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin     = margin(10, 10, 10, 10)
  )

# Factor all Method columns
for (nm in c("sparse_data", "sparse_prop",
             "complex1_data", "complex1_prop",
             "complex2_data", "complex2_prop",
             "complex3_data", "complex3_prop")) {
  d[[nm]]$Method <- factor(d[[nm]]$Method, levels = method_levels)
}

# =============================================================================
# Helper: Complex R² boxplot panel
# =============================================================================
make_complex_r2 <- function(dat) {
  ggplot(dat, aes(x = Method, y = R2, fill = Method)) +
    geom_boxplot(outlier.size = 0.8) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = " ", y = expression(bold("Pearson's")~bold(R^2))) +
    theme_s3 +
    theme(axis.text.x = element_text(color = "transparent"),
          axis.ticks.x = element_line(color = "transparent"))
}

# =============================================================================
# Helper: Proportion of replicates won bar plot
# =============================================================================
# `scenario_label` is shown as the panel's x-axis title (descriptive scenario
# name, replacing the boxed column titles previously rendered above row 1).
make_prop_bar <- function(dat, scenario_label) {
  ggplot(dat, aes(x = Method, y = Proportion, fill = Method)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = scenario_label, y = "Proportion of Replicates") +
    theme_s3 +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(face = "bold", size = 13, lineheight = 0.9,
                                  margin = margin(t = 6))
    )
}

# =============================================================================
# Panel A: Sparse R² boxplots
# =============================================================================
cat("Creating Panel A: Sparse R² boxplots...\n")

d$sparse_data$PVE_label <- factor(d$sparse_data$Total_PVE,
                                   levels = c(0.03, 0.06, 0.09, 0.12, 0.15))

panel_A <- ggplot(d$sparse_data,
                  aes(x = PVE_label, y = R2, fill = Method)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Total PVE", y = expression(bold("Pearson's")~bold(R^2))) +
  theme_s3

# =============================================================================
# Panels B, C, D: Complex 1, 2, 3 R² boxplots
# =============================================================================
cat("Creating Panels B-D: Complex R² boxplots...\n")

panel_B <- make_complex_r2(d$complex1_data)
panel_C <- make_complex_r2(d$complex2_data)
panel_D <- make_complex_r2(d$complex3_data)

# =============================================================================
# Panel E: Sparse Proportion of Replicates Won (pooled)
# =============================================================================
cat("Creating Panel E: Sparse proportion of replicates won...\n")

panel_E <- make_prop_bar(d$sparse_prop, "Sparse\n ")  # trailing blank line keeps height
                                                      # equal to the 2-line complex titles

# =============================================================================
# Panels F, G, H: Complex 1, 2, 3 Proportion of Replicates Won
# =============================================================================
cat("Creating Panels F-H: Complex proportion of replicates won...\n")

panel_F <- make_prop_bar(d$complex1_prop,
                          "Oligogenic Effects on a\nPolygenic Background")
panel_G <- make_prop_bar(d$complex2_prop,
                          "Oligogenic Effects on a Moderate\nInfinitesimal Background")
panel_H <- make_prop_bar(d$complex3_prop,
                          "Oligogenic Effects on an Extensive\nInfinitesimal Background")

# =============================================================================
# Shared legend
# =============================================================================
legend_plot <- ggplot(d$sparse_prop,
                      aes(x = Method, y = Proportion, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = method_colors) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))

legend_grob <- cowplot::get_legend(legend_plot)

# =============================================================================
# Combine: 2x4 grid + legend (scenario labels carried by row-2 x-axis titles)
# =============================================================================
cat("Combining panels...\n")

row1 <- plot_grid(panel_A, panel_B, panel_C, panel_D, nrow = 1,
                  labels = c("A", "B", "C", "D"),
                  label_size = 20, label_fontface = "bold")

row2 <- plot_grid(panel_E, panel_F, panel_G, panel_H, nrow = 1,
                  labels = c("E", "F", "G", "H"),
                  label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(
  row1, row2, legend_grob,
  ncol = 1,
  rel_heights = c(1, 1.12, 0.06)
)

# Save
output_pdf <- file.path(script_dir, "S3.pdf")
ggsave(output_pdf, final_figure,
       width = 16, height = 10, units = "in", bg = "white")

output_png <- file.path(script_dir, "S3.png")
ggsave(output_png, final_figure,
       width = 16, height = 10, units = "in", dpi = 300, bg = "white")

cat(sprintf("\nFigure saved to:\n  - %s\n  - %s\n", output_pdf, output_png))
