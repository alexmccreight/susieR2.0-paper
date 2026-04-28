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
make_prop_bar <- function(dat) {
  ggplot(dat, aes(x = Method, y = Proportion, fill = Method)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Proportion of Replicates") +
    theme_s3 +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
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

panel_E <- make_prop_bar(d$sparse_prop)

# =============================================================================
# Panels F, G, H: Complex 1, 2, 3 Proportion of Replicates Won
# =============================================================================
cat("Creating Panels F-H: Complex proportion of replicates won...\n")

panel_F <- make_prop_bar(d$complex1_prop)
panel_G <- make_prop_bar(d$complex2_prop)
panel_H <- make_prop_bar(d$complex3_prop)

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
# Column titles
# =============================================================================
make_header <- function(label, fill = "grey90", border = "grey40") {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = fill, color = border, linewidth = 0.8) +
    annotate("text", x = 0.5, y = 0.5, label = label,
             fontface = "bold", size = 4.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void() +
    theme(plot.margin = margin(2, 4, 2, 4))
}

title_sparse   <- make_header("Sparse")
title_complex1 <- make_header("Oligogenic Effects on a\nPolygenic Background")
title_complex2 <- make_header("Oligogenic Effects on a Moderate\nInfinitesimal Background")
title_complex3 <- make_header("Oligogenic Effects on an Extensive\nInfinitesimal Background")

title_row <- plot_grid(title_sparse, title_complex1, title_complex2, title_complex3,
                       nrow = 1)

# =============================================================================
# Combine: title row + 2x4 grid + legend
# =============================================================================
cat("Combining panels...\n")

row1 <- plot_grid(panel_A, panel_B, panel_C, panel_D, nrow = 1,
                  labels = c("A", "B", "C", "D"),
                  label_size = 20, label_fontface = "bold")

row2 <- plot_grid(panel_E, panel_F, panel_G, panel_H, nrow = 1,
                  labels = c("E", "F", "G", "H"),
                  label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(
  title_row, row1, row2, legend_grob,
  ncol = 1,
  rel_heights = c(0.14, 1, 1, 0.06)
)

# Save
output_pdf <- file.path(script_dir, "S3.pdf")
ggsave(output_pdf, final_figure,
       width = 16, height = 10, units = "in", bg = "white")

output_png <- file.path(script_dir, "S3.png")
ggsave(output_png, final_figure,
       width = 16, height = 10, units = "in", dpi = 300, bg = "white")

cat(sprintf("\nFigure saved to:\n  - %s\n  - %s\n", output_pdf, output_png))
