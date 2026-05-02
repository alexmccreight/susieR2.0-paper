#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S5: Credible Set Size & Purity
# =============================================================================
# Panels A-D (2x2): per-scenario barplot of CS counts by size bin
#   A = Sparse, B = Complex, C = Complex S1, D = Complex S2
# Panel E: per-scenario median purity dotplot (wide row beneath A-D)
#
# Input:  data/s5_cs_data.rds  (built by extract_S5_data.R)
# Output: S5.pdf, S5.png
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

# =============================================================================
# Paths
# =============================================================================

s5_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S5"
data_dir <- file.path(s5_dir, "data")
data_path <- file.path(data_dir, "s5_cs_data.rds")

if (!file.exists(data_path)) {
  stop("Run extract_S5_data.R first to produce ", data_path)
}

# =============================================================================
# Load
# =============================================================================

cat("Loading data...\n")
dat <- readRDS(data_path)
all_cs          <- dat$cs_data
method_levels   <- dat$meta$method_levels      # SuSiE, SuSiE-ash, SuSiE-inf
scenario_levels <- dat$meta$scenario_levels    # Sparse, Complex, Complex S1, Complex S2
size_bin_levels <- dat$meta$size_bin_levels    # 1, 2, 3-5, 6-10, >10

# Defensive: re-apply factor levels in case the saved rds was modified.
all_cs$method   <- factor(all_cs$method,   levels = method_levels)
all_cs$scenario <- factor(all_cs$scenario, levels = scenario_levels)
all_cs$cs_size_bin <- factor(all_cs$cs_size_bin, levels = size_bin_levels)

cat(sprintf("Total CSs: %d\n", nrow(all_cs)))
print(table(all_cs$scenario, all_cs$method))

# =============================================================================
# Aesthetics
# =============================================================================

method_colors <- c(
  "SuSiE"     = "#4A90E2",
  "SuSiE-ash" = "#E53935",
  "SuSiE-inf" = "#7CB342"
)

# Display labels for panel-E x-axis (descriptive simulation names)
scenario_labels <- c(
  "Sparse"     = "Sparse",
  "Complex"    = "Oligogenic Effects on a\nPolygenic Background",
  "Complex S1" = "Oligogenic Effects on a Moderate\nInfinitesimal Background",
  "Complex S2" = "Oligogenic Effects on an Extensive\nInfinitesimal Background"
)

theme_s5 <- theme_classic(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold", size = 18),
    axis.text        = element_text(color = "black", size = 15),
    axis.text.x      = element_text(face = "bold"),
    legend.position  = "none",
    plot.margin      = margin(10, 10, 10, 10)
  )

# =============================================================================
# Helper: build CS-size barplot for one scenario
# =============================================================================

make_cs_barplot <- function(cs_df) {
  # Total counts per method
  total_df <- as.data.frame(table(cs_df$method), stringsAsFactors = FALSE)
  colnames(total_df) <- c("method", "count")
  total_df$category <- "Total"

  # Size-bin counts per method
  size_dist <- as.data.frame(
    table(cs_df$method, cs_df$cs_size_bin),
    stringsAsFactors = FALSE
  )
  colnames(size_dist) <- c("method", "category", "count")

  grouped_df <- rbind(total_df[, c("method", "category", "count")],
                      size_dist[, c("method", "category", "count")])
  grouped_df$method   <- factor(grouped_df$method, levels = method_levels)
  grouped_df$category <- factor(grouped_df$category,
                                levels = c("Total", size_bin_levels))

  ggplot(grouped_df, aes(x = category, y = count, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey60",
               linewidth = 0.4) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "CS Size", y = "Number of CSs") +
    theme_s5
}

# =============================================================================
# Panels A-D: per-scenario CS-size barplots
# =============================================================================

cat("Creating Panels A-D...\n")
panel_A <- make_cs_barplot(all_cs[all_cs$scenario == "Sparse",     ])
panel_B <- make_cs_barplot(all_cs[all_cs$scenario == "Complex",    ])
panel_C <- make_cs_barplot(all_cs[all_cs$scenario == "Complex S1", ])
panel_D <- make_cs_barplot(all_cs[all_cs$scenario == "Complex S2", ])

# =============================================================================
# Panel E: per-scenario median purity dotplot
# =============================================================================

cat("Creating Panel E (median purity)...\n")

purity_summary <- all_cs %>%
  group_by(scenario, method) %>%
  summarise(median_purity = median(min_abs_corr, na.rm = TRUE),
            .groups = "drop")

panel_E <- ggplot(purity_summary,
                  aes(x = scenario, y = median_purity, color = method)) +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = method_colors) +
  scale_x_discrete(labels = scenario_labels) +
  scale_y_continuous(limits = c(0.9, 1.0),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Scenario", y = "Median Purity") +
  theme_s5 +
  theme(
    axis.text.x = element_text(face = "bold", size = 11, lineheight = 0.9),
    legend.position = "none"
  )

# =============================================================================
# Shared legend (extracted once from a dummy plot)
# =============================================================================

legend_plot <- ggplot(purity_summary,
                      aes(x = scenario, y = median_purity, fill = method)) +
  geom_col() +
  scale_fill_manual(values = method_colors) +
  theme(legend.position      = "bottom",
        legend.justification = "center",
        legend.title         = element_blank(),
        legend.text          = element_text(size = 14, face = "bold"))

legend_grob <- cowplot::get_legend(legend_plot)

# =============================================================================
# Compose: 2x2 (A-D) + wide bottom (E) + legend
# =============================================================================

cat("Combining panels...\n")

top_row <- plot_grid(panel_A, panel_B, nrow = 1,
                     labels = c("A", "B"),
                     label_size = 20, label_fontface = "bold")

mid_row <- plot_grid(panel_C, panel_D, nrow = 1,
                     labels = c("C", "D"),
                     label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(
  top_row, mid_row,
  panel_E, legend_grob,
  ncol        = 1,
  rel_heights = c(1, 1, 0.9, 0.06),
  labels      = c("", "", "E", ""),
  label_size  = 20,
  label_fontface = "bold"
)

# =============================================================================
# Save
# =============================================================================

output_pdf <- file.path(s5_dir, "S5.pdf")
ggsave(output_pdf, final_figure,
       width = 14, height = 16, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(s5_dir, "S5.png")
ggsave(output_png, final_figure,
       width = 14, height = 16, units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

cat("\nDone!\n")
