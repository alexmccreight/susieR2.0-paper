#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel A — Plotting
# =============================================================================
# Grouped bar chart: total CS counts + CS size distribution per method
#
# Input:  panel_A_data.rds (from prepare_panel_A_data.R)
# Output: panel_A.pdf, panel_A.png, panel_A_plots.rds
# =============================================================================

library(ggplot2)

# =============================================================================
# Paths
# =============================================================================

script_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2/figure_2_panel_A"
data_path  <- file.path(script_dir, "panel_A_data.rds")

# =============================================================================
# Shared aesthetics
# =============================================================================

method_colors <- c(
  "SuSiE"     = "#4A90E2",
  "SuSiE-ash" = "#E53935",
  "SuSiE-inf" = "#7CB342"
)

size_bin_levels <- c("1", "2", "3-5", "6-10", ">10")

# =============================================================================
# Load data
# =============================================================================

cat("Loading panel A data...\n")
cs_data <- readRDS(data_path)

# =============================================================================
# Build grouped bar data: Total + size bins
# =============================================================================

cat("Creating Panel A (grouped bar chart)...\n")

# Total counts per method
total_df <- as.data.frame(table(cs_data$method), stringsAsFactors = FALSE)
colnames(total_df) <- c("method", "count")
total_df$category <- "Total"

# Size-bin counts per method
size_dist <- as.data.frame(
  table(cs_data$method, cs_data$cs_size_bin),
  stringsAsFactors = FALSE
)
colnames(size_dist) <- c("method", "category", "count")

# Combine
grouped_df <- rbind(total_df[, c("method", "category", "count")],
                    size_dist[, c("method", "category", "count")])
grouped_df$method <- factor(grouped_df$method,
                            levels = c("SuSiE", "SuSiE-ash", "SuSiE-inf"))
grouped_df$category <- factor(grouped_df$category,
                              levels = c("Total", size_bin_levels))

# =============================================================================
# Plot
# =============================================================================

plot_A <- ggplot(grouped_df, aes(x = category, y = count, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey60",
             linewidth = 0.4) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "CS size", y = "Number of credible sets") +
  theme_classic(base_size = 20) +
  theme(
    axis.title       = element_text(face = "bold", size = 20),
    axis.text        = element_text(color = "black", size = 18),
    axis.text.x      = element_text(face = "bold"),
    legend.position  = "none",
    plot.margin      = margin(5, 10, 5, 5)
  )

# =============================================================================
# Save
# =============================================================================

cat("Saving panel A...\n")

output_pdf <- file.path(script_dir, "panel_A.pdf")
ggsave(output_pdf, plot_A, width = 8, height = 5, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(script_dir, "panel_A.png")
ggsave(output_png, plot_A, width = 8, height = 5, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

# Save plot object for use by create_figure2.R
panel_A_plots <- list(A = plot_A)
plots_path <- file.path(script_dir, "panel_A_plots.rds")
saveRDS(panel_A_plots, plots_path)
cat(sprintf("Saved: %s\n", plots_path))

cat("\nDone!\n")
