#!/usr/bin/env Rscript

# Figure 2, panel A: number of 95% credible sets for SuSiE and SuSiE-inf, in
# total and by credible-set size.
#
# Usage:  plot_panel_A.R [--include-ash]
#         --include-ash draws the three-method panel of Supplementary Figure
#         S10 from Supplementary_Figures/S10/data/panel_A_data.rds and writes
#         only panel_A_plots.rds, in that folder.
# Input:  panel_A_data.rds  (prepare_panel_A_data.R)
# Output: panel_A.pdf, panel_A.png,
#         panel_A_plots.rds  (list with one element, A; read by create_figure2.R)

library(ggplot2)

source("R/paths.R")
source("R/aesthetics.R")

args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--include-ash")) {
  stop("Unknown argument(s): ", paste(args[!args %in% "--include-ash"], collapse = ", "))
}
include_ash <- "--include-ash" %in% args

output_dir <- if (include_ash) {
  file.path(supp_dir(10), "data")
} else {
  file.path(fig_dir(2), "figure_2_panel_A")
}
cs_data <- readRDS(file.path(output_dir, "panel_A_data.rds"))

# Counts per method: all credible sets, then each size bin.
total_df <- as.data.frame(table(cs_data$method), stringsAsFactors = FALSE)
colnames(total_df) <- c("method", "count")
total_df$category <- "Total"

size_dist <- as.data.frame(
  table(cs_data$method, cs_data$cs_size_bin),
  stringsAsFactors = FALSE
)
colnames(size_dist) <- c("method", "category", "count")

grouped_df <- rbind(total_df[, c("method", "category", "count")], size_dist)
grouped_df$method <- factor(grouped_df$method, levels = levels(cs_data$method))
grouped_df$category <- factor(grouped_df$category,
                              levels = c("Total", levels(cs_data$cs_size_bin)))

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

if (!include_ash) {
  save_figure(plot_A, output_dir, "panel_A", width = 8, height = 5)
}

plots_path <- file.path(output_dir, "panel_A_plots.rds")
saveRDS(list(A = plot_A), plots_path)
cat("Saved:", plots_path, "\n")
