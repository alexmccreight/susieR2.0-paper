#!/usr/bin/env Rscript

# Supplementary Figure S3: credible set size and purity in the simulations.
#   A-D  number of 95% credible sets, in total and by size, for the sparse (A),
#        oligogenic on a polygenic background (B), moderate infinitesimal (C)
#        and extensive infinitesimal (D) scenarios
#   E    median purity of the credible sets in each scenario
#
# Input:  data/s3_cs_data.rds  (extract_S3_data.R)
# Output: S3.pdf, S3.png

library(ggplot2)
library(cowplot)
suppressPackageStartupMessages(library(dplyr))

source("R/paths.R")
source("R/aesthetics.R")

s3_dir <- supp_dir(3)

dat <- readRDS(file.path(s3_dir, "data", "s3_cs_data.rds"))
all_cs          <- dat$cs_data
method_levels   <- dat$meta$method_levels      # SuSiE, SuSiE-ash, SuSiE-inf
scenario_levels <- dat$meta$scenario_levels    # Sparse, Complex, Complex S1, Complex S2
size_bin_levels <- dat$meta$size_bin_levels    # 1, 2, 3-5, 6-10, >10

# Scenario names shown on the x-axis of panel E.
scenario_labels <- c(
  "Sparse"     = "Sparse",
  "Complex"    = "Oligogenic Effects on a\nPolygenic Background",
  "Complex S1" = "Oligogenic Effects on a Moderate\nInfinitesimal Background",
  "Complex S2" = "Oligogenic Effects on an Extensive\nInfinitesimal Background"
)

theme_base <- theme_classic(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold", size = 18),
    axis.text        = element_text(color = "black", size = 15),
    axis.text.x      = element_text(face = "bold"),
    legend.position  = "none",
    plot.margin      = margin(10, 10, 10, 10)
  )

# ---- A-D: credible set counts by size ----

make_cs_barplot <- function(cs_df) {
  total_df <- as.data.frame(table(cs_df$method), stringsAsFactors = FALSE)
  colnames(total_df) <- c("method", "count")
  total_df$category <- "Total"

  size_dist <- as.data.frame(
    table(cs_df$method, cs_df$cs_size_bin),
    stringsAsFactors = FALSE
  )
  colnames(size_dist) <- c("method", "category", "count")

  grouped_df <- rbind(total_df[, c("method", "category", "count")], size_dist)
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
    theme_base
}

size_panels <- lapply(scenario_levels, function(s) {
  make_cs_barplot(all_cs[all_cs$scenario == s, ])
})

# ---- E: median purity ----

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
  theme_base +
  theme(axis.text.x = element_text(face = "bold", size = 11, lineheight = 0.9))

# ---- Assemble ----

legend <- get_legend(
  ggplot(purity_summary, aes(x = scenario, y = median_purity, fill = method)) +
    geom_col() +
    scale_fill_manual(values = method_colors) +
    theme(legend.position      = "bottom",
          legend.justification = "center",
          legend.title         = element_blank(),
          legend.text          = element_text(size = 14, face = "bold"))
)

size_grid <- plot_grid(plotlist = size_panels, nrow = 2,
                       labels = c("A", "B", "C", "D"),
                       label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(
  size_grid, panel_E, legend,
  ncol        = 1,
  rel_heights = c(2, 0.9, 0.06),
  labels      = c("", "E", ""),
  label_size  = 20,
  label_fontface = "bold"
)

save_figure(final_figure, s3_dir, "S3", width = 14, height = 16)
