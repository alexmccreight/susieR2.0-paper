#!/usr/bin/env Rscript

# Supplementary Figure S4: gene expression prediction in the simulations.
#   A-D  cross-validated R2 of each replicate
#   E-H  proportion of replicates in which each method had the highest R2
# One column per scenario: sparse (A, E; R2 by total PVE, proportions pooled
# over K = 1-5), then oligogenic effects on a polygenic (B, F), a moderate
# infinitesimal (C, G) and an extensive infinitesimal (D, H) background.
#
# Input:  data/prediction_data.rds  (extract_prediction_data.R)
# Output: S4.pdf, S4.png

library(ggplot2)
library(cowplot)

source("R/paths.R")
source("R/aesthetics.R")

s4_dir <- supp_dir(4)

d <- readRDS(file.path(s4_dir, "data", "prediction_data.rds"))
for (nm in names(d)) {
  d[[nm]]$Method <- factor(d[[nm]]$Method, levels = ALL_METHODS)
}

theme_base <- theme_classic(base_size = 14) +
  theme(
    axis.title      = element_text(face = "bold", size = 18),
    axis.text       = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin     = margin(10, 10, 10, 10)
  )

r2_label <- expression(bold("Pearson's")~bold(R^2))

# ---- A-D: R2 boxplots ----

panel_A <- ggplot(d$sparse_data,
                  aes(x = factor(Total_PVE), y = R2, fill = Method)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Total PVE", y = r2_label) +
  theme_base

# The x-axis text and title are drawn invisibly rather than removed, so that
# these panels keep the height of panel A.
make_complex_r2 <- function(dat) {
  ggplot(dat, aes(x = Method, y = R2, fill = Method)) +
    geom_boxplot(outlier.size = 0.8) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = " ", y = r2_label) +
    theme_base +
    theme(axis.text.x = element_text(color = "transparent"),
          axis.ticks.x = element_line(color = "transparent"))
}

r2_panels <- c(
  list(panel_A),
  lapply(c("complex1", "complex2", "complex3"), function(s) {
    make_complex_r2(d[[paste0(s, "_data")]])
  })
)

# ---- E-H: proportion of replicates in which each method was best ----

# The scenario name is the x-axis title. "Sparse" has a blank second line so
# that it takes the same height as the two-line names.
scenario_labels <- c(
  sparse   = "Sparse\n ",
  complex1 = "Oligogenic Effects on a\nPolygenic Background",
  complex2 = "Oligogenic Effects on a Moderate\nInfinitesimal Background",
  complex3 = "Oligogenic Effects on an Extensive\nInfinitesimal Background"
)

make_prop_bar <- function(dat, scenario_label) {
  ggplot(dat, aes(x = Method, y = Proportion, fill = Method)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = scenario_label, y = "Proportion of Replicates") +
    theme_base +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(face = "bold", size = 13, lineheight = 0.9,
                                  margin = margin(t = 6))
    )
}

prop_panels <- lapply(names(scenario_labels), function(s) {
  make_prop_bar(d[[paste0(s, "_prop")]], scenario_labels[[s]])
})

# ---- Assemble ----

legend <- get_legend(
  ggplot(d$sparse_prop, aes(x = Method, y = Proportion, fill = Method)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = method_colors) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14))
)

row1 <- plot_grid(plotlist = r2_panels, nrow = 1,
                  labels = c("A", "B", "C", "D"),
                  label_size = 20, label_fontface = "bold")

row2 <- plot_grid(plotlist = prop_panels, nrow = 1,
                  labels = c("E", "F", "G", "H"),
                  label_size = 20, label_fontface = "bold")

# Row 2 is taller to leave room for the scenario names.
final_figure <- plot_grid(
  row1, row2, legend,
  ncol = 1,
  rel_heights = c(1, 1.12, 0.06)
)

save_figure(final_figure, s4_dir, "S4", width = 16, height = 10)
