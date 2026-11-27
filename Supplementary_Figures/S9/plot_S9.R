#!/usr/bin/env Rscript

# Supplementary Figure S9: correlation-based credible set extension in the
# complex simulations. One row per scenario, 95% CS power on the left and FDR
# on the right, against the number N of top variants treated as causal:
#   A, B   extensive infinitesimal background ("Complex S2")
#   C, D   moderate infinitesimal background  ("Complex S1")
#   E, F   polygenic background              ("Complex")
# Each panel has six lines: SuSiE, SuSiE-ash and SuSiE-inf, each without the
# extension (method colour, solid) and with it (lighter tint, dashed).
#
# Input:  data/s9_finemapping_data.rds  (extract_S9_data.R)
# Output: S9.pdf, S9.png

library(ggplot2)
library(cowplot)

source("R/paths.R")
source("R/aesthetics.R")

s9_dir <- supp_dir(9)
agg    <- readRDS(file.path(s9_dir, "data", "s9_finemapping_data.rds"))

# The extended sets take the method colour mixed about half-way to white.
series_colors <- c(
  method_colors,
  "SuSiE w/ Extension"     = "#A6C8F0",
  "SuSiE-inf w/ Extension" = "#BCD9A0",
  "SuSiE-ash w/ Extension" = "#F29C9A"
)
# The series levels alternate between a method and its extension.
series_linetypes <- setNames(rep(c("solid", "dashed"), 3), levels(agg$series))

theme_topn <- theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text  = element_text(color = "black", size = 11),
    legend.position = "none",
    plot.margin = margin(6, 6, 6, 6)
  )

# ---- Panels ----

# Power panels run N = 3 -> 25 and FDR panels 25 -> 3, so that performance
# declines from left to right in every panel (as in Figure 1C-F).
make_panel <- function(sc, metric) {
  N_order <- c(3, 5, 10, 15, 20, 25)
  if (metric == "FDR") N_order <- rev(N_order)
  y_mean <- paste0(metric, "_mean")
  y_se   <- paste0(metric, "_se")

  df <- agg[agg$scenario == sc, ]
  df$N_position <- match(df$N, N_order)

  p <- ggplot(df, aes(x = N_position, y = .data[[y_mean]],
                      color = series, linetype = series, group = series)) +
    geom_line(linewidth = 1.0) +
    geom_errorbar(
      aes(ymin = .data[[y_mean]] - .data[[y_se]], ymax = .data[[y_mean]] + .data[[y_se]]),
      width = 0.25, linewidth = 0.5, linetype = "solid"
    ) +
    geom_point(size = 2.6) +
    scale_color_manual(values = series_colors) +
    scale_linetype_manual(values = series_linetypes) +
    scale_x_continuous(breaks = 1:6, labels = N_order) +
    labs(x = "Top N as Causal Variants", y = paste("95% CS", metric)) +
    theme_topn +
    coord_cartesian(ylim = if (metric == "Power") c(0, 0.9) else c(0, 0.35),
                    xlim = c(0.5, 6.5))

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted", color = "red", linewidth = 0.8)
  }
  p
}

build_row <- function(sc, labels) {
  plot_grid(make_panel(sc, "Power"), make_panel(sc, "FDR"), nrow = 1,
            labels = labels, label_size = 18, label_fontface = "bold")
}

# ---- Assemble ----

# The legend shows filled squares in two rows, one column per method: without
# the extension above, with it below.
legend <- get_legend(
  ggplot(agg[agg$scenario == "Complex", ], aes(x = N, y = Power_mean, fill = series)) +
    geom_col() +
    scale_fill_manual(values = series_colors) +
    guides(fill = guide_legend(nrow = 2)) +
    theme(legend.position = "bottom",
          legend.title    = element_blank(),
          legend.text     = element_text(size = 12),
          legend.key.size = unit(0.6, "cm"))
)

final_figure <- plot_grid(
  build_row("Complex S2", c("A", "B")),
  build_row("Complex S1", c("C", "D")),
  build_row("Complex",    c("E", "F")),
  legend,
  ncol = 1,
  rel_heights = c(1, 1, 1, 0.20)
)

save_figure(final_figure, s9_dir, "S9", width = 12, height = 15)
