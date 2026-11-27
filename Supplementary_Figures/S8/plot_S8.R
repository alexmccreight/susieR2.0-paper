#!/usr/bin/env Rscript

# Supplementary Figure S8: correlation-based credible set extension in the
# sparse simulations.
#   A, B   95% CS power and FDR against total PVE (K = 1-5 causal variants) for
#          SuSiE, SuSiE-ash and SuSiE-inf, each without the extension (method
#          colour) and with it (lighter tint)
#
# Input:  data/s8_finemapping_data.rds  (extract_S8_data.R)
# Output: S8.pdf, S8.png

library(ggplot2)
library(cowplot)

source("R/paths.R")
source("R/aesthetics.R")

s8_dir <- supp_dir(8)
dat    <- readRDS(file.path(s8_dir, "data", "s8_finemapping_data.rds"))

# The extended sets take the method colour mixed about half-way to white.
series_colors <- c(
  method_colors,
  "SuSiE w/ Extension"     = "#A6C8F0",
  "SuSiE-inf w/ Extension" = "#BCD9A0",
  "SuSiE-ash w/ Extension" = "#F29C9A"
)

theme_bar <- theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold", size = 18),
    axis.text    = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin  = margin(10, 10, 10, 10)
  )

# ---- A, B: power and FDR ----

create_bar_panel <- function(data, metric = "Power") {
  if (metric == "Power") {
    y_label  <- "95% CS Power"
    y_mean   <- "Power_mean"
    y_se     <- "Power_se"
    y_limits <- c(0.3, 0.85)
  } else {
    y_label  <- "95% CS FDR"
    y_mean   <- "FDR_mean"
    y_se     <- "FDR_se"
    y_limits <- c(0, 0.1)
  }

  p <- ggplot(data, aes(x = factor(Total_PVE), y = .data[[y_mean]], fill = series)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.8) +
    geom_errorbar(
      aes(ymin = .data[[y_mean]] - .data[[y_se]],
          ymax = .data[[y_mean]] + .data[[y_se]]),
      position = position_dodge(width = 0.85),
      width = 0.25, linewidth = 0.4
    ) +
    scale_fill_manual(values = series_colors) +
    scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Total PVE", y = y_label) +
    theme_bar +
    coord_cartesian(ylim = y_limits)

  for (K in 1:5) {
    p <- p + annotate("text", x = K, y = y_limits[2] * 0.99,
                      label = paste0("K = ", K), size = 4.5, fontface = "bold")
  }

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted",
                        color = "red", linewidth = 0.8)
  }

  p
}

panel_A <- create_bar_panel(dat, "Power")
panel_B <- create_bar_panel(dat, "FDR")

# ---- Assemble ----

# Two legend rows, one column per method: without the extension above, with it
# below.
legend <- get_legend(
  ggplot(dat, aes(x = factor(Total_PVE), y = Power_mean, fill = series)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = series_colors) +
    guides(fill = guide_legend(nrow = 2)) +
    theme(legend.position = "bottom",
          legend.title    = element_blank(),
          legend.text     = element_text(size = 13))
)

row <- plot_grid(panel_A, panel_B, nrow = 1,
                 labels = c("A", "B"),
                 label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(row, legend, ncol = 1, rel_heights = c(1, 0.16))

save_figure(final_figure, s8_dir, "S8", width = 16, height = 7)
