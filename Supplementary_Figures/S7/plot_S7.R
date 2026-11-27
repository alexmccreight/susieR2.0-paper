#!/usr/bin/env Rscript

# Supplementary Figure S7: credible set construction strategies in the sparse
# simulations. Five variants: SuSiE and SuSiE+BLiP with a purity threshold of
# 0.5 or 0.8, and SuSiE with attainable coverage.
#   A, B   95% CS power and FDR against total PVE (K = 1-5 causal variants)
#   C      purity of the reported sets, pooled over K and replicates
#   D      size of the reported sets, pooled over K and replicates
#   E      cumulative distribution of the full posterior mass in the sets that
#          susieR and SparsePro report from the same fitted posteriors
#
# Input:  data/s7_metrics.rds             (extract_S7_data.R)
#         data/paired-posterior-mass.csv  (simulate_S7E.R)
# Output: S7.pdf, S7.png

library(ggplot2)
library(cowplot)
suppressPackageStartupMessages(library(dplyr))

source("R/paths.R")

s7_dir   <- supp_dir(7)
data_dir <- file.path(s7_dir, "data")

dat    <- readRDS(file.path(data_dir, "s7_metrics.rds"))
paired <- read.csv(file.path(data_dir, "paired-posterior-mass.csv"))

# Display labels, in plotting order.
method_labels <- c(
  "SuSiE(0.5)"       = "SuSiE (purity = 0.5)",
  "SuSiE+BLiP(0.5)"  = "SuSiE+BLiP (purity = 0.5)",
  "SuSiE(0.8)"       = "SuSiE (purity = 0.8)",
  "SuSiE+BLiP(0.8)"  = "SuSiE+BLiP (purity = 0.8)",
  "SuSiE+attainable" = "SuSiE+attainable"
)

relabel <- function(df) {
  df$Method <- factor(method_labels[as.character(df$Method)],
                      levels = unname(method_labels))
  df
}
rep_metrics <- relabel(dat$rep_metrics)
cs_purity   <- relabel(dat$cs_purity)
cs_size     <- relabel(dat$cs_size)

method_colors <- c(
  "SuSiE (purity = 0.5)"      = "#4A90E2",   # SuSiE blue (default purity)
  "SuSiE (purity = 0.8)"      = "#1E5A9B",   # darker blue (purity-filtered)
  "SuSiE+attainable"          = "#E67E22",   # orange (attainable)
  "SuSiE+BLiP (purity = 0.5)" = "#9B59B6",   # purple (BLiP)
  "SuSiE+BLiP (purity = 0.8)" = "#5B2C6F"    # darker purple (BLiP, purity-filtered)
)

base_theme <- theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(color = "black", size = 11),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

# ---- A, B: power and FDR ----

agg <- rep_metrics %>%
  group_by(Method, K, Total_PVE) %>%
  summarise(
    Power_mean = mean(Power, na.rm = TRUE),
    Power_se   = sd(Power,  na.rm = TRUE) / sqrt(sum(!is.na(Power))),
    FDR_mean   = mean(FDR,   na.rm = TRUE),
    FDR_se     = sd(FDR,    na.rm = TRUE) / sqrt(sum(!is.na(FDR))),
    .groups = "drop"
  )

create_bar_panel <- function(data, metric = "Power") {
  if (metric == "Power") {
    y_label  <- "95% CS Power"
    y_mean   <- "Power_mean"
    y_se     <- "Power_se"
    y_limits <- c(0.3, 1.0)
  } else {
    y_label  <- "95% CS FDR"
    y_mean   <- "FDR_mean"
    y_se     <- "FDR_se"
    y_limits <- c(0, 0.1)
  }

  # The lower end of an error bar stops at the axis limit.
  p <- ggplot(data, aes(x = factor(Total_PVE), y = .data[[y_mean]], fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
    geom_errorbar(
      aes(ymin = pmax(y_limits[1], .data[[y_mean]] - .data[[y_se]]),
          ymax = .data[[y_mean]] + .data[[y_se]]),
      position = position_dodge(width = 0.8),
      width = 0.15, linewidth = 0.5
    ) +
    scale_fill_manual(values = method_colors) +
    scale_x_discrete(expand = expansion(add = c(0.3, 0.3))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Total PVE", y = y_label) +
    base_theme +
    coord_cartesian(ylim = y_limits)

  for (K in 1:5) {
    p <- p + annotate("text", x = K + 0.05, y = y_limits[2] * 0.99,
                      label = paste0("K = ", K), size = 4.5, fontface = "bold")
  }

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted",
                        color = "red", linewidth = 0.8)
  }

  p
}

panel_A <- create_bar_panel(agg, metric = "Power")
panel_B <- create_bar_panel(agg, metric = "FDR")

# ---- C: purity ----

panel_C <- ggplot(cs_purity, aes(x = Method, y = purity, fill = Method)) +
  geom_violin(scale = "width", alpha = 0.85,
              color = "grey20", linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "grey20",
               fill = "white", alpha = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1.02),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Purity") +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ---- D: credible set size ----

panel_D <- ggplot(cs_size, aes(x = Method, y = size, fill = Method)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4,
               color = "grey20", linewidth = 0.4) +
  scale_fill_manual(values = method_colors) +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "CS Size (log scale)") +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ---- E: full posterior mass in the reported sets ----
# Both reporting rules receive the same fitted single-effect posteriors. The
# dashed line marks 95%; the key counts the sets that fall below it.

reporting_methods <- c("susieR", "SparsePro (2023)")
reporting_colors <- c("susieR" = "#E67E22", "SparsePro (2023)" = "#16758A")
mass_columns <- c("mass_susier", "mass_sparsepro")
curves <- do.call(rbind, lapply(seq_along(reporting_methods), function(i) {
  mass <- sort(paired[[mass_columns[i]]] * 100)
  data.frame(Method = reporting_methods[i], mass = c(90, mass, 100),
             cumulative = c(0, seq_along(mass) / length(mass) * 100, 100))
}))
# The tolerance keeps a mass of 0.95 up to rounding error out of the count.
counts <- vapply(paired[mass_columns], function(x) sum(x < 0.95 - 1e-12), integer(1))

panel_E <- ggplot(curves, aes(mass, cumulative, color = Method)) +
  annotate("rect", xmin = 90, xmax = 95, ymin = -Inf, ymax = Inf,
           fill = "#F4F5F6") +
  geom_vline(xintercept = 95, linetype = "dashed", color = "#555555", linewidth = 0.6) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = reporting_colors) +
  scale_x_continuous(breaks = c(90, 92.5, 95, 97.5, 100),
                     labels = c("90", "92.5", "95", "97.5", "100"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(90, 100), ylim = c(0, 102)) +
  labs(x = "Full posterior mass in credible set (%)", y = "Cumulative sets (%)") +
  base_theme + theme(plot.margin = margin(18, 8, 8, 18))

reporting_key <- ggdraw() +
  draw_text("Credible sets below 95%", x = 0.02, y = 0.83,
            hjust = 0, size = 14, fontface = "bold")
for (i in seq_along(reporting_methods)) {
  y <- c(0.55, 0.29)[i]
  color <- reporting_colors[i]
  reporting_key <- reporting_key +
    draw_line(x = c(0.02, 0.09), y = c(y, y), color = color, linewidth = 1.1) +
    draw_text(reporting_methods[i], x = 0.12, y = y, hjust = 0, size = 12) +
    draw_text(sprintf("%d/%d (%.1f%%)", counts[i], nrow(paired),
                      100 * counts[i] / nrow(paired)),
              x = 0.98, y = y, hjust = 1, size = 12, color = color, fontface = "bold")
}

# ---- Assemble ----

# Legend for panels A-D, placed under panels C and D.
legend <- get_legend(
  panel_A + theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 11),
    legend.key.width = unit(18, "pt")
  )
)

top_row <- plot_grid(panel_A, panel_B, nrow = 1,
                     labels = c("A", "B"), label_size = 16, label_fontface = "bold")
bottom_row <- plot_grid(panel_C, panel_D, nrow = 1,
                        labels = c("C", "D"), label_size = 16, label_fontface = "bold")
coverage_row <- plot_grid(panel_E, reporting_key, nrow = 1, rel_widths = c(0.63, 0.37),
                          labels = c("E", ""), label_size = 16, label_fontface = "bold")

# Panels A-D and their legend fill a 14 x 10 in block; Panel E takes a further
# 3.25 in underneath (13.25 in = 954 pt, a whole number of points, so the PDF
# device does not shift the page content).
AD_REL_H <- c(1, 1, 0.12)
AD_IN    <- 10
E_IN     <- 3.25
E_REL_H  <- sum(AD_REL_H) * E_IN / AD_IN

final_figure <- plot_grid(top_row, bottom_row, legend, coverage_row,
                          ncol = 1, rel_heights = c(AD_REL_H, E_REL_H))

save_figure(final_figure, s7_dir, "S7", width = 14, height = AD_IN + E_IN)
