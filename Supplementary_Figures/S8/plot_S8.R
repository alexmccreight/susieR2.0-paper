#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S8: Complex Power / FDR -- no extension vs LD extension (0.99)
# =============================================================================
# 3 rows (scenario) x 2 cols (Power, FDR). Row order top->bottom:
#   Complex S2, Complex S1, Complex.
# Each panel has 6 lines = 3 methods x {no extension, LD ext 0.99}, vs Top-N.
# Encoding: color = series (base color + lighter "w/ Extension" tint),
#           linetype = version (solid = no extension, dashed = w/ extension).
# Style otherwise matches Figure 1C-1F (Top-N axis, error bars, points,
# Power ascending / FDR descending, FDR 0.05 reference line).
#
# Input:  data/s8_finemapping_data.rds  (built by extract_S8_data.R)
# Output: S8.pdf, S8.png
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

source("R/paths.R")

s8_dir    <- supp_dir(8)
data_path <- file.path(s8_dir, "data", "s8_finemapping_data.rds")
if (!file.exists(data_path)) stop("Run extract_S8_data.R first.")

agg <- readRDS(data_path)

# =============================================================================
# Aesthetics
# =============================================================================

series_levels <- c("SuSiE", "SuSiE w/ Extension",
                   "SuSiE-ash", "SuSiE-ash w/ Extension",
                   "SuSiE-inf", "SuSiE-inf w/ Extension")
agg$series <- factor(agg$series, levels = series_levels)
agg$scenario <- factor(agg$scenario, levels = c("Complex S2", "Complex S1", "Complex"))

series_colors <- c(
  "SuSiE"                  = "#4A90E2",
  "SuSiE w/ Extension"     = "#A6C8F0",
  "SuSiE-inf"              = "#7CB342",
  "SuSiE-inf w/ Extension" = "#BCD9A0",
  "SuSiE-ash"              = "#E53935",
  "SuSiE-ash w/ Extension" = "#F29C9A"
)
series_linetypes <- c(
  "SuSiE"                  = "solid",
  "SuSiE w/ Extension"     = "dashed",
  "SuSiE-inf"              = "solid",
  "SuSiE-inf w/ Extension" = "dashed",
  "SuSiE-ash"              = "solid",
  "SuSiE-ash w/ Extension" = "dashed"
)

theme_panel <- theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text  = element_text(color = "black", size = 11),
    legend.position = "none",
    plot.margin = margin(6, 6, 6, 6)
  )

# =============================================================================
# Panel builder (one scenario, one metric) -- mirrors Fig 1C-1F
# =============================================================================

make_panel <- function(sc, metric) {
  df <- agg[agg$scenario == sc, ]
  if (metric == "Power") {
    df$N_position <- with(df, dplyr::case_when(
      N == 3 ~ 1, N == 5 ~ 2, N == 10 ~ 3, N == 15 ~ 4, N == 20 ~ 5, N == 25 ~ 6))
    y_mean <- "Power_mean"; y_se <- "Power_se"
    y_lab  <- "95% CS Power"; y_lim <- c(0, 0.9)
    x_labs <- c("3", "5", "10", "15", "20", "25")
  } else {
    df$N_position <- with(df, dplyr::case_when(
      N == 25 ~ 1, N == 20 ~ 2, N == 15 ~ 3, N == 10 ~ 4, N == 5 ~ 5, N == 3 ~ 6))
    y_mean <- "FDR_mean"; y_se <- "FDR_se"
    y_lab  <- "95% CS FDR"; y_lim <- c(0, 0.35)
    x_labs <- c("25", "20", "15", "10", "5", "3")
  }

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
    scale_x_continuous(breaks = 1:6, labels = x_labs) +
    labs(x = "Top N as Causal Variants", y = y_lab) +
    theme_panel +
    coord_cartesian(ylim = y_lim, xlim = c(0.5, 6.5))

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted", color = "red", linewidth = 0.8)
  }
  p
}

# =============================================================================
# Build the 3 scenario rows (each: title + [Power | FDR])
# =============================================================================

cat("Building S8 panels...\n")

build_row <- function(sc, lab_pw, lab_fd) {
  pw <- make_panel(sc, "Power")
  fd <- make_panel(sc, "FDR")
  plot_grid(pw, fd, nrow = 1,
            labels = c(lab_pw, lab_fd),
            label_size = 18, label_fontface = "bold")
}

row_S2      <- build_row("Complex S2", "A", "B")
row_S1      <- build_row("Complex S1", "C", "D")
row_complex <- build_row("Complex",    "E", "F")

# =============================================================================
# Shared legend (color + linetype, one row)
# =============================================================================

# Legend uses filled SQUARES (fill aesthetic) in a 2-row x 3-col, column-major
# layout: row 1 = no extension, row 2 = w/ extension; columns = SuSiE, ash, inf.
legend_plot <- ggplot(agg[agg$scenario == "Complex", ],
                      aes(x = N, y = Power_mean, fill = series)) +
  geom_col() +
  scale_fill_manual(values = series_colors) +
  guides(fill = guide_legend(nrow = 2, byrow = FALSE)) +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"))

legend <- get_legend(legend_plot)

# =============================================================================
# Assemble + save
# =============================================================================

final_figure <- plot_grid(
  row_S2, row_S1, row_complex, legend,
  ncol = 1,
  rel_heights = c(1, 1, 1, 0.20)
)

out_pdf <- file.path(s8_dir, "S8.pdf")
ggsave(out_pdf, final_figure, width = 12, height = 15, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s8_dir, "S8.png")
ggsave(out_png, final_figure, width = 12, height = 15, units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
