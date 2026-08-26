#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S7: Sparse Power / FDR -- no extension vs LD extension (0.99)
# =============================================================================
# Same style as S1A/S1B, but two bars per method per K (6 bars per K):
#   base color  = "SuSiE" / "SuSiE-inf" / "SuSiE-ash"   (Xcorr, no extension)
#   light tint  = "... w/ Extension"                    (Xcorr + LD ext 0.99)
#
# Panel A = 95% CS Power, Panel B = 95% CS FDR (+ 0.05 reference line).
#
# Input:  data/s7_finemapping_data.rds  (built by extract_S7_data.R)
# Output: S7.pdf, S7.png
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(cowplot)
})

source("R/paths.R")

s7_dir    <- supp_dir(7)
data_path <- file.path(s7_dir, "data", "s7_finemapping_data.rds")
if (!file.exists(data_path)) stop("Run extract_S7_data.R first.")

dat <- readRDS(data_path)

# =============================================================================
# Aesthetics
# =============================================================================

series_levels <- c("SuSiE", "SuSiE w/ Extension",
                   "SuSiE-ash", "SuSiE-ash w/ Extension",
                   "SuSiE-inf", "SuSiE-inf w/ Extension")
dat$series <- factor(dat$series, levels = series_levels)

# base method colors + lighter "w/ Extension" tints (~50% toward white)
series_colors <- c(
  "SuSiE"                  = "#4A90E2",
  "SuSiE w/ Extension"     = "#A6C8F0",
  "SuSiE-inf"              = "#7CB342",
  "SuSiE-inf w/ Extension" = "#BCD9A0",
  "SuSiE-ash"              = "#E53935",
  "SuSiE-ash w/ Extension" = "#F29C9A"
)

theme_bar <- theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold", size = 18),
    axis.text    = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin  = margin(10, 10, 10, 10)
  )

# =============================================================================
# Bar panel (mirrors S1 create_bar_panel, fill = 6-level series)
# =============================================================================

create_bar_panel <- function(data, metric = "Power") {
  if (metric == "Power") {
    y_label  <- "95% CS Power"; y_mean <- "Power_mean"; y_se <- "Power_se"
    y_limits <- c(0.3, 0.85)
  } else {
    y_label  <- "95% CS FDR";   y_mean <- "FDR_mean";   y_se <- "FDR_se"
    y_limits <- c(0, 0.1)
  }

  k_labels <- data.frame(Total_PVE = c(0.03, 0.06, 0.09, 0.12, 0.15), K = 1:5)

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

  for (i in seq_len(nrow(k_labels))) {
    p <- p + annotate("text", x = i, y = y_limits[2] * 0.99,
                       label = paste0("K = ", k_labels$K[i]),
                       size = 4.5, fontface = "bold")
  }

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted",
                         color = "red", linewidth = 0.8)
  }

  p
}

cat("Creating S7 bar panels...\n")
panel_A <- create_bar_panel(dat, "Power")
panel_B <- create_bar_panel(dat, "FDR")

# =============================================================================
# Shared bottom legend (6 entries in one row)
# =============================================================================

legend_plot <- ggplot(dat, aes(x = factor(Total_PVE), y = Power_mean, fill = series)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = series_colors) +
  guides(fill = guide_legend(nrow = 2, byrow = FALSE)) +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 13))

legend <- get_legend(legend_plot)

# =============================================================================
# Assemble: A + B side by side, legend beneath
# =============================================================================

row <- plot_grid(panel_A, panel_B, nrow = 1,
                 labels = c("A", "B"),
                 label_size = 20, label_fontface = "bold")

final_figure <- plot_grid(row, legend, ncol = 1, rel_heights = c(1, 0.16))

out_pdf <- file.path(s7_dir, "S7.pdf")
ggsave(out_pdf, final_figure, width = 16, height = 7, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s7_dir, "S7.png")
ggsave(out_png, final_figure, width = 16, height = 7, units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
