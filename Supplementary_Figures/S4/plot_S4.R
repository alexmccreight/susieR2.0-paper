#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S4: 2x2 figure
# =============================================================================
# Panel A (top-left):     Power dot+errorbar across K, 5 method variants
# Panel B (top-right):    FDR   dot+errorbar across K, 5 method variants
# Panel C (bottom-left):  Purity violin (5 variants, pooled across K and reps)
# Panel D (bottom-right): CS size box plot (5 variants, pooled)
#
# Input:  data/s4_metrics.rds  (built by extract_S4_data.R)
# Output: S4.pdf, S4.png
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

# =============================================================================
# Paths
# =============================================================================

s4_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S4"
data_dir <- file.path(s4_dir, "data")
data_path <- file.path(data_dir, "s4_metrics.rds")

if (!file.exists(data_path)) {
  stop("Run extract_S4_data.R first to produce ", data_path)
}

# =============================================================================
# Load
# =============================================================================

cat("Loading data...\n")
dat <- readRDS(data_path)

rep_metrics <- dat$rep_metrics
cs_purity   <- dat$cs_purity
cs_size     <- dat$cs_size

# =============================================================================
# Filter out SuSiE+BLiP(q=0.10) and relabel methods for display.
# =============================================================================

# Mapping: source method label  ->  display label
# (SuSiE+BLiP(q0.10) intentionally absent so it gets dropped from all panels.)
rename_map <- c(
  "SuSiE(0.5)"        = "SuSiE (purity = 0.5)",
  "SuSiE(0.8)"        = "SuSiE (purity = 0.8)",
  "SuSiE+BLiP(q0.05)" = "SuSiE+BLiP",
  "SuSiE+attainable"  = "SuSiE+attainable"
)
keep_methods <- names(rename_map)
method_levels <- unname(rename_map)   # left-to-right display order

apply_rename <- function(df) {
  df <- df[as.character(df$Method) %in% keep_methods, , drop = FALSE]
  df$Method <- factor(rename_map[as.character(df$Method)], levels = method_levels)
  df
}
rep_metrics <- apply_rename(rep_metrics)
cs_purity   <- apply_rename(cs_purity)
cs_size     <- apply_rename(cs_size)

# =============================================================================
# Aesthetics
# =============================================================================

method_colors <- c(
  "SuSiE (purity = 0.5)" = "#4A90E2",   # SuSiE blue (default purity)
  "SuSiE (purity = 0.8)" = "#1E5A9B",   # darker blue (purity-filtered)
  "SuSiE+BLiP"           = "#9B59B6",   # purple (BLiP)
  "SuSiE+attainable"     = "#E67E22"    # orange (attainable)
)

base_theme <- theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(color = "black", size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(size = 10),
    plot.margin = margin(8, 8, 8, 8)
  )

# =============================================================================
# Aggregate per (Method, K) for Panels A and B
# =============================================================================

agg <- rep_metrics %>%
  group_by(Method, K, Total_PVE) %>%
  summarise(
    Power_mean = mean(Power, na.rm = TRUE),
    Power_se   = sd(Power,  na.rm = TRUE) / sqrt(sum(!is.na(Power))),
    FDR_mean   = mean(FDR,   na.rm = TRUE),
    FDR_se     = sd(FDR,    na.rm = TRUE) / sqrt(sum(!is.na(FDR))),
    .groups = "drop"
  )

# =============================================================================
# Panels A & B: Power / FDR bar plots (matches S1 sparse style)
# =============================================================================

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
    y_limits <- c(0,   0.1)
  }

  k_labels <- data.frame(Total_PVE = c(0.03, 0.06, 0.09, 0.12, 0.15), K = 1:5)

  p <- ggplot(data, aes(x = factor(Total_PVE), y = .data[[y_mean]],
                         fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.75) +
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

  for (i in seq_len(nrow(k_labels))) {
    p <- p + annotate("text", x = i + 0.05, y = y_limits[2] * 0.99,
                       label = paste0("K = ", k_labels$K[i]),
                       size = 4.5, fontface = "bold")
  }

  if (metric == "FDR") {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dotted",
                         color = "red", linewidth = 0.8)
  }

  p
}

panel_A <- create_bar_panel(agg, metric = "Power")
panel_B <- create_bar_panel(agg, metric = "FDR")

# =============================================================================
# Panel C: Purity violin
# =============================================================================

panel_C <- ggplot(cs_purity, aes(x = Method, y = purity, fill = Method)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85,
              color = "grey20", linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "grey20",
               fill = "white", alpha = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1.02),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Purity") +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

# =============================================================================
# Panel D: CS size box plot (log y to handle the range)
# =============================================================================

panel_D <- ggplot(cs_size, aes(x = Method, y = size, fill = Method)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4,
               color = "grey20", linewidth = 0.4) +
  scale_fill_manual(values = method_colors) +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "CS Size (log scale)") +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

# =============================================================================
# Compose 2x2 with shared legend below panels C/D, centered
# =============================================================================

# Build a centered horizontal legend (placed under the bottom row)
shared_legend <- get_legend(
  panel_A + theme(
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.justification = "center",
    legend.title         = element_blank(),
    legend.text          = element_text(size = 11),
    legend.key.width     = unit(18, "pt"),
    legend.box.margin    = margin(0, 0, 0, 0)
  )
)

panel_A_nolegend <- panel_A + theme(legend.position = "none")
panel_B_nolegend <- panel_B + theme(legend.position = "none")

top_row <- plot_grid(panel_A_nolegend, panel_B_nolegend,
                     nrow = 1, rel_widths = c(1, 1),
                     labels = c("A", "B"), label_size = 16,
                     label_fontface = "bold")

bottom_row <- plot_grid(panel_C, panel_D,
                        nrow = 1, rel_widths = c(1, 1),
                        labels = c("C", "D"), label_size = 16,
                        label_fontface = "bold")

final_figure <- plot_grid(top_row, bottom_row, shared_legend, ncol = 1,
                          rel_heights = c(1, 1, 0.12))

# =============================================================================
# Save
# =============================================================================

fig_width  <- 14
fig_height <- 10

out_pdf <- file.path(s4_dir, "S4.pdf")
ggsave(out_pdf, final_figure, width = fig_width, height = fig_height,
       units = "in", bg = "white")
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s4_dir, "S4.png")
ggsave(out_png, final_figure, width = fig_width, height = fig_height,
       units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
