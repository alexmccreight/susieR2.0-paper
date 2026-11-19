#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S4: 2x2 figure with Panel E (manuscript Figure S7)
# =============================================================================
# Panel A (top-left):     Power bars and standard errors across K
# Panel B (top-right):    FDR bars and standard errors across K
# Panel C (bottom-left):  Purity violin (pooled across K and replicates)
# Panel D (bottom-right): CS size box plot (pooled)
# Panel E (third row):    Full posterior mass under two CS reporting rules
#
# Inputs:
#   - data/s4_metrics.rds  (built by extract_S4_data.R)
#   - data/paired-posterior-mass.csv  (built by simulate_S4E.R)
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

source("R/paths.R")

s4_dir   <- supp_dir(4)
data_dir <- file.path(s4_dir, "data")
data_path <- file.path(data_dir, "s4_metrics.rds")
paired_path <- file.path(data_dir, "paired-posterior-mass.csv")

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
# Panel E: Cumulative distribution of full posterior mass in reported CSs
# =============================================================================
# Both methods receive the same fitted single-effect posterior probabilities.
# The dashed line marks 95%; the key reports how many sets fall below it.
# Read the saved table here. Simulation is run separately in simulate_S4E.R.

paired <- read.csv(paired_path)
reporting_methods <- c("susieR", "SparsePro (2023)")
reporting_colors <- c("susieR" = "#E67E22", "SparsePro (2023)" = "#16758A")
mass_columns <- c("mass_susier", "mass_sparsepro")
curves <- do.call(rbind, lapply(seq_along(reporting_methods), function(i) {
  mass <- sort(paired[[mass_columns[i]]] * 100)
  data.frame(Method = reporting_methods[i], mass = c(90, mass, 100),
             cumulative = c(0, seq_along(mass) / length(mass) * 100, 100))
}))
counts <- vapply(paired[mass_columns], function(x) sum(x < 0.95 - 1e-12), integer(1))

panel_E <- ggplot(curves, aes(mass, cumulative, color = Method)) +
  annotate("rect", xmin = 90, xmax = 95, ymin = -Inf, ymax = Inf,
           fill = "#F4F5F6") +
  geom_vline(xintercept = 95, linetype = "dashed", color = "#555555", linewidth = 0.6) +
  geom_step(linewidth = 1.1, direction = "hv") +
  scale_color_manual(values = reporting_colors) +
  scale_x_continuous(breaks = c(90, 92.5, 95, 97.5, 100),
                     labels = c("90", "92.5", "95", "97.5", "100"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(90, 100), ylim = c(0, 102)) +
  labs(x = "Full posterior mass in credible set (%)", y = "Cumulative sets (%)") +
  base_theme + theme(legend.position = "none", plot.margin = margin(18, 8, 8, 18))

# Build the Panel E key with the number and percentage of CSs below 95%.
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

# =============================================================================
# Compose 2x2 with shared legend below panels C/D, then Panel E
# =============================================================================

# Build a centered horizontal legend (placed under the second row).
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
panels <- align_plots(panel_A + theme(legend.position = "none"),
                      panel_B + theme(legend.position = "none"),
                      panel_C, panel_D, align = "hv", axis = "tblr")
top_rows <- plot_grid(plotlist = panels, ncol = 2, labels = LETTERS[1:4],
                       label_size = 16, label_fontface = "bold")
coverage_row <- plot_grid(panel_E, reporting_key, nrow = 1,
                          rel_widths = c(0.63, 0.37), labels = c("E", ""),
                          label_size = 16, label_fontface = "bold")
final_figure <- plot_grid(top_rows, shared_legend, coverage_row, ncol = 1,
                          rel_heights = c(2, 0.13, 0.68))

# =============================================================================
# Save
# =============================================================================

fig_width  <- 12
fig_height <- 11.2

out_pdf <- file.path(s4_dir, "S4.pdf")
ggsave(out_pdf, final_figure, width = fig_width, height = fig_height,
       units = "in", bg = "white")
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s4_dir, "S4.png")
ggsave(out_png, final_figure, width = fig_width, height = fig_height,
       units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
