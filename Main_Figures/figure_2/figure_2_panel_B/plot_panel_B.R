#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel B — Plotting (UpSet-style)
# =============================================================================
# Creates an UpSet-style plot showing cross-method concordance of credible sets:
#   Top:    Bar chart of signal counts per concordance group
#   Bottom: Dot matrix showing method membership for each group
#
# Built as a single ggplot to guarantee bar-dot alignment.
#
# Input:  panel_B_data.rds (from prepare_panel_B_data.R)
# Output: panel_B.pdf, panel_B.png, panel_B_plot.rds
# =============================================================================

library(ggplot2)

# =============================================================================
# Paths
# =============================================================================

script_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2/figure_2_panel_B"
data_path  <- file.path(script_dir, "panel_B_data.rds")

# =============================================================================
# Shared aesthetics
# =============================================================================

concordance_colors <- c(
  "Consensus"   = "#888888",
  "SuSiE&ash"   = "#9C27B0",   # Blue + Red -> Purple
  "SuSiE&inf"   = "#2E9D8F",   # Blue + Green -> Teal
  "ash&inf"     = "#FF9800",   # Red + Green -> Orange
  "SuSiE"       = "#4A90E2",   # Blue (matches Panel A SuSiE)
  "SuSiE-ash"   = "#E53935",   # Red (matches Panel A SuSiE-ash)
  "SuSiE-inf"   = "#7CB342"    # Green (matches Panel A SuSiE-inf)
)

# =============================================================================
# Load data
# =============================================================================

cat("Loading panel B data...\n")
panel_data    <- readRDS(data_path)
concordance   <- panel_data$concordance
method_totals <- panel_data$method_totals

# Remap group names for consistent method-name capitalization
group_remap <- c(
  "Consensus"   = "Consensus",
  "SuSiE + ASH" = "SuSiE&ash",
  "SuSiE + INF" = "SuSiE&inf",
  "ASH + INF"   = "ash&inf",
  "SuSiE only"  = "SuSiE",
  "ASH only"    = "SuSiE-ash",
  "INF only"    = "SuSiE-inf"
)
levels(concordance$group) <- group_remap[levels(concordance$group)]

# Order groups structurally: Consensus -> Pairwise -> Unique
concordance <- concordance[order(concordance$group), ]
concordance$x_pos <- seq_len(nrow(concordance))

# =============================================================================
# Build dot matrix data using the same numeric x positions
# =============================================================================

cat("Building dot matrix...\n")

# Dot y positions below the bars (negative y space)
y_max <- max(concordance$n_signals)
dot_spacing <- y_max * 0.10
dot_y <- c("SuSiE" = -dot_spacing, "SuSiE-ash" = -2 * dot_spacing,
           "SuSiE-inf" = -3 * dot_spacing)

# Build dot data
dot_data <- expand.grid(
  x_pos  = concordance$x_pos,
  method = c("SuSiE", "SuSiE-ash", "SuSiE-inf"),
  stringsAsFactors = FALSE
)
dot_data$y_pos <- dot_y[dot_data$method]

# Determine active/inactive for each dot
dot_data$active <- mapply(function(xp, m) {
  row <- concordance[concordance$x_pos == xp, ]
  if (m == "SuSiE")     return(row$has_susie)
  if (m == "SuSiE-ash") return(row$has_ash)
  if (m == "SuSiE-inf") return(row$has_inf)
  FALSE
}, dot_data$x_pos, dot_data$method)

# Build connecting line segments for active dots within each group
line_data <- do.call(rbind, lapply(concordance$x_pos, function(xp) {
  active_y <- dot_data$y_pos[dot_data$x_pos == xp & dot_data$active]
  if (length(active_y) < 2) return(NULL)
  data.frame(x_pos = xp, ymin = min(active_y), ymax = max(active_y),
             stringsAsFactors = FALSE)
}))

# =============================================================================
# Build single ggplot
# =============================================================================

cat("Creating plot...\n")

panel_B <- ggplot() +
  # Bars
  geom_col(data = concordance,
           aes(x = x_pos, y = n_signals, fill = group),
           width = 0.6) +
  # Count labels above bars
  geom_text(data = concordance,
            aes(x = x_pos, y = n_signals, label = n_signals),
            vjust = -0.5, size = 6, fontface = "bold") +
  # Connecting lines between active dots
  { if (!is.null(line_data))
      geom_segment(data = line_data,
                   aes(x = x_pos, xend = x_pos, y = ymin, yend = ymax),
                   linewidth = 1.2, color = "black") } +
  # Inactive dots
  geom_point(data = dot_data[!dot_data$active, ],
             aes(x = x_pos, y = y_pos),
             size = 4, color = "gray80") +
  # Active dots
  geom_point(data = dot_data[dot_data$active, ],
             aes(x = x_pos, y = y_pos),
             size = 4, color = "black") +
  # Y-axis line from 0 upward only (avoid cutting into dot labels below)
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = Inf,
           color = "black", linewidth = 0.5) +
  # X-axis line along y = 0
  annotate("segment", x = -Inf, xend = nrow(concordance) + 0.5, y = 0, yend = 0,
           color = "black", linewidth = 0.5) +
  # Method labels on left side (outside plot area)
  annotate("text", x = 0.35, y = dot_y,
           label = names(dot_y), hjust = 1, size = 7, fontface = "bold") +
  # Scales
  scale_fill_manual(values = concordance_colors) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_cartesian(xlim = c(0.5, nrow(concordance) + 0.5), clip = "off") +
  labs(y = "Number of Signals") +
  theme_classic(base_size = 20) +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.x   = element_blank(),
    axis.line.y   = element_blank(),
    axis.title.y  = element_text(face = "bold", size = 20),
    axis.text.y   = element_text(color = "black", size = 18),
    legend.position = "none",
    plot.margin   = margin(5, 10, 80, 100)
  )

# =============================================================================
# Save
# =============================================================================

output_pdf <- file.path(script_dir, "panel_B.pdf")
ggsave(output_pdf, panel_B, width = 6, height = 5, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(script_dir, "panel_B.png")
ggsave(output_png, panel_B, width = 6, height = 5, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

# Save plot object for use by create_figure2.R
plot_path <- file.path(script_dir, "panel_B_plot.rds")
saveRDS(panel_B, plot_path)
cat(sprintf("Saved: %s\n", plot_path))

cat("\nDone!\n")
