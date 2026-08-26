#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, cross-method concordance panel — Plotting
# Labelled "C" in main Figure 2 and "B" in Supplementary S11; the directory is
# figure_2_panel_B in both cases. Renders as a plain labelled bar chart for two
# methods and as an UpSet-style bars + dot matrix for three.
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

source("R/paths.R")
source("R/aesthetics.R")

script_dir <- file.path(fig_dir(2), "figure_2_panel_B")

# Optional overrides so the three-method (S11) version can be rebuilt without
# touching the main-figure outputs:
#   Rscript plot_panel_B.R <panel_B_data.rds> <output basename>
.cli <- commandArgs(trailingOnly = TRUE)
.arg <- function(i, default) if (length(.cli) >= i && nzchar(.cli[i])) .cli[i] else default
data_path <- .arg(1, file.path(script_dir, "panel_B_data.rds"))
out_base  <- .arg(2, file.path(script_dir, "panel_B"))

# =============================================================================
# Shared aesthetics
# =============================================================================

# concordance_colors comes from R/aesthetics.R. Group labels emitted by
# prepare_panel_B_data.R already match its keys.

# =============================================================================
# Load data
# =============================================================================

cat("Loading panel B data...\n")
panel_data    <- readRDS(data_path)
concordance   <- panel_data$concordance
method_totals <- panel_data$method_totals
membership    <- panel_data$group_membership   # groups x methods, logical
plot_methods  <- panel_data$methods            # display names, plot order

cat(sprintf("  %d groups across %d methods: %s\n",
            nrow(concordance), length(plot_methods),
            paste(plot_methods, collapse = ", ")))

# Groups already arrive ordered Consensus -> pairwise -> single-method.
concordance$x_pos <- seq_len(nrow(concordance))

# =============================================================================
# Build dot matrix data using the same numeric x positions
# =============================================================================

# An UpSet dot matrix earns its keep when the combinations are hard to name.
# With two methods there are only three groups, and the AlphaGenome table in the
# adjacent panel already lists them by name in the same three colours, so the
# matrix is pure redundancy and the bars are identified by colour alone. Three
# methods (Supplementary S11) keeps the matrix, where seven combinations do need
# it. Dropping the matrix also frees the wide left margin it required, which is
# what created the gap between this panel and its neighbour.
use_dots <- length(plot_methods) > 2

if (use_dots) {
  cat("Building dot matrix...\n")

  y_max <- max(concordance$n_signals)
  dot_spacing <- y_max * 0.10
  dot_y <- setNames(-dot_spacing * seq_along(plot_methods), plot_methods)

  dot_data <- expand.grid(
    x_pos  = concordance$x_pos,
    method = plot_methods,
    stringsAsFactors = FALSE
  )
  dot_data$y_pos <- dot_y[dot_data$method]

  # A dot is active when that method belongs to the group at this x position.
  dot_data$active <- membership[cbind(
    match(as.character(concordance$group)[dot_data$x_pos], rownames(membership)),
    match(dot_data$method, plot_methods)
  )]

  # Connecting line segments for active dots within each group
  line_data <- do.call(rbind, lapply(concordance$x_pos, function(xp) {
    active_y <- dot_data$y_pos[dot_data$x_pos == xp & dot_data$active]
    if (length(active_y) < 2) return(NULL)
    data.frame(x_pos = xp, ymin = min(active_y), ymax = max(active_y),
               stringsAsFactors = FALSE)
  }))
} else {
  cat("Bars identified by colour; dot matrix not needed for ", length(plot_methods),
      " methods\n", sep = "")
  dot_data <- NULL
  line_data <- NULL
  dot_y <- numeric(0)
}

# No x-axis labels: the bars are identified by colour, and the Group column of
# the AlphaGenome table in the adjacent panel names those same three colours,
# so it serves as the shared legend for both panels. Keeping labels here as
# well just repeated it, and they were long enough to need angling.

# =============================================================================
# Build single ggplot
# =============================================================================

cat("Creating plot...\n")

# The theme axis lines are blanked below and redrawn as annotations (so the
# y-axis can stop at 0 rather than running down through the dot rows). Match
# theme_classic's own weight exactly -- it derives from base_size/22 -- otherwise
# this panel's axes read visibly thinner than every other panel in the figure.
BASE_SIZE  <- 20
AXIS_LWD   <- BASE_SIZE / 22

panel_B <- ggplot() +
  # Bars
  geom_col(data = concordance,
           aes(x = x_pos, y = n_signals, fill = group),
           width = 0.6) +
  # Count labels above bars
  geom_text(data = concordance,
            aes(x = x_pos, y = n_signals, label = n_signals),
            vjust = -0.5, size = 6, fontface = "bold") +
  # Dot matrix (three-method layout only)
  { if (use_dots)
      list(
        geom_segment(data = line_data,
                     aes(x = x_pos, xend = x_pos, y = ymin, yend = ymax),
                     linewidth = 1.2, color = "black"),
        geom_point(data = dot_data[!dot_data$active, ],
                   aes(x = x_pos, y = y_pos), size = 4, color = "gray80"),
        geom_point(data = dot_data[dot_data$active, ],
                   aes(x = x_pos, y = y_pos), size = 4, color = "black"),
        annotate("text", x = 0.35, y = dot_y,
                 label = names(dot_y), hjust = 1, size = 7, fontface = "bold")
      ) } +
  # Y-axis line from 0 upward only (avoid cutting into dot labels below)
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = Inf,
           color = "black", linewidth = AXIS_LWD) +
  # X-axis line along y = 0
  annotate("segment", x = -Inf, xend = nrow(concordance) + 0.5, y = 0, yend = 0,
           color = "black", linewidth = AXIS_LWD) +
  # Scales
  scale_fill_manual(values = concordance_colors) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_x_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(0.5, nrow(concordance) + 0.5), clip = "off") +
  labs(y = "Number of Signals") +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.x   = element_blank(),
    axis.line.y   = element_blank(),
    axis.title.y  = element_text(face = "bold", size = 20),
    axis.text.y   = element_text(color = "black", size = 18),
    legend.position = "none",
    # With the dot matrix the bottom margin holds the dot rows and the left
    # margin holds their method labels; without it neither is needed.
    plot.margin   = if (use_dots) margin(5, 10, 25 + 27 * length(plot_methods), 100)
                    else margin(5, 10, 5, 10)
  )

# =============================================================================
# Save
# =============================================================================

output_pdf <- paste0(out_base, ".pdf")
ggsave(output_pdf, panel_B, width = 6, height = 5, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- paste0(out_base, ".png")
ggsave(output_png, panel_B, width = 6, height = 5, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

# Save plot object for use by create_figure2.R
plot_path <- paste0(out_base, "_plot.rds")
saveRDS(panel_B, plot_path)
cat(sprintf("Saved: %s\n", plot_path))

cat("\nDone!\n")
