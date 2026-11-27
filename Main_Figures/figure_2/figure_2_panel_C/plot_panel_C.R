#!/usr/bin/env Rscript

# Figure 2, panel C: number of signals recovered by both methods or by one
# method only.
#
# With two methods the three bars are identified by colour alone: the table in
# panel D names the same groups in the same colours. With three methods
# (Supplementary Figure S10) an UpSet-style dot matrix under the bars shows
# which methods make up each of the seven groups.
#
# Usage:  plot_panel_C.R [panel_C_data.rds] [output basename]
# Input:  panel_C_data.rds  (prepare_panel_C_data.R)
# Output: panel_C.pdf, panel_C.png, panel_C_plot.rds

library(ggplot2)

source("R/paths.R")
source("R/aesthetics.R")

script_dir <- file.path(fig_dir(2), "figure_2_panel_C")
data_path  <- cli_arg(1, file.path(script_dir, "panel_C_data.rds"))
out_base   <- cli_arg(2, file.path(script_dir, "panel_C"))

panel_data   <- readRDS(data_path)
concordance  <- panel_data$concordance        # ordered Consensus, pairwise, single-method
membership   <- panel_data$group_membership   # groups x methods, logical
plot_methods <- panel_data$methods            # display names, plot order

concordance$x_pos <- seq_len(nrow(concordance))

# ---- Dot matrix (three methods only) ----

use_dots <- length(plot_methods) > 2

if (use_dots) {
  # One row of dots per method, below y = 0.
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

  # A line joins the active dots of each group.
  line_data <- do.call(rbind, lapply(concordance$x_pos, function(xp) {
    active_y <- dot_data$y_pos[dot_data$x_pos == xp & dot_data$active]
    if (length(active_y) < 2) return(NULL)
    data.frame(x_pos = xp, ymin = min(active_y), ymax = max(active_y))
  }))
}

# ---- Plot ----

# The axis lines are drawn as segments so that the y-axis stops at 0 instead of
# running down through the dot rows. In Figure 2 their width is the one
# theme_classic gives its own axis lines (base_size / 22), to match the other
# panels; the three-method panel of Figure S10 uses 0.5.
BASE_SIZE  <- 20
AXIS_LWD   <- if (use_dots) 0.5 else BASE_SIZE / 22

panel_C <- ggplot() +
  geom_col(data = concordance,
           aes(x = x_pos, y = n_signals, fill = group),
           width = 0.6) +
  geom_text(data = concordance,
            aes(x = x_pos, y = n_signals, label = n_signals),
            vjust = -0.5, size = 6, fontface = "bold") +
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
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = Inf,
           color = "black", linewidth = AXIS_LWD) +
  annotate("segment", x = -Inf, xend = nrow(concordance) + 0.5, y = 0, yend = 0,
           color = "black", linewidth = AXIS_LWD) +
  scale_fill_manual(values = concordance_colors) +
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
    # With the dot matrix the bottom margin holds the three dot rows and the
    # left margin holds their method labels; without it neither is needed.
    plot.margin   = if (use_dots) margin(5, 10, 80, 100) else margin(5, 10, 5, 10)
  )

save_figure(panel_C, dirname(out_base), basename(out_base), width = 6, height = 5)

plot_path <- paste0(out_base, "_plot.rds")
saveRDS(panel_C, plot_path)
cat("Saved:", plot_path, "\n")
