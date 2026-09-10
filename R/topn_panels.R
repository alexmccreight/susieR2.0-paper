# =============================================================================
# Shared Top-N power / FDR panel for the complex-architecture benchmarks
# =============================================================================
# Single source of truth for Figure 1 panels C-F and their three-method
# counterparts in Supplementary Figure S10. Both call plot_topn_panel() on the
# same summary tables (figure_1/data/complex_S{1,2}_finemapping_summary.rds),
# differing only in which methods are drawn.
#
# Axis convention (deliberate, stated in the manuscript legends): power panels
# read N = 3 -> 25 from left to right and FDR panels read 25 -> 3, so that
# performance declines from left to right in every panel.
#
# Usage:
#     source("R/aesthetics.R")
#     source("R/topn_panels.R")
#     p <- plot_topn_panel(summary_df, metric = "Power", methods = MAIN_METHODS)
# =============================================================================

plot_topn_panel <- function(summary_df,
                            metric   = c("Power", "FDR"),
                            methods  = MAIN_METHODS,
                            y_limits = NULL) {
  metric <- match.arg(metric)

  df <- summary_df[summary_df$Method %in% methods, , drop = FALSE]
  df$Method <- factor(df$Method, levels = methods)

  if (metric == "Power") {
    x_pos    <- c(`3` = 1, `5` = 2, `10` = 3, `15` = 4, `20` = 5, `25` = 6)
    x_labels <- c("3", "5", "10", "15", "20", "25")
    if (is.null(y_limits)) y_limits <- c(0, 0.9)
  } else {
    x_pos    <- c(`25` = 1, `20` = 2, `15` = 3, `10` = 4, `5` = 5, `3` = 6)
    x_labels <- c("25", "20", "15", "10", "5", "3")
    if (is.null(y_limits)) y_limits <- c(0, 0.35)
  }
  df$N_position <- unname(x_pos[as.character(df$N)])

  y_col  <- paste0(metric, "_mean")
  se_col <- paste0(metric, "_se")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = N_position, y = .data[[y_col]],
                                        color = Method, group = Method)) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[y_col]] - .data[[se_col]],
                   ymax = .data[[y_col]] + .data[[se_col]]),
      width = 0.25, linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 3.5)

  if (metric == "FDR") {
    p <- p + ggplot2::geom_hline(yintercept = 0.05, linetype = "dotted",
                                 color = "red", linewidth = 0.8)
  }

  p +
    ggplot2::scale_color_manual(values = method_colors, breaks = methods) +
    ggplot2::scale_x_continuous(breaks = 1:6, labels = x_labels) +
    ggplot2::labs(x = "Top N as Causal Variants", y = paste0("95% CS ", metric)) +
    theme_panel +
    ggplot2::coord_cartesian(ylim = y_limits, xlim = c(0.5, 6.5))
}
