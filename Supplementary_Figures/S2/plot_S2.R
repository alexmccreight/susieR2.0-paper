#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S2: Unmappable Effects Figure
# =============================================================================
# Row 1: Ground truth (2 panels, 50% each)
#   A: True effect sizes    B: Z-scores (-log10 p)
# Row 2: PIP plots (4 panels, 25% each)
#   C: SuSiE (L=10)   D: SuSiE (L=40)   E: SuSiE-inf   F: SuSiE-ash
# Row 3: CS on true effects (4 panels, 25% each)
#   G: SuSiE (L=10)   H: SuSiE (L=40)   I: SuSiE-inf   J: SuSiE-ash
# =============================================================================

library(ggplot2)
library(cowplot)

script_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S2"
data_dir   <- file.path(script_dir, "data")

cat("Creating supplementary figure S2: unmappable effects example...\n")

# =============================================================================
# Load cached results
# =============================================================================
cat("Loading cached results...\n")

data_stats <- readRDS(file.path(data_dir, "data_and_sumstats.rds"))
X <- data_stats$X
y <- data_stats$y
b <- data_stats$b
z_scores <- data_stats$z_scores
pvals <- data_stats$pvals

fit_susie    <- readRDS(file.path(data_dir, "fit_susie.rds"))
fit_inf      <- readRDS(file.path(data_dir, "fit_inf.rds"))
fit_ash      <- readRDS(file.path(data_dir, "fit_ash.rds"))
fit_susie_L40 <- readRDS(file.path(data_dir, "fit_susie_L40.rds"))

# =============================================================================
# susie_plot_gg: ggplot2 version of susie_plot
# =============================================================================
susie_plot_gg <- function(model, y, add_bar = FALSE, pos = NULL, b = NULL,
                          max_cs = 400, add_legend = FALSE,
                          point_size = 1.5, cs_point_size = 3.5,
                          causal_point_size = 3.5,
                          title = NULL) {

  is_susie <- inherits(model, "susie")

  cs_colors <- c(
    "dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "gold1",
    "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
    "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
    "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )

  if (y == "z") {
    if (is_susie) {
      if (is.null(model$z)) stop("z-scores are not available from SuSiE fit")
      zneg <- -abs(model$z)
    } else {
      zneg <- -abs(model)
    }
    p_vals <- -log10(2 * pnorm(zneg))
    ylab <- expression(-log[10](p))
  } else if (y == "PIP") {
    if (is_susie) {
      p_vals <- model$pip
    } else {
      p_vals <- model
    }
    ylab <- "PIP"
  } else {
    if (is_susie) {
      stop("Need to specify z or PIP for SuSiE fits")
    }
    p_vals <- model
    ylab <- y
  }

  if (is.null(b)) b <- rep(0, length(p_vals))
  if (is.null(pos)) pos <- 1:length(p_vals)

  df <- data.frame(
    position = pos,
    value = p_vals[pos],
    is_causal = b[pos] != 0
  )

  df$in_cs <- FALSE
  df$cs_id <- NA
  df$cs_color <- NA
  legend_labels <- NULL
  legend_colors <- NULL

  if (is_susie && !is.null(model$sets$cs)) {
    color_idx <- 1
    for (cs_idx in seq_along(model$sets$cs)) {
      cs_vars <- model$sets$cs[[cs_idx]]
      purity <- if (!is.null(model$sets$purity)) model$sets$purity[cs_idx, 1] else 1

      include_cs <- FALSE
      if (!is.null(model$sets$purity) && max_cs < 1 && purity >= max_cs) {
        include_cs <- TRUE
      } else if (length(cs_vars) < max_cs) {
        include_cs <- TRUE
      }

      if (include_cs) {
        cs_in_pos <- intersect(pos, cs_vars)
        if (length(cs_in_pos) > 0) {
          mask <- df$position %in% cs_in_pos
          df$in_cs[mask] <- TRUE
          df$cs_id[mask] <- cs_idx
          df$cs_color[mask] <- cs_colors[color_idx]

          effect_label <- if (!is.null(model$sets$cs_index)) {
            model$sets$cs_index[cs_idx]
          } else {
            cs_idx
          }

          if (length(cs_in_pos) == 1) {
            label_text <- paste0("L", effect_label, ": C=1")
          } else {
            label_text <- paste0("L", effect_label, ": C=", length(cs_in_pos),
                                "/R=", round(purity, 4))
          }

          legend_labels <- c(legend_labels, label_text)
          legend_colors <- c(legend_colors, cs_colors[color_idx])

          color_idx <- color_idx + 1
          if (color_idx > length(cs_colors)) color_idx <- 1
        }
      }
    }
  }

  p <- ggplot(df, aes(x = .data$position, y = .data$value))

  if (add_bar && any(df$in_cs)) {
    cs_df <- df[df$in_cs, ]
    p <- p + geom_segment(
      data = cs_df,
      aes(x = .data$position, xend = .data$position, y = 0, yend = .data$value),
      color = "gray", linewidth = 0.8
    )
  }

  df_noncausal <- df[!df$is_causal, ]
  df_causal <- df[df$is_causal, ]

  if (nrow(df_noncausal) > 0) {
    p <- p + geom_point(data = df_noncausal, color = "black",
                        size = point_size, alpha = 0.7)
  }
  if (nrow(df_causal) > 0) {
    p <- p + geom_point(data = df_causal, color = "red",
                        size = causal_point_size, alpha = 0.9)
  }

  if (any(df$in_cs)) {
    cs_df <- df[df$in_cs, ]
    p <- p + geom_point(
      data = cs_df,
      color = cs_df$cs_color,
      size = cs_point_size,
      shape = 21,
      stroke = 2,
      fill = NA
    )
  }

  p <- p +
    labs(x = "Index", y = ylab, title = title) +
    theme_cowplot(font_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = if (add_legend) "right" else "none"
    )

  return(p)
}

# =============================================================================
# Helper: CS on true effects plot
# =============================================================================
make_cs_effects_plot <- function(fit, b, title) {
  colors <- c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "gold1",
              "skyblue2", "#FB9A99", "palegreen2")

  df <- data.frame(
    pos = seq_along(b),
    effect = abs(b),
    causal = b != 0,
    in_cs = FALSE,
    cs_id = NA
  )

  if (!is.null(fit$sets$cs)) {
    for (i in seq_along(fit$sets$cs)) {
      cs_idx <- fit$sets$cs[[i]]
      df$in_cs[cs_idx] <- TRUE
      df$cs_id[cs_idx] <- i
    }
  }

  p <- ggplot(df, aes(x = pos, y = effect)) +
    geom_point(color = "black", size = 1.5, alpha = 0.7) +
    labs(x = "Index", y = "Absolute Effect Size", title = title) +
    theme_cowplot(font_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))

  if (any(df$in_cs)) {
    cs_df <- df[df$in_cs, ]
    p <- p + geom_point(data = cs_df, aes(color = factor(cs_id)), size = 3.5) +
      scale_color_manual(values = colors)
  }

  p
}

# =============================================================================
# ROW 1: Ground Truth (2 panels)
# =============================================================================
cat("Creating Row 1: Ground truth panels...\n")

df_effects <- data.frame(
  pos = seq_along(b),
  effect = abs(b),
  causal = b != 0
)

panel_A <- ggplot(df_effects, aes(x = pos, y = effect)) +
  geom_point(data = df_effects[!df_effects$causal, ], color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = df_effects[df_effects$causal, ], color = "red", size = 3.5, alpha = 0.9) +
  labs(x = "Index", y = "Absolute Effect Size") +
  theme_cowplot(font_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(5, 10, 5, 5))

df_zscore <- data.frame(
  pos = seq_along(z_scores),
  neglog10p = -log10(pvals),
  causal = b != 0
)

panel_B <- ggplot(df_zscore, aes(x = pos, y = neglog10p)) +
  geom_point(data = df_zscore[!df_zscore$causal, ], color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = df_zscore[df_zscore$causal, ], color = "red", size = 3.5, alpha = 0.9) +
  labs(x = "Index", y = expression(-log[10](p))) +
  theme_cowplot(font_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 10))

row1 <- cowplot::plot_grid(panel_A, panel_B, nrow = 1, rel_widths = c(1, 1),
                           labels = c("A", "B"), label_size = 14, label_fontface = "bold")

# =============================================================================
# ROW 2: PIP Plots (4 panels)
# =============================================================================
cat("Creating Row 2: PIP plots...\n")

panel_C <- susie_plot_gg(fit_susie, y = "PIP", b = b, title = "SuSiE") +
  theme(plot.margin = margin(5, 5, 5, 5))
panel_D <- susie_plot_gg(fit_susie_L40, y = "PIP", b = b, title = "SuSiE (L=40)") +
  theme(plot.margin = margin(5, 5, 5, 5))
panel_E <- susie_plot_gg(fit_inf, y = "PIP", b = b, title = "SuSiE-inf") +
  theme(plot.margin = margin(5, 5, 5, 5))
panel_F <- susie_plot_gg(fit_ash, y = "PIP", b = b, title = "SuSiE-ash") +
  theme(plot.margin = margin(5, 5, 5, 5))

row2 <- cowplot::plot_grid(panel_C, panel_D, panel_E, panel_F, nrow = 1,
                           rel_widths = c(1, 1, 1, 1),
                           labels = c("C", "D", "E", "F"),
                           label_size = 14, label_fontface = "bold")

# =============================================================================
# ROW 3: CS on True Effects (4 panels)
# =============================================================================
cat("Creating Row 3: CS on true effects plots...\n")

panel_G <- make_cs_effects_plot(fit_susie, b, NULL) +
  theme(plot.margin = margin(5, 5, 5, 5))
panel_H <- make_cs_effects_plot(fit_susie_L40, b, NULL) +
  theme(plot.margin = margin(5, 5, 5, 5))
panel_I <- make_cs_effects_plot(fit_inf, b, NULL) +
  theme(plot.margin = margin(5, 5, 5, 5))
panel_J <- make_cs_effects_plot(fit_ash, b, NULL) +
  theme(plot.margin = margin(5, 5, 5, 5))

row3 <- cowplot::plot_grid(panel_G, panel_H, panel_I, panel_J, nrow = 1,
                           rel_widths = c(1, 1, 1, 1),
                           labels = c("G", "H", "I", "J"),
                           label_size = 14, label_fontface = "bold")

# =============================================================================
# Combine all rows
# =============================================================================
cat("Combining all rows into final figure...\n")

final_figure <- cowplot::plot_grid(
  row1, row2, row3,
  ncol = 1,
  rel_heights = c(1, 1, 1)
)

output_pdf <- file.path(script_dir, "S2.pdf")
ggsave(output_pdf, final_figure,
       width = 13, height = 9, units = "in", bg = "white")

output_png <- file.path(script_dir, "S2.png")
ggsave(output_png, final_figure,
       width = 13, height = 9, units = "in", dpi = 300, bg = "white")

cat(sprintf("\nFigure saved to:\n  - %s\n  - %s\n", output_pdf, output_png))

# Print summary
cat("\n=== Summary ===\n")
cat(sprintf("SuSiE (L=10):  %d credible sets\n", length(fit_susie$sets$cs)))
cat(sprintf("SuSiE (L=40):  %d credible sets\n", length(fit_susie_L40$sets$cs)))
cat(sprintf("SuSiE-inf:     %d credible sets\n", length(fit_inf$sets$cs)))
cat(sprintf("SuSiE-ash:     %d credible sets\n", length(fit_ash$sets$cs)))
