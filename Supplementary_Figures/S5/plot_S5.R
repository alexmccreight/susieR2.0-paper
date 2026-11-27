#!/usr/bin/env Rscript

# Supplementary Figure S5: illustrative locus under a polygenic architecture
# (susieR's `unmappable_data` example: n = 1,000, p = 5,000, with 3 sparse,
# 5 oligogenic and 15 polygenic effects).
#   A, B   true absolute effect sizes and marginal association (-log10 p),
#          causal variants in red
#   C-E    PIPs for SuSiE, SuSiE-inf and SuSiE-ash, credible sets circled
#   F-H    the same credible sets marked on the true effect sizes
#
# Input:  data/data_and_sumstats.rds, data/fit_{susie,inf,ash}.rds  (fit_models.R)
# Output: S5.pdf, S5.png

library(ggplot2)
library(cowplot)

source("R/paths.R")

s5_dir   <- supp_dir(5)
data_dir <- file.path(s5_dir, "data")

sumstats <- readRDS(file.path(data_dir, "data_and_sumstats.rds"))
b        <- sumstats$b
causal   <- b != 0

fits <- list(
  "SuSiE"     = readRDS(file.path(data_dir, "fit_susie.rds")),
  "SuSiE-inf" = readRDS(file.path(data_dir, "fit_inf.rds")),
  "SuSiE-ash" = readRDS(file.path(data_dir, "fit_ash.rds"))
)

cs_colors <- c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "gold1",
               "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F")

# ---- A, B: ground truth ----

truth_panel <- function(value, y_label, margins) {
  df <- data.frame(pos = seq_along(value), value = value, causal = causal)
  ggplot(df, aes(x = pos, y = value)) +
    geom_point(data = df[!df$causal, ], color = "black", size = 1.5, alpha = 0.7) +
    geom_point(data = df[df$causal, ], color = "red", size = 3.5, alpha = 0.9) +
    labs(x = "Index", y = y_label) +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", plot.margin = margins)
}

panel_A <- truth_panel(abs(b), "Absolute Effect Size", margin(5, 10, 5, 5))
panel_B <- truth_panel(-log10(sumstats$pvals), expression(-log[10](p)), margin(5, 5, 5, 10))

# ---- C-E: PIPs with credible sets circled ----

pip_panel <- function(fit, title) {
  df <- data.frame(pos = seq_along(fit$pip), pip = fit$pip, causal = causal, cs_color = NA)
  # As in susieR's susie_plot(), sets of 400 or more variants are not circled.
  shown <- fit$sets$cs[lengths(fit$sets$cs) < 400]
  for (i in seq_along(shown)) df$cs_color[shown[[i]]] <- cs_colors[i]
  in_cs <- df[!is.na(df$cs_color), ]

  ggplot(df, aes(x = pos, y = pip)) +
    geom_point(data = df[!df$causal, ], color = "black", size = 1.5, alpha = 0.7) +
    geom_point(data = df[df$causal, ], color = "red", size = 3.5, alpha = 0.9) +
    geom_point(data = in_cs, color = in_cs$cs_color, size = 3.5,
               shape = 21, stroke = 2, fill = NA) +
    labs(x = "Index", y = "PIP", title = title) +
    theme_cowplot(font_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "none")
}

pip_panels <- Map(pip_panel, fits, names(fits))

# ---- F-H: credible sets on the true effects ----

cs_panel <- function(fit) {
  df <- data.frame(pos = seq_along(b), effect = abs(b), cs_id = NA)
  for (i in seq_along(fit$sets$cs)) df$cs_id[fit$sets$cs[[i]]] <- i

  ggplot(df, aes(x = pos, y = effect)) +
    geom_point(color = "black", size = 1.5, alpha = 0.7) +
    geom_point(data = df[!is.na(df$cs_id), ], aes(color = factor(cs_id)), size = 3.5) +
    scale_color_manual(values = cs_colors) +
    labs(x = "Index", y = "Absolute Effect Size") +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none")
}

cs_panels <- lapply(fits, cs_panel)

# ---- Assemble ----

row1 <- plot_grid(panel_A, panel_B, nrow = 1, align = "hv", axis = "tblr",
                  labels = c("A", "B"), label_size = 15, label_fontface = "bold")

# Align the two lower rows so that each method occupies one column.
lower <- align_plots(plotlist = c(pip_panels, cs_panels), align = "v", axis = "lr")
row2 <- plot_grid(plotlist = lower[1:3], nrow = 1,
                  labels = c("C", "D", "E"), label_size = 15, label_fontface = "bold")
row3 <- plot_grid(plotlist = lower[4:6], nrow = 1,
                  labels = c("F", "G", "H"), label_size = 15, label_fontface = "bold")

final_figure <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(0.9, 1, 1))

save_figure(final_figure, s5_dir, "S5", width = 11, height = 8.5)
