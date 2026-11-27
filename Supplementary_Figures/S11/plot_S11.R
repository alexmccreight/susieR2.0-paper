#!/usr/bin/env Rscript

# Supplementary Figure S11: mean runtime per method. One bar each for SuSiE,
# SuSiE-ash and SuSiE-inf: the mean wall-clock time of a fit over the 500
# replicates of the runtime simulation, with standard-error bars.
#
# Input:  data/s11_runtime_data.rds  (extract_S11_data.R)
# Output: S11.pdf, S11.png

library(ggplot2)

source("R/paths.R")
source("R/aesthetics.R")

s11_dir <- supp_dir(11)
dat     <- readRDS(file.path(s11_dir, "data", "s11_runtime_data.rds"))

theme_bar <- theme_classic(base_size = 14) +
  theme(
    axis.title      = element_text(face = "bold", size = 18),
    axis.text       = element_text(color = "black", size = 15),
    axis.text.x     = element_text(face = "bold", size = 16),
    legend.position = "none",
    plot.margin     = margin(12, 14, 12, 12)
  )

# Each value is printed above its bar, raised by the largest standard error so
# that it clears the error bar; y_top leaves headroom for it.
y_top <- max(dat$mean_time + dat$se_time) * 1.12

p <- ggplot(dat, aes(x = Method, y = mean_time, fill = Method)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_errorbar(
    aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
    width = 0.18, linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.1f s", mean_time)),
    vjust = -0.6, nudge_y = max(dat$se_time) * 1.1,
    fontface = "bold", size = 5
  ) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, y_top)) +
  labs(x = NULL, y = "Runtime (seconds)") +
  theme_bar

save_figure(p, s11_dir, "S11", width = 7, height = 6)
