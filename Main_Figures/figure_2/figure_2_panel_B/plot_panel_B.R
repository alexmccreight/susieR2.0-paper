#!/usr/bin/env Rscript

# Figure 2, panel B: mean cross-validated R2 (+/- 1 SE) of six prediction
# methods, over all gene-tissue pairs of the three ROSMAP brain tissues (DLPFC,
# AC, PCC).
#
# Input:  panel_B_data.rds  (prepare_panel_B_data.R)
# Output: panel_B.pdf, panel_B.png, panel_B_plot.rds

library(ggplot2)

source("R/paths.R")
source("R/aesthetics.R")

script_dir <- file.path(fig_dir(2), "figure_2_panel_B")
d          <- readRDS(file.path(script_dir, "panel_B_data.rds"))

# Mean and standard error per method; aggregate() drops the missing R2 values.
mean_df <- aggregate(rsq ~ method, data = d, FUN = function(x) {
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
})
mean_df <- data.frame(
  method   = mean_df$method,
  mean_rsq = mean_df$rsq[, "mean"],
  se       = mean_df$rsq[, "se"]
)

panel_B <- ggplot(mean_df, aes(x = method, y = mean_rsq, color = method)) +
  geom_point(size = 8) +
  geom_errorbar(aes(ymin = mean_rsq - se, ymax = mean_rsq + se),
                width = 0.2, linewidth = 1.2) +
  scale_color_manual(values = c(method_colors[MAIN_METHODS], prediction_colors)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = NULL, y = expression(bold("Mean CV" ~ R^2))) +
  theme_classic(base_size = 20) +
  theme(
    axis.title.y  = element_text(face = "bold", size = 20),
    axis.text.x   = element_text(color = "black", size = 17, angle = 45,
                                 hjust = 1, face = "bold"),
    axis.text.y   = element_text(color = "black", size = 18),
    legend.position = "none",
    plot.margin   = margin(10, 10, 5, 5)
  )

save_figure(panel_B, script_dir, "panel_B", width = 6, height = 5)

plot_path <- file.path(script_dir, "panel_B_plot.rds")
saveRDS(panel_B, plot_path)
cat("Saved:", plot_path, "\n")
