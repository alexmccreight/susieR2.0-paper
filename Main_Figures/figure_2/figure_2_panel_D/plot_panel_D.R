#!/usr/bin/env Rscript

# Figure 2, panel D: mean PIP-weighted AlphaGenome score (+/- 1 SE) of the
# signals in each cross-method group, for RNA-seq (circles) and DNase
# (triangles). Stars mark means above 0.5 (one-sided t-test, BH-adjusted:
# * < 0.05, ** < 0.01, *** < 0.001). The table underneath gives, per group and
# modality, the mean difference from 0.5 and the share of signals above 0.5;
# its coloured group names also serve as the colour key of panel C.
#
# The number of groups comes from the input: three for Figure 2, seven for the
# three-method version in Supplementary Figure S10 (see run_all.R).
#
# Usage:  plot_panel_D.R [panel_D_data.rds] [output directory]
# Input:  panel_D_data.rds  (prepare_panel_D_data.R)
# Output: panel_D.pdf, panel_D.png, panel_D_plot.rds

library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)

source("R/paths.R")
source("R/aesthetics.R")

script_dir <- file.path(fig_dir(2), "figure_2_panel_D")
data_path  <- cli_arg(1, file.path(script_dir, "panel_D_data.rds"))
out_dir    <- cli_arg(2, script_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

panel_data <- readRDS(data_path)
scores     <- panel_data$scores
summary_df <- panel_data$summary
summary_df$modality <- factor(summary_df$modality, levels = c("RNA-seq", "DNase"))

# ---- Dot plot ----

# Modality on the x-axis, with the groups dodged within each modality.
dodge_w <- 0.75

# Stars sit just above the upper error bar.
summary_df$star_y <- summary_df$mean + summary_df$se + 0.012

panel_D <- ggplot(summary_df, aes(x = modality, y = mean, color = group,
                                  group = group, shape = modality)) +
  # Null expectation
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50",
             linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 1.2,
                position = position_dodge(width = dodge_w)) +
  geom_point(size = 8, position = position_dodge(width = dodge_w)) +
  geom_text(aes(y = star_y, label = ifelse(sig == "ns", "", sig)),
            size = 7, color = "black", fontface = "bold",
            position = position_dodge(width = dodge_w)) +
  scale_color_manual(values = concordance_colors) +
  scale_shape_manual(values = c("RNA-seq" = 16, "DNase" = 17)) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = NULL, y = expression(bold("Mean CS Score"))) +
  theme_classic(base_size = 20) +
  theme(
    axis.title.y    = element_text(face = "bold", size = 20),
    axis.text.x     = element_text(color = "black", size = 18, face = "bold"),
    axis.text.y     = element_text(color = "black", size = 18),
    legend.position = "none",
    plot.margin     = margin(5, 10, 5, 5)
  )

# ---- Table ----

group_order <- levels(summary_df$group)

# Mean difference from 0.5 and share of signals above 0.5, as text; "NA" for a
# group with no signal in that modality.
score_cells <- function(g, modality) {
  x <- scores$cs_score[scores$group == g & scores$modality == modality]
  if (length(x) == 0) x <- NA
  c(sprintf("%+.3f", mean(x) - 0.5), sprintf("%.0f%%", 100 * mean(x > 0.5)))
}

table_df <- do.call(rbind, lapply(group_order, function(g) {
  data.frame(g, t(score_cells(g, "RNA-seq")), t(score_cells(g, "DNase")))
}))
colnames(table_df) <- c("Group", "Mean Diff", "%>0.5", "Mean Diff", "%>0.5")

tt <- ttheme_minimal(
  base_size = 20,
  core    = list(fg_params = list(hjust = 0.5, x = 0.5, fontsize = 16),
                 padding   = unit(c(4, 8), "pt")),
  colhead = list(fg_params = list(fontface = "bold", fontsize = 17, hjust = 0.5, x = 0.5),
                 padding   = unit(c(4, 8), "pt"))
)

tbl_grob <- tableGrob(table_df, rows = NULL, theme = tt)

# Column widths follow the widest string of each column (header at 17 pt, body
# at 16 pt, both measured in bold) plus 0.16 in of padding, so that no cell is
# clipped whatever the number of groups.
grDevices::pdf(NULL)
graphics::par(ps = 12)
nat <- vapply(seq_along(table_df), function(j) {
  max(graphics::strwidth(colnames(table_df)[j], units = "inches", cex = 17 / 12, font = 2),
      graphics::strwidth(table_df[[j]], units = "inches", cex = 16 / 12, font = 2))
}, numeric(1))
invisible(grDevices::dev.off())

nat <- nat + 0.16
tbl_grob$widths <- unit(nat / sum(nat), "npc")

# Group names in the colours of the dot plot.
group_colors_vec <- concordance_colors[table_df$Group]
for (i in seq_along(group_colors_vec)) {
  row_idx <- i + 1
  tbl_grob$grobs[[which(tbl_grob$layout$t == row_idx &
                         tbl_grob$layout$l == 1 &
                         tbl_grob$layout$name == "core-fg")]]$gp <-
    gpar(fontsize = 16, col = group_colors_vec[i], fontface = "bold")
}

# A top row with "RNA" and "DNase" spanning their two columns, each underlined.
tbl_grob <- gtable::gtable_add_rows(tbl_grob, heights = unit(1.2, "lines"), pos = 0)
for (span in list(list("RNA", 2, 3), list("DNase", 4, 5))) {
  tbl_grob <- gtable::gtable_add_grob(tbl_grob,
    textGrob(span[[1]], gp = gpar(fontsize = 17, fontface = "bold")),
    t = 1, l = span[[2]], r = span[[3]])
  tbl_grob <- gtable::gtable_add_grob(tbl_grob,
    segmentsGrob(x0 = 0.1, x1 = 0.9, y0 = 0, y1 = 0,
                 gp = gpar(lwd = 2)),
    t = 1, l = span[[2]], r = span[[3]], z = Inf)
}

# ---- Combine ----

# The table's share of the height follows its number of rows: 1.2 relative
# units for the 9 rows of the seven-group table, against 1.5 for the dot plot.
# The figure height follows too, so that the dot plot keeps the size it has in
# the seven-group version, which is 6.5 in tall.
DOT_REL_H     <- 1.5
TBL_ROW_H     <- 1.2 / 9
n_tbl_rows    <- nrow(table_df) + 2            # + column headers + RNA/DNase row
tbl_rel_h     <- TBL_ROW_H * n_tbl_rows

panel_D_combined <- plot_grid(
  panel_D, tbl_grob,
  ncol = 1,
  rel_heights = c(DOT_REL_H, tbl_rel_h)
)

DOT_IN     <- 6.5 * (DOT_REL_H / (DOT_REL_H + 1.2))
fig_height <- DOT_IN * (DOT_REL_H + tbl_rel_h) / DOT_REL_H

save_figure(panel_D_combined, out_dir, "panel_D", width = 6, height = fig_height)

plot_path <- file.path(out_dir, "panel_D_plot.rds")
saveRDS(panel_D_combined, plot_path)
cat("Saved:", plot_path, "\n")
