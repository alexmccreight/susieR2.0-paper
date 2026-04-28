#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel D — Plotting (Dot Plot with Error Bars)
# =============================================================================
# Creates a dot plot of mean PIP-weighted AlphaGenome CS score for
# RNA-seq (circles) and DNase (triangles) across 7 concordance groups,
# with SE error bars and significance stars.
#
# Input:  panel_D_data.rds (from prepare_panel_D_data.R)
# Output: panel_D.pdf, panel_D.png, panel_D_plot.rds
# =============================================================================

library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)

# =============================================================================
# Paths
# =============================================================================

script_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2/figure_2_panel_D"
data_path  <- file.path(script_dir, "panel_D_data.rds")

# =============================================================================
# Shared aesthetics (match Panel B concordance colors)
# =============================================================================

concordance_colors <- c(
  "Consensus"   = "#888888",
  "SuSiE&ash"   = "#9C27B0",
  "SuSiE&inf"   = "#2E9D8F",
  "ash&inf"     = "#FF9800",
  "SuSiE"       = "#4A90E2",
  "SuSiE-ash"   = "#E53935",
  "SuSiE-inf"   = "#7CB342"
)

modality_shapes <- c("RNA-seq" = 16, "DNase" = 17)  # circle, triangle

# =============================================================================
# Load data
# =============================================================================

cat("Loading panel D data...\n")
panel_data <- readRDS(data_path)
summary_df <- panel_data$summary

# Set modality order: RNA-seq (circle) first, DNase (triangle) second
summary_df$modality <- factor(summary_df$modality, levels = c("RNA-seq", "DNase"))
panel_data$scores$modality <- factor(panel_data$scores$modality, levels = c("RNA-seq", "DNase"))

# Remap group names for consistent method-name capitalization
rename_groups <- c(
  "SuSiE + ASH" = "SuSiE&ash",
  "SuSiE + INF" = "SuSiE&inf",
  "ASH + INF"   = "ash&inf",
  "SuSiE only"  = "SuSiE",
  "ASH only"    = "SuSiE-ash",
  "INF only"    = "SuSiE-inf"
)
levels(summary_df$group) <- ifelse(levels(summary_df$group) %in% names(rename_groups),
                                   rename_groups[levels(summary_df$group)],
                                   levels(summary_df$group))
levels(panel_data$scores$group) <- ifelse(levels(panel_data$scores$group) %in% names(rename_groups),
                                          rename_groups[levels(panel_data$scores$group)],
                                          levels(panel_data$scores$group))

# =============================================================================
# Build dot plot
# =============================================================================

cat("Creating dot plot...\n")

# Dodge width for separating RNA-seq and DNase within each group
dodge_w <- 0.4

# Position for significance stars (above upper error bar)
summary_df$star_y <- summary_df$mean + summary_df$se + 0.012

# Sample size labels
summary_df$n_label <- paste0("n=", summary_df$n)

panel_D <- ggplot(summary_df, aes(x = group, y = mean, color = group,
                                  shape = modality)) +
  # Reference line at null expectation (0.5)
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50",
             linewidth = 0.5) +
  # Error bars (± 1 SE)
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 1.2,
                position = position_dodge(width = dodge_w)) +
  # Points
  geom_point(size = 8, position = position_dodge(width = dodge_w)) +
  # Significance stars
  geom_text(aes(y = star_y, label = ifelse(sig == "ns", "", sig)),
            size = 7, color = "black", fontface = "bold",
            position = position_dodge(width = dodge_w),
            show.legend = FALSE) +
  # Scales
  scale_color_manual(values = concordance_colors, guide = "none") +
  scale_shape_manual(values = modality_shapes,
                     guide = guide_legend(override.aes = list(size = 5, color = "black"))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = NULL,
       y = expression(bold("Mean CS Score")),
       shape = NULL) +
  theme_classic(base_size = 20) +
  theme(
    axis.title.y    = element_text(face = "bold", size = 20),
    axis.text.x     = element_blank(),
    axis.text.y     = element_text(color = "black", size = 18),
    legend.position = "bottom",
    legend.text     = element_text(size = 16),
    legend.margin   = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0),
    plot.margin     = margin(5, 10, 5, 5)
  )

# =============================================================================
# Build summary table: mean diff from 0.5 and % > 0.5
# =============================================================================

cat("Building summary table...\n")

scores <- panel_data$scores
group_order <- levels(summary_df$group)

# Compute metrics per group × modality
table_rows <- lapply(group_order, function(g) {
  rna <- scores$cs_score[scores$group == g & scores$modality == "RNA-seq" &
                          !is.na(scores$cs_score)]
  dna <- scores$cs_score[scores$group == g & scores$modality == "DNase" &
                          !is.na(scores$cs_score)]

  data.frame(
    Group     = g,
    `Mean Diff`  = sprintf("%+.3f", mean(rna) - 0.5),
    `%>0.5`      = sprintf("%.0f%%", 100 * mean(rna > 0.5)),
    `Mean Diff`  = sprintf("%+.3f", if (length(dna) > 0) mean(dna) - 0.5 else NA),
    `%>0.5`      = sprintf("%.0f%%", if (length(dna) > 0) 100 * mean(dna > 0.5) else NA),
    check.names = FALSE, stringsAsFactors = FALSE
  )
})
table_df <- do.call(rbind, table_rows)

# Color each group name in the table
group_colors_vec <- concordance_colors[table_df$Group]

# Build tableGrob with compact styling
tt <- ttheme_minimal(
  base_size = 20,
  core    = list(fg_params = list(hjust = 0.5, x = 0.5, fontsize = 19),
                 padding   = unit(c(4, 14), "pt")),
  colhead = list(fg_params = list(fontface = "bold", fontsize = 20, hjust = 0.5, x = 0.5),
                 padding   = unit(c(4, 8), "pt"))
)

tbl_grob <- tableGrob(table_df, rows = NULL, theme = tt)

# Force table to span full width; give Group column a bit more space
ncols <- ncol(table_df)
col_widths <- c(0.28, rep((1 - 0.28) / (ncols - 1), ncols - 1))
tbl_grob$widths <- unit(col_widths, "npc")

# Color group name cells to match dot plot
for (i in seq_along(group_colors_vec)) {
  row_idx <- i + 1
  tbl_grob$grobs[[which(tbl_grob$layout$t == row_idx &
                         tbl_grob$layout$l == 1 &
                         tbl_grob$layout$name == "core-fg")]]$gp <-
    gpar(fontsize = 19, col = group_colors_vec[i], fontface = "bold")
}

# Add spanning "RNA" and "DNase" header row with lines above column headers
# Insert a new row at the top for the spanning labels
tbl_grob <- gtable::gtable_add_rows(tbl_grob, heights = unit(1.2, "lines"), pos = 0)

# "RNA" label spanning columns 2-3
tbl_grob <- gtable::gtable_add_grob(tbl_grob,
  textGrob("RNA", gp = gpar(fontsize = 20, fontface = "bold")),
  t = 1, l = 2, r = 3)

# "DNase" label spanning columns 4-5
tbl_grob <- gtable::gtable_add_grob(tbl_grob,
  textGrob("DNase", gp = gpar(fontsize = 20, fontface = "bold")),
  t = 1, l = 4, r = 5)

# Line under "RNA" (above its sub-headers)
tbl_grob <- gtable::gtable_add_grob(tbl_grob,
  segmentsGrob(x0 = 0.1, x1 = 0.9, y0 = 0, y1 = 0,
               gp = gpar(lwd = 2)),
  t = 1, l = 2, r = 3, z = Inf)

# Line under "DNase" (above its sub-headers)
tbl_grob <- gtable::gtable_add_grob(tbl_grob,
  segmentsGrob(x0 = 0.1, x1 = 0.9, y0 = 0, y1 = 0,
               gp = gpar(lwd = 2)),
  t = 1, l = 4, r = 5, z = Inf)

# Combine dot plot and table
panel_D_combined <- plot_grid(
  panel_D, tbl_grob,
  ncol = 1,
  rel_heights = c(1.5, 1.2)
)

# =============================================================================
# Save
# =============================================================================

output_pdf <- file.path(script_dir, "panel_D.pdf")
ggsave(output_pdf, panel_D_combined, width = 6, height = 6.5, units = "in",
       bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(script_dir, "panel_D.png")
ggsave(output_png, panel_D_combined, width = 6, height = 6.5, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

# Save plot object for use by create_figure2.R
plot_path <- file.path(script_dir, "panel_D_plot.rds")
saveRDS(panel_D_combined, plot_path)
cat(sprintf("Saved: %s\n", plot_path))

cat("\nDone!\n")
