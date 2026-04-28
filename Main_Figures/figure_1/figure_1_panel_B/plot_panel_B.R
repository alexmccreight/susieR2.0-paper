#!/usr/bin/env Rscript
# =============================================================================
# Fold-Change Speed Benchmark: susieR 1.0 vs susieR 2.0
# =============================================================================
# Grouped bar chart showing fold speedup over susieR 1.0 (R)
# for susieR 2.0 (R) and susieR 2.0 (X), across p/n ratios (p = 5,000 fixed).
#
# Two stories in one plot:
#   - Algorithm improvement: consistent ~7-8x (orange bars)
#   - X-input advantage:     grows with p/n  (gap between teal and orange)
#
# Input:  rss_lambda_comparison_data/*.rds
# Output: panel_B.{pdf,png}, panel_B_plot.rds
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# Paths
# =============================================================================

fig1_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_1"
script_dir <- file.path(fig1_dir, "figure_1_panel_B")
data_dir   <- file.path(fig1_dir, "data", "panel_B", "rss_lambda_comparison_data")

# =============================================================================
# Load and summarize
# =============================================================================

cat("Loading benchmark results...\n")

files <- list.files(data_dir, pattern = "[.]rds$", full.names = TRUE)

summary_list <- lapply(files, function(f) {
  d <- readRDS(f)
  data.frame(
    n = d$parameters$n,
    p = d$parameters$p,
    median_1.0_R = median(sapply(d$replicates, `[[`, "time_1.0_R")),
    median_2.0_R = median(sapply(d$replicates, `[[`, "time_2.0_R")),
    median_2.0_X = median(sapply(d$replicates, `[[`, "time_2.0_X"))
  )
})

df_all <- bind_rows(summary_list)

# =============================================================================
# Fixed p = 5,000 — compute fold speedups over susieR 1.0
# =============================================================================

df <- df_all %>%
  filter(p == 5000) %>%
  mutate(
    pn_ratio = as.integer(p / n),
    fold_R   = median_1.0_R / median_2.0_R,
    fold_X   = median_1.0_R / median_2.0_X
  ) %>%
  arrange(pn_ratio)

cat(sprintf("  %d conditions (p = 5000, p/n = %s)\n",
            nrow(df), paste(df$pn_ratio, collapse = ", ")))

# Pivot for ggplot
df_fold <- df %>%
  select(pn_ratio, fold_R, fold_X) %>%
  pivot_longer(
    cols      = c(fold_R, fold_X),
    names_to  = "method",
    values_to = "fold"
  ) %>%
  mutate(
    method = factor(method,
      levels = c("fold_R", "fold_X"),
      labels = c("LD Matrix", "Genotype Matrix")
    ),
    pn_label   = factor(pn_ratio, levels = sort(unique(pn_ratio))),
    fold_label = paste0(sprintf("%.1f", fold), "\u00d7")
  )

# =============================================================================
# Plot
# =============================================================================

cat("Creating fold-speedup plot...\n")

dodge_w <- 0.7

method_colors <- c(
  "LD Matrix"       = "#DF8A2C",
  "Genotype Matrix" = "#2E6B9E"
)

p_fold <- ggplot(df_fold, aes(x = pn_label, y = fold, fill = method)) +
  geom_col(position = position_dodge(width = dodge_w), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red",
             linewidth = 0.8) +
  scale_fill_manual(values = method_colors, name = NULL) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    breaks = seq(0, 18, by = 2)
  ) +
  labs(
    x = "Features / Samples",
    y = "Speedup (fold)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title   = element_text(face = "bold", size = 10),
    axis.text    = element_text(color = "black", size = 9),
    legend.position  = c(0.5, 0.97),
    legend.justification = c(0.5, 1),
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    plot.margin      = margin(5, 10, 5, 5)
  )

# =============================================================================
# Save
# =============================================================================

cat("Saving outputs...\n")

ggsave(file.path(script_dir, "panel_B.pdf"), p_fold,
       width = 4.5, height = 4, units = "in", bg = "white")
ggsave(file.path(script_dir, "panel_B.png"), p_fold,
       width = 4.5, height = 4, units = "in", dpi = 300, bg = "white")
saveRDS(p_fold, file.path(script_dir, "panel_B_plot.rds"))

cat("Done!\n")
