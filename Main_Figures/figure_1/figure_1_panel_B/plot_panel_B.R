#!/usr/bin/env Rscript
# =============================================================================
# Panel B (SuSiE-inf): Fold-speedup over the OLD reference at fixed p = 5000
# =============================================================================
# Grouped bar chart of fold speedup over OLD-LD (the ORIGINAL RSS + eigen
# method), across p/n ratios (p = 5,000 fixed):
#   - "Genotype Matrix" = median(time_old_LD) / median(time_new_X)   [susie(X,y), contribution ii]
#   - "LD Matrix"       = median(time_old_LD) / median(time_new_LD)  [susie_rss(R), contribution i]
# Uses ratio-of-medians. Same grouped-bar style as the RSS-lambda panel that
# now lives in Supplementary Figure S6.
#
# Input:  ../data/panel_B/susie_inf_comparison_data/benchmark_susieinf_*_nrep*.rds
# Output: panel_B.{pdf,png}, panel_B_plot.rds
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# =============================================================================
# Paths
# =============================================================================

source("R/paths.R")
source("R/aesthetics.R")

fig1_dir    <- fig_dir(1)
script_dir  <- file.path(fig1_dir, "figure_1_panel_B")
results_dir <- file.path(fig1_dir, "data", "panel_B", "susie_inf_comparison_data")

# ---- Load combined results (exclude per-job _rep files) ----
files <- list.files(results_dir, pattern = "^benchmark_susieinf_.*nrep[0-9]+\\.rds$",
                    full.names = TRUE)
files <- files[!grepl("_rep[0-9]+-[0-9]+\\.rds$", files)]
if (length(files) == 0) {
  stop("No combined result files in ", results_dir,
       " — run the sweep + combine step first.")
}

summary_list <- lapply(files, function(f) {
  d <- readRDS(f)
  reps <- d$replicates[!sapply(d$replicates, is.null)]
  med <- function(field) median(sapply(reps, `[[`, field), na.rm = TRUE)
  data.frame(
    n = d$parameters$n,
    p = d$parameters$p,
    median_old_LD = med("time_old_LD"),
    median_new_X  = med("time_new_X"),
    median_new_LD = med("time_new_LD")
  )
})
df_all <- bind_rows(summary_list)

# ---- Fixed p = 5000: fold speedups over OLD-LD (original RSS) ----
df <- df_all %>%
  filter(p == 5000) %>%
  mutate(
    pn_ratio = as.integer(round(p / n)),
    fold_X   = median_old_LD / median_new_X,    # Genotype Matrix (NEW-X individual, contribution ii)
    fold_LD  = median_old_LD / median_new_LD    # LD Matrix (NEW-LD susie_rss, contribution i)
  ) %>%
  arrange(pn_ratio)

cat(sprintf("  %d conditions (p = 5000, p/n = %s)\n",
            nrow(df), paste(df$pn_ratio, collapse = ", ")))

df_fold <- df %>%
  select(pn_ratio, fold_LD, fold_X) %>%
  pivot_longer(cols = c(fold_LD, fold_X), names_to = "method", values_to = "fold") %>%
  mutate(
    method   = factor(method, levels = c("fold_LD", "fold_X"),
                      labels = c("LD Matrix", "Genotype Matrix")),
    pn_label = factor(pn_ratio, levels = sort(unique(pn_ratio))),
    fold_label = paste0(sprintf("%.1f", fold), "×")
  )

# ---- Plot (same style as the S6 RSS-lambda panel) ----
dodge_w <- 0.7
method_colors <- c("LD Matrix" = "#DF8A2C", "Genotype Matrix" = "#2E6B9E")

p_fold <- ggplot(df_fold, aes(x = pn_label, y = fold, fill = method)) +
  geom_col(position = position_dodge(width = dodge_w), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = method_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Feature-Sample Ratio", y = "Speedup (fold)") +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold", size = 10),
    axis.text  = element_text(color = "black", size = 9),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    plot.margin = margin(5, 10, 5, 5)
  )

# =============================================================================
# Save
# =============================================================================

ggsave(file.path(script_dir, "panel_B.pdf"), p_fold,
       width = 4.5, height = 4, units = "in", bg = "white")
ggsave(file.path(script_dir, "panel_B.png"), p_fold,
       width = 4.5, height = 4, units = "in", dpi = 300, bg = "white")
saveRDS(p_fold, file.path(script_dir, "panel_B_plot.rds"))

cat("Done. Fold speedups (p = 5000):\n")
print(df %>% select(n, pn_ratio, fold_X, fold_LD))
