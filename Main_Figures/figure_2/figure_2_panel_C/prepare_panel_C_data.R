#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel C — Data Preparation
# =============================================================================
# Loads cross-validated TWAS prediction metrics, filters to the 6 methods
# of interest (excluding SuSiE-ash and Mr.ASH), renames methods for
# display, and saves a clean data frame for plotting.
#
# Input:  twas_cv_all_metrics_combined.rds
# Output: panel_C_data.rds
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

fig2_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2"
script_dir <- file.path(fig2_dir, "figure_2_panel_C")
input_file <- file.path(fig2_dir, "data", "twas_cv_all_metrics_combined.rds")
output_file <- file.path(script_dir, "panel_C_data.rds")

# =============================================================================
# Load data
# =============================================================================

cat("Loading TWAS CV metrics...\n")
d <- readRDS(input_file)
cat(sprintf("  %d rows, %d methods, %d tissues\n",
            nrow(d), length(unique(d$method)), length(unique(d$tissue))))

# =============================================================================
# Filter to 6 methods (exclude susie_ash, susie_inf)
# =============================================================================

keep_methods <- c("susie", "susie_inf", "enet", "lasso", "bayes_r", "bayes_l")
d <- d[d$method %in% keep_methods, ]
cat(sprintf("  After filtering to 6 methods: %d rows\n", nrow(d)))

# =============================================================================
# Rename methods for display
# =============================================================================

method_labels <- c(
  "susie"     = "SuSiE",
  "susie_inf" = "SuSiE-inf",
  "enet"      = "Elastic Net",
  "lasso"     = "LASSO",
  "bayes_r"   = "BayesR",
  "bayes_l"   = "BayesL"
)

d$method <- method_labels[d$method]

# Order by descending mean R²
mean_rsq <- tapply(d$rsq, d$method, mean, na.rm = TRUE)
method_order <- names(sort(mean_rsq, decreasing = TRUE))
d$method <- factor(d$method, levels = method_order)

# =============================================================================
# Summary
# =============================================================================

cat("\n--- Mean R² by method (descending) ---\n")
for (m in method_order) {
  mn <- mean(d$rsq[d$method == m], na.rm = TRUE)
  n  <- sum(!is.na(d$rsq[d$method == m]))
  cat(sprintf("  %-16s  mean R² = %.6f  (n = %d)\n", m, mn, n))
}

cat(sprintf("\n  Total observations: %d\n", nrow(d)))
cat(sprintf("  NAs in rsq: %d\n", sum(is.na(d$rsq))))

# =============================================================================
# Save
# =============================================================================

saveRDS(d, output_file)
cat(sprintf("\nSaved: %s\n", output_file))
