#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel B — Data Preparation
# =============================================================================
# Loads the CS group assignments (from Jaccard >= 0.75 concordance analysis),
# computes signal-level concordance counts for the 7 intersection groups,
# and outputs a summary data frame for UpSet-style plotting.
#
# Input:  alphagenome_cs_group_assignments.rds (from alphagenome_cs_group_comparison.R)
# Output: panel_B_data.rds
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

fig2_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2"
input_file <- file.path(fig2_dir, "data", "alphagenome_cs_group_assignments.rds")
output_dir <- file.path(fig2_dir, "figure_2_panel_B")

# =============================================================================
# Load data
# =============================================================================

cat("Loading group assignments...\n")
signals <- readRDS(input_file)
cat(sprintf("  %d CS entries across %d signals\n",
            nrow(signals), length(unique(signals$signal_id))))

# =============================================================================
# Classify signals into 7 concordance groups
# =============================================================================

cat("Classifying signals into concordance groups...\n")

signal_methods <- tapply(signals$method, signals$signal_id, function(m) {
  sort(unique(m))
}, simplify = FALSE)

classify_group <- function(methods) {
  has_std <- "standard" %in% methods
  has_ash <- "ash" %in% methods
  has_inf <- "inf" %in% methods

  if (has_std && has_ash && has_inf) return("Consensus")
  if (has_std && has_ash && !has_inf) return("SuSiE + ASH")
  if (has_std && !has_ash && has_inf) return("SuSiE + INF")
  if (!has_std && has_ash && has_inf) return("ASH + INF")
  if (has_std && !has_ash && !has_inf) return("SuSiE only")
  if (!has_std && has_ash && !has_inf) return("ASH only")
  if (!has_std && !has_ash && has_inf) return("INF only")
  return("Unknown")
}

signal_groups <- sapply(signal_methods, classify_group)

# =============================================================================
# Build signal-level concordance counts
# =============================================================================

group_labels <- c("Consensus", "SuSiE + ASH", "SuSiE + INF", "ASH + INF",
                  "SuSiE only", "ASH only", "INF only")

signal_counts <- table(factor(signal_groups, levels = group_labels))

concordance_df <- data.frame(
  group       = factor(group_labels, levels = group_labels),
  n_signals   = as.integer(signal_counts),
  stringsAsFactors = FALSE
)

# Add method membership indicators (for the UpSet dot matrix)
concordance_df$has_susie <- concordance_df$group %in%
  c("Consensus", "SuSiE + ASH", "SuSiE + INF", "SuSiE only")
concordance_df$has_ash   <- concordance_df$group %in%
  c("Consensus", "SuSiE + ASH", "ASH + INF", "ASH only")
concordance_df$has_inf   <- concordance_df$group %in%
  c("Consensus", "SuSiE + INF", "ASH + INF", "INF only")

# Compute total CSs per method
method_totals <- data.frame(
  method = c("SuSiE", "SuSiE-ash", "SuSiE-inf"),
  total  = c(
    sum(signals$method == "standard"),
    sum(signals$method == "ash"),
    sum(signals$method == "inf")
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# Verification
# =============================================================================

cat("\n--- Signal counts per concordance group ---\n")
for (i in seq_len(nrow(concordance_df))) {
  cat(sprintf("  %-15s %5d signals\n",
              concordance_df$group[i], concordance_df$n_signals[i]))
}
cat(sprintf("  %-15s %5d\n", "TOTAL", sum(concordance_df$n_signals)))

cat("\n--- Total CSs per method ---\n")
for (i in seq_len(nrow(method_totals))) {
  cat(sprintf("  %-10s %5d CSs\n", method_totals$method[i], method_totals$total[i]))
}

# =============================================================================
# Save
# =============================================================================

output <- list(
  concordance   = concordance_df,
  method_totals = method_totals
)

output_path <- file.path(output_dir, "panel_B_data.rds")
saveRDS(output, output_path)
cat(sprintf("\nSaved: %s\n", output_path))
