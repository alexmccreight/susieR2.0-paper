#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel A — Data Preparation
# =============================================================================
# Loads ROSMAP eQTL fine-mapping results, restricts to gene-tissue pairs common
# across all three methods, filters to 95% coverage credible sets, and outputs a
# combined long-format data frame for plotting.
#
# The gene-tissue universe is still defined by ALL three methods, so that main
# Figure 2 and Supplementary S11 share an identical denominator. Only the
# main-figure methods (MAIN_METHODS: SuSiE, SuSiE-inf) contribute credible sets;
# the three-method version is Supplementary Figure S11.
#
# Input:  {method}_credible_sets.rds, {method}_gene_tissue_summary.rds
# Output: panel_A_data.rds
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

source("R/paths.R")
source("R/aesthetics.R")

fig2_dir   <- fig_dir(2)
data_dir   <- file.path(fig2_dir, "data", "panel_A")
output_dir <- file.path(fig2_dir, "figure_2_panel_A")

# =============================================================================
# Load data
# =============================================================================

cat("Loading data...\n")

summary_files <- c(
  SuSiE       = "susie_standard_gene_tissue_summary.rds",
  `SuSiE-ash` = "susie_ash_gene_tissue_summary.rds",
  `SuSiE-inf` = "susie_inf_gene_tissue_summary.rds"
)

cs_files <- c(
  SuSiE       = "susie_standard_credible_sets.rds",
  `SuSiE-ash` = "susie_ash_credible_sets.rds",
  `SuSiE-inf` = "susie_inf_credible_sets.rds"
)

# All three summaries define the universe; only MAIN_METHODS contribute CSs.
summaries <- lapply(summary_files, function(f) readRDS(file.path(data_dir, f)))
cs_list   <- lapply(cs_files[MAIN_METHODS], function(f) readRDS(file.path(data_dir, f)))

# =============================================================================
# Identify common gene x tissue pairs across all three methods
# =============================================================================

cat("Identifying common gene x tissue pairs...\n")

make_pair <- function(df) paste(df$gene_id, df$tissue, sep = ":")

pairs_per_method <- lapply(summaries, make_pair)
common_pairs     <- Reduce(intersect, pairs_per_method)

cat(sprintf("  Common gene x tissue pairs: %d\n", length(common_pairs)))
cat(sprintf("  Unique genes: %d\n", length(unique(sub(":.*$", "", common_pairs)))))
cat("  Pairs by tissue:\n")
print(table(sub("^.*:", "", common_pairs)))

# =============================================================================
# Filter credible sets to common pairs and 95% coverage
# =============================================================================

cat("\nFiltering credible sets to common pairs (95% coverage)...\n")

filter_cs <- function(cs_df, method_name) {
  cs_df$pair <- make_pair(cs_df)
  cs_filtered <- cs_df[cs_df$pair %in% common_pairs &
                        cs_df$coverage_level == "0.95", ]
  cs_filtered$method <- method_name
  cs_filtered$pair <- NULL
  cs_filtered
}

cs_combined <- do.call(rbind, mapply(filter_cs, cs_list, names(cs_list),
                                     SIMPLIFY = FALSE))
rownames(cs_combined) <- NULL

# =============================================================================
# Derive CS size bins (matching panel specification)
# =============================================================================

cs_combined$cs_size_bin <- cut(
  cs_combined$n_variants,
  breaks = c(0, 1, 2, 5, 10, Inf),
  labels = c("1", "2", "3-5", "6-10", ">10"),
  right  = TRUE
)

# Set method factor levels for consistent ordering
cs_combined$method <- factor(cs_combined$method, levels = MAIN_METHODS)

# =============================================================================
# Verification summary
# =============================================================================

cat("\n--- Verification ---\n")
cat(sprintf("Total rows: %d\n", nrow(cs_combined)))
cat("\nCredible sets per method:\n")
print(table(cs_combined$method))

cat("\nMedian CS size per method:\n")
for (m in levels(cs_combined$method)) {
  cat(sprintf("  %s: %.1f\n", m, median(cs_combined$n_variants[cs_combined$method == m])))
}

cat("\nMedian CS purity (min_abs_corr) per method:\n")
for (m in levels(cs_combined$method)) {
  cat(sprintf("  %s: %.4f\n", m,
              median(cs_combined$min_abs_corr[cs_combined$method == m], na.rm = TRUE)))
}

cat("\nCS size bin distribution:\n")
print(table(cs_combined$method, cs_combined$cs_size_bin))

# =============================================================================
# Save
# =============================================================================

output_path <- file.path(output_dir, "panel_A_data.rds")
saveRDS(cs_combined, output_path)
cat(sprintf("\nSaved: %s\n", output_path))
