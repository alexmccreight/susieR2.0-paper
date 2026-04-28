#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel D — Data Preparation
# =============================================================================
# Loads AlphaGenome CS scores annotated with concordance group labels,
# restricts to "complex regions" (gene-tissue pairs where standard SuSiE
# identifies >= 2 credible sets), filters to RNA-seq and DNase modalities,
# renames groups for display, computes per-group-per-modality summary
# statistics and one-sided t-tests (mean > 0.5), and saves for plotting.
#
# Input:  alphagenome_cs_group_scores.csv
#         alphagenome_cs_group_assignments.rds
# Output: panel_D_data.rds
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

fig2_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2"
script_dir <- file.path(fig2_dir, "figure_2_panel_D")
data_dir   <- file.path(fig2_dir, "data")
scores_file      <- file.path(data_dir, "alphagenome_cs_group_scores.csv")
assignments_file <- file.path(data_dir, "alphagenome_cs_group_assignments.rds")

# =============================================================================
# Load data
# =============================================================================

cat("Loading AlphaGenome CS group scores...\n")
d <- read.csv(scores_file, stringsAsFactors = FALSE)
cat(sprintf("  %d rows, %d modalities, %d groups\n",
            nrow(d), length(unique(d$modality)), length(unique(d$group))))

cat("Loading group assignments...\n")
assignments <- readRDS(assignments_file)

# =============================================================================
# Identify complex regions: gene-tissue pairs with >= 2 standard SuSiE CSs
# =============================================================================

cat("\nIdentifying complex regions (standard SuSiE >= 2 CSs)...\n")

std_assignments <- assignments[assignments$method == "standard", ]
std_assignments$pair <- paste(std_assignments$gene_id, std_assignments$tissue, sep = ":")
cs_per_pair <- table(std_assignments$pair)
complex_pairs <- names(cs_per_pair[cs_per_pair >= 2])

cat(sprintf("  Total gene-tissue pairs with standard SuSiE CSs: %d\n",
            length(cs_per_pair)))
cat(sprintf("  Complex regions (>= 2 CSs): %d\n", length(complex_pairs)))

# Filter scores to complex regions
d$pair <- paste(d$gene_id, d$tissue, sep = ":")
d <- d[d$pair %in% complex_pairs, ]
cat(sprintf("  Score rows after filtering to complex regions: %d\n", nrow(d)))

# =============================================================================
# Filter to RNA-seq and DNase modalities
# =============================================================================

d <- d[d$modality %in% c("RNA_SEQ", "DNASE"), ]
cat(sprintf("  After filtering to RNA_SEQ + DNASE: %d rows\n", nrow(d)))

# Rename modalities for display
d$modality <- ifelse(d$modality == "RNA_SEQ", "RNA-seq", "DNase")

# =============================================================================
# Rename groups for display (match Panel B labels)
# =============================================================================

group_labels <- c(
  "Consensus"          = "Consensus",
  "Standard+ASH"       = "SuSiE + ASH",
  "Standard+Inf"       = "SuSiE + INF",
  "ASH+Inf"            = "ASH + INF",
  "Standard-specific"  = "SuSiE only",
  "ASH-specific"       = "ASH only",
  "Inf-specific"       = "INF only"
)

d$group <- group_labels[d$group]

# Order groups to match Panel B (descending by n_signals)
group_order <- c("Consensus", "SuSiE + ASH", "SuSiE + INF", "ASH + INF",
                 "SuSiE only", "ASH only", "INF only")
d$group <- factor(d$group, levels = group_order)

# =============================================================================
# Compute per-group-per-modality summary statistics
# =============================================================================

cat("\nComputing per-group-per-modality statistics...\n")

modalities <- c("RNA-seq", "DNase")

summary_df <- do.call(rbind, lapply(modalities, function(mod) {
  do.call(rbind, lapply(group_order, function(g) {
    x <- d$cs_score[d$group == g & d$modality == mod & !is.na(d$cs_score)]
    n <- length(x)
    if (n < 2) {
      return(data.frame(group = g, modality = mod, n = n,
                        mean = ifelse(n == 1, x, NA),
                        se = NA, pval_raw = NA,
                        stringsAsFactors = FALSE))
    }
    mn <- mean(x)
    se <- sd(x) / sqrt(n)

    # One-sided t-test: H1: mu > 0.5 (null expectation)
    tt <- t.test(x, mu = 0.5, alternative = "greater")

    data.frame(
      group    = g,
      modality = mod,
      n        = n,
      mean     = mn,
      se       = se,
      pval_raw = tt$p.value,
      stringsAsFactors = FALSE
    )
  }))
}))

# FDR correction across all tests (Benjamini-Hochberg), excluding NAs
non_na <- !is.na(summary_df$pval_raw)
summary_df$pval_fdr <- NA
summary_df$pval_fdr[non_na] <- p.adjust(summary_df$pval_raw[non_na], method = "BH")

# Significance stars
summary_df$sig <- ifelse(is.na(summary_df$pval_fdr), "ns",
                  ifelse(summary_df$pval_fdr < 0.001, "***",
                  ifelse(summary_df$pval_fdr < 0.01,  "**",
                  ifelse(summary_df$pval_fdr < 0.05,  "*", "ns"))))

summary_df$group <- factor(summary_df$group, levels = group_order)

# =============================================================================
# Print summary
# =============================================================================

for (mod in modalities) {
  cat(sprintf("\n--- %s CS scores by concordance group (complex regions only) ---\n", mod))
  cat(sprintf("  %-15s %5s %8s %8s %10s %10s %4s\n",
              "Group", "N", "Mean", "SE", "p (raw)", "p (FDR)", "Sig"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  sub <- summary_df[summary_df$modality == mod, ]
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-15s %5d %8.4f %8.4f %10.2e %10.2e %4s\n",
                sub$group[i], sub$n[i],
                ifelse(is.na(sub$mean[i]), 0, sub$mean[i]),
                ifelse(is.na(sub$se[i]), 0, sub$se[i]),
                ifelse(is.na(sub$pval_raw[i]), 1, sub$pval_raw[i]),
                ifelse(is.na(sub$pval_fdr[i]), 1, sub$pval_fdr[i]),
                sub$sig[i]))
  }
}

# =============================================================================
# Save
# =============================================================================

output <- list(
  scores  = d,
  summary = summary_df
)

output_path <- file.path(script_dir, "panel_D_data.rds")
saveRDS(output, output_path)
cat(sprintf("\nSaved: %s\n", output_path))
