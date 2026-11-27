#!/usr/bin/env Rscript

# Figure 2, panel D (data): PIP-weighted AlphaGenome scores of the signals in
# each cross-method group, for RNA-seq and DNase.
#
# Restricted to complex regions: gene-tissue pairs in which SuSiE reports at
# least two credible sets. The unit of analysis is the signal. A signal matched
# across methods (Jaccard >= 0.75) has one credible set per member method, with
# nearly identical scores; these are averaged so that the signal counts once.
# Each group and modality is tested against the null expectation of 0.5 with a
# one-sided t-test (mean > 0.5), Benjamini-Hochberg adjusted over all tests.
#
# The groups are read from the input, so the same script prepares the
# three-method version for Supplementary Figure S10 (see run_all.R).
#
# Usage:  prepare_panel_D_data.R [scores.csv] [assignments.rds] [output.rds]
# Input:  data/alphagenome_cs_group_scores__standard_inf.csv
#         data/alphagenome_cs_group_assignments__standard_inf.rds
# Output: figure_2_panel_D/panel_D_data.rds  (list: $scores, $summary)

source("R/paths.R")
source("R/aesthetics.R")

fig2_dir <- fig_dir(2)
data_dir <- file.path(fig2_dir, "data")

scores_file      <- cli_arg(1, file.path(data_dir, "alphagenome_cs_group_scores__standard_inf.csv"))
assignments_file <- cli_arg(2, file.path(data_dir, "alphagenome_cs_group_assignments__standard_inf.rds"))
output_path      <- cli_arg(3, file.path(fig2_dir, "figure_2_panel_D", "panel_D_data.rds"))

d           <- read.csv(scores_file)
assignments <- readRDS(assignments_file)

# ---- Complex regions, RNA-seq and DNase ----

std_assignments <- assignments[assignments$method == "standard", ]
cs_per_pair     <- table(paste(std_assignments$gene_id, std_assignments$tissue, sep = ":"))
complex_pairs   <- names(cs_per_pair[cs_per_pair >= 2])

d$pair <- paste(d$gene_id, d$tissue, sep = ":")
d <- d[d$pair %in% complex_pairs & d$modality %in% c("RNA_SEQ", "DNASE"), ]
d$modality <- ifelse(d$modality == "RNA_SEQ", "RNA-seq", "DNase")

# ---- One observation per signal and modality ----

# Each input row is one credible set; the credible sets of a signal are averaged.
d <- d[!is.na(d$cs_score), ]
d <- aggregate(cs_score ~ signal_id + gene_id + tissue + pair + modality + group,
               data = d, FUN = mean)

# GROUP_DISPLAY (R/aesthetics.R) maps the group labels of the input to display
# labels; its order is the plotting order.
unknown <- setdiff(d$group, names(GROUP_DISPLAY))
if (length(unknown)) stop("Unrecognised group label(s): ", paste(unknown, collapse = ", "))

group_order <- unname(GROUP_DISPLAY[names(GROUP_DISPLAY) %in% d$group])
d$group     <- factor(unname(GROUP_DISPLAY[d$group]), levels = group_order)

# ---- Summary per group and modality ----

modalities <- c("RNA-seq", "DNase")

summary_df <- do.call(rbind, lapply(modalities, function(mod) {
  do.call(rbind, lapply(group_order, function(g) {
    x <- d$cs_score[d$group == g & d$modality == mod]
    n <- length(x)
    if (n < 2) {
      return(data.frame(group = g, modality = mod, n = n,
                        mean = ifelse(n == 1, x, NA),
                        se = NA, pval_raw = NA))
    }
    data.frame(
      group    = g,
      modality = mod,
      n        = n,
      mean     = mean(x),
      se       = sd(x) / sqrt(n),
      pval_raw = t.test(x, mu = 0.5, alternative = "greater")$p.value
    )
  }))
}))

non_na <- !is.na(summary_df$pval_raw)
summary_df$pval_fdr <- NA
summary_df$pval_fdr[non_na] <- p.adjust(summary_df$pval_raw[non_na], method = "BH")

summary_df$sig <- ifelse(is.na(summary_df$pval_fdr), "ns",
                  ifelse(summary_df$pval_fdr < 0.001, "***",
                  ifelse(summary_df$pval_fdr < 0.01,  "**",
                  ifelse(summary_df$pval_fdr < 0.05,  "*", "ns"))))

summary_df$group <- factor(summary_df$group, levels = group_order)

saveRDS(list(scores = d, summary = summary_df), output_path)
cat("Saved:", output_path, "\n")
