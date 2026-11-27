#!/usr/bin/env Rscript

# Figure 2, panels C and D, and Supplementary Figure S10: cross-method
# concordance groups of the 95% credible sets, with their AlphaGenome scores.
# Step 3 of the AlphaGenome scoring pipeline.
#
# Within a gene-tissue pair, two credible sets from different methods match
# when the Jaccard index of their variant sets is >= 0.75. A signal is a
# connected component of the match graph (union-find) and is labelled by the
# methods it contains: "Consensus" (all of them), "<Method>-specific" (one) or
# the methods joined by "+". Each AlphaGenome CS score row then takes the
# signal and group of its credible set.
#
# Only gene-tissue pairs listed in the gene_tissue_summary of all three methods
# are analysed, whichever methods are compared.
#
# Input:  data/panel_A/susie_{standard,ash,inf}_gene_tissue_summary.rds,
#         data/panel_A/susie_<method>_credible_sets.rds  (extract_bvsr_results.R)
#         data/alphagenome_cs_scores.csv                 (alphagenome_cs_scoring.py)
# Output: data/alphagenome_cs_group_assignments<suffix>.rds  one row per CS
#         data/alphagenome_cs_group_scores<suffix>.csv       one row per CS and modality
#         data/alphagenome_cs_group_summary<suffix>.txt      text summary of the scores
#         <suffix> is empty when all three methods are compared, otherwise "__"
#         and the methods joined by "_" (__standard_inf by default).
# Usage:  Rscript .../alphagenome_cs_group_comparison.R [methods]
#         methods: comma-separated, from standard, ash, inf. The default,
#         standard,inf, is main Figure 2; standard,ash,inf is Supplementary
#         Figure S10.

source("R/paths.R")
source("R/aesthetics.R")

data_dir   <- file.path(fig_dir(2), "data", "panel_A")
output_dir <- file.path(fig_dir(2), "data")

# ---- Methods ----

# PAIR_METHODS fixes the universe: a gene-tissue pair is analysed only if it is
# in the gene_tissue_summary of every one of these methods. Holding it at all
# three gives main Figure 2 and Supplementary Figure S10 the same 7,941 pairs.
PAIR_METHODS <- c("standard", "ash", "inf")

# CS_METHODS are the methods whose credible sets are matched. A method left out
# is removed before the union-find, not after it: otherwise one of its credible
# sets could still bridge two sets of the remaining methods into one signal and
# overstate their agreement.
CS_METHODS <- trimws(strsplit(cli_arg(1, "standard,inf"), ",")[[1]])
stopifnot(all(CS_METHODS %in% PAIR_METHODS), length(CS_METHODS) >= 2)

OUT_SUFFIX <- if (setequal(CS_METHODS, PAIR_METHODS)) "" else
  paste0("__", paste(CS_METHODS, collapse = "_"))

GROUP_LABEL <- c(standard = "Standard", ash = "ASH", inf = "Inf")
JACCARD_THRESHOLD <- 0.75

# ---- Load the 95% credible sets ----

common_pairs <- Reduce(intersect, lapply(PAIR_METHODS, function(m) {
  s <- readRDS(file.path(data_dir, paste0("susie_", m, "_gene_tissue_summary.rds")))
  paste(s$gene_id, s$tissue, sep = ":")
}))

load_cs95 <- function(method) {
  cs <- readRDS(file.path(data_dir, paste0("susie_", method, "_credible_sets.rds")))
  cs$pair <- paste(cs$gene_id, cs$tissue, sep = ":")
  cs <- cs[cs$pair %in% common_pairs & cs$coverage_level == "0.95", ]
  cs$method <- method
  cs$variant_set <- strsplit(cs$variants, ",\\s*")
  cs[, c("gene_id", "tissue", "cs_id", "method", "pair",
         "n_variants", "top_variant", "top_pip", "variant_set")]
}

cs_list   <- sapply(CS_METHODS, load_cs95, simplify = FALSE)
cs_counts <- vapply(cs_list, nrow, integer(1))
all_cs    <- do.call(rbind, cs_list)

# ---- Match credible sets across methods ----

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

find_root <- function(parent, i) {
  while (parent[i] != i) i <- parent[i]
  i
}

# Signals are numbered pair by pair in order of first appearance and, within a
# pair, by the position of their root credible set.
all_cs$signal_id <- NA_integer_
n_signals <- 0L

for (p in unique(all_cs$pair)) {
  rows   <- which(all_cs$pair == p)
  n      <- length(rows)
  method <- all_cs$method[rows]
  sets   <- all_cs$variant_set[rows]

  # Union-find over the credible sets of this pair.
  parent <- seq_len(n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (method[i] != method[j] &&
          jaccard(sets[[i]], sets[[j]]) >= JACCARD_THRESHOLD) {
        parent[find_root(parent, i)] <- find_root(parent, j)
      }
    }
  }

  roots     <- vapply(seq_len(n), function(i) find_root(parent, i), integer(1))
  component <- match(roots, sort(unique(roots)))
  all_cs$signal_id[rows] <- n_signals + component
  n_signals <- n_signals + max(component)
}

signals <- all_cs[order(all_cs$signal_id),
                  c("signal_id", "gene_id", "tissue", "method", "cs_id",
                    "top_variant", "top_pip", "n_variants")]
rownames(signals) <- NULL

# ---- Label each signal by the methods it contains ----

classify_group <- function(methods) {
  present <- CS_METHODS[CS_METHODS %in% methods]   # in CS_METHODS order
  if (length(present) == length(CS_METHODS)) return("Consensus")
  if (length(present) == 1) return(paste0(GROUP_LABEL[[present]], "-specific"))
  paste(GROUP_LABEL[present], collapse = "+")
}

signals$group <- ave(signals$method, signals$signal_id, FUN = classify_group)

# ---- Attach the AlphaGenome scores ----

cs_scores <- read.csv(file.path(output_dir, "alphagenome_cs_scores.csv"))
annotated <- merge(
  cs_scores,
  signals[, c("gene_id", "tissue", "method", "cs_id", "signal_id", "group")],
  by = c("gene_id", "tissue", "method", "cs_id")
)

# Score rows of methods outside CS_METHODS drop out of the merge by design. Any
# other unmatched row means the scores and the credible sets are out of step.
n_unmatched <- sum(cs_scores$method %in% CS_METHODS) - nrow(annotated)
if (n_unmatched > 0) {
  warning(n_unmatched, " CS score rows could not be matched to a credible set")
}

# ---- Text summary of the scores by group ----

valid <- annotated[!is.na(annotated$cs_score), ]

GROUP_ORDER <- names(GROUP_DISPLAY)
MODALITIES  <- c("RNA_SEQ", "DNASE", "CHIP_HISTONE", "CHIP_TF")

# Percentage of scores above the null expectation of 0.5.
pct_above_null <- function(x) 100 * sum(x > 0.5) / length(x)

lines <- c(
  strrep("=", 80),
  "AlphaGenome CS Score Comparison by Method-Agreement Group",
  "CS_Score = sum(PIP * |quantile_score|) / sum(PIP); Null = 0.5",
  sprintf("Jaccard threshold for matching: %.2f", JACCARD_THRESHOLD),
  strrep("=", 80),
  "\n--- Overall per Group ---",
  sprintf("%-20s %6s %8s %8s %8s %7s",
          "Group", "N_CS", "Mean", "Median", "Std", "%>0.5"),
  strrep("-", 65)
)
for (g in GROUP_ORDER) {
  x <- valid$cs_score[valid$group == g]
  if (length(x) == 0) next
  lines <- c(lines, sprintf("%-20s %6d %8.4f %8.4f %8.4f %6.1f%%",
                            g, length(x), mean(x), median(x), sd(x), pct_above_null(x)))
}

lines <- c(
  lines,
  "\n\n--- Per Group x Modality ---",
  sprintf("%-20s %-15s %6s %8s %8s %7s",
          "Group", "Modality", "N_CS", "Mean", "Median", "%>0.5"),
  strrep("-", 75)
)
for (g in GROUP_ORDER) {
  for (mod in MODALITIES) {
    x <- valid$cs_score[valid$group == g & valid$modality == mod]
    if (length(x) == 0) next
    lines <- c(lines, sprintf("%-20s %-15s %6d %8.4f %8.4f %6.1f%%",
                              g, mod, length(x), mean(x), median(x), pct_above_null(x)))
  }
}

lines <- c(
  lines,
  "\n\n--- Wilcoxon Rank-Sum: Consensus vs Other Groups ---",
  sprintf("%-20s %-15s %8s %12s %10s",
          "Group", "Modality", "N_group", "W-stat", "P-value"),
  strrep("-", 70)
)
for (g in setdiff(GROUP_ORDER, "Consensus")) {
  for (mod in MODALITIES) {
    consensus <- valid$cs_score[valid$group == "Consensus" & valid$modality == mod]
    x         <- valid$cs_score[valid$group == g & valid$modality == mod]
    if (length(consensus) < 5 || length(x) < 5) next

    wt  <- wilcox.test(consensus, x)
    sig <- if (wt$p.value < 0.001) "***" else if (wt$p.value < 0.01) "**" else
      if (wt$p.value < 0.05) "*" else ""
    lines <- c(lines, sprintf("%-20s %-15s %8d %12.1f %10.2e %s",
                              g, mod, length(x), wt$statistic, wt$p.value, sig))
  }
}

lines <- c(lines, "\n\n--- Per Tissue x Group ---")
for (tis in c("AC", "DLPFC", "PCC")) {
  if (!any(valid$tissue == tis)) next
  lines <- c(
    lines,
    sprintf("\n  %s:", tis),
    sprintf("  %-20s %6s %8s %8s %7s", "Group", "N_CS", "Mean", "Median", "%>0.5")
  )
  for (g in GROUP_ORDER) {
    x <- valid$cs_score[valid$tissue == tis & valid$group == g]
    if (length(x) == 0) next
    lines <- c(lines, sprintf("  %-20s %6d %8.4f %8.4f %6.1f%%",
                              g, length(x), mean(x), median(x), pct_above_null(x)))
  }
}

# Every credible set belongs to exactly one signal, so the counts always agree.
lines <- c(lines, "\n\n--- Verification ---")
for (m in CS_METHODS) {
  n_assigned <- sum(signals$method == m)
  lines <- c(lines, sprintf("  %s: %d assigned, %d expected — %s",
                            m, n_assigned, cs_counts[[m]],
                            if (n_assigned == cs_counts[[m]]) "OK" else "MISMATCH"))
}

# ---- Save ----

out_path <- function(name, ext) file.path(output_dir, paste0(name, OUT_SUFFIX, ext))
assignments_file <- out_path("alphagenome_cs_group_assignments", ".rds")
scores_file      <- out_path("alphagenome_cs_group_scores", ".csv")
summary_file     <- out_path("alphagenome_cs_group_summary", ".txt")

saveRDS(signals, assignments_file)
write.csv(annotated, scores_file, row.names = FALSE)
writeLines(lines, summary_file)
cat("Saved:", assignments_file, "\n")
cat("Saved:", scores_file, "\n")
cat("Saved:", summary_file, "\n")
