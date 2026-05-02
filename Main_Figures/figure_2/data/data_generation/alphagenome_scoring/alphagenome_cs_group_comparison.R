# =============================================================================
# AlphaGenome Group Comparison — CS Classification + Score Annotation
# =============================================================================
# Step 3 of the AlphaGenome scoring pipeline.
#
# Reads:
#   - susie_{standard,ash,inf}_credible_sets.rds   (from Step 1, extract_bvsr_results.R)
#   - susie_{standard,ash,inf}_gene_tissue_summary.rds
#   - alphagenome_cs_scores.csv                    (from Step 2, alphagenome_cs_scoring.py)
#
# Pipeline:
#   1. Load 95% credible sets across the 3 methods, restrict to gene-tissue
#      pairs common to all 3.
#   2. Match cross-method CSs by Jaccard >= 0.75 over their variant sets;
#      assign every CS a signal_id via union-find on the resulting graph.
#   3. Classify each signal into one of 7 cross-method groups:
#        Consensus           Standard + ASH + Inf
#        Standard-specific   Standard only
#        ASH-specific        ASH only
#        Inf-specific        Inf only
#        Standard+ASH        Standard + ASH (no Inf)
#        Standard+Inf        Standard + Inf (no ASH)
#        ASH+Inf             ASH + Inf (no Standard)
#   4. Merge AlphaGenome CS-level scores with the group labels.
#   5. Compare score distributions (overall, per modality, per tissue) and
#      run Wilcoxon rank-sum tests of Consensus vs each other group.
#
# Outputs:
#   alphagenome_cs_group_assignments.rds   (Figure 2 panels B & D input)
#   alphagenome_cs_group_scores.csv        (Figure 2 panel D input)
#   alphagenome_cs_group_summary.txt       (text summary; not a figure input)
#
# Usage: Rscript alphagenome_cs_group_comparison.R
# =============================================================================

data_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/alphagenome/data"
output_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark"

# =============================================================================
# Step 1: Load data and build variant sets per CS
# =============================================================================
cat("============================================================\n")
cat("STEP 1: Loading data and building variant sets\n")
cat("============================================================\n")

# Load summaries to find common gene-tissue pairs
std_summary <- readRDS(file.path(data_dir, "susie_standard_gene_tissue_summary.rds"))
ash_summary <- readRDS(file.path(data_dir, "susie_ash_gene_tissue_summary.rds"))
inf_summary <- readRDS(file.path(data_dir, "susie_inf_gene_tissue_summary.rds"))

std_pairs <- paste(std_summary$gene_id, std_summary$tissue, sep = ":")
ash_pairs <- paste(ash_summary$gene_id, ash_summary$tissue, sep = ":")
inf_pairs <- paste(inf_summary$gene_id, inf_summary$tissue, sep = ":")
common_pairs <- Reduce(intersect, list(std_pairs, ash_pairs, inf_pairs))
cat("Common gene-tissue pairs:", length(common_pairs), "\n")

# Load 95% credible sets for all 3 methods
load_cs95 <- function(method_name) {
  cs <- readRDS(file.path(data_dir, paste0("susie_", method_name, "_credible_sets.rds")))
  cs$pair <- paste(cs$gene_id, cs$tissue, sep = ":")
  cs <- cs[cs$pair %in% common_pairs & cs$coverage_level == "0.95", ]
  cs$method <- method_name

  # Parse variants column into a list of character vectors
  cs$variant_set <- strsplit(as.character(cs$variants), ",\\s*")

  cat(sprintf("  %s: %d credible sets at 95%% coverage\n", method_name, nrow(cs)))
  cs
}

std_cs <- load_cs95("standard")
ash_cs <- load_cs95("ash")
inf_cs <- load_cs95("inf")

all_cs <- rbind(
  std_cs[, c("gene_id", "tissue", "cs_id", "method", "pair",
             "n_variants", "top_variant", "top_pip", "variant_set")],
  ash_cs[, c("gene_id", "tissue", "cs_id", "method", "pair",
             "n_variants", "top_variant", "top_pip", "variant_set")],
  inf_cs[, c("gene_id", "tissue", "cs_id", "method", "pair",
             "n_variants", "top_variant", "top_pip", "variant_set")]
)
cat(sprintf("\nTotal CS across all methods: %d\n", nrow(all_cs)))

# =============================================================================
# Step 2: Match CS across methods using Jaccard >= 0.75
# =============================================================================
cat("\n============================================================\n")
cat("STEP 2: Matching CS across methods (Jaccard >= 0.75)\n")
cat("============================================================\n")

JACCARD_THRESHOLD <- 0.75

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union_size <- length(union(a, b))
  if (union_size == 0) return(0)
  inter / union_size
}

# Simple union-find implementation
uf_make <- function(n) seq_len(n)

uf_find <- function(parent, i) {
  while (parent[i] != i) {
    parent[i] <- parent[parent[i]]  # path compression
    i <- parent[i]
  }
  i
}

uf_union <- function(parent, i, j) {
  ri <- uf_find(parent, i)
  rj <- uf_find(parent, j)
  if (ri != rj) parent[ri] <- rj
  parent
}

# Process each gene-tissue pair
unique_pairs <- unique(all_cs$pair)
cat(sprintf("Processing %d gene-tissue pairs...\n", length(unique_pairs)))

signal_list <- list()
signal_counter <- 0L

for (p in unique_pairs) {
  pair_cs <- all_cs[all_cs$pair == p, ]
  n <- nrow(pair_cs)

  if (n == 0) next

  if (n == 1) {
    # Single CS — trivially its own signal
    signal_counter <- signal_counter + 1L
    signal_list[[length(signal_list) + 1]] <- data.frame(
      signal_id = signal_counter,
      gene_id = pair_cs$gene_id[1],
      tissue = pair_cs$tissue[1],
      method = pair_cs$method[1],
      cs_id = pair_cs$cs_id[1],
      top_variant = pair_cs$top_variant[1],
      top_pip = pair_cs$top_pip[1],
      n_variants = pair_cs$n_variants[1],
      stringsAsFactors = FALSE
    )
    next
  }

  # Build union-find over CS indices within this pair
  parent <- uf_make(n)

  # Compare all pairs from different methods
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (pair_cs$method[i] == pair_cs$method[j]) next  # same method, skip
      jac <- jaccard(pair_cs$variant_set[[i]], pair_cs$variant_set[[j]])
      if (jac >= JACCARD_THRESHOLD) {
        parent <- uf_union(parent, i, j)
      }
    }
  }

  # Extract connected components
  roots <- sapply(1:n, function(i) uf_find(parent, i))
  components <- split(1:n, roots)

  for (comp in components) {
    signal_counter <- signal_counter + 1L
    comp_rows <- pair_cs[comp, ]
    for (r in seq_len(nrow(comp_rows))) {
      signal_list[[length(signal_list) + 1]] <- data.frame(
        signal_id = signal_counter,
        gene_id = comp_rows$gene_id[r],
        tissue = comp_rows$tissue[r],
        method = comp_rows$method[r],
        cs_id = comp_rows$cs_id[r],
        top_variant = comp_rows$top_variant[r],
        top_pip = comp_rows$top_pip[r],
        n_variants = comp_rows$n_variants[r],
        stringsAsFactors = FALSE
      )
    }
  }
}

signals <- do.call(rbind, signal_list)
cat(sprintf("Total signals (connected components): %d\n", max(signals$signal_id)))
cat(sprintf("Total CS entries in signals table: %d\n", nrow(signals)))

# =============================================================================
# Step 3: Assign group labels
# =============================================================================
cat("\n============================================================\n")
cat("STEP 3: Assigning group labels\n")
cat("============================================================\n")

# For each signal, determine which methods are present
signal_methods <- tapply(signals$method, signals$signal_id, function(m) {
  sort(unique(m))
})

classify_group <- function(methods) {
  has_std <- "standard" %in% methods
  has_ash <- "ash" %in% methods
  has_inf <- "inf" %in% methods

  if (has_std && has_ash && has_inf) return("Consensus")
  if (has_std && !has_ash && !has_inf) return("Standard-specific")
  if (!has_std && has_ash && !has_inf) return("ASH-specific")
  if (!has_std && !has_ash && has_inf) return("Inf-specific")
  if (has_std && has_ash && !has_inf) return("Standard+ASH")
  if (has_std && !has_ash && has_inf) return("Standard+Inf")
  if (!has_std && has_ash && has_inf) return("ASH+Inf")
  return("Unknown")
}

signal_groups <- sapply(signal_methods, classify_group)
group_df <- data.frame(
  signal_id = as.integer(names(signal_groups)),
  group = unname(signal_groups),
  stringsAsFactors = FALSE
)

signals <- merge(signals, group_df, by = "signal_id")

# Print group counts
cat("\n--- Signal counts per group ---\n")
signal_group_counts <- table(group_df$group)
for (g in names(sort(signal_group_counts, decreasing = TRUE))) {
  cat(sprintf("  %-20s %5d signals\n", g, signal_group_counts[g]))
}
cat(sprintf("  %-20s %5d\n", "TOTAL", sum(signal_group_counts)))

cat("\n--- CS entries per group ---\n")
cs_group_counts <- table(signals$group)
for (g in names(sort(cs_group_counts, decreasing = TRUE))) {
  cat(sprintf("  %-20s %5d CS entries\n", g, cs_group_counts[g]))
}

# Per-method breakdown
cat("\n--- CS per method per group ---\n")
for (m in c("standard", "ash", "inf")) {
  cat(sprintf("\n  %s:\n", m))
  m_signals <- signals[signals$method == m, ]
  m_group_counts <- table(m_signals$group)
  for (g in names(sort(m_group_counts, decreasing = TRUE))) {
    cat(sprintf("    %-20s %5d CS\n", g, m_group_counts[g]))
  }
  cat(sprintf("    %-20s %5d\n", "TOTAL", sum(m_group_counts)))
}

# Save group assignments
assignments_file <- file.path(output_dir, "alphagenome_cs_group_assignments.rds")
saveRDS(signals, assignments_file)
cat(sprintf("\nSaved group assignments to %s\n", assignments_file))

# =============================================================================
# Step 4: Merge with AlphaGenome CS scores
# =============================================================================
cat("\n============================================================\n")
cat("STEP 4: Merging with AlphaGenome CS scores\n")
cat("============================================================\n")

cs_scores <- read.csv(file.path(output_dir, "alphagenome_cs_scores.csv"),
                       stringsAsFactors = FALSE)
cat(sprintf("Loaded %d CS score rows\n", nrow(cs_scores)))

# Join: attach signal_id and group to each CS score row
score_key <- signals[, c("gene_id", "tissue", "method", "cs_id", "signal_id", "group")]
annotated <- merge(cs_scores, score_key,
                   by = c("gene_id", "tissue", "method", "cs_id"),
                   all.x = FALSE)

cat(sprintf("Matched %d / %d CS score rows to groups\n", nrow(annotated), nrow(cs_scores)))

# Check for unmatched
unmatched <- nrow(cs_scores) - nrow(annotated)
if (unmatched > 0) {
  cat(sprintf("WARNING: %d CS score rows could not be matched (likely non-common pairs)\n", unmatched))
}

# Save annotated scores
annotated_file <- file.path(output_dir, "alphagenome_cs_group_scores.csv")
write.csv(annotated, annotated_file, row.names = FALSE)
cat(sprintf("Saved annotated scores to %s\n", annotated_file))

# =============================================================================
# Step 5: Compare AlphaGenome score distributions across groups
# =============================================================================
cat("\n============================================================\n")
cat("STEP 5: Comparing AlphaGenome CS scores across groups\n")
cat("============================================================\n")

valid <- annotated[!is.na(annotated$cs_score), ]

lines <- character()
lines <- c(lines,
  "================================================================================",
  "AlphaGenome CS Score Comparison by Method-Agreement Group",
  "CS_Score = sum(PIP * |quantile_score|) / sum(PIP); Null = 0.5",
  sprintf("Jaccard threshold for matching: %.2f", JACCARD_THRESHOLD),
  "================================================================================")

# --- Overall per group ---
lines <- c(lines, "\n--- Overall per Group ---")
lines <- c(lines, sprintf("%-20s %6s %8s %8s %8s %7s",
                           "Group", "N_CS", "Mean", "Median", "Std", "%>0.5"))
lines <- c(lines, paste(rep("-", 65), collapse = ""))

group_order <- c("Consensus", "Standard+ASH", "Standard+Inf", "ASH+Inf",
                 "Standard-specific", "ASH-specific", "Inf-specific")

for (g in group_order) {
  sub <- valid[valid$group == g, ]
  if (nrow(sub) == 0) next
  n <- nrow(sub)
  m <- mean(sub$cs_score)
  med <- median(sub$cs_score)
  s <- sd(sub$cs_score)
  pct <- 100 * sum(sub$cs_score > 0.5) / n
  lines <- c(lines, sprintf("%-20s %6d %8.4f %8.4f %8.4f %6.1f%%", g, n, m, med, s, pct))
}

# --- Per group x modality ---
lines <- c(lines, "\n\n--- Per Group x Modality ---")
lines <- c(lines, sprintf("%-20s %-15s %6s %8s %8s %7s",
                           "Group", "Modality", "N_CS", "Mean", "Median", "%>0.5"))
lines <- c(lines, paste(rep("-", 75), collapse = ""))

for (g in group_order) {
  for (mod in c("RNA_SEQ", "DNASE", "CHIP_HISTONE", "CHIP_TF")) {
    sub <- valid[valid$group == g & valid$modality == mod, ]
    if (nrow(sub) == 0) next
    n <- nrow(sub)
    m <- mean(sub$cs_score)
    med <- median(sub$cs_score)
    pct <- 100 * sum(sub$cs_score > 0.5) / n
    lines <- c(lines, sprintf("%-20s %-15s %6d %8.4f %8.4f %6.1f%%",
                               g, mod, n, m, med, pct))
  }
}

# --- Statistical tests: Consensus vs each other group ---
lines <- c(lines, "\n\n--- Wilcoxon Rank-Sum: Consensus vs Other Groups ---")
lines <- c(lines, sprintf("%-20s %-15s %8s %12s %10s",
                           "Group", "Modality", "N_group", "W-stat", "P-value"))
lines <- c(lines, paste(rep("-", 70), collapse = ""))

consensus_scores <- valid[valid$group == "Consensus", ]

for (g in group_order[-1]) {
  group_scores <- valid[valid$group == g, ]
  for (mod in c("RNA_SEQ", "DNASE", "CHIP_HISTONE", "CHIP_TF")) {
    con_mod <- consensus_scores[consensus_scores$modality == mod, "cs_score"]
    grp_mod <- group_scores[group_scores$modality == mod, "cs_score"]

    if (length(con_mod) < 5 || length(grp_mod) < 5) next

    wt <- tryCatch(
      wilcox.test(con_mod, grp_mod),
      error = function(e) NULL
    )
    if (!is.null(wt)) {
      sig <- ""
      if (wt$p.value < 0.001) sig <- "***"
      else if (wt$p.value < 0.01) sig <- "**"
      else if (wt$p.value < 0.05) sig <- "*"

      lines <- c(lines, sprintf("%-20s %-15s %8d %12.1f %10.2e %s",
                                 g, mod, length(grp_mod), wt$statistic, wt$p.value, sig))
    }
  }
}

# --- Per tissue breakdown ---
lines <- c(lines, "\n\n--- Per Tissue x Group ---")
for (tis in c("AC", "DLPFC", "PCC")) {
  tis_df <- valid[valid$tissue == tis, ]
  if (nrow(tis_df) == 0) next
  lines <- c(lines, sprintf("\n  %s:", tis))
  lines <- c(lines, sprintf("  %-20s %6s %8s %8s %7s",
                             "Group", "N_CS", "Mean", "Median", "%>0.5"))
  for (g in group_order) {
    sub <- tis_df[tis_df$group == g, ]
    if (nrow(sub) == 0) next
    n <- nrow(sub)
    m <- mean(sub$cs_score)
    med <- median(sub$cs_score)
    pct <- 100 * sum(sub$cs_score > 0.5) / n
    lines <- c(lines, sprintf("  %-20s %6d %8.4f %8.4f %6.1f%%", g, n, m, med, pct))
  }
}

# --- Verification ---
lines <- c(lines, "\n\n--- Verification ---")
for (m in c("standard", "ash", "inf")) {
  m_assigned <- sum(signals$method == m)
  if (m == "standard") m_expected <- nrow(std_cs)
  else if (m == "ash") m_expected <- nrow(ash_cs)
  else m_expected <- nrow(inf_cs)
  status <- if (m_assigned == m_expected) "OK" else "MISMATCH"
  lines <- c(lines, sprintf("  %s: %d assigned, %d expected — %s",
                             m, m_assigned, m_expected, status))
}

text <- paste(lines, collapse = "\n")
cat(text)

summary_file <- file.path(output_dir, "alphagenome_cs_group_summary.txt")
writeLines(text, summary_file)
cat(sprintf("\n\nSaved summary to %s\n", summary_file))

cat("\nDone!\n")
