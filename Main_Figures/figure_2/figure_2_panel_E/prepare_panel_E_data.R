#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel E — Data Preparation
# =============================================================================
# Loads SuSiE fit objects for ENSG00000163431 (DLPFC tissue) from all three
# methods (standard, ash, inf), extracts PIPs, variant positions, and
# credible set membership, and saves a clean plotting-ready data frame.
#
# Gene:   ENSG00000163431 (chr1)
# Tissue: DLPFC (dorsolateral prefrontal cortex)
# Story:  SuSiE-ash finds a novel 6-variant CS at ~201.806 Mb (CS1, lead
#         PIP 0.368) that standard and inf do not detect. All three methods
#         share a signal at ~201.917 Mb. Standard also has a diffuse
#         secondary CS2 spanning ~201.788-201.935 Mb.
#
# Input:  STANDARD/ASH/INF fit RDS files (in this directory)
# Output: panel_E_data.rds
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

fig2_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2"
script_dir <- file.path(fig2_dir, "figure_2_panel_E")
data_dir   <- file.path(fig2_dir, "data", "panel_E")

fit_files <- list(
  standard = file.path(data_dir, "STANDARD-ROSMAP_DeJager_eQTL_twas.chr1_ENSG00000163431.univariate_bvsr.rds"),
  ash      = file.path(data_dir, "ASH-ROSMAP_DeJager_eQTL_twas.chr1_ENSG00000163431.univariate_bvsr.rds"),
  inf      = file.path(data_dir, "INF-ROSMAP_DeJager_eQTL_twas.chr1_ENSG00000163431.univariate_bvsr.rds")
)

gene_id    <- "ENSG00000163431"
tissue_key <- "DLPFC_DeJager_eQTL_ENSG00000163431"

method_labels <- c(standard = "SuSiE", ash = "SuSiE-ash", inf = "SuSiE-inf")
method_order  <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")

# =============================================================================
# Extract data from each method
# =============================================================================

cat("Extracting PIP data for ENSG00000163431 in DLPFC...\n\n")

all_variants <- list()
cs_summaries <- list()

for (method_key in names(fit_files)) {
  cat(sprintf("  Loading %s...\n", method_labels[method_key]))

  raw <- readRDS(fit_files[[method_key]])
  tissue_data <- raw[[gene_id]][[tissue_key]]
  fit <- tissue_data$susie_result_trimmed
  var_names <- tissue_data$variant_names

  # Parse genomic positions from variant names (format: chr:pos:ref:alt)
  positions <- as.numeric(sapply(strsplit(var_names, ":"), `[`, 2))

  # Build per-variant data frame
  n_vars <- length(var_names)
  df <- data.frame(
    variant_name = var_names,
    position     = positions,
    pip          = fit$pip,
    in_cs        = rep(FALSE, n_vars),
    cs_id        = rep(NA_integer_, n_vars),
    method       = method_labels[method_key],
    stringsAsFactors = FALSE
  )

  # Mark CS membership
  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
    for (cs_idx in seq_along(fit$sets$cs)) {
      cs_vars <- fit$sets$cs[[cs_idx]]
      df$in_cs[cs_vars] <- TRUE
      df$cs_id[cs_vars] <- cs_idx
    }
  }

  all_variants[[method_key]] <- df

  # CS summary
  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
    for (cs_idx in seq_along(fit$sets$cs)) {
      cs_vars <- fit$sets$cs[[cs_idx]]
      cs_pips <- fit$pip[cs_vars]
      lead_idx <- cs_vars[which.max(cs_pips)]

      cs_summaries[[paste(method_key, cs_idx, sep = "_")]] <- data.frame(
        method       = method_labels[method_key],
        cs_id        = cs_idx,
        n_variants   = length(cs_vars),
        lead_variant = var_names[lead_idx],
        lead_pip     = max(cs_pips),
        coverage     = fit$sets$coverage[cs_idx],
        min_abs_corr = fit$sets$purity[cs_idx, 1],
        stringsAsFactors = FALSE
      )
    }
  } else {
    cat(sprintf("    No credible sets found for %s\n", method_labels[method_key]))
  }
}

# Combine
pip_data <- do.call(rbind, all_variants)
pip_data$method <- factor(pip_data$method, levels = method_order)
rownames(pip_data) <- NULL

cs_summary <- do.call(rbind, cs_summaries)
cs_summary$method <- factor(cs_summary$method, levels = method_order)
rownames(cs_summary) <- NULL

# Get region info from one of the fits
region_info <- readRDS(fit_files[["standard"]])[[gene_id]][[tissue_key]]$region_info

# =============================================================================
# Print summary
# =============================================================================

cat("\n--- PIP data summary ---\n")
cat(sprintf("  Total variant-method rows: %d\n", nrow(pip_data)))
cat(sprintf("  Variants per method: %d\n", nrow(pip_data) / length(method_order)))
cat(sprintf("  Genomic range: chr1:%s-%s\n",
            format(min(pip_data$position), big.mark = ","),
            format(max(pip_data$position), big.mark = ",")))

cat("\n--- Credible set summary ---\n")
for (i in seq_len(nrow(cs_summary))) {
  s <- cs_summary[i, ]
  cat(sprintf("  %-10s CS%d: %d variants, lead=%s (PIP=%.4f), purity=%.4f\n",
              s$method, s$cs_id, s$n_variants, s$lead_variant,
              s$lead_pip, s$min_abs_corr))
}

cat(sprintf("\n--- Variants in CS per method ---\n"))
for (m in method_order) {
  n_in_cs <- sum(pip_data$in_cs[pip_data$method == m])
  max_pip <- max(pip_data$pip[pip_data$method == m])
  cat(sprintf("  %-10s %d variants in CS, max PIP = %.4f\n", m, n_in_cs, max_pip))
}

# =============================================================================
# Save
# =============================================================================

output <- list(
  pip_data    = pip_data,
  cs_summary  = cs_summary,
  region_info = region_info,
  gene_id     = gene_id,
  tissue      = "DLPFC"
)

output_path <- file.path(script_dir, "panel_E_data.rds")
saveRDS(output, output_path)
cat(sprintf("\nSaved: %s\n", output_path))
