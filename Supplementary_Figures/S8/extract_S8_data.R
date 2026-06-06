#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S8: Extract complex Top-N Power/FDR for the
#                  no-extension vs LD-extension (0.99) comparison
# =============================================================================
# For the three complex scenarios, recompute Top-N Power/FDR (same Top-N
# definition as Figure 1C-1F) for BOTH credible-set constructions:
#   xcorr_noext -> "SuSiE" / "SuSiE-inf" / "SuSiE-ash"   (no extension)
#   xcorr_ext   -> "... w/ Extension"                    (LD ext 0.99)
# so S8 can show 2 lines per method per panel (6 lines).
#
# For each replicate and each N in {3,5,10,15,20,25}:
#   Top-N causal set = N variants with largest |beta|
#   Power = fraction of Top-N captured by >=1 CS; FDR = CS-with-no-Top-N / n_cs
#
# Inputs (read-only):  500_rep_results/{complex,complex_S1,complex_S2}/sim_*_ldext.rds
# Output:              data/s8_finemapping_data.rds
# =============================================================================

ldext_root <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results"
s8_dir     <- "/Users/alexmccreight/StatFunGen/susieR2.0-paper/Supplementary_Figures/S8"
output_dir <- file.path(s8_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

N_values     <- c(3, 5, 10, 15, 20, 25)
versions     <- c("xcorr_noext", "xcorr_ext")
method_names <- c("susie" = "SuSiE", "susie_inf" = "SuSiE-inf", "susie_ash" = "SuSiE-ash")

scenarios <- list(
  list(label = "Complex",
       file  = file.path(ldext_root, "complex",
         "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInf15_n1000_ldext.rds")),
  list(label = "Complex S1",
       file  = file.path(ldext_root, "complex_S1",
         "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_ldext.rds")),
  list(label = "Complex S2",
       file  = file.path(ldext_root, "complex_S2",
         "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_ldext.rds"))
)

# -----------------------------------------------------------------------------
# Top-N power / FDR for a credible-set list (identical to Fig 1 extraction)
# -----------------------------------------------------------------------------
compute_metrics <- function(cs_list, causal_set) {
  if (is.null(cs_list) || length(cs_list) == 0) {
    return(list(Power = 0, FDR = NA))
  }
  n_causal <- length(causal_set)
  n_cs     <- length(cs_list)
  causal_covered <- sapply(causal_set, function(v) any(sapply(cs_list, function(cs) v %in% cs)))
  power <- sum(causal_covered) / n_causal
  cs_has_causal <- sapply(cs_list, function(cs) any(cs %in% causal_set))
  fdr <- sum(!cs_has_causal) / n_cs
  list(Power = power, FDR = fdr)
}

# -----------------------------------------------------------------------------
# Process each scenario
# -----------------------------------------------------------------------------
all_rows <- list()

for (sc in scenarios) {
  cat(sprintf("\n=== %s ===\n", sc$label))
  if (!file.exists(sc$file)) stop("Missing ldext file: ", sc$file)
  sim  <- readRDS(sc$file)
  reps <- sim$replicates

  for (r in seq_along(reps)) {
    rep <- reps[[r]]
    if (is.null(rep) || !identical(rep$status, "ok")) next

    beta       <- rep$beta
    ranked_idx <- order(abs(beta), decreasing = TRUE)

    for (N in N_values) {
      causal_set <- ranked_idx[1:N]
      for (mkey in names(method_names)) {
        if (is.null(rep$credible_sets[[mkey]])) next
        for (v in versions) {
          m <- compute_metrics(rep$credible_sets[[mkey]][[v]]$cs, causal_set)
          all_rows[[length(all_rows) + 1]] <- data.frame(
            scenario  = sc$label,
            Method    = method_names[[mkey]],
            version   = v,
            N         = N,
            Power     = m$Power,
            FDR       = m$FDR,
            replicate = r
          )
        }
      }
    }
    if (r %% 100 == 0) cat(sprintf("  processed %d / %d replicates\n", r, length(reps)))
  }
  rm(sim, reps); gc(verbose = FALSE)
}

rep_data <- do.call(rbind, all_rows)

# -----------------------------------------------------------------------------
# Aggregate scenario x Method x version x N  (mean +/- SE across reps)
# -----------------------------------------------------------------------------
agg <- do.call(rbind, lapply(
  split(rep_data, list(rep_data$scenario, rep_data$Method, rep_data$version, rep_data$N), drop = TRUE),
  function(df) {
    data.frame(
      scenario   = df$scenario[1],
      Method     = df$Method[1],
      version    = df$version[1],
      N          = df$N[1],
      Power_mean = mean(df$Power, na.rm = TRUE),
      Power_se   = sd(df$Power, na.rm = TRUE) / sqrt(sum(!is.na(df$Power))),
      FDR_mean   = mean(df$FDR, na.rm = TRUE),
      FDR_se     = sd(df$FDR, na.rm = TRUE) / sqrt(sum(!is.na(df$FDR)))
    )
  }))
rownames(agg) <- NULL

# 6-level series factor (method-major, version-minor) -- same labels as S7
series_levels <- c("SuSiE", "SuSiE w/ Extension",
                   "SuSiE-ash", "SuSiE-ash w/ Extension",
                   "SuSiE-inf", "SuSiE-inf w/ Extension")
agg$series <- ifelse(agg$version == "xcorr_ext",
                     paste0(agg$Method, " w/ Extension"), agg$Method)
agg$series <- factor(agg$series, levels = series_levels)

# scenario factor ordered top -> bottom for the figure (S2, S1, Complex)
agg$scenario <- factor(agg$scenario, levels = c("Complex S2", "Complex S1", "Complex"))

output_file <- file.path(output_dir, "s8_finemapping_data.rds")
saveRDS(agg, output_file)

cat(sprintf("\nSaved: %s\n", output_file))
cat(sprintf("Dimensions: %d rows  (expect 108 = 3 scenarios x 6 series x 6 N)\n", nrow(agg)))
cat("\n=== N=3 spot check (Power_mean / FDR_mean) ===\n")
print(agg[agg$N == 3, c("scenario", "series", "Power_mean", "FDR_mean")])
