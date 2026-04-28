#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S2: Unmappable Effects Example
# =============================================================================
# Fits all models from the susie_unmappable_effects vignette and saves as RDS.
#
# Output:
#   - data_and_sumstats.rds  (X, y, beta, z-scores, p-values)
#   - fit_susie.rds          (SuSiE, L = 10)
#   - fit_inf.rds            (SuSiE-inf, L = 10)
#   - fit_ash.rds            (SuSiE-ash, L = 10)
#   - fit_susie_L40.rds      (SuSiE, L = 40)
# =============================================================================

pkgload::load_all("/Users/alexmccreight/StatFunGen/susieR", quiet = TRUE)

s2_dir     <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S2"
output_dir <- file.path(s2_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load data
# =============================================================================
cat("Loading data...\n")
data(unmappable_data)
X <- unmappable_data$X
y <- as.vector(unmappable_data$y)
b <- unmappable_data$beta

# Compute summary statistics
cat("Computing summary statistics...\n")
sumstats <- univariate_regression(X, y)
z_scores <- sumstats$betahat / sumstats$sebetahat
pvals <- 2 * pnorm(-abs(z_scores))

saveRDS(list(X = X, y = y, b = b, z_scores = z_scores, pvals = pvals),
        file.path(output_dir, "data_and_sumstats.rds"))

# =============================================================================
# Fit SuSiE (L = 10)
# =============================================================================
cat("Fitting SuSiE (L = 10)...\n")
t0 <- proc.time()
fit_susie <- susie(X, y, L = 10)
t1 <- proc.time()
cat(sprintf("  Completed in %.1f seconds\n", (t1 - t0)[3]))
saveRDS(fit_susie, file.path(output_dir, "fit_susie.rds"))

# =============================================================================
# Fit SuSiE-inf (L = 10)
# =============================================================================
cat("Fitting SuSiE-inf (this may take several minutes)...\n")
t0 <- proc.time()
fit_inf <- susie(X, y, L = 10, unmappable_effects = "inf")
t1 <- proc.time()
cat(sprintf("  Completed in %.1f seconds\n", (t1 - t0)[3]))
saveRDS(fit_inf, file.path(output_dir, "fit_inf.rds"))

# =============================================================================
# Fit SuSiE-ash (L = 10)
# =============================================================================
cat("Fitting SuSiE-ash (this may take several minutes)...\n")
t0 <- proc.time()
fit_ash <- susie(X, y, L = 10, unmappable_effects = "ash")
t1 <- proc.time()
cat(sprintf("  Completed in %.1f seconds\n", (t1 - t0)[3]))
saveRDS(fit_ash, file.path(output_dir, "fit_ash.rds"))

# =============================================================================
# Fit SuSiE (L = 40)
# =============================================================================
cat("Fitting SuSiE (L = 40)...\n")
t0 <- proc.time()
fit_susie_L40 <- susie(X, y, L = 40)
t1 <- proc.time()
cat(sprintf("  Completed in %.1f seconds\n", (t1 - t0)[3]))
saveRDS(fit_susie_L40, file.path(output_dir, "fit_susie_L40.rds"))

# =============================================================================
# Summary
# =============================================================================
cat("\nAll analyses complete! Saved files:\n")
cat("  - data_and_sumstats.rds\n")
cat("  - fit_susie.rds\n")
cat("  - fit_inf.rds\n")
cat("  - fit_ash.rds\n")
cat("  - fit_susie_L40.rds\n")

cat("\n=== Summary ===\n")
cat(sprintf("SuSiE (L=10):  %d credible sets\n", length(fit_susie$sets$cs)))
cat(sprintf("SuSiE-inf:     %d credible sets\n", length(fit_inf$sets$cs)))
cat(sprintf("SuSiE-ash:     %d credible sets\n", length(fit_ash$sets$cs)))
cat(sprintf("SuSiE (L=40):  %d credible sets\n", length(fit_susie_L40$sets$cs)))
