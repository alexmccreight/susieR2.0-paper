#!/usr/bin/env Rscript

# Supplementary Figure S5: fit SuSiE, SuSiE-inf and SuSiE-ash (L = 10) to
# susieR's `unmappable_data` example and save the fits for plot_S5.R.
#
# The susieR source tree is loaded with pkgload (see susier_root() in
# R/paths.R; set SUSIER_ROOT if it is not at ../susieR).
#
# Output: data/data_and_sumstats.rds  (X, y, true effects, z-scores, p-values)
#         data/fit_susie.rds, data/fit_inf.rds, data/fit_ash.rds

source("R/paths.R")

pkgload::load_all(susier_root(), quiet = TRUE)

output_dir <- file.path(supp_dir(5), "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data(unmappable_data)
X <- unmappable_data$X
y <- as.vector(unmappable_data$y)
b <- unmappable_data$beta

sumstats <- univariate_regression(X, y)
z_scores <- sumstats$betahat / sumstats$sebetahat
pvals    <- 2 * pnorm(-abs(z_scores))

fits <- list(
  susie = susie(X, y, L = 10),
  inf   = susie(X, y, L = 10, unmappable_effects = "inf"),
  ash   = susie(X, y, L = 10, unmappable_effects = "ash")
)

outputs <- c(
  list(data_and_sumstats = list(X = X, y = y, b = b, z_scores = z_scores, pvals = pvals)),
  setNames(fits, paste0("fit_", names(fits)))
)
for (name in names(outputs)) {
  out <- file.path(output_dir, paste0(name, ".rds"))
  saveRDS(outputs[[name]], out)
  cat("Saved:", out, "\n")
}
