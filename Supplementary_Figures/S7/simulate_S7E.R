#!/usr/bin/env Rscript

# Supplementary Figure S7, panel E: full posterior mass of the credible sets
# that two reporting rules build from the same fitted SuSiE single-effect
# posteriors.
#   SparsePro   original 2023 get_effect() implementation
#   susieR      attainable screening followed by 95% CS construction
#
# Simulation:
#   - Public 121-variant ADSP sketch example from susieR.
#   - K = 1 or 2 causal variants; total signal strength nh^2 = 24 or 80.
#   - 150 replicates per setting; seed = 83011.
#   - Retain fitted effects for which both rules report a CS.
#
# Requirements:
#   - susieR 0.16.6 installed from commit
#     ef213feed2cb82419677661a8c986e1504df2c73 (the fitting build).
#   - Python 3 with numpy, pandas, and scipy (the PYTHON environment variable
#     selects the interpreter).
#   - Network access: the SparsePro source, the ADSP example and the susieR
#     reporting code are downloaded at the commits pinned below.
#
# Usage:  Rscript Supplementary_Figures/S7/simulate_S7E.R [output_dir]
# Output: data/paired-posterior-mass.csv, or the same file in output_dir

source("R/paths.R")

# ---- Paths and pinned source versions ----

s7_dir <- supp_dir(7)
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1] else file.path(s7_dir, "data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
python <- Sys.getenv("PYTHON", "python3")
if (system2(python, c("-c", shQuote("import numpy, pandas, scipy")),
            stdout = FALSE, stderr = FALSE) != 0L) {
  stop("Use a Python with numpy, pandas and scipy installed; set PYTHON if needed.")
}
work <- tempfile("s7e-")
dir.create(work)

if (as.character(packageVersion("susieR")) != "0.16.6") {
  stop("Install the pinned susieR 0.16.6 build listed in Supplementary_Figures/index.md.")
}
sparsepro_commit <- "e13ac87d6e650d7682ff1f33f6f9898d66932d64"
susier_commit <- "ed6e491e3d0235ea5fa2e8545ea206b5d759eb6e"
sources <- c(
  sparsepro_zld.py = sprintf("zhwm/SparsePro/%s/sparsepro_zld.py", sparsepro_commit),
  rss_mismatch_cases.RData = sprintf("stephenslab/susieR/%s/data/rss_mismatch_cases.RData", susier_commit),
  susie_get_functions.R = sprintf("stephenslab/susieR/%s/R/susie_get_functions.R", susier_commit)
)
for (name in names(sources)) {
  download.file(paste0("https://raw.githubusercontent.com/", sources[name]),
                file.path(work, name), mode = "wb", quiet = TRUE)
}

# ---- Load the susieR reporting rule and ADSP LD ----

reporting <- new.env(parent = asNamespace("susieR"))
sys.source(file.path(work, "susie_get_functions.R"), envir = reporting)
data <- new.env()
load(file.path(work, "rss_mismatch_cases.RData"), envir = data)
R <- unname(cor(data$rss_mismatch_cases$ld_mismatch_sketch$X))
stopifnot(ncol(R) == 121L)
eig <- eigen(R, symmetric = TRUE)
root <- sweep(eig$vectors, 2, sqrt(pmax(eig$values, 0)), "*")

# ---- Simulate z-scores, fit SuSiE and record the susieR sets ----

settings <- expand.grid(strength = c(24L, 80L), K = 1:2, rep = seq_len(150))
posteriors <- sets <- vector("list", nrow(settings))
set.seed(83011)
cat("Running 600 simulations for Panel E...\n")

for (i in seq_len(nrow(settings))) {
  k <- settings$K[i]
  strength <- settings$strength[i]
  causal <- sample.int(ncol(R), k)
  effect <- sample(c(-1, 1), k, replace = TRUE) * sqrt(strength / k)
  z <- drop(root %*% rnorm(ncol(R)) + R[, causal, drop = FALSE] %*% effect)
  fit <- suppressMessages(suppressWarnings(
    susieR::susie_rss(z, R, n = 50000, L = 10, z_method = "score")
  ))
  if (!isTRUE(fit$converged)) next
  alpha <- fit$alpha
  colnames(alpha) <- paste0("alpha_", seq_len(ncol(alpha)))
  truth <- paste(causal, collapse = " ")
  posteriors[[i]] <- data.frame(case = i, K = k, strength, truth,
                                 effect = seq_len(nrow(alpha)), alpha)
  cs <- suppressMessages(reporting$susie_get_cs_attainable(fit))
  if (!length(cs$cs)) next

  sets[[i]] <- data.frame(
    case = i, K = k, strength,
    effect = as.integer(sub("^L", "", names(cs$cs))),
    mass = cs$coverage, size = lengths(cs$cs),
    contains_truth = vapply(cs$cs, function(s) any(causal %in% s), logical(1))
  )
}

# ---- Apply the original SparsePro reporting rule to the same posteriors ----
# Python imports the original class without running its command-line interface.
# Each reported CS is measured against the full single-effect posterior, then
# paired with the susieR CS for the same fitted effect.

write.csv(do.call(rbind, posteriors), file.path(work, "posteriors.csv"), row.names = FALSE)
write.csv(do.call(rbind, sets), file.path(work, "susier-sets.csv"), row.names = FALSE)
python_code <- r"(
import ast
from contextlib import redirect_stdout
import io
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def sparsepro_class(path):
    # Load the original class without running the file's top-level CLI.
    tree = ast.parse(path.read_text())
    tree.body = [
        node for node in tree.body
        if isinstance(node, (ast.Import, ast.ImportFrom))
        or (isinstance(node, ast.ClassDef) and node.name == "SparsePro")
    ]
    module = {}
    exec(compile(tree, str(path), "exec"), module)
    return module["SparsePro"]


def pair_sets(work):
    keys = ["case", "K", "strength", "effect"]
    posteriors = pd.read_csv(work / "posteriors.csv")
    sparsepro = sparsepro_class(work / "sparsepro_zld.py")
    reported = []
    for case, fit in posteriors.groupby("case"):
        alpha = fit.filter(regex=r"^alpha_").to_numpy()
        effects, p = alpha.shape
        model = sparsepro(p, effects, np.ones(p), 0, 1)
        model.gamma = alpha.T
        with redirect_stdout(io.StringIO()):
            sets, *_ = model.get_effect()
        truth = np.array(str(fit["truth"].iloc[0]).split(), dtype=int) - 1
        for effect, members in sets.items():
            if len(members):
                reported.append((
                    case, fit["K"].iloc[0], fit["strength"].iloc[0], effect + 1,
                    alpha[effect, members].sum(), len(members), np.isin(members, truth).any(),
                ))
    columns = keys + ["mass", "size", "contains_truth"]
    paired = pd.DataFrame(reported, columns=columns).merge(
        pd.read_csv(work / "susier-sets.csv"), on=keys,
        suffixes=("_sparsepro", "_susier"), validate="one_to_one",
    )
    print(f"{posteriors['case'].nunique()} converged fits, {len(paired)} paired sets")
    return paired.sort_values(keys).reset_index(drop=True)


work = Path(sys.argv[1])
paired = pair_sets(work)
paired.to_csv(work / "paired-posterior-mass.csv", index=False)
)"
status <- system2(python, c("-W", "ignore::SyntaxWarning", "-", shQuote(work)),
                   input = python_code)
if (status != 0L) stop("SparsePro reporting failed.")
paired <- read.csv(file.path(work, "paired-posterior-mass.csv"))

# ---- Save results and print summary ----

output_file <- file.path(output_dir, "paired-posterior-mass.csv")
if (!file.copy(file.path(work, "paired-posterior-mass.csv"), output_file, overwrite = TRUE)) {
  stop("Could not write ", output_file)
}
unlink(work, recursive = TRUE)

cat(sprintf("Saved: %s\n", output_file))
cat("\n=== Summary ===\n")
cat(sprintf("SparsePro: %d/%d CSs below 95%% posterior mass\n",
            sum(paired$mass_sparsepro < 0.95 - 1e-12), nrow(paired)))
cat(sprintf("susieR:    %d/%d CSs below 95%% posterior mass\n",
            sum(paired$mass_susier < 0.95 - 1e-12), nrow(paired)))
