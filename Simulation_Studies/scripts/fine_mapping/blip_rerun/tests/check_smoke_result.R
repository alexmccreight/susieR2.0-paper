#!/usr/bin/env Rscript

# Supplementary Figure S7, SuSiE+BLiP re-run: checks run by the smoke test.
#   1. The local and HPC routes wrote the same replicate folders: alpha.npy and
#      corr_sub.npy byte for byte, meta.json apart from the timestamp,
#      blip_arms.json apart from timings and solver diagnostics.
#   2. At least one replicate has a SuSiE(0.5) CS above the size cap, and on
#      those replicates BLiP A_capped has lower power than SuSiE(0.5) and
#      B_uncapped has higher power than A_capped.
# Also prints the pooled table at level q stratified by the largest SuSiE(0.5)
# CS, and the per-replicate power of the key methods where that CS exceeds the
# cap. Exits with status 1 when a check fails.
#
# Input: <local_dir>/metrics/blip_metrics.rds and the replicate folders under
#        <local_dir> and <hpc_dir>
#
# Usage: Rscript check_smoke_result.R local_dir=<dir> [hpc_dir=<dir>] [q=0.05]

suppressMessages({
  library(jsonlite)
  library(dplyr)
})

.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.script_dir, "..", "blip_common.R"))

strip <- function(x, drop) x[setdiff(names(x), drop)]

norm_meta <- function(path) {
  m <- fromJSON(path, simplifyVector = FALSE)
  m$export <- strip(m$export, "timestamp")
  m
}

norm_arms <- function(path) {
  a <- fromJSON(path, simplifyVector = FALSE)
  a <- strip(a, c("elapsed_seconds", "timestamp"))
  a$arms <- lapply(a$arms, function(by_q) {
    lapply(by_q, function(r) strip(r, c("elapsed_seconds", "problem_status")))
  })
  a
}

compare_routes <- function(local_dir, hpc_dir) {
  rep_dirs <- list.dirs(local_dir, recursive = TRUE)
  rep_dirs <- rep_dirs[grepl("^rep_[0-9]+$", basename(rep_dirs))]
  n_bad <- 0L
  for (d in rep_dirs) {
    rel <- sub(paste0("^", normalizePath(local_dir), "/"), "", normalizePath(d))
    h <- file.path(hpc_dir, rel)
    problems <- character(0)
    if (!dir.exists(h)) {
      problems <- "missing in the HPC route"
    } else {
      for (f in c("alpha.npy", "corr_sub.npy")) {
        if (!identical(unname(tools::md5sum(file.path(d, f))), unname(tools::md5sum(file.path(h, f))))) {
          problems <- c(problems, paste(f, "differs"))
        }
      }
      if (!identical(norm_meta(file.path(d, "meta.json")), norm_meta(file.path(h, "meta.json")))) {
        problems <- c(problems, "meta.json differs")
      }
      if (!identical(norm_arms(file.path(d, "blip_arms.json")), norm_arms(file.path(h, "blip_arms.json")))) {
        problems <- c(problems, "blip_arms.json differs")
      }
    }
    if (length(problems) > 0) {
      n_bad <- n_bad + 1L
      cat(sprintf("  %s: %s\n", rel, paste(problems, collapse = "; ")))
    }
  }
  cat(sprintf("Local vs HPC route: %d replicate folders compared, %d differ\n", length(rep_dirs), n_bad))
  n_bad == 0L && length(rep_dirs) > 0
}

main <- function() {
  a <- parse_key_value_args(list(local_dir = NULL, hpc_dir = "", q = 0.05))
  if (is.null(a$local_dir)) stop("local_dir must be specified", call. = FALSE)
  checks <- c()

  if (nzchar(a$hpc_dir)) {
    checks["local and HPC routes give identical per-replicate results"] <- compare_routes(a$local_dir, a$hpc_dir)
  }

  m   <- readRDS(file.path(a$local_dir, "metrics", "blip_metrics.rds"))
  cap <- m$meta$size_cap
  gt  <- sprintf("max_cs>%d", cap)
  le  <- sprintf("max_cs<=%d", cap)
  s   <- m$summary %>% filter(q == a$q, K_label == "pooled")

  show <- s %>%
    filter(stratum_label %in% c(le, gt, "all")) %>%
    mutate(stratum_label = factor(stratum_label, levels = c(le, gt, "all"))) %>%
    arrange(stratum_label, method) %>%
    select(stratum_label, method, n_reps, power_mean, res_power_mean, fdr_mean,
           n_cs_mean, size_median, size_max, frac_size_gt_cap, purity_median)
  cat(sprintf("\nq = %g, pooled, stratified by the largest SuSiE(0.5) CS:\n", a$q))
  print(as.data.frame(show), digits = 3, row.names = FALSE)

  key_methods <- c("SuSiE(0.5)", sprintf("SuSiE capped(%d)", cap), "BLiP A_capped", "BLiP B_uncapped")
  per_rep <- m$rep_metrics %>%
    filter(q == a$q, stratum == gt, method %in% key_methods) %>%
    select(replicate_id, max_susie05_cs, method, n_cs, power, size_mean) %>%
    arrange(replicate_id, method)
  cat(sprintf("\nReplicates whose largest SuSiE(0.5) CS exceeds %d:\n", cap))
  print(as.data.frame(per_rep), digits = 3, row.names = FALSE)

  pw <- function(method, stratum) {
    v <- s$power_mean[s$method == method & s$stratum_label == stratum]
    if (length(v) == 0) NA_real_ else v
  }
  n_gt <- s$n_reps[s$method == "SuSiE(0.5)" & s$stratum_label == gt]
  checks[sprintf("some replicates have a SuSiE(0.5) CS above %d", cap)] <- length(n_gt) == 1 && n_gt >= 1
  checks["A_capped has lower power than SuSiE(0.5) on those replicates"] <-
    isTRUE(pw("BLiP A_capped", gt) < pw("SuSiE(0.5)", gt))
  checks["B_uncapped has higher power than A_capped on those replicates"] <-
    isTRUE(pw("BLiP B_uncapped", gt) > pw("BLiP A_capped", gt))

  cat(sprintf("\nPower on replicates with largest CS <= %d: SuSiE(0.5) %.3f, A_capped %.3f, B_uncapped %.3f\n",
              cap, pw("SuSiE(0.5)", le), pw("BLiP A_capped", le), pw("BLiP B_uncapped", le)))
  cat(sprintf("Power on replicates with largest CS  > %d: SuSiE(0.5) %.3f, A_capped %.3f, B_uncapped %.3f\n",
              cap, pw("SuSiE(0.5)", gt), pw("BLiP A_capped", gt), pw("BLiP B_uncapped", gt)))

  cat("\n")
  for (nm in names(checks)) cat(sprintf("  [%s] %s\n", if (isTRUE(checks[[nm]])) "PASS" else "FAIL", nm))
  if (!all(unlist(checks))) {
    cat("\nSMOKE TEST FAILED\n")
    quit(status = 1)
  }
  cat("\nSMOKE TEST PASSED\n")
}

if (sys.nframe() == 0L) main()
