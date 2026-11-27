#!/usr/bin/env Rscript

# Supplementary Figure S7, SuSiE+BLiP re-run, step 3: per-replicate metrics and
# summaries for the SuSiE credible sets and every BLiP arm.
#
# BLiP at level q is paired with SuSiE credible sets of coverage 1 - q
# (q = 0.05 with 95%, q = 0.10 with 90%). Methods compared at each q:
#   SuSiE(0.5), SuSiE(0.8)  sets with purity (min |r|) >= 0.5 / 0.8
#   SuSiE(no purity)        every set reaching the coverage
#   SuSiE(0.5) size<=cap    SuSiE(0.5) sets of at most size_cap variants
#   SuSiE capped(cap)       a component is kept only when its top size_cap
#                           alphas exceed 1 - q, no purity filter (the SuSiE
#                           comparator in the BLiP authors' UK Biobank analysis)
#   BLiP <arm>              every arm present in blip_arms.json
# Power and FDR are computed as in extract_S7_data.R. Resolution-adjusted power
# is the sum of 1/|G| over the sets G containing a causal SNP, divided by K.
# Replicates are stratified by whether their largest SuSiE(0.5) set at the
# paired coverage exceeds size_cap.
#
# Input:  <input_dir>/K<k>/rep_<id>/{meta.json, alpha.npy, corr_sub.npy}
#           (export_blip_inputs.R) and blip_arms.json (run_blip_arms.py);
#           replicates without both JSON files are skipped
# Output: <output_dir>/blip_metrics.rds, a list with
#           rep_metrics  one row per replicate x q x method: power, FDR,
#                        resolution-adjusted power, FDP given >= 1 set, number
#                        and mean size of sets, stratum
#           group_rows   one row per reported set: size, purity (min |r| from
#                        corr_sub.npy), whether it contains a causal SNP
#           summary      both tables aggregated per q x method x K (each K and
#                        "pooled") x stratum (each and "all")
#           cum_true     per summary block, the number of true sets of size
#                        <= k for k = 1..p, and that number over the causal SNPs
#           meta         settings, versions, skipped and failed folders
#         <output_dir>/blip_metrics_summary.csv, blip_rep_metrics.csv
#
# Usage: Rscript blip_metrics.R input_dir=<export root> \
#          [output_dir=<input_dir>/metrics] [k_values=<all found>] \
#          [q_values=0.05,0.1] [size_cap=25] [n_cores=1]

suppressMessages({
  library(jsonlite)
  library(dplyr)
  library(parallel)
})

.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.script_dir, "blip_common.R"))
locate_paper_paths(.script_dir)

read_json_lists <- function(path) {
  fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
}

# ---- Method lists for one replicate at one q ----

rep_methods <- function(meta, arms, alpha, q, size_cap) {
  cov_key <- sprintf("%g", 1 - q)
  entries <- meta$susie_cs[[cov_key]]
  if (is.null(entries)) stop("meta.json has no SuSiE CS list for coverage ", cov_key)
  cs    <- lapply(entries, function(e) as.integer(e$cs))
  pur   <- vapply(entries, function(e) as.numeric(e$min_abs_corr), numeric(1))
  sizes <- lengths(cs)

  methods <- list()
  methods[["SuSiE(0.5)"]]       <- cs[pur >= 0.5]
  methods[["SuSiE(0.8)"]]       <- cs[pur >= 0.8]
  methods[["SuSiE(no purity)"]] <- cs
  methods[[sprintf("SuSiE(0.5) size<=%d", size_cap)]] <- cs[pur >= 0.5 & sizes <= size_cap]
  methods[[sprintf("SuSiE capped(%d)", size_cap)]]    <- susie_capped_sets(alpha, q, size_cap)

  qk <- sprintf("q%g", q)
  for (arm in names(arms$arms)) {
    res <- arms$arms[[arm]][[qk]]
    if (is.null(res)) next
    methods[[paste0("BLiP ", arm)]] <- lapply(res$groups, as.integer)
  }
  methods
}

# ---- Per-replicate rows ----

process_rep <- function(rep_dir, q_values, size_cap) {
  meta  <- read_json_lists(file.path(rep_dir, "meta.json"))
  arms  <- read_json_lists(file.path(rep_dir, "blip_arms.json"))
  alpha <- read_npy(file.path(rep_dir, "alpha.npy"))
  corr_sub   <- read_npy(file.path(rep_dir, "corr_sub.npy"))
  corr_index <- as.integer(meta$corr_index)
  causal     <- as.integer(meta$causal)
  k      <- as.integer(meta$K)
  rep_id <- as.integer(meta$replicate_id)

  # corr_sub.npy omits members of impure SuSiE sets; a reported group that is
  # exactly such a set takes the purity susie_get_cs() exported for it.
  exported_purity <- list()
  for (entries in meta$susie_cs) for (e in entries) {
    exported_purity[[paste(sort(as.integer(e$cs)), collapse = ",")]] <- as.numeric(e$min_abs_corr)
  }
  purity_of <- function(groups) {
    if (length(groups) == 0) return(numeric(0))
    vapply(groups, function(g) {
      m <- match(g, corr_index)
      if (!anyNA(m)) return(compute_purity_per_cs(list(m), corr_sub))
      p <- exported_purity[[paste(sort(g), collapse = ",")]]
      if (is.null(p)) NA_real_ else p
    }, numeric(1))
  }

  rep_rows <- list()
  grp_rows <- list()
  for (q in q_values) {
    methods  <- rep_methods(meta, arms, alpha, q, size_cap)
    susie05  <- methods[["SuSiE(0.5)"]]
    max_cs05 <- if (length(susie05) > 0) max(lengths(susie05)) else 0L
    stratum  <- if (max_cs05 > size_cap) sprintf("max_cs>%d", size_cap) else sprintf("max_cs<=%d", size_cap)

    for (m in names(methods)) {
      cs      <- methods[[m]]
      met     <- calculate_cs_metrics(cs, causal)
      is_true <- vapply(cs, function(g) any(g %in% causal), logical(1))
      sizes   <- lengths(cs)

      rep_rows[[length(rep_rows) + 1]] <- data.frame(
        K = k, replicate_id = rep_id, q = q, coverage = 1 - q, method = m,
        n_cs      = length(cs),
        power     = met$power,
        fdr       = met$fdr,
        res_power = if (length(cs) > 0) sum(1 / sizes[is_true]) / length(causal) else 0,
        fdp_cond  = if (length(cs) > 0) met$fdr else NA_real_,
        size_mean = if (length(cs) > 0) met$size else NA_real_,
        max_susie05_cs = max_cs05,
        stratum   = stratum,
        n_retained = as.integer(arms$n_retained),
        stringsAsFactors = FALSE)

      if (length(cs) > 0) {
        grp_rows[[length(grp_rows) + 1]] <- data.frame(
          K = k, replicate_id = rep_id, q = q, method = m,
          size = sizes, purity = purity_of(cs), is_true = is_true,
          stratum = stratum, stringsAsFactors = FALSE)
      }
    }
  }

  list(rep_rows = do.call(rbind, rep_rows),
       grp_rows = if (length(grp_rows) > 0) do.call(rbind, grp_rows) else NULL,
       p = as.integer(meta$p), settings = arms$settings, versions = arms$versions)
}

# ---- Summaries ----

expand_blocks <- function(df) {
  df <- bind_rows(mutate(df, K_label = as.character(K)), mutate(df, K_label = "pooled"))
  bind_rows(mutate(df, stratum_label = stratum), mutate(df, stratum_label = "all"))
}

se <- function(x) if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_

build_summary <- function(rep_metrics, group_rows, size_cap) {
  by <- c("q", "method", "K_label", "stratum_label")
  summary_rep <- expand_blocks(rep_metrics) %>%
    group_by(across(all_of(by))) %>%
    summarise(
      n_reps         = n(),
      power_mean     = mean(power),
      power_se       = se(power),
      fdr_mean       = mean(fdr),
      fdr_se         = se(fdr),
      res_power_mean = mean(res_power),
      res_power_se   = se(res_power),
      fdp_cond_mean  = if (any(!is.na(fdp_cond))) mean(fdp_cond, na.rm = TRUE) else NA_real_,
      n_cs_mean      = mean(n_cs),
      zero_disc_frac = mean(n_cs == 0),
      .groups = "drop")

  if (nrow(group_rows) > 0) {
    summary_grp <- expand_blocks(group_rows) %>%
      group_by(across(all_of(by))) %>%
      summarise(
        n_groups          = n(),
        size_median       = median(size),
        size_max          = max(size),
        frac_size_gt_cap  = mean(size > size_cap),
        purity_median     = if (any(!is.na(purity))) median(purity, na.rm = TRUE) else NA_real_,
        frac_purity_lt05  = if (any(!is.na(purity))) mean(purity < 0.5, na.rm = TRUE) else NA_real_,
        .groups = "drop")
    summary_rep <- left_join(summary_rep, summary_grp, by = by)
  }
  summary_rep
}

build_cum_true <- function(rep_metrics, group_rows, summary, p_max) {
  key <- function(d) paste(d$q, d$method, d$K_label, d$stratum_label, sep = "|")
  blocks <- distinct(summary, q, method, K_label, stratum_label)
  gt <- group_rows[group_rows$is_true, , drop = FALSE]
  gt_split <- if (nrow(gt) > 0) {
    gt_exp <- expand_blocks(gt)
    split(gt_exp$size, key(gt_exp))
  } else list()
  causal_by_block <- expand_blocks(distinct(rep_metrics, K, replicate_id, q, stratum)) %>%
    group_by(q, K_label, stratum_label) %>%
    summarise(n_causal = sum(K), .groups = "drop")

  curves <- lapply(seq_len(nrow(blocks)), function(i) {
    b <- blocks[i, ]
    sizes <- gt_split[[key(b)]]
    if (is.null(sizes)) sizes <- integer(0)
    data.frame(q = b$q, method = b$method, K_label = b$K_label,
               stratum_label = b$stratum_label, k = seq_len(p_max),
               cum_true = cumsum(tabulate(pmin(sizes, p_max), nbins = p_max)),
               stringsAsFactors = FALSE)
  })
  bind_rows(curves) %>%
    left_join(causal_by_block, by = c("q", "K_label", "stratum_label")) %>%
    mutate(cum_power = cum_true / n_causal)
}

# ---- Main ----

main <- function() {
  a <- parse_key_value_args(list(
    input_dir  = NULL,
    output_dir = "",
    k_values   = "",
    q_values   = "0.05,0.1",
    size_cap   = 25L,
    n_cores    = 1L
  ))
  if (is.null(a$input_dir)) stop("input_dir must be specified", call. = FALSE)
  input_dir  <- path.expand(a$input_dir)
  output_dir <- if (nzchar(a$output_dir)) path.expand(a$output_dir) else file.path(input_dir, "metrics")
  q_values   <- as.numeric(strsplit(a$q_values, ",")[[1]])
  size_cap   <- a$size_cap
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================\n")
  cat("BLiP re-run: metrics\n")
  cat("========================================\n")
  cat("Input dir:  ", input_dir, "\n")
  cat("Output dir: ", output_dir, "\n")
  cat("q values:   ", paste(q_values, collapse = ", "), "   size cap:", size_cap, "\n")

  k_dirs <- list.files(input_dir, pattern = "^K[0-9]+$", full.names = TRUE)
  if (nzchar(a$k_values)) {
    keep <- sprintf("K%s", trimws(strsplit(a$k_values, ",")[[1]]))
    k_dirs <- k_dirs[basename(k_dirs) %in% keep]
  }
  rep_dirs <- unlist(lapply(k_dirs, function(d) list.files(d, pattern = "^rep_[0-9]+$", full.names = TRUE)))
  if (length(rep_dirs) == 0) stop("No replicate folders under ", input_dir, call. = FALSE)
  has_both <- file.exists(file.path(rep_dirs, "meta.json")) &
              file.exists(file.path(rep_dirs, "blip_arms.json"))
  skipped  <- rep_dirs[!has_both]
  rep_dirs <- rep_dirs[has_both]
  cat(sprintf("Replicates with arms output: %d (skipping %d without)\n",
              length(rep_dirs), length(skipped)))

  job_start <- Sys.time()
  results <- mclapply(rep_dirs, function(d) {
    tryCatch(process_rep(d, q_values, size_cap),
             error = function(e) list(error = conditionMessage(e), rep_dir = d))
  }, mc.cores = a$n_cores)

  failed <- Filter(function(r) !is.null(r$error), results)
  for (f in failed) cat(sprintf("  ERROR %s: %s\n", f$rep_dir, f$error))
  ok <- Filter(function(r) is.null(r$error), results)
  if (length(ok) == 0) stop("No replicate could be processed", call. = FALSE)

  rep_metrics <- bind_rows(lapply(ok, `[[`, "rep_rows"))
  group_rows  <- bind_rows(lapply(ok, `[[`, "grp_rows"))
  if (nrow(group_rows) == 0) {
    group_rows <- data.frame(K = integer(0), replicate_id = integer(0), q = numeric(0),
                             method = character(0), size = integer(0), purity = numeric(0),
                             is_true = logical(0), stratum = character(0))
  }
  p_max <- max(vapply(ok, function(r) r$p, integer(1)))

  summary  <- build_summary(rep_metrics, group_rows, size_cap)
  cum_true <- build_cum_true(rep_metrics, group_rows, summary, p_max)

  method_levels <- unique(rep_metrics$method)
  rep_metrics$method <- factor(rep_metrics$method, levels = method_levels)
  group_rows$method  <- factor(group_rows$method,  levels = method_levels)
  summary$method     <- factor(summary$method,     levels = method_levels)
  cum_true$method    <- factor(cum_true$method,    levels = method_levels)

  out <- list(
    rep_metrics = as.data.frame(rep_metrics),
    group_rows  = as.data.frame(group_rows),
    summary     = as.data.frame(summary),
    cum_true    = as.data.frame(cum_true),
    meta = list(
      input_dir       = input_dir,
      size_cap        = size_cap,
      q_values        = q_values,
      K_values        = sort(unique(rep_metrics$K)),
      n_replicates    = nrow(distinct(rep_metrics, K, replicate_id)),
      p_max           = p_max,
      method_levels   = method_levels,
      blip_settings   = ok[[1]]$settings,
      blip_versions   = ok[[1]]$versions,
      skipped_dirs    = skipped,
      failed_dirs     = vapply(failed, `[[`, "", "rep_dir"),
      elapsed_seconds = as.numeric(difftime(Sys.time(), job_start, units = "secs")),
      timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  )
  saveRDS(out, file.path(output_dir, "blip_metrics.rds"))
  write.csv(out$summary,     file.path(output_dir, "blip_metrics_summary.csv"), row.names = FALSE)
  write.csv(out$rep_metrics, file.path(output_dir, "blip_rep_metrics.csv"),     row.names = FALSE)

  show <- summary %>%
    filter(K_label == "pooled") %>%
    arrange(q, stratum_label, method) %>%
    select(q, stratum_label, method, n_reps, power_mean, res_power_mean, fdr_mean,
           n_cs_mean, zero_disc_frac, size_median, size_max, frac_size_gt_cap,
           purity_median, frac_purity_lt05)
  cat("\nPooled over K (power/FDR as in Figure S7; res_power = sum 1/|G| over true groups / K):\n")
  print(as.data.frame(show), digits = 3, row.names = FALSE)

  cat(sprintf("\nSaved %s (%d replicates, %d methods) in %.1f s.\n",
              file.path(output_dir, "blip_metrics.rds"), out$meta$n_replicates,
              length(method_levels), out$meta$elapsed_seconds))
  if (length(failed) > 0) quit(status = 1)
}

if (sys.nframe() == 0L) main()
