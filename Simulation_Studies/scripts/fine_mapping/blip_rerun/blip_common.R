# Supplementary Figure S7, SuSiE+BLiP re-run: helpers shared by
# export_blip_inputs.R, blip_metrics.R, combine_blip_to_s7.R,
# command_generator_blip.R and the test scripts. Sourced, never run on its own.
#
# reconstruct_genotype_matrix(), calculate_cs_metrics() and
# compute_purity_per_cs() are copies of the functions of the same name in
# Supplementary_Figures/S7/extract_S7_data.R and must stay identical to them,
# so that the re-run computes its numbers exactly as the figure does.

# ---- Paths and arguments ----

# Walk upwards from `start` to the repository root and source R/paths.R there.
# Returns the root, or NULL when the scripts are run from a copy outside the
# repository (then only absolute paths work).
locate_paper_paths <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, ".susieR2-paper-root"))) {
      source(file.path(d, "R", "paths.R"))
      return(invisible(d))
    }
    parent <- dirname(d)
    if (identical(parent, d)) return(invisible(NULL))
    d <- parent
  }
}

# Parse `key=value` arguments. `defaults` names the accepted keys and gives
# each one's type (character, integer, numeric or logical); NULL means a
# required character argument. Quotes around values are stripped, and the
# value may itself contain '='.
parse_key_value_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  out <- defaults
  for (arg in args) {
    kv <- regmatches(arg, regexpr("=", arg), invert = TRUE)[[1]]
    if (length(kv) != 2) stop("Arguments must be key=value, got: ", arg, call. = FALSE)
    key   <- kv[1]
    value <- gsub("^'|'$", "", gsub('^"|"$', "", kv[2]))
    if (!key %in% names(defaults)) {
      stop("Unknown argument '", key, "'. Accepted: ",
           paste(names(defaults), collapse = ", "), call. = FALSE)
    }
    d <- defaults[[key]]
    out[[key]] <- if (is.null(d) || is.character(d)) {
      value
    } else if (is.logical(d)) {
      as.logical(value)
    } else if (is.integer(d)) {
      as.integer(value)
    } else {
      as.numeric(value)
    }
  }
  out
}

# Resolve a path given on the command line: absolute (or ~) paths and paths
# that exist relative to the working directory are used as given; anything
# else is taken relative to <benchmark_root()>/final_scripts/500_rep_results,
# the directory extract_S7_data.R reads the raw simulation output from.
resolve_benchmark_path <- function(p) {
  if (grepl("^(/|~)", p)) return(path.expand(p))
  if (file.exists(p)) return(normalizePath(p))
  if (!exists("benchmark_root", mode = "function")) {
    stop("Relative path '", p, "' given but R/paths.R could not be located; ",
         "use an absolute path.", call. = FALSE)
  }
  file.path(benchmark_root(), "final_scripts", "500_rep_results", p)
}

# ---- NumPy .npy files (format 1.0, C order, little-endian) ----

write_npy <- function(x, path, dtype = c("float64", "float32", "int32")) {
  dtype <- match.arg(dtype)
  descr <- switch(dtype, float64 = "<f8", float32 = "<f4", int32 = "<i4")
  shape <- if (is.matrix(x)) dim(x) else length(x)
  shape_str <- if (length(shape) == 1) paste0(shape, ",") else paste(shape, collapse = ", ")
  header <- sprintf("{'descr': '%s', 'fortran_order': False, 'shape': (%s), }",
                    descr, shape_str)
  pad <- (64L - ((10L + nchar(header) + 1L) %% 64L)) %% 64L
  header <- paste0(header, strrep(" ", pad), "\n")

  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(as.raw(c(0x93, charToRaw("NUMPY"), 0x01, 0x00)), con)
  writeBin(nchar(header), con, size = 2, endian = "little")
  writeBin(charToRaw(header), con)
  vals <- if (is.matrix(x)) as.vector(t(x)) else as.vector(x)
  if (dtype == "int32") {
    writeBin(as.integer(vals), con, size = 4, endian = "little")
  } else {
    writeBin(as.numeric(vals), con, size = if (dtype == "float64") 8 else 4,
             endian = "little")
  }
  invisible(path)
}

read_npy <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 6)
  if (!identical(magic, as.raw(c(0x93, charToRaw("NUMPY"))))) stop("Not a .npy file: ", path)
  version <- readBin(con, "raw", 2)
  hlen <- if (as.integer(version[1]) == 1L) {
    readBin(con, "integer", 1, size = 2, endian = "little", signed = FALSE)
  } else {
    readBin(con, "integer", 1, size = 4, endian = "little")
  }
  header  <- rawToChar(readBin(con, "raw", hlen))
  descr   <- sub(".*'descr': *'([^']+)'.*", "\\1", header)
  fortran <- grepl("'fortran_order': *True", header)
  shape   <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", header)
  shape   <- as.integer(strsplit(gsub("[^0-9,]", "", shape), ",")[[1]])
  spec <- switch(descr,
                 "<f8" = list(what = "double",  size = 8L),
                 "<f4" = list(what = "double",  size = 4L),
                 "<i4" = list(what = "integer", size = 4L),
                 stop("Unsupported .npy dtype '", descr, "' in ", path))
  v <- readBin(con, spec$what, prod(shape), size = spec$size, endian = "little")
  if (length(shape) == 2) {
    if (fortran) matrix(v, shape[1], shape[2]) else matrix(v, shape[1], shape[2], byrow = TRUE)
  } else {
    v
  }
}

# ---- Copies of functions in Supplementary_Figures/S7/extract_S7_data.R ----

reconstruct_genotype_matrix <- function(ld_block_path, n = NULL) {
  ld_block <- readRDS(ld_block_path)
  X <- ld_block$genotypes

  if (!is.null(n) && n < nrow(X)) {
    X <- X[1:n, , drop = FALSE]
  }

  col_means <- colMeans(X, na.rm = TRUE)
  na_idx <- which(is.na(X), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    X[na_idx] <- col_means[na_idx[, 2]]
  }

  maf <- colMeans(X) / 2
  maf <- pmin(maf, 1 - maf)
  X <- X[, maf >= 0.01, drop = FALSE]

  col_vars <- apply(X, 2, var)
  X <- X[, col_vars > 0, drop = FALSE]

  if (ncol(X) > 5000) {
    X <- X[, 1:5000, drop = FALSE]
  }

  X <- scale(X)
  return(X)
}

calculate_cs_metrics <- function(test.cs, causal_indices) {
  cs_size <- 0
  cs_fdr <- 0
  cs_power <- 0

  if (length(test.cs) > 0) {
    cs_size <- length(unlist(test.cs)) / length(test.cs)
    TP_fdr <- sum(sapply(test.cs, function(cs) any(cs %in% causal_indices)))
    FP_fdr <- length(test.cs) - TP_fdr
    cs_fdr <- if ((TP_fdr + FP_fdr) > 0) FP_fdr / (TP_fdr + FP_fdr) else 0
    TP_power <- sum(causal_indices %in% unlist(test.cs))
    FN_power <- length(causal_indices) - TP_power
    cs_power <- TP_power / (TP_power + FN_power)
  }

  return(list(power = cs_power, fdr = cs_fdr, size = cs_size))
}

compute_purity_per_cs <- function(cs_list, R) {
  if (length(cs_list) == 0) return(numeric(0))
  vapply(cs_list, function(idx) {
    if (length(idx) <= 1) return(1.0)
    sub <- abs(R[idx, idx, drop = FALSE])
    diag(sub) <- NA_real_
    min(sub, na.rm = TRUE)
  }, numeric(1))
}

# ---- SuSiE sets under BLiP's size cap ----

# SuSiE's own per-component sets, as formed in the BLiP authors' UK Biobank
# wrapper: sort a component's alphas, cumulate the top `max_size`, keep the
# component when that sum exceeds 1 - q, and report the smallest top-k set
# reaching it. No purity filter; identical sets are collapsed as in
# susie_get_cs(dedup = TRUE). Indices are 1-based.
susie_capped_sets <- function(alpha, q, max_size) {
  sets <- list()
  for (l in seq_len(nrow(alpha))) {
    o <- order(alpha[l, ], decreasing = TRUE)[seq_len(min(max_size, ncol(alpha)))]
    cums <- cumsum(alpha[l, o])
    hit <- which(cums > 1 - q)
    if (length(hit) == 0) next
    sets[[length(sets) + 1]] <- sort(o[seq_len(hit[1])])
  }
  unique(sets)
}

# The arm labels run_blip_arms.py writes.
blip_arm_names <- function() {
  c("A_capped", "B_uncapped",
    "C_hier_capped", "C_hier_uncapped",
    "D_purity_capped", "D_purity_uncapped", "D_purity08_uncapped")
}
