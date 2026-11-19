work <- commandArgs(trailingOnly = TRUE)[1]
if (as.character(packageVersion("susieR")) != "0.16.6") {
  stop("Install the pinned susieR 0.16.6 build described in README.md.")
}

data <- new.env()
load(file.path(work, "rss_mismatch_cases.RData"), envir = data)
R <- unname(cor(data$rss_mismatch_cases$ld_mismatch_sketch$X))
stopifnot(ncol(R) == 121L)

reporting <- new.env(parent = asNamespace("susieR"))
sys.source(file.path(work, "susie_get_functions.R"), envir = reporting)

eig <- eigen(R, symmetric = TRUE)
root <- sweep(eig$vectors, 2, sqrt(pmax(eig$values, 0)), "*")
settings <- expand.grid(strength = c(24L, 80L), K = 1:2, rep = seq_len(150))
posteriors <- sets <- vector("list", nrow(settings))
set.seed(83011)

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

write.csv(do.call(rbind, posteriors), file.path(work, "posteriors.csv"), row.names = FALSE)
write.csv(do.call(rbind, sets), file.path(work, "susier-sets.csv"), row.names = FALSE)
