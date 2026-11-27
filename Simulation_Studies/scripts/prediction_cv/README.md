# Prediction CV benchmark (Supplementary Figure S4)

Five-fold cross-validated gene expression prediction with SuSiE, SuSiE-inf and
SuSiE-ash on simulated architectures: eight parameter sets (five sparse, three
complex), 300 replicates each, n = 1000. The complex scenarios have h2g = 0.25
and K = 3 sparse effects.

| Scenario | S4 panels | Architecture |
|---|---|---|
| `sparse` | A, E | K = 1-5 causal variants (one parameter set per K), h2 = 0.03 per variant |
| `complex` | B, F | polygenic background: h2g shares 0.5 sparse / 0.35 oligogenic / 0.15 infinitesimal, 5 oligogenic effects, infinitesimal effects on 15 variants |
| `complex_s1` | C, G | moderate infinitesimal background: as `complex`, with infinitesimal effects on all remaining variants |
| `complex_s2` | D, H | extensive infinitesimal background: shares 0.5 / 0.15 / 0.35, 10 oligogenic effects, infinitesimal effects on all remaining variants |

## Files

| File | Role |
|---|---|
| `command_generator_parallel_cv.R` | writes one worker command per replicate chunk and one combine command per parameter set |
| `prediction_cv_script_parallel.R` | worker: simulates, splits, fits and scores the replicates of one chunk |
| `combine_parallel_cv_results.R` | merges the per-chunk files of one parameter set |
| `local/run_prediction_cv_local.sh` | runs the same scripts on one machine with `xargs -P`, resumable, then rebuilds the figure |
| `local/run_one_job.sh` | runs one worker command for the two scripts beside it and logs the outcome |
| `local/run_smoke_test.sh` | one replicate of two scenarios, end to end, in about 5 minutes |

The figure is built from the combined files by `extract_prediction_data.R`
(per-replicate mean R² and per-replicate winner) and `plot_S4.R` in
`Supplementary_Figures/S4/`.

## Design, per replicate i

1. LD block: the i-th file (recycled) of the block list sorted in byte order
   (`sort(..., method = "radix")`, the C-locale order of the cluster; a
   locale-dependent sort puts `LD_chr1_*` before `LD_chr10_*` and changes every
   replicate). First 1000 samples; mean-impute missing genotypes; keep variants
   with MAF ≥ 0.01 and non-zero variance; keep the first 5000 of them.
2. Phenotype: `generate_data()` with seed `i + 1000` on the standardized
   matrix. Sparse: K causal variants with pairwise |r| ≤ 0.1, unit effects,
   residual variance set for total h2 = 0.03 K. Complex: sparse, oligogenic and
   infinitesimal components as described in the supplement.
3. Folds: `set.seed(i + 1000); sample(rep(1:5, length.out = n))`.
4. For each fold, on the training rows: `susie(X, y, L = 10, max_iter = 100)`;
   `susie(..., unmappable_effects = "inf", convergence_method = "pip")`;
   `susie(..., unmappable_effects = "ash", convergence_method = "pip")`.
   Columns constant within the training fold are dropped for fitting and get a
   zero coefficient at prediction.
5. Prediction on the held-out rows from `coef(fit)` (intercept plus per-variant
   effects, including the unmappable component). Pearson R² = cor(y, ŷ)² and
   RMSE per fold; the replicate's value is the mean over folds.

Steps 1 and 2 are those of the fine-mapping benchmark
(`../fine_mapping/simulation_script_parallel.R`), so replicate i of a scenario
has the same LD block, causal indices and effect sizes in both.

The models never see test-fold statistics. The worker hands the unscaled
filtered genotypes to the CV loop; `susie()` centres and standardizes the
training fold internally (for the unmappable methods standardization is forced)
and `coef()` returns effects on the input scale. Standardizing the whole matrix
first gives the same predictions, because a linear predictor with an intercept
is invariant to per-column affine transforms applied identically to both folds;
fitting on the unscaled training fold makes the separation explicit. The MAF
filter, the imputation and the 5000-variant cap are applied to all 1000 samples
before the phenotype is simulated: they define the simulated genotype matrix,
use no phenotype information, and are not redone within folds.

SuSiE-inf and SuSiE-ash run with `convergence_method = "pip"` (the package
default is `"elbo"`); pass `convergence_method=elbo` to the generator to change
it. The worker records the susieR version and the convergence setting in every
output file.

## Running on the cluster

```
Rscript command_generator_parallel_cv.R [key=value ...]
```

This writes `commands_cv_to_submit.txt` and `commands_cv_combine.txt`. Submit
every line of the first as a job, then run the lines of the second. A worker
job needs one core; a 2-replicate job takes about 7 minutes with the current
susieR. The defaults are the cluster layout
(`scripts_dir=/home/apm2217/data/benchmark`, `ld_blocks_dir=oligo_LD_blocks`,
`output_dir=/home/apm2217/output`) with `total_replicates=300` and
`replicates_per_job=2`; the script header lists every key.

A worker job writes `<output_dir>/<base>_rep<start>-<end>.rds` and the combiner
`<output_dir>/<base>.rds`, where `<base>` starts with
`pred_cv_nrep300_sparse_h2persnp0.03_K<K>` or `pred_cv_nrep300_oligo_h2g0.25_K3`
and ends with `_n1000_nfolds5`. `extract_prediction_data.R` reads the five
combined sparse files from
`$SUSIER2_BENCHMARK_ROOT/final_scripts/prediction_cv_results/sparse` and the
three complex ones from `.../complex`.

## Running locally

```
bash ../fine_mapping/blip_rerun/local/setup_local_env.sh   # once; builds the development susieR
export SUSIER2_BENCHMARK_ROOT=/path/to/susieR2.0-benchmark   # holds final_scripts/500_rep_results/ld_blocks
bash local/run_smoke_test.sh                                 # optional
caffeinate -ims bash local/run_prediction_cv_local.sh
```

Output goes to `$SUSIER2_BENCHMARK_ROOT/final_scripts/prediction_cv_results/{sparse,complex}`,
the layout `extract_prediction_data.R` reads, with the command files in
`commands/` and the logs in `logs/` beside them. The study costs about 150
CPU-hours, and a worker peaks at about 4.3 GB resident (SuSiE-ash on
p = 5000), so the driver bounds the number of parallel jobs by memory as well
as cores: on a 48 GB laptop that is 8 jobs and about 19 hours. `N_JOBS`
overrides it; `STEPS="combine figure"` rebuilds only the figure, from existing
results and the command files of an earlier run. The driver's header lists
every environment variable.

What protects an unattended run:

- A worker job handles one replicate (local default) and writes its result
  file through a temporary name and a rename, so a crash never leaves a
  truncated file. A job whose output exists exits at once, so re-running the
  same command resumes and loses at most one replicate per worker.
- Every fit is wrapped: a failing `susie()` call or prediction gives NA for
  that method and fold, with the error message stored in the fold result and
  `n_fit_errors` in the replicate; the job continues. The combine step prints,
  per combined file, how many replicates are present, how many fits failed and
  how many mean R² values are non-finite.
- `logs/progress.log` gets one line per finished job, so `wc -l` counts
  progress (a re-run adds a line for every job it skips). Failed jobs are
  listed in `logs/failed_commands.txt` and their logs kept as `logs/job.*`; the
  driver stops before combining if any job failed.
- `caffeinate -ims` keeps the Mac awake on power with the lid open; closing
  the lid still sleeps it (the jobs pause and resume, they do not fail).
