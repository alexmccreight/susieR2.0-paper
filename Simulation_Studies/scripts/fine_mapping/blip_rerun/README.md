# SuSiE+BLiP re-run for Supplementary Figure S7

## What this is

SuSiE credible sets compared with SuSiE+BLiP on the sparse benchmark
(n = 1000, p = 5000 variants from real LD blocks, K = 1..5 causal variants,
500 replicates per K). BLiP runs in the BLiP authors' Python package `pyblip`
0.4.2 and follows their UK Biobank recipe for turning a SuSiE fit into
candidate groups (`blip_wrapper.py` in github.com/amspector100/ukbb_blip).

On this benchmark 36-46% of the SuSiE credible sets exceed 25 variants, the
authors' group size cap. The size cap, the purity filter and the
candidate-group family are therefore explicit arms, the SuSiE side is reported
under the same restrictions, and the authors' resolution-adjusted power is
given next to the unweighted power of Figure S7.

| Step | Script | Role |
|---|---|---|
| 1 | `export_blip_inputs.R` | rebuild X and the SuSiE credible sets of each replicate |
| 2 | `run_blip_arms.py` | run the BLiP arms on each replicate |
| 3 | `blip_metrics.R` | per-replicate metrics, stratified summaries, cumulative resolution curves |
| 4 | `combine_blip_to_s7.R` | one arm in the format `Supplementary_Figures/S7/extract_S7_data.R` reads |

`blip_common.R` holds the shared helpers. Two routes run the four steps:
`local/run_blip_local.sh` on one machine and `command_generator_blip.R` on the
HPC (generator -> workers -> combiner, as `command_generator_parallel.R`).
Every replicate is seeded, so both write the same per-replicate files.

## What you need

Environments, built once (needs network access):

```
[SUSIER_ROOT=/path/to/susieR] bash local/setup_local_env.sh [env dir]   # default $HOME/.susieR2_blip_env
```

- `<env dir>/venv`: Python with `pyblip==0.4.2`, cvxopt (the GLPK and GLPK_MI
  solvers for cvxpy), cvxpy, numpy, scipy, networkx and tqdm. `uv` creates it
  with Python `PY_VERSION` (default 3.12); without `uv`, `python3 -m venv` does.
- `<env dir>/Rlib`: susieR built from the development checkout at
  `SUSIER_ROOT` (default `../susieR`, as in `R/paths.R`), never from CRAN; its
  dependencies, jsonlite and dplyr come from CRAN when missing.

Raw inputs, under `<benchmark>/final_scripts/500_rep_results/`:

- `sparse/sim_nrep500_sparse_h2persnp0.03_K{1..5}_n1000.rds`; each replicate
  needs `fits$susie$fit$alpha`, `fits$susie$fit$V`, `causal_indices` and
  `ld_block_name`
- `ld_blocks/<ld_block_name>.rds` for every block a replicate names
  (`ld_block$genotypes` is the raw matrix)

`<benchmark>` is `SUSIER2_BENCHMARK_ROOT` (default `../susieR2.0-benchmark`,
see `R/paths.R`). The data are not in the repository; copies are on the HPC
under `/home/apm2217/data/benchmark` and on the old machine under
`/Users/alexmccreight/StatFunGen/susieR2.0-benchmark`.

## How to run

### Local route

```
export SUSIER2_BENCHMARK_ROOT=/path/to/susieR2.0-benchmark
caffeinate -i bash local/run_blip_local.sh      # caffeinate keeps a Mac awake
```

For each K the script runs step 1 (simulation file loaded once, replicates
forked over `N_CORES`) and step 2 (process pool of `N_CORES`), then step 3,
then step 4 for `A_capped` and `B_uncapped`. Output goes to
`<benchmark>/final_scripts/500_rep_results/blip_rerun`, logs to `logs/` inside
it. Finished replicates are skipped, so an interrupted run resumes where it
stopped. `STEPS`, `K_VALUES`, `N_CORES`, `OUT_DIR`, `ARMS_ARGS`, `CONVERT_ARMS`
and the other variables in the script header override the defaults: for
example, `STEPS="arms metrics" ARMS_ARGS="--prenarrow"` reruns the arms with
prenarrowing, overwriting every `blip_arms.json` unless `OUT_DIR` points at a
separate export.

Rough cost on 16 workers: seconds per replicate for the export (dominated by
`cor(X)` on 1000 x 5000) and per replicate and arm for the linear program.
Peak memory is one simulation file per K plus about 0.3 GB (X and R) per
worker.

### HPC route

1. Copy this folder next to the data, e.g.
   `/home/apm2217/data/benchmark/blip_rerun`. `R/paths.R` is not found from
   there, which is fine because the command files use absolute paths.
2. Once, on the login node:
   `SUSIER_ROOT=/path/to/susieR bash local/setup_local_env.sh /home/apm2217/blip_env`.
   The generated commands call plain `Rscript`, so jobs need
   `R_LIBS_USER=/home/apm2217/blip_env/Rlib` in their environment.
3. `Rscript command_generator_blip.R`, with `key=value` arguments for anything
   that differs from the defaults in its header (`python=`, `output_dir=`,
   ...; 500 replicates per K in jobs of 25). It writes
   `commands_blip_export.txt`, `commands_blip_arms.txt` and
   `commands_blip_combine.txt`.
4. Submit every line of the export file; when all have finished, every line of
   the arms file; when all have finished, run the combine file in order. Each
   export job loads its simulation file (a few GB) once; each arms job needs
   one core and modest memory.
5. Bring back `metrics/` and the `s7_<arm>/` folders (small).

### Smoke test

```
PYTHON=<venv python> R_LIBS_USER=<private R lib> bash tests/run_smoke_test.sh [work dir]
```

Simulates 20 replicates (n = 800, p = 1500, K = 2, h2 = 0.03 per SNP; a
Gaussian AR(1) block with rho = 0.99 stands in for genotypes, so that several
SuSiE credible sets exceed 25 variants), runs both routes (7 replicates per
HPC job) and checks that they wrote the same replicate folders and that, on
the replicates whose largest set exceeds 25, `A_capped` has lower power than
`SuSiE(0.5)` and `B_uncapped` higher power than `A_capped`. Both variables
default to the environment in `$HOME/.susieR2_blip_env`.

### Single steps

Each script also runs by itself; its header gives the usage line with every
argument and default (`python run_blip_arms.py --help` for step 2). Step 1
skips a replicate that has `meta.json`, whatever its settings, unless
`overwrite=TRUE`; step 2 skips one whose `blip_arms.json` carries the same
settings. A relative `sim_file` or `ld_blocks_dir` that does not exist from
the working directory resolves against
`<benchmark>/final_scripts/500_rep_results`.

## What comes out

```
<output_dir>/
  K<k>/
    sim_parameters.json    sim$parameters of the raw file, plus its file name
    rep_<id>/              <id> zero-padded to four digits
      alpha.npy            L x p SuSiE alpha (float64)
      corr_sub.npy         float32 block of cor(X) over meta$corr_index
      meta.json            replicate_id, K, n, p, L, causal, pip, V, ld_block_name,
                           corr_index, susie_cs (sets with purity), export
                           settings, susieR version and path
      blip_arms.json       groups and PEPs per arm and q, settings, versions,
                           number of retained SNPs, SuSiE sets used,
                           candidate-group counts, pyblip's problem status
  metrics/                 blip_metrics.rds, blip_metrics_summary.csv, blip_rep_metrics.csv
  s7_<arm>/                <sim file base name>_blip.rds per K
  logs/                    local route only
```

All indices are 1-based (`index_base` in `meta.json`); `run_blip_arms.py`
converts internally and writes 1-based indices back.

`blip_metrics.rds` holds four tables: `rep_metrics` (per replicate x q x
method), `group_rows` (per set), `summary` (per q x method x K x stratum, with
`pooled` and `all` margins) and `cum_true` (for the same blocks, the number of
true sets of size <= k for k = 1..p, and that number over the causal SNPs).
The CSV files are `summary` and `rep_metrics`.

Each `s7_<arm>/` file holds `replicates[[i]]$credible_sets$susie_blip_q0.05$cs`
and `$susie_blip_q0.1$cs` (slot i is replicate i), the simulation parameters
and the arm's settings and versions. `extract_S7_data.R` uses the q = 0.05
groups of two arms: `D_purity_uncapped`, read from
`<benchmark>/final_scripts/500_rep_results/blip_purity0.5/`, and
`D_purity08_uncapped`, read from `blip_purity0.8/` beside it. Neither is
converted by default, so convert them (`convert_arms=` on the HPC route):

```
STEPS=combine CONVERT_ARMS="D_purity_uncapped D_purity08_uncapped" bash local/run_blip_local.sh
```

Then copy or symlink the files of `s7_D_purity_uncapped/` into
`blip_purity0.5/` and those of `s7_D_purity08_uncapped/` into
`blip_purity0.8/` (or run step 4 with `output_dir=` set to those folders), and
run `extract_S7_data.R` and `plot_S7.R` as usual.

## Settings

| Setting | Value |
|---|---|
| genotypes (`n`) | first 1000 rows of the LD block, mean-imputed, MAF >= 0.01, zero-variance columns dropped, first 5000 columns, scaled, as in `extract_S7_data.R` |
| SuSiE credible sets | `susie_get_cs(fit, Xcorr = cor(X), coverage, min_abs_corr = 0)` at coverage 0.95 and 0.90, each exported with its purity; those with purity >= 0.5 are checked against the Figure S7 call (`min_abs_corr = 0.5`) |
| purity seed | the replicate id, set before each `susie_get_cs` call (`purity_seed`): given `Xcorr` the development susieR computes purity over the whole set, whereas CRAN susieR subsamples 100 SNPs |
| PIPs | `susie_get_pip(fit)` |
| correlation block (`pip_floor`, `corr_set_purity`) | SNPs with PIP >= 0.01, members of exported sets with purity >= 0.5, causal SNPs |
| FDR levels (`--q`) | 0.05, paired with 95% credible sets, and 0.10, paired with 90% sets |
| PIP prefilter (`--pip-threshold`) | PIP >= 0.01 (the authors' wrapper keeps `1 - PIP <= 0.99`) |
| maximum PEP (`--max-pep`) | 0.25, for candidate groups and inside BLiP |
| size cap, capped arms (`--max-size`) | 25 |
| size cap, uncapped arms (`--uncapped-max-size`) | `auto`: max(100, largest SuSiE set at the paired coverage, any purity); an integer; or `0` for none (sequential windows up to the number of retained SNPs, SuSiE sets up to p), which can make the linear program large on diffuse replicates |
| prenarrowing (`--prenarrow`) | off |
| BLiP call | `weight_fn = 'inverse_size'`, `error = 'fdr'`, `deterministic = True`; `np.random`, with which BLiP perturbs its weights, seeded with the replicate id |
| solvers (`--solver`, `--bin-solver`) | GLPK (from cvxopt) for the linear program, passed explicitly because `pyblip`'s default falls back to ECOS, which cvxpy >= 1.7 does not ship; GLPK_MI for the integer step in `binarize_selections` |
| strata and capped SuSiE comparators (`size_cap`) | 25 variants |

- **Correlation block.** A full float32 `cor(X)` is ~100 MB per replicate
  (~250 GB for the 2500 replicates); the block makes the export a few hundred
  MB. Members of impure sets are left out: near-null components produce sets
  spanning most of the region, and no purity-matched arm can report such a
  set. A `--pip-threshold` below `pip_floor`, or a purity filter below
  `corr_set_purity`, requires a new export.
- **Candidate groups.** Over the SNPs passing the PIP prefilter,
  `pyblip.create_groups.sequential_groups` gives every contiguous window up to
  the size cap whose PEP (the product over components of 1 minus the window's
  alpha mass) is below the maximum. SuSiE's own sets are added from the full
  alphas: when a component's top alphas (as many as the size cap) sum to more
  than 1 - q, the smallest top-k set exceeding 1 - q, with PEP = 1 - its sum.
  Its SNPs outside the retained ones are appended to the retained index; a set
  that is already a sequential group is not added twice.
  `pyblip.create_groups.susie_groups` is not called directly: on the filtered
  alphas it returns nothing whenever a component's retained mass is below
  1 - q.

| Arm | Size cap | Extra candidate groups | Purity filter on candidate groups |
|---|---|---|---|
| `A_capped` | capped | none | none |
| `B_uncapped` | uncapped | none | none |
| `C_hier_capped`, `C_hier_uncapped` | as A, B | hierarchical clustering | none |
| `D_purity_capped`, `D_purity_uncapped` | as A, B | none | min abs(r) >= 0.5 |
| `D_purity08_uncapped` | as B | none | min abs(r) >= 0.8, the analogue of SuSiE's purity-0.8 sets |

Hierarchical groups come from single, average and complete linkage on the
distance `abs(1 - r)` over the retained SNPs, as `pyblip.susie_groups` builds
them when given X; a non-contiguous group within the size cap is added when
its PEP, computed from the alphas, is below the maximum. The purity filter is
`pyblip.create_groups.filter_by_purity(min_purity, corr_matrix)` on every
candidate group before BLiP; SuSiE's own sets whose exported purity is below
it are dropped first (`n_own_dropped_impure`): the filter would remove them
anyway, and the correlation block need not cover their members.

| Method in the metrics | Definition |
|---|---|
| `SuSiE(0.5)`, `SuSiE(0.8)` | sets with min abs(r) >= 0.5 / 0.8, as in `extract_S7_data.R` |
| `SuSiE(no purity)` | every set reaching the coverage |
| `SuSiE(0.5) size<=25` | `SuSiE(0.5)` sets of at most `size_cap` variants |
| `SuSiE capped(25)` | a component is kept only when its top `size_cap` alphas exceed 1 - q, no purity filter: the SuSiE comparator in the authors' UK Biobank analysis, and the sets step 2 feeds to BLiP in the capped arms |
| `BLiP <arm>` | every arm present in `blip_arms.json` |

Unweighted power and FDR come from `calculate_cs_metrics()` as
`extract_S7_data.R` defines them: power is the fraction of causal SNPs in the
union of sets; FDR is the fraction of sets without a causal SNP, 0 when there
are no sets. Resolution-adjusted power, the BLiP authors' metric, is the sum
of 1/|G| over the sets G containing a causal SNP, divided by K. Also given:
the FDP conditional on at least one set, the number, size and purity (min
abs(r) from `corr_sub.npy`) of the sets, and every summary within strata
defined by whether the replicate's largest `SuSiE(0.5)` set exceeds `size_cap`.

A figure caption for an arm states the settings above that apply to it,
whether SuSiE sets are shown as in the main analysis (purity >= 0.5) or under
BLiP's restriction (same size cap, no purity filter), both power definitions
and any stratification.
