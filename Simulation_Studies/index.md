# Simulation Studies

The simulation pipelines behind Figure 1C–F and Supplementary Figures S2–S4,
S6–S9 and S11. The designs are described in the manuscript's Methods and
Supplementary Notes; the scripts here are the ones that were run. They run on a
cluster and keep its paths (`/home/apm2217/...`) as defaults.

| Pipeline | Folder | Feeds |
|---|---|---|
| Fine-mapping | `scripts/fine_mapping/` | Figure 1C–F; Figures S2, S3, S6–S9 and S11 |
| SuSiE+BLiP re-run | `scripts/fine_mapping/blip_rerun/` | Figure S7 |
| Prediction | `scripts/prediction_cv/` | Figure S4 |

Each pipeline has a command generator, a per-replicate worker and a combiner.
The combined outputs are the raw simulation output that the `extract_*` scripts
under **Main Figures** and **Supplementary Figures** read.

## Fine-mapping (`scripts/fine_mapping/`)

- `command_generator_parallel.R` writes the per-job commands for the parameter grid.
- `simulation_script_parallel.R` is the worker: it simulates a phenotype and fits
  SuSiE, SuSiE-ash and SuSiE-inf.
- `combine_parallel_results.R` merges the per-job outputs into one file per scenario.

## SuSiE+BLiP re-run (`scripts/fine_mapping/blip_rerun/`)

Applies BLiP to the SuSiE fits of the sparse simulations with the BLiP authors'
Python package (`pyblip` 0.4.2), with the group-size cap, PEP cap, PIP prefilter
and purity filter as explicit settings.

- `command_generator_blip.R` writes the per-job commands.
- `export_blip_inputs.R` and `run_blip_arms.py` are the workers: they rebuild the
  genotypes and the SuSiE credible sets, then run BLiP under each setting.
- `blip_metrics.R` and `combine_blip_to_s7.R` compute the metrics and write the
  combined file that `Supplementary_Figures/S7/extract_S7_data.R` reads.
- `local/run_blip_local.sh` runs the same steps on one machine.

The `README.md` in that folder documents every setting.

## Prediction (`scripts/prediction_cv/`)

- `command_generator_parallel_cv.R` writes the per-job commands for all four scenarios.
- `prediction_cv_script_parallel.R` is the worker: 5-fold cross-validated
  prediction R² for SuSiE, SuSiE-ash and SuSiE-inf, with each model fit on the
  training fold only.
- `combine_parallel_cv_results.R` merges the per-job outputs.
- `local/run_prediction_cv_local.sh` runs the same steps on one machine.

The `README.md` in that folder documents the design.

## Elsewhere

The runtime benchmark behind Supplementary Figure S1 is in
`Supplementary_Figures/S1/data_generation/`.
