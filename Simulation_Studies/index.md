# Simulation Studies

This directory holds the code used to benchmark susieR 2.0 across two task settings. Methodology and parameter choices are described in the paper's Methods section; the scripts here are the exact pipeline used to generate the results shown in **Figure 1 (panels B–F)** and **Supplementary Figures S1, S3, S4, S5**.

## Pipelines

Both pipelines follow a 3-step HPC pattern (generator → parallel worker → combiner).

### Fine-mapping (`scripts/fine_mapping/`)

- `command_generator_parallel.R` — emits per-job HPC commands for the parameter grid
- `simulation_script_parallel.R` — per-rep worker: simulates phenotype + fits SuSiE, SuSiE-ash, SuSiE-inf
- `combine_parallel_results.R` — merges per-job outputs into one rds per scenario

Outputs feed: Figure 1 (panels C–F), Figures S1, S4, S5.

### Prediction CV (`scripts/prediction_cv/`)

- `command_generator_parallel_cv.R` — emits per-job HPC commands
- `prediction_cv_script_parallel.R` — per-rep worker: 5-fold cross-validated prediction R² for SuSiE, SuSiE-ash, SuSiE-inf
- `combine_parallel_cv_results.R` — merges per-job outputs

Outputs feed: Figure S3.
