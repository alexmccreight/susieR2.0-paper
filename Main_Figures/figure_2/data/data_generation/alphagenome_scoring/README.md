# AlphaGenome scoring pipeline

Builds the cross-method concordance groups and the AlphaGenome credible-set
(CS) scores shown in Figure 2, panels C (concordance) and D (AlphaGenome
scores), and in the same panels of Supplementary Figure S10. The panel scripts
are in `Main_Figures/figure_2/figure_2_panel_C/` (concordance) and
`Main_Figures/figure_2/figure_2_panel_D/` (scores).

| Step | Script | Needs |
|---|---|---|
| 1 | `extract_bvsr_results.R` | the per-gene fine-mapping results (mounted bucket) |
| 2 | `alphagenome_cs_scoring.py` | an AlphaGenome API key |
| 3 | `alphagenome_cs_group_comparison.R` | only files in this repository |

Steps 1 and 2 read and write outside the repository; the files that step 3
reads are stored in `Main_Figures/figure_2/data/`. The methods are keyed
`standard` (SuSiE), `ash` (SuSiE-ash) and `inf` (SuSiE-inf).

## Step 1: `extract_bvsr_results.R`

Flattens the per-gene fine-mapping results (`*univariate_bvsr.rds`, written by
the jobs in `Data_Applications/scripts/`) into tables. Needs the R packages
`dplyr` and `tidyr`.

```bash
Rscript extract_bvsr_results.R                           # all three methods
Rscript extract_bvsr_results.R method='susie_standard'   # or susie_ash, susie_inf
Rscript extract_bvsr_results.R input_dir='/path' output_dir='/path'
```

**Input.** Every `.rds` file in these folders under `input_dir` (default
`/home/apm2217/data`, the mount point of
`statfungen/ftp_fgc_xqtl/interactive_sessions/apm2217`):

| Method | Folder |
|---|---|
| `susie_standard` | `real_data_analysis/susie_standard/fine_mapping` |
| `susie_ash` | `real_data_analysis/updated_ash_inf_02012026/susie_ash/susie_twas/fine_mapping` |
| `susie_inf` | `real_data_analysis/updated_ash_inf_02012026/susie_inf/susie_twas/fine_mapping` |

**Output.** In `output_dir` (default `/home/apm2217/output`), for each method:

| File | Content | Read by |
|---|---|---|
| `<method>_credible_sets.rds` | one row per credible set (coverage 0.95, 0.7 and 0.5) | Figure 2 panel A, step 3 |
| `<method>_gene_tissue_summary.rds` | one row per gene-tissue pair | Figure 2 panel A, steps 2 and 3 |
| `<method>_variant_cs_membership.rds` | one row per variant of each credible set | step 2 |
| `<method>_top_loci.rds` | top-loci variants with PIP and CS membership | |
| `<method>_failed_files.txt` | files that could not be processed, if any | |

## Step 2: `alphagenome_cs_scoring.py`

Scores the variants of the 95% credible sets with AlphaGenome and reduces the
scores to one value per credible set and modality. Needs an AlphaGenome API
key, the Python packages `alphagenome`, `pandas` and `pyarrow` (`scipy` for
the tests in the summary), and `Rscript` on the `PATH`.

```bash
python alphagenome_cs_scoring.py <API_KEY>              # full run
python alphagenome_cs_scoring.py <API_KEY> --resume     # continue from the checkpoint
python alphagenome_cs_scoring.py <API_KEY> --test       # score 5 variants only
python alphagenome_cs_scoring.py --step1-only           # only list the CS variants
python alphagenome_cs_scoring.py --summarize-only       # no API calls: rebuild the CS scores
                                                        # and summary from the saved raw scores
```

**Paths.** `DATA_DIR` and `OUTPUT_DIR` at the top of the script are hard-coded
absolute paths (`/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/alphagenome/data`
and `/Users/alexmccreight/StatFunGen/susieR2.0-benchmark`). Edit them before
running.

**Input.** In `DATA_DIR`, from step 1:
`susie_{standard,ash,inf}_gene_tissue_summary.rds` and
`susie_{standard,ash,inf}_variant_cs_membership.rds`.

**Method.** The variants are those of the 95% credible sets in gene-tissue
pairs present for all three methods. Each unique autosomal variant is scored
once, in a 1 Mb window, for the modalities of its tissues only (RNA_SEQ, DNASE,
CHIP_HISTONE, CHIP_TF), and only the tracks with the ontology term set for that
tissue and modality are kept (`TISSUE_TRACK_MAP` in the script). The score of a
variant for a modality is the mean `quantile_score` of those tracks (for
RNA_SEQ, of the credible set's own gene). The CS score is
sum(PIP x |quantile_score|) / sum(PIP) over the scored variants of the
credible set; its null expectation is 0.5.

**Output.** In `OUTPUT_DIR`:

| File | Content |
|---|---|
| `alphagenome_cs_variants.csv` | all CS variants with gene, tissue, CS, PIP and method |
| `alphagenome_cs_raw_scores.parquet` | the AlphaGenome score rows, merged from `alphagenome_cs_raw_scores_parts/` |
| `alphagenome_cs_checkpoint.json` | variants already scored, for `--resume` |
| `alphagenome_cs_variant_scores.csv` | one score per variant, credible set and modality |
| `alphagenome_cs_scores.csv` | one score per credible set and modality: the input of step 3 |
| `alphagenome_cs_summary.txt` | score distributions by method, modality and tissue |

The AlphaGenome API is rate-limited. Rows and checkpoint are saved every 100
variants, and the run stops after 10 consecutive errors, so a multi-hour run
can be continued with `--resume`.

## Step 3: `alphagenome_cs_group_comparison.R`

Matches the 95% credible sets across methods, labels each signal by the
methods that found it, and attaches the AlphaGenome CS scores. Run it from the
repository root, once per figure (`run_all.R` runs both commands):

```bash
# SuSiE and SuSiE-inf: main Figure 2
Rscript Main_Figures/figure_2/data/data_generation/alphagenome_scoring/alphagenome_cs_group_comparison.R

# all three methods: Supplementary Figure S10
Rscript Main_Figures/figure_2/data/data_generation/alphagenome_scoring/alphagenome_cs_group_comparison.R standard,ash,inf
```

The argument is a comma-separated list of at least two of `standard`, `ash`,
`inf`. Group labels and the output suffix follow the order given, so keep the
order `standard`, `ash`, `inf`.

**Input.** In `Main_Figures/figure_2/data/`:
`panel_A/susie_{standard,ash,inf}_gene_tissue_summary.rds` (always all three),
`panel_A/susie_<method>_credible_sets.rds` for the methods compared, and
`alphagenome_cs_scores.csv`.

**Method.** Only gene-tissue pairs listed in the `gene_tissue_summary` of all
three methods are analysed (7,941 pairs), whichever methods are compared.
Within a pair, two credible sets from different methods match when the Jaccard
index of their variant sets is >= 0.75. A signal is a connected component of
the match graph (union-find), labelled by the methods it contains:

| Run | Groups |
|---|---|
| default | `Consensus`, `Standard-specific`, `Inf-specific` |
| `standard,ash,inf` | `Consensus`, `Standard+ASH`, `Standard+Inf`, `ASH+Inf`, `Standard-specific`, `ASH-specific`, `Inf-specific` |

**Output.** In `Main_Figures/figure_2/data/`:

| Default run | `standard,ash,inf` run | Content |
|---|---|---|
| `alphagenome_cs_group_assignments__standard_inf.rds` | `alphagenome_cs_group_assignments.rds` | one row per credible set, with `signal_id` and `group` |
| `alphagenome_cs_group_scores__standard_inf.csv` | `alphagenome_cs_group_scores.csv` | the CS scores of the methods compared, with `signal_id` and `group` |
| `alphagenome_cs_group_summary__standard_inf.txt` | `alphagenome_cs_group_summary.txt` | scores by group, modality and tissue; Wilcoxon rank-sum tests of Consensus against each other group |

`figure_2_panel_C/prepare_panel_C_data.R` reads the assignments;
`figure_2_panel_D/prepare_panel_D_data.R` reads the scores and the assignments.
