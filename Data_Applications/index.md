# Data Applications

This section contains the real-data application behind **Figure 2**: fine-mapping
of ROSMAP bulk-brain eQTLs with SuSiE, SuSiE-inf and SuSiE-ash, followed by
functional scoring of the resulting credible sets.

## Cohort and design

| | |
|---|---|
| Cohort | ROSMAP (De Jager) bulk-brain eQTL |
| Tissues | AC, DLPFC, PCC |
| Genes attempted | 3,000 (`scripts/gene_list_3000.txt`) |
| Genes with results | 2,765 (SuSiE), 2,744 (SuSiE-inf), 2,681 (SuSiE-ash) |
| Gene–tissue pairs analysed | 7,941, covering 2,654 genes (the pairs present in the output of all three methods) |
| Association windows | TADB enhanced cis-coding regions |
| LD reference | EUR |
| MAF / variant caps | `--min-twas-maf 0.05`, `--max-cv-variants 5000` |

Fine-mapping ran through the `mnm_regression.ipynb` `susie_twas` workflow of the
[xqtl-protocol](https://github.com/StatFunGen/xqtl-protocol). The three method
arms differ in the `unmappable_effects` option. The SuSiE-inf and SuSiE-ash arms
also set `convergence_method='pip'` and pass `--skip-twas-weights`, and they use
larger jobs (30 GB, 5 genes per task) than the SuSiE arm (14 GB, 8 genes per task).

| Method | Option | Script |
|---|---|---|
| SuSiE | *(none)* | `scripts/run_array_susie_standard.sh` |
| SuSiE-inf | `unmappable_effects='inf'` | `scripts/run_array_susie_inf.slurm` |
| SuSiE-ash | `unmappable_effects='ash'` | `scripts/run_array_susie_ash.slurm` |

All three are SLURM array jobs that shard `gene_list_3000.txt` across tasks and
write per-gene `*.univariate_bvsr.rds` files, prefixed `STANDARD-`, `INF-` and
`ASH-` respectively. Those files are the input to every Figure 2 panel.

`scripts/panel_E_candidate_loci.txt` is the shortlist of candidate loci
considered for the Figure 2E locus example; `ENSG00000163431` (*LMOD1*) was the one
selected. Its three result files are kept in `Main_Figures/figure_2/data/panel_E/`
without their per-sample fields (`sample_names` and `preset_variants_result`).

```{note}
These scripts are kept as they were run on the lab cluster, with its absolute
paths (`/mnt/lustre/...`) and SLURM settings. They document the analysis and
would need their paths changed to run elsewhere. The figure scripts in
**Main Figures** run from the data in this repository.
```

## Downstream steps

Three further steps produce the Figure 2 inputs. They live under
`Main_Figures/figure_2/data/data_generation/alphagenome_scoring/`, and the
README there has the details:

1. **`extract_bvsr_results.R`** — collects the per-gene fits into per-method
   credible-set and gene–tissue summary tables.
2. **`alphagenome_cs_scoring.py`** — scores credible-set variants with
   [AlphaGenome](https://github.com/google-deepmind/alphagenome) across RNA-seq,
   DNase, ChIP-histone and ChIP-TF modalities.
3. **`alphagenome_cs_group_comparison.R`** — matches credible sets across
   methods by Jaccard ≥ 0.75 (union-find over the resulting graph), assigns each
   a cross-method concordance group, and joins the AlphaGenome scores onto it.

Step 3 takes the methods to compare as an argument, which is how the main figure
(SuSiE and SuSiE-inf) and Supplementary Figure S10 (all three methods) come from
the same code.
