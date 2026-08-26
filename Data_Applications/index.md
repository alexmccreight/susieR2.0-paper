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
| Genes with results | 2,765 |
| Gene–tissue pairs analysed | 7,941 (pairs where all three methods returned a fit) |
| Association windows | TADB enhanced cis-coding regions |
| LD reference | EUR |
| MAF / variant caps | `--min-twas-maf 0.05`, `--max-cv-variants 5000` |

Fine-mapping ran through the `mnm_regression.ipynb` `susie_twas` workflow of the
[xqtl-protocol](https://github.com/StatFunGen/xqtl-protocol). The three method
arms differ only in the `unmappable_effects` option:

| Method | Option | Script |
|---|---|---|
| SuSiE | *(none)* | `scripts/run_array_susie_standard.sh` |
| SuSiE-inf | `unmappable_effects='inf'` | `scripts/run_array_susie_inf.slurm` |
| SuSiE-ash | `unmappable_effects='ash'` | `scripts/run_array_susie_ash.slurm` |

All three are SLURM array jobs that shard `gene_list_3000.txt` across tasks and
write per-gene `*.univariate_bvsr.rds` files, prefixed `STANDARD-`, `INF-` and
`ASH-` respectively. Those files are the input to every Figure 2 panel.

`scripts/panel_E_candidate_loci.txt` is the shortlist of candidate loci
considered for the Figure 2 Panel E locus example; `ENSG00000163431` was the one
selected.

```{note}
These scripts are recorded as-run on the lab HPC and contain cluster-specific
absolute paths (`/mnt/lustre/...`) and SLURM directives. They are provided to
document the exact analysis rather than as a turnkey pipeline, and would need
their paths repointed to run elsewhere. The figure-building code in
**Main Figures** is fully portable and runs from vendored data.
```

## Downstream steps

Once fine-mapping completed, three further steps produced the Figure 2 inputs.
All live under `Main_Figures/figure_2/data/data_generation/alphagenome_scoring/`:

1. **`extract_bvsr_results.R`** — collects the per-gene fits into per-method
   credible-set and gene–tissue summary tables.
2. **`alphagenome_cs_scoring.py`** — scores credible-set variants with
   [AlphaGenome](https://github.com/google-deepmind/alphagenome) across RNA-seq,
   DNase, ChIP-histone and ChIP-TF modalities.
3. **`alphagenome_cs_group_comparison.R`** — matches credible sets across
   methods by Jaccard ≥ 0.75 (union-find over the resulting graph), assigns each
   a cross-method concordance group, and joins the AlphaGenome scores onto it.

Step 3 is parameterised by which methods contribute credible sets, which is how
the main figure (SuSiE + SuSiE-inf) and Supplementary Figure S11 (all three
methods) are generated from the same code.
