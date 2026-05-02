# AlphaGenome Scoring Pipeline

End-to-end pipeline that produces the AlphaGenome-derived inputs to Figure 2
panels B (UpSet) and D (functional enrichment):

- `data/alphagenome_cs_group_assignments.rds`
- `data/alphagenome_cs_group_scores.csv`

Plus the upstream `panel_A/` files used by panel A and re-used by every
downstream alphagenome step.

## Pipeline

```
*univariate_bvsr.rds                   (per-gene SuSiE fine-mapping; produced
                                        by an upstream pipeline not in this folder)
        │
        ▼
1. extract_bvsr_results.R              (parses raw fits into structured tables)
        │
        ├──> data/panel_A/                                ← Figure 2, Panel A input
        │      susie_{standard,ash,inf}_credible_sets.rds
        │      susie_{standard,ash,inf}_gene_tissue_summary.rds
        │
        ├──> susie_{...}_top_loci.rds                     (per-method top-loci tables)
        └──> susie_{...}_variant_cs_membership.rds        (per-CS variant lists)
        │
        ▼
2. alphagenome_cs_scoring.py           (PIP-weighted per-CS AlphaGenome scoring)
        │
        ▼
   alphagenome_cs_scores.csv            (per-CS scores; intermediate)
        │
        ▼
3. alphagenome_cs_group_comparison.R   (Jaccard ≥ 0.75 cross-method
                                        classification + score merge)
        │
        ├──> data/alphagenome_cs_group_assignments.rds   ← Fig 2 panels B, D
        └──> data/alphagenome_cs_group_scores.csv        ← Fig 2 panel D
```

## Step 1 — `extract_bvsr_results.R`

Parses raw per-gene `*univariate_bvsr.rds` fine-mapping outputs into
structured tables. For each of the three methods (standard, ash, inf) it
emits four files: `top_loci`, `credible_sets`, `gene_tissue_summary`, and
`variant_cs_membership`.

```bash
Rscript extract_bvsr_results.R
```

**Inputs:** raw `*univariate_bvsr.rds` files (the upstream fine-mapping
pipeline that generates these is **not yet** included in this repo).

**Outputs:** populates `data/panel_A/` with the 6 files used by Figure 2
panel A, plus `*_top_loci.rds` and `*_variant_cs_membership.rds` (used as
intermediates by Step 2).

## Step 2 — `alphagenome_cs_scoring.py`

Scores every variant in every 95% credible set across the three SuSiE methods
using the AlphaGenome API. Each variant is queried only with the modalities
required by its tissue(s); results are filtered to ontology-matched tracks
and PIP-weighted to a single CS-level score per modality.

Requires an AlphaGenome API key.

```bash
python alphagenome_cs_scoring.py <API_KEY>
python alphagenome_cs_scoring.py <API_KEY> --resume        # resume from checkpoint
python alphagenome_cs_scoring.py --step1-only              # extract variants only
python alphagenome_cs_scoring.py --summarize-only          # skip API; just resummarize
```

**Inputs:** `data/panel_A/susie_{...}_credible_sets.rds` and
`susie_{...}_gene_tissue_summary.rds` (from Step 1).

**Outputs** (written to the script's working directory):
- `alphagenome_cs_variants.csv` — all CS variants with metadata
- `alphagenome_cs_raw_scores.parquet` — all tidy_scores rows from the API
- `alphagenome_cs_variant_scores.csv` — per-variant filtered scores
- `alphagenome_cs_scores.csv` — PIP-weighted CS-level scores  ← **input to Step 3**
- `alphagenome_cs_summary.txt` — distribution comparison

## Step 3 — `alphagenome_cs_group_comparison.R`

Classifies each 95% credible set into one of seven cross-method groups based
on Jaccard ≥ 0.75 matching, then merges the AlphaGenome CS-level scores from
Step 2 with the group labels.

Groups:
1. **Consensus** — Standard + ASH + Inf
2. **Standard-specific** — Standard only
3. **ASH-specific** — ASH only
4. **Inf-specific** — Inf only
5. **Standard + ASH** — Standard + ASH (not Inf)
6. **Standard + Inf** — Standard + Inf (not ASH)
7. **ASH + Inf** — ASH + Inf (not Standard)

```bash
Rscript alphagenome_cs_group_comparison.R
```

**Inputs**
- `data/panel_A/susie_{...}_credible_sets.rds` (from Step 1)
- `alphagenome_cs_scores.csv` (from Step 2)

**Outputs** (the two files consumed by Figure 2):
- `alphagenome_cs_group_assignments.rds`
- `alphagenome_cs_group_scores.csv`

## Notes

- All three scripts hard-code paths into the benchmark workspace
  (`/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/...`). When porting,
  update each script's `data_dir` / `output_dir` to point at this repo's
  `figure_2/data/`.
- The AlphaGenome API is rate-limited; Step 2 supports `--resume` and writes
  checkpoint files so a multi-hour run can be resumed cleanly.
