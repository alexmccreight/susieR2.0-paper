"""
AlphaGenome PIP-weighted credible set scoring for SuSiE method comparison.

Step 2 of the AlphaGenome scoring pipeline (extract_bvsr_results.R produces
the inputs; alphagenome_cs_group_comparison.R consumes the main output).

Scores variants in each 95% credible set with AlphaGenome using only the
tissue-matched modalities and tracks needed per variant. Each variant is
queried with only the modalities required by its tissue(s); results are
filtered to ontology-matched tracks before storage and PIP-weighted to a
single CS-level score per modality.

Inputs (from extract_bvsr_results.R):
    susie_{standard,ash,inf}_credible_sets.rds
    susie_{standard,ash,inf}_gene_tissue_summary.rds

Usage:
    python alphagenome_cs_scoring.py <API_KEY>
    python alphagenome_cs_scoring.py <API_KEY> --resume
    python alphagenome_cs_scoring.py <API_KEY> --test          # 5 variants only
    python alphagenome_cs_scoring.py --step1-only              # extract variants only
    python alphagenome_cs_scoring.py --summarize-only          # skip scoring

Output (numbered sub-steps within this script, unrelated to pipeline steps):
    alphagenome_cs_variants.csv        - all CS variants with metadata (sub-step 1)
    alphagenome_cs_raw_scores.parquet  - all tidy_scores rows (sub-step 2)
    alphagenome_cs_variant_scores.csv  - per-variant filtered scores (sub-step 3)
    alphagenome_cs_scores.csv          - PIP-weighted CS-level scores (sub-step 4)
                                       <-- input to alphagenome_cs_group_comparison.R
    alphagenome_cs_summary.txt         - distribution comparison (sub-step 5)
"""

import sys
import os
import json
import time
import argparse
import subprocess
import statistics
from pathlib import Path

import pandas as pd

# ===========================================================================
# Paths
# ===========================================================================
DATA_DIR = "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/alphagenome/data"
OUTPUT_DIR = "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark"

VARIANTS_FILE = os.path.join(OUTPUT_DIR, "alphagenome_cs_variants.csv")
RAW_SCORES_FILE = os.path.join(OUTPUT_DIR, "alphagenome_cs_raw_scores.parquet")
CHECKPOINT_FILE = os.path.join(OUTPUT_DIR, "alphagenome_cs_checkpoint.json")
VARIANT_SCORES_FILE = os.path.join(OUTPUT_DIR, "alphagenome_cs_variant_scores.csv")
CS_SCORES_FILE = os.path.join(OUTPUT_DIR, "alphagenome_cs_scores.csv")
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "alphagenome_cs_summary.txt")

# Tissue-to-track mapping for tissue-matched scoring
# Each eQTL tissue maps to specific AlphaGenome tracks per modality
TISSUE_TRACK_MAP = {
    "DLPFC": {
        "RNA_SEQ": "UBERON:0009834",
        "DNASE": "UBERON:0009834",
        "CHIP_HISTONE": "UBERON:0009834",
        "CHIP_TF": "UBERON:0009834",
    },
    "AC": {
        "RNA_SEQ": "UBERON:0009835",
        "CHIP_HISTONE": "UBERON:0002967",
    },
    "PCC": {
        "DNASE": "UBERON:0002740",
        "CHIP_HISTONE": "UBERON:0002967",
    },
}


# ===========================================================================
# Step 1: Extract all CS variants via R subprocess
# ===========================================================================
def extract_cs_variants():
    """Extract all variants in 95% credible sets for common gene-tissue pairs."""
    print("=" * 60)
    print("STEP 1: Extracting all CS variants")
    print("=" * 60)

    r_script = f"""
    data_dir <- "{DATA_DIR}"

    # Load summaries to find common gene-tissue pairs
    std_summary <- readRDS(file.path(data_dir, "susie_standard_gene_tissue_summary.rds"))
    ash_summary <- readRDS(file.path(data_dir, "susie_ash_gene_tissue_summary.rds"))
    inf_summary <- readRDS(file.path(data_dir, "susie_inf_gene_tissue_summary.rds"))

    common_pairs <- Reduce(intersect, list(
        paste(std_summary$gene_id, std_summary$tissue, sep=":"),
        paste(ash_summary$gene_id, ash_summary$tissue, sep=":"),
        paste(inf_summary$gene_id, inf_summary$tissue, sep=":")
    ))
    cat("Common gene-tissue pairs:", length(common_pairs), "\\n")

    extract_all_cs_variants <- function(method_name) {{
        mem <- readRDS(file.path(data_dir, paste0("susie_", method_name, "_variant_cs_membership.rds")))
        mem$pair <- paste(mem$gene_id, mem$tissue, sep=":")
        mem <- mem[mem$pair %in% common_pairs & mem$coverage_level == "0.95", ]
        mem$method <- method_name
        cat(method_name, ": ", nrow(mem), " variant entries\\n")
        mem[, c("variant_id", "gene_id", "tissue", "cs_id", "pip", "method")]
    }}

    std_vars <- extract_all_cs_variants("standard")
    ash_vars <- extract_all_cs_variants("ash")
    inf_vars <- extract_all_cs_variants("inf")

    all_vars <- rbind(std_vars, ash_vars, inf_vars)
    cat("Total variant entries:", nrow(all_vars), "\\n")
    cat("Unique variant_ids:", length(unique(all_vars$variant_id)), "\\n")

    write.csv(all_vars, "{VARIANTS_FILE}", row.names=FALSE)
    cat("Saved to {VARIANTS_FILE}\\n")
    """

    result = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True,
    )
    print(result.stdout)
    if result.returncode != 0:
        print("R stderr:", result.stderr)
        sys.exit(1)

    return load_cs_variants()


def load_cs_variants():
    """Load CS variants CSV and report summary stats."""
    df = pd.read_csv(VARIANTS_FILE)
    print(f"\nLoaded {len(df)} variant entries from {VARIANTS_FILE}")

    for method in ["standard", "ash", "inf"]:
        mdf = df[df["method"] == method]
        n_variants = mdf["variant_id"].nunique()
        n_cs = mdf.groupby(["gene_id", "tissue", "cs_id"]).ngroups
        print(f"  {method}: {len(mdf)} entries, {n_variants} unique variants, {n_cs} credible sets")

    unique_vids = df["variant_id"].nunique()
    print(f"\nTotal unique variant_ids to score: {unique_vids}")

    # Count by tissue
    for tissue in sorted(df["tissue"].unique()):
        n = df[df["tissue"] == tissue]["variant_id"].nunique()
        print(f"  {tissue}: {n} unique variants")

    return df


# ===========================================================================
# Step 2: Score unique variants with AlphaGenome (tissue-targeted)
# ===========================================================================

def build_variant_scoring_plan(df):
    """Build per-variant scoring requirements from the CS variants DataFrame.

    Returns:
        variant_modalities: dict mapping variant_id -> set of modality names needed
        needed_ontology_curies: flat set of all ontology CURIEs we care about
    """
    variant_tissues = df.groupby("variant_id")["tissue"].apply(set).to_dict()

    variant_modalities = {}
    all_curies = set()
    for vid, tissues in variant_tissues.items():
        modalities = set()
        for tissue in tissues:
            tissue_map = TISSUE_TRACK_MAP.get(tissue, {})
            modalities.update(tissue_map.keys())
            all_curies.update(tissue_map.values())
        variant_modalities[vid] = modalities

    return variant_modalities, all_curies


def filter_scores_by_ontology(scores, needed_curies):
    """Filter a list of AnnData score objects to only keep target ontology tracks.

    Args:
        scores: list of AnnData objects from score_variant()
        needed_curies: set of ontology CURIE strings to keep

    Returns:
        list of filtered AnnData objects (empty AnnDatas are dropped)
    """
    filtered = []
    for adata in scores:
        if "ontology_curie" in adata.var.columns:
            mask = adata.var["ontology_curie"].isin(needed_curies)
            if mask.any():
                filtered.append(adata[:, mask.values])
        else:
            filtered.append(adata)
    return filtered


def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r") as f:
            return set(json.load(f).get("scored_variants", []))
    return set()


def save_checkpoint(scored):
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump({"scored_variants": list(scored), "timestamp": time.time()}, f)


PARTS_DIR = os.path.join(OUTPUT_DIR, "alphagenome_cs_raw_scores_parts")


def _flush_parquet_part(rows, columns, part_num):
    """Write a batch of rows as a numbered parquet part file."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    os.makedirs(PARTS_DIR, exist_ok=True)
    part_path = os.path.join(PARTS_DIR, f"part_{part_num:05d}.parquet")
    df = pd.DataFrame(rows, columns=columns)
    table = pa.Table.from_pandas(df)
    pq.write_table(table, part_path, compression="snappy")


def merge_parquet_parts():
    """Merge all part files into the final parquet file."""
    import pyarrow.parquet as pq

    if not os.path.exists(PARTS_DIR):
        print("No parts directory found.")
        return

    part_files = sorted(
        f for f in os.listdir(PARTS_DIR) if f.endswith(".parquet")
    )
    if not part_files:
        print("No part files found.")
        return

    print(f"Merging {len(part_files)} part files...")
    dataset = pq.ParquetDataset(PARTS_DIR)
    table = dataset.read()
    pq.write_table(table, RAW_SCORES_FILE, compression="snappy")
    print(f"Merged into {RAW_SCORES_FILE} ({len(table)} rows)")


PARQUET_COLUMNS = [
    "variant_id", "gene_id", "gene_name", "gene_type",
    "output_type", "track_name", "gtex_tissue",
    "biosample_name", "ontology_curie", "raw_score", "quantile_score",
]


def score_variants(cs_df, api_key, resume=False, test=False):
    """Score unique variants with AlphaGenome using only tissue-needed modalities."""
    print("\n" + "=" * 60)
    print("STEP 2: Scoring variants with AlphaGenome (tissue-targeted)")
    print("=" * 60)

    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome

    dna_model = dna_client.create(api_key)
    organism = dna_client.Organism.HOMO_SAPIENS
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']

    # Map modality names to scorer objects
    modality_to_scorer = {
        "RNA_SEQ": variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'],
        "DNASE": variant_scorers.RECOMMENDED_VARIANT_SCORERS['DNASE'],
        "CHIP_HISTONE": variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_HISTONE'],
        "CHIP_TF": variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_TF'],
    }

    # Build per-variant scoring plan
    variant_modalities, needed_curies = build_variant_scoring_plan(cs_df)
    print(f"Target ontology CURIEs: {sorted(needed_curies)}")
    print(f"Sequence length: 1MB")

    # Report modality distribution
    from collections import Counter
    mod_counts = Counter()
    for mods in variant_modalities.values():
        for m in mods:
            mod_counts[m] += 1
    for mod, count in sorted(mod_counts.items()):
        print(f"  {mod}: {count} variants need this modality")

    # Sort for deterministic order, filter to autosomes
    variant_ids = sorted(variant_modalities.keys())
    auto_vids = []
    for vid in variant_ids:
        chrom = vid.split(":")[0].replace("chr", "")
        try:
            c = int(chrom)
            if 1 <= c <= 22:
                auto_vids.append(vid)
        except ValueError:
            pass
    print(f"Autosomal variants: {len(auto_vids)}")

    if test:
        auto_vids = auto_vids[:5]
        print(f"TEST MODE: scoring only {len(auto_vids)} variants")

    # Resume handling
    scored_ids = set()
    if resume:
        scored_ids = load_checkpoint()
        print(f"Resuming: {len(scored_ids)} already scored")

    remaining = [v for v in auto_vids if v not in scored_ids]
    total = len(remaining)
    print(f"Variants to score: {total}")

    batch_rows = []
    # Count existing part files so new parts get unique numbers
    os.makedirs(PARTS_DIR, exist_ok=True)
    part_num = len([f for f in os.listdir(PARTS_DIR) if f.endswith(".parquet")])
    start_time = time.time()
    consecutive_errors = 0

    for idx, vid in enumerate(remaining):
        parts = vid.split(":")
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]

        # Per-variant scorer selection based on tissue membership
        vid_mods = variant_modalities.get(vid, set())
        selected_scorers = [modality_to_scorer[m] for m in sorted(vid_mods)
                            if m in modality_to_scorer]
        if not selected_scorers:
            print(f"  SKIP {vid}: no modalities needed")
            scored_ids.add(vid)
            continue

        try:
            variant = genome.Variant(
                chromosome=chrom, position=pos,
                reference_bases=ref, alternate_bases=alt, name=vid,
            )
            interval = variant.reference_interval.resize(sequence_length)

            scores = dna_model.score_variant(
                interval=interval, variant=variant,
                variant_scorers=selected_scorers, organism=organism,
            )

            # Filter to only tissue-matched ontology tracks
            filtered_scores = filter_scores_by_ontology(scores, needed_curies)

            df = variant_scorers.tidy_scores([filtered_scores])

            if df is not None and len(df) > 0:
                for _, row in df.iterrows():
                    batch_rows.append([
                        vid,
                        row.get("gene_id", ""),
                        row.get("gene_name", ""),
                        row.get("gene_type", ""),
                        row.get("output_type", ""),
                        row.get("track_name", ""),
                        row.get("gtex_tissue", ""),
                        row.get("biosample_name", ""),
                        row.get("ontology_curie", ""),
                        row.get("raw_score", ""),
                        row.get("quantile_score", ""),
                    ])

                # In test mode, print detailed info for the first variant
                if test and idx == 0:
                    print(f"\n--- Scorers used: {sorted(vid_mods)} ---")
                    print(f"--- tidy_scores columns: {list(df.columns)} ---")
                    print(f"--- output_type values: {df['output_type'].unique().tolist()} ---")
                    print(f"--- Total rows for 1 variant (filtered): {len(df)} ---")
                    if "ontology_curie" in df.columns:
                        print(f"--- ontology_curie values: {df['ontology_curie'].dropna().unique().tolist()} ---")
                    for ot in df['output_type'].unique():
                        sub = df[df['output_type'] == ot]
                        bnames = sub['biosample_name'].dropna().unique()[:5].tolist()
                        print(f"--- {ot}: {len(sub)} rows, biosample_names: {bnames} ---")

            scored_ids.add(vid)
            consecutive_errors = 0

        except Exception as e:
            consecutive_errors += 1
            print(f"  ERROR {vid}: {e}")
            if consecutive_errors >= 10:
                print("  10 consecutive errors — saving checkpoint and stopping.")
                if batch_rows:
                    _flush_parquet_part(batch_rows, PARQUET_COLUMNS, part_num)
                    part_num += 1
                save_checkpoint(scored_ids)
                sys.exit(1)
            continue

        # Progress
        if (idx + 1) % 25 == 0 or (idx + 1) == total:
            elapsed = time.time() - start_time
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            eta = (total - idx - 1) / rate / 60 if rate > 0 else 0
            mods_str = "+".join(sorted(vid_mods))
            print(f"  [{idx+1}/{total}] {rate:.2f} var/sec, "
                  f"ETA: {eta:.0f} min, scored: {len(scored_ids)}, "
                  f"last: {mods_str}")

        # Flush to parquet part file every 100 variants
        if (idx + 1) % 100 == 0 and batch_rows:
            _flush_parquet_part(batch_rows, PARQUET_COLUMNS, part_num)
            part_num += 1
            batch_rows = []
            save_checkpoint(scored_ids)

    # Flush remaining rows
    if batch_rows:
        _flush_parquet_part(batch_rows, PARQUET_COLUMNS, part_num)

    save_checkpoint(scored_ids)
    print(f"\nScoring complete! {len(scored_ids)} variants scored.")

    # Merge part files into final parquet
    merge_parquet_parts()


# ===========================================================================
# Step 3: Per-variant tissue-matched scores
# ===========================================================================
def compute_variant_scores(cs_df):
    """Join CS variants with raw AlphaGenome scores, filtered to tissue-matched tracks.

    For each CS variant × modality, find the matching raw score row(s) and
    compute the mean quantile_score (averaging over strands/replicates).

    For RNA_SEQ (gene-centric), the join also matches on gene_id so each CS
    variant gets only the score for its specific gene.
    """
    print("\n" + "=" * 60)
    print("STEP 3: Computing per-variant tissue-matched scores")
    print("=" * 60)

    raw = pd.read_parquet(RAW_SCORES_FILE)
    print(f"Loaded {len(raw)} raw score rows from parquet")

    # Ensure numeric quantile_score
    raw["quantile_score"] = pd.to_numeric(raw["quantile_score"], errors="coerce")

    # Strip ENSEMBL version suffix for gene_id matching (e.g. ENSG...12 -> ENSG...)
    raw["gene_id_base"] = raw["gene_id"].fillna("").str.split(".").str[0]

    all_var_scores = []

    for tissue, track_map in TISSUE_TRACK_MAP.items():
        tissue_cs = cs_df[cs_df["tissue"] == tissue]
        if tissue_cs.empty:
            print(f"  {tissue}: no CS variants")
            continue

        n_variants = tissue_cs["variant_id"].nunique()
        print(f"\n  {tissue}: {len(tissue_cs)} CS entries, {n_variants} unique variants")

        for modality, curie in sorted(track_map.items()):
            # Filter raw scores to this modality + ontology track
            mod_scores = raw[
                (raw["output_type"] == modality) & (raw["ontology_curie"] == curie)
            ]

            if mod_scores.empty:
                print(f"    {modality} ({curie}): 0 raw rows — skipping")
                continue

            if modality == "RNA_SEQ":
                # Gene-centric: join on variant_id AND gene_id
                merged = tissue_cs.merge(
                    mod_scores[["variant_id", "gene_id_base", "quantile_score"]],
                    left_on=["variant_id", "gene_id"],
                    right_on=["variant_id", "gene_id_base"],
                    how="inner",
                    suffixes=("", "_score"),
                )
            else:
                # Variant-level: join on variant_id only
                merged = tissue_cs.merge(
                    mod_scores[["variant_id", "quantile_score"]],
                    on="variant_id",
                    how="inner",
                )

            if merged.empty:
                print(f"    {modality} ({curie}): 0 matched variants")
                continue

            # Average quantile_score per unique (variant, gene, tissue, method, cs_id)
            group_cols = ["variant_id", "gene_id", "tissue", "method", "cs_id", "pip"]
            var_agg = (
                merged.groupby(group_cols)["quantile_score"].mean().reset_index()
            )
            var_agg["modality"] = modality
            var_agg["ontology_curie"] = curie

            n_matched = var_agg["variant_id"].nunique()
            print(
                f"    {modality} ({curie}): {len(mod_scores)} raw rows, "
                f"{n_matched}/{n_variants} variants matched"
            )

            all_var_scores.append(var_agg)

    if not all_var_scores:
        print("No variant scores computed!")
        return pd.DataFrame()

    result = pd.concat(all_var_scores, ignore_index=True)
    result.to_csv(VARIANT_SCORES_FILE, index=False)
    print(f"\nTotal: {len(result)} variant-modality score entries")
    print(f"Saved to {VARIANT_SCORES_FILE}")
    return result


# ===========================================================================
# Step 4: PIP-weighted CS-level scores
# ===========================================================================
def compute_cs_scores(cs_df, var_scores_df):
    """Compute PIP-weighted CS-level scores per modality.

    CS_Score = sum(PIP_j * |quantile_score_j|) / sum(PIP_j)

    Null expectation for this metric is 0.5 since |Uniform[-1,1]| = Uniform[0,1].
    """
    print("\n" + "=" * 60)
    print("STEP 4: Computing PIP-weighted CS-level scores")
    print("=" * 60)

    # Get total CS sizes from the original CS data
    cs_sizes = (
        cs_df.groupby(["gene_id", "tissue", "method", "cs_id"])
        .size()
        .reset_index(name="cs_size")
    )

    group_cols = ["gene_id", "tissue", "method", "cs_id", "modality"]

    results = []
    for name, group in var_scores_df.groupby(group_cols):
        gene_id, tissue, method, cs_id, modality = name
        pip = group["pip"].values
        qs = group["quantile_score"].values
        valid = ~pd.isna(qs)

        n_scored = int(valid.sum())
        if n_scored == 0:
            cs_score = None
        else:
            pip_v = pip[valid]
            qs_v = qs[valid]
            cs_score = float((pip_v * abs(qs_v)).sum() / pip_v.sum())

        results.append({
            "gene_id": gene_id,
            "tissue": tissue,
            "method": method,
            "cs_id": cs_id,
            "modality": modality,
            "n_scored": n_scored,
            "cs_score": cs_score,
        })

    result = pd.DataFrame(results)

    # Attach total CS size
    result = result.merge(cs_sizes, on=["gene_id", "tissue", "method", "cs_id"], how="left")

    result.to_csv(CS_SCORES_FILE, index=False)
    print(f"Computed {len(result)} CS-modality scores")
    print(f"Saved to {CS_SCORES_FILE}")

    for method in ["standard", "ash", "inf"]:
        m_rows = result[result["method"] == method]
        n_cs = m_rows.groupby(["gene_id", "tissue", "cs_id"]).ngroups
        valid_scores = m_rows["cs_score"].dropna()
        print(
            f"  {method}: {len(m_rows)} scores across {n_cs} CS, "
            f"mean={valid_scores.mean():.4f}, median={valid_scores.median():.4f}"
        )

    return result


# ===========================================================================
# Step 5: Distribution comparison across methods
# ===========================================================================
def summarize_results(cs_scores_df):
    """Compare CS score distributions across methods with summary stats and tests."""
    print("\n" + "=" * 60)
    print("STEP 5: Distribution comparison across methods")
    print("=" * 60)

    lines = []
    lines.append("=" * 80)
    lines.append("AlphaGenome CS Score Distribution Comparison")
    lines.append("CS_Score = sum(PIP * |quantile_score|) / sum(PIP)")
    lines.append("Null expectation: 0.5 (|quantile| ~ Uniform[0,1])")
    lines.append("=" * 80)

    valid = cs_scores_df[cs_scores_df["cs_score"].notna()].copy()

    # ---- Overall summary per method × modality ----
    lines.append("\n--- Per Method x Modality Summary ---")
    lines.append(
        f"{'Method':<12} {'Modality':<15} {'N_CS':>6} {'Mean':>8} {'Median':>8} "
        f"{'Std':>8} {'>0.5':>6} {'>0.7':>6} {'%>0.5':>7}"
    )
    lines.append("-" * 85)

    for method in ["standard", "ash", "inf"]:
        for modality in ["RNA_SEQ", "DNASE", "CHIP_HISTONE", "CHIP_TF"]:
            mask = (valid["method"] == method) & (valid["modality"] == modality)
            sub = valid[mask]
            if sub.empty:
                continue

            scores = sub["cs_score"]
            n = len(scores)
            mean_s = scores.mean()
            med_s = scores.median()
            std_s = scores.std()
            above_05 = int((scores > 0.5).sum())
            above_07 = int((scores > 0.7).sum())
            pct_05 = 100.0 * above_05 / n

            lines.append(
                f"{method:<12} {modality:<15} {n:>6} {mean_s:>8.4f} {med_s:>8.4f} "
                f"{std_s:>8.4f} {above_05:>6} {above_07:>6} {pct_05:>6.1f}%"
            )

    # ---- Paired comparison across methods ----
    lines.append("\n\n--- Paired Comparison (same gene-tissue-cs_id scored by both methods) ---")

    for modality in ["RNA_SEQ", "DNASE", "CHIP_HISTONE", "CHIP_TF"]:
        mod_df = valid[valid["modality"] == modality]
        if mod_df.empty:
            continue

        lines.append(f"\n  {modality}:")

        pivot = mod_df.pivot_table(
            index=["gene_id", "tissue", "cs_id"],
            columns="method",
            values="cs_score",
        )

        for m1, m2 in [("ash", "standard"), ("ash", "inf"), ("standard", "inf")]:
            if m1 not in pivot.columns or m2 not in pivot.columns:
                continue
            paired = pivot[[m1, m2]].dropna()
            if len(paired) < 5:
                continue

            diff = paired[m1] - paired[m2]
            n_higher = int((diff > 0).sum())
            n_lower = int((diff < 0).sum())
            n_tied = int((diff == 0).sum())
            mean_diff = diff.mean()
            med_diff = diff.median()

            lines.append(
                f"    {m1} vs {m2}: n={len(paired)}, "
                f"mean_diff={mean_diff:+.4f}, median_diff={med_diff:+.4f}, "
                f"{m1}_higher={n_higher}, {m2}_higher={n_lower}, tied={n_tied}"
            )

            try:
                from scipy.stats import wilcoxon
                stat, pval = wilcoxon(paired[m1], paired[m2])
                lines.append(f"      Wilcoxon signed-rank: stat={stat:.1f}, p={pval:.2e}")
            except Exception:
                pass

    # ---- Per-tissue breakdown ----
    lines.append("\n\n--- Per Tissue x Method Summary ---")
    for tissue in ["AC", "DLPFC", "PCC"]:
        tissue_df = valid[valid["tissue"] == tissue]
        if tissue_df.empty:
            continue
        lines.append(f"\n  {tissue}:")
        for method in ["standard", "ash", "inf"]:
            sub = tissue_df[tissue_df["method"] == method]
            if sub.empty:
                continue
            n_cs = sub.groupby(["gene_id", "cs_id"]).ngroups
            mean_s = sub["cs_score"].mean()
            med_s = sub["cs_score"].median()
            above_05 = int((sub["cs_score"] > 0.5).sum())
            pct = 100.0 * above_05 / len(sub)
            lines.append(
                f"    {method:12s}: {n_cs:>5} CS, mean={mean_s:.4f}, "
                f"median={med_s:.4f}, >0.5: {above_05}/{len(sub)} ({pct:.1f}%)"
            )

    # ---- Method-unique CS analysis ----
    lines.append("\n\n--- Method-Unique CS (found by one method but not others) ---")
    all_cs_keys = valid.groupby(["gene_id", "tissue", "cs_id", "method"]).first().reset_index()
    cs_method_sets = (
        all_cs_keys.groupby(["gene_id", "tissue", "cs_id"])["method"]
        .apply(set)
        .reset_index()
    )

    for method in ["standard", "ash", "inf"]:
        others = {"standard", "ash", "inf"} - {method}
        unique_mask = cs_method_sets["method"].apply(
            lambda s: method in s and not s.intersection(others)
        )
        unique_cs = cs_method_sets[unique_mask]
        if unique_cs.empty:
            lines.append(f"\n  {method}-unique CS: 0")
            continue

        unique_keys = unique_cs[["gene_id", "tissue", "cs_id"]]
        unique_scores = valid.merge(unique_keys, on=["gene_id", "tissue", "cs_id"])
        unique_scores = unique_scores[unique_scores["method"] == method]

        lines.append(
            f"\n  {method}-unique CS: {len(unique_cs)} CS, "
            f"{len(unique_scores)} scores"
        )
        if not unique_scores.empty:
            mean_s = unique_scores["cs_score"].mean()
            med_s = unique_scores["cs_score"].median()
            above_05 = int((unique_scores["cs_score"] > 0.5).sum())
            pct = 100.0 * above_05 / len(unique_scores)
            lines.append(
                f"    mean={mean_s:.4f}, median={med_s:.4f}, "
                f">0.5: {above_05}/{len(unique_scores)} ({pct:.1f}%)"
            )

    text = "\n".join(lines)
    with open(SUMMARY_FILE, "w") as f:
        f.write(text)

    print(text)
    print(f"\nSaved to {SUMMARY_FILE}")

    return cs_scores_df


# ===========================================================================
# Main
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(description="AlphaGenome CS scoring (v2)")
    parser.add_argument("api_key", nargs="?", help="AlphaGenome API key")
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    parser.add_argument("--test", action="store_true", help="Score only 5 variants")
    parser.add_argument("--step1-only", action="store_true",
                        help="Only extract CS variants (Step 1)")
    parser.add_argument("--summarize-only", action="store_true",
                        help="Skip scoring, re-summarize existing results")
    args = parser.parse_args()

    # Step 1: Extract CS variants
    if args.step1_only:
        extract_cs_variants()
        return

    if args.summarize_only:
        if not os.path.exists(VARIANTS_FILE):
            print("No variants file found. Run Step 1 first.")
            sys.exit(1)
        if not os.path.exists(RAW_SCORES_FILE):
            print("No raw scores file found. Run Step 2 first.")
            sys.exit(1)

        df = load_cs_variants()

        # Step 3: Per-variant tissue-matched scores
        if os.path.exists(VARIANT_SCORES_FILE):
            print(f"\nLoading existing {VARIANT_SCORES_FILE}")
            var_scores = pd.read_csv(VARIANT_SCORES_FILE)
            print(f"  {len(var_scores)} variant-modality score entries")
        else:
            var_scores = compute_variant_scores(df)

        # Step 4: PIP-weighted CS-level scores
        cs_scores = compute_cs_scores(df, var_scores)

        # Step 5: Distribution comparison
        summarize_results(cs_scores)
        return

    if not args.api_key:
        print("Usage: python alphagenome_cs_scoring.py <API_KEY>")
        print("       python alphagenome_cs_scoring.py --step1-only")
        sys.exit(1)

    # Step 1
    if not args.resume or not os.path.exists(VARIANTS_FILE):
        df = extract_cs_variants()
    else:
        df = load_cs_variants()

    # Step 2: Score unique variants (tissue-targeted)
    score_variants(df, args.api_key, resume=args.resume, test=args.test)

    # Step 3: Per-variant tissue-matched scores
    var_scores = compute_variant_scores(df)

    # Step 4: PIP-weighted CS-level scores
    cs_scores = compute_cs_scores(df, var_scores)

    # Step 5: Distribution comparison
    summarize_results(cs_scores)


if __name__ == "__main__":
    main()
