#!/bin/bash
# Figure 2 and Supplementary Figure S10: SuSiE-inf fine-mapping of the ROSMAP
# eQTL data (AC, DLPFC, PCC), as run on the lab cluster.
#
# SLURM array over the gene list: each task runs the xqtl-protocol susie_twas
# workflow (mnm_regression.ipynb) for LINES_PER_JOB genes, one at a time.
# This arm sets unmappable_effects='inf' and passes --skip-twas-weights.
# All paths are specific to the cluster.
#
# Input:  random_3000_overlap_genes.txt (gene list); ROSMAP genotype,
#         expression and covariate files
# Output: <home_dir>/susie_inf/susie_twas/ (the workflow's --cwd)
# Usage:  sbatch run_array_susie_inf.slurm

#SBATCH --job-name=susie_twas_array
#SBATCH --output=/mnt/lustre/lab/gwang/home/rl3328/array/logs/mnm_array_%A_%a.out
#SBATCH --error=/mnt/lustre/lab/gwang/home/rl3328/array/logs/mnm_array_%A_%a.err
#SBATCH --array=0-599
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --constraint="cpu4mem32a|cpu4mem32b|cpu8mem64a|cpu8mem64b|cpu16mem128a|cpu16mem128b"


export PATH="/mnt/lustre/lab/gwang/home/rl3328/.pixi/bin:$PATH"
source ~/.bashrc

TMPDIR=/mnt/lustre/lab/gwang/home/rl3328/tmp/slurm_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

home_dir=/mnt/lustre/home/rl3328/rl3328
cd ${home_dir}

trap "rm -rf $TMPDIR" EXIT

# Gene list: each line holds one Ensembl gene ID (gene_list_3000.txt in this folder).
OOM_FILE=/mnt/lustre/home/rl3328/rl3328/susie_ash/random_3000_overlap_genes.txt

# Task i takes lines i * LINES_PER_JOB + 1 to (i + 1) * LINES_PER_JOB of the list:
# 600 tasks x 5 genes = 3,000 genes.
LINES_PER_JOB=5  # 500/46 ≈ 10.87, round up
START_LINE=$((SLURM_ARRAY_TASK_ID * LINES_PER_JOB + 1))
END_LINE=$(( (SLURM_ARRAY_TASK_ID + 1) * LINES_PER_JOB ))

for LINE_NUM in $(seq $START_LINE $END_LINE); do
    FILENAME=$(awk -v line=$LINE_NUM 'NR==line {print}' $OOM_FILE)
    [ -z "$FILENAME" ] && break  # Stop if past end of file
    
    REGION_NAME=$(echo $FILENAME | sed -n 's/.*\(ENSG[0-9]\+\).*/\1/p')
    echo "[$(date)] Processing region: $REGION_NAME"
    sos run /mnt/lustre/home/rl3328/rl3328/xqtl-protocol/code/mnm_analysis/mnm_methods/mnm_regression.ipynb susie_twas \
        --name ROSMAP_DeJager_eQTL_twas \
        --genoFile /mnt/lustre/home/rl3328/rl3328/gpQTL_project/TensorQTL/genotype_meta/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.genotype_by_chrom_files.txt \
        --phenoFile \
            /mnt/lustre/data/gaowang/ROSMAP/eQTL/AC/phenotype/analysis_ready/AC_samples_list.rnaseqc.gene_tpm.low_expression_filtered.outlier_removed.tmm.expression.region_list.txt \
            /mnt/lustre/data/gaowang/ROSMAP/eQTL/DLPFC/phenotype/analysis_ready/DLPFC_samples_list.rnaseqc.gene_tpm.low_expression_filtered.outlier_removed.tmm.expression.remove_duplicates.region_list.txt \
            /mnt/lustre/data/gaowang/ROSMAP/eQTL/PCC/phenotype/analysis_ready/PCC_samples_list.rnaseqc.gene_tpm.low_expression_filtered.outlier_removed.tmm.expression.region_list.txt \
        --covFile  \
            /mnt/lustre/data/gaowang/ROSMAP/eQTL/AC/covariate/analysis_ready/AC_samples_list.rnaseqc.gene_tpm.low_expression_filtered.outlier_removed.tmm.expression.ROSMAP_xqtl_complete_samples_covariates_sex_death_pmi_study_transpose.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.unrelated.plink_qc.prune.pca.Marchenko_PC.gz \
            /mnt/lustre/data/gaowang/ROSMAP/eQTL/DLPFC/covariate/analysis_ready/DLPFC_samples_list.rnaseqc.gene_tpm.low_expression_filtered.outlier_removed.tmm.expression.ROSMAP_xqtl_complete_samples_covariates_sex_death_pmi_study_transpose.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.unrelated.plink_qc.prune.pca.Marchenko_PC.gz \
            /mnt/lustre/data/gaowang/ROSMAP/eQTL/PCC/covariate/analysis_ready/PCC_samples_list.rnaseqc.gene_tpm.low_expression_filtered.outlier_removed.tmm.expression.ROSMAP_xqtl_complete_samples_covariates_sex_death_pmi_study_transpose.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.unrelated.plink_qc.prune.pca.Marchenko_PC.gz \
        --customized_association_windows /mnt/lustre/data/gaowang/resource/TADB_enhanced_cis.coding.bed \
        --phenotype-names AC_DeJager_eQTL DLPFC_DeJager_eQTL PCC_DeJager_eQTL \
        --keep_variants /mnt/lustre/data/gaowang/GWAS/AD/AD_GWAS_NO_QC_NO_imputed_concatenate_sumstat/AD_Bellenguez_2022.sumstats.tsv.gz \
        --max-cv-variants 5000 \
        --min-twas-maf 0.05 \
        --ld_reference_meta_file /mnt/lustre/data/gaowang/resource/LD_EUR/ld_meta_file.tsv --cwd ./susie_inf/susie_twas/ -s build \
        --finemapping_extra_opts "list(unmappable_effects='inf',convergence_method='pip')" \
        --region-name $REGION_NAME --skip-twas-weights
done
