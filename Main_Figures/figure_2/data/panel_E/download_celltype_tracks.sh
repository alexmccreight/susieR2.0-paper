#!/usr/bin/env bash
# =============================================================================
# Download FILER cell-type-specific annotation tracks for ENSG00000163431 locus
# Region: chr1:201750000-201970000
# Genome build: hg38
#
# Cell types: Excitatory neurons, Glutamatergic neurons, Astrocytes,
#             Bipolar neurons, generic Neurons
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TABIX="/Users/alexmccreight/.pixi/bin/tabix"
REGION="chr1:201750000-201970000"
OUTDIR="${SCRIPT_DIR}/celltype_tracks"

mkdir -p "${OUTDIR}"

echo "Downloading cell-type-specific FILER tracks for ${REGION}..."
echo ""

FILER="https://tf.lisanwanglab.org/GADB/Annotationtracks"

# =============================================================================
# EXCITATORY NEURONS — Mint-ChIP-seq (histone marks)
# =============================================================================

echo "=== EXCITATORY NEURONS ==="

echo "  [1] Excitatory neuron H3K27ac (ENCFF126UEK) — active enhancers"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF126UEK.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_H3K27ac.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_H3K27ac.bed") intervals"

echo "  [2] Excitatory neuron H3K4me3 (ENCFF185EQB) — active promoters"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF185EQB.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_H3K4me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_H3K4me3.bed") intervals"

echo "  [3] Excitatory neuron H3K36me3 (ENCFF395ETQ) — gene body"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF395ETQ.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_H3K36me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_H3K36me3.bed") intervals"

echo "  [4] Excitatory neuron H3K27me3 (ENCFF431GVT) — polycomb repressed"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF431GVT.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_H3K27me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_H3K27me3.bed") intervals"

# --- Excitatory neuron DNase-seq (narrowpeak = open chromatin peaks) ---

echo "  [5] Excitatory neuron DNase narrowpeak rep1 (ENCFF108EIS)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF108EIS.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep1.bed") intervals"

echo "  [6] Excitatory neuron DNase narrowpeak rep2 (ENCFF446BGZ)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF446BGZ.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep2.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep2.bed") intervals"

echo "  [7] Excitatory neuron DNase narrowpeak rep3 (ENCFF521ZBJ)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF521ZBJ.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep3.bed") intervals"

echo "  [8] Excitatory neuron DNase narrowpeak rep4 (ENCFF579UIZ)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF579UIZ.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep4.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_DNase_narrowpeak_rep4.bed") intervals"

# --- Excitatory neuron DNase-seq (footprints = TF binding sites) ---

echo "  [9] Excitatory neuron DNase footprints rep1 (ENCFF030WVA)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF030WVA.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_DNase_footprints_rep1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_DNase_footprints_rep1.bed") intervals"

echo "  [10] Excitatory neuron DNase footprints rep2 (ENCFF481NUD)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF481NUD.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_DNase_footprints_rep2.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_DNase_footprints_rep2.bed") intervals"

# --- Excitatory neuron TF ChIP-seq ---

echo "  [11] Excitatory neuron CTCF (ENCFF816BTR)"
${TABIX} "${FILER}/ENCODE/data/TF-ChIP-seq/narrowpeak/hg38/2/ENCFF816BTR.bed.gz" ${REGION} \
  > "${OUTDIR}/ExcNeuron_CTCF.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/ExcNeuron_CTCF.bed") intervals"

echo ""

# =============================================================================
# GLUTAMATERGIC NEURONS — ATAC-seq (open chromatin)
# =============================================================================

echo "=== GLUTAMATERGIC NEURONS ==="

echo "  [12] Glutamatergic neuron ATAC-seq rep1 (ENCFF319MVQ)"
${TABIX} "${FILER}/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF319MVQ.bed.gz" ${REGION} \
  > "${OUTDIR}/GlutNeuron_ATAC_rep1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/GlutNeuron_ATAC_rep1.bed") intervals"

echo "  [13] Glutamatergic neuron ATAC-seq rep2 (ENCFF425MGV)"
${TABIX} "${FILER}/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF425MGV.bed.gz" ${REGION} \
  > "${OUTDIR}/GlutNeuron_ATAC_rep2.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/GlutNeuron_ATAC_rep2.bed") intervals"

echo "  [14] Glutamatergic neuron ATAC-seq rep3 (ENCFF444WLW)"
${TABIX} "${FILER}/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF444WLW.bed.gz" ${REGION} \
  > "${OUTDIR}/GlutNeuron_ATAC_rep3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/GlutNeuron_ATAC_rep3.bed") intervals"

echo "  [15] Glutamatergic neuron ATAC-seq rep4 (ENCFF857ZIX)"
${TABIX} "${FILER}/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF857ZIX.bed.gz" ${REGION} \
  > "${OUTDIR}/GlutNeuron_ATAC_rep4.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/GlutNeuron_ATAC_rep4.bed") intervals"

echo ""

# =============================================================================
# ASTROCYTES — Mint-ChIP-seq (histone marks)
# =============================================================================

echo "=== ASTROCYTES ==="

echo "  [16] Astrocyte H3K27ac (ENCFF321HNK) — active enhancers"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF321HNK.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_H3K27ac.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_H3K27ac.bed") intervals"

echo "  [17] Astrocyte H3K4me1 (ENCFF823AHN) — enhancers"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF823AHN.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_H3K4me1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_H3K4me1.bed") intervals"

echo "  [18] Astrocyte H3K4me3 (ENCFF611BIT) — active promoters"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF611BIT.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_H3K4me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_H3K4me3.bed") intervals"

echo "  [19] Astrocyte H3K27me3 (ENCFF313EKJ) — polycomb repressed"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF313EKJ.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_H3K27me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_H3K27me3.bed") intervals"

echo "  [20] Astrocyte H3K9me3 (ENCFF321TVF) — heterochromatin"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF321TVF.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_H3K9me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_H3K9me3.bed") intervals"

echo "  [21] Astrocyte H3K36me3 (ENCFF365YMX) — gene body"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF365YMX.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_H3K36me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_H3K36me3.bed") intervals"

# --- Astrocyte DNase-seq (narrowpeak) — pick 3 representative replicates ---

echo "  [22] Astrocyte DNase narrowpeak rep1 (ENCFF149LZX)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF149LZX.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_DNase_narrowpeak_rep1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_DNase_narrowpeak_rep1.bed") intervals"

echo "  [23] Astrocyte DNase narrowpeak rep2 (ENCFF391KKV)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF391KKV.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_DNase_narrowpeak_rep2.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_DNase_narrowpeak_rep2.bed") intervals"

echo "  [24] Astrocyte DNase narrowpeak rep3 (ENCFF419HUP)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF419HUP.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_DNase_narrowpeak_rep3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_DNase_narrowpeak_rep3.bed") intervals"

# --- Astrocyte DNase-seq (footprints) ---

echo "  [25] Astrocyte DNase footprints rep1 (ENCFF077UPK)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF077UPK.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_DNase_footprints_rep1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_DNase_footprints_rep1.bed") intervals"

echo "  [26] Astrocyte DNase footprints rep2 (ENCFF164JOB)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF164JOB.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_DNase_footprints_rep2.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_DNase_footprints_rep2.bed") intervals"

echo "  [27] Astrocyte DNase footprints rep3 (ENCFF334FGD)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF334FGD.bed.gz" ${REGION} \
  > "${OUTDIR}/Astrocyte_DNase_footprints_rep3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Astrocyte_DNase_footprints_rep3.bed") intervals"

echo ""

# =============================================================================
# BIPOLAR NEURONS — DNase-seq
# =============================================================================

echo "=== BIPOLAR NEURONS ==="

echo "  [28] Bipolar neuron DNase narrowpeak rep1 (ENCFF094CBM)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF094CBM.bed.gz" ${REGION} \
  > "${OUTDIR}/BipolarNeuron_DNase_narrowpeak_rep1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/BipolarNeuron_DNase_narrowpeak_rep1.bed") intervals"

echo "  [29] Bipolar neuron DNase narrowpeak rep2 (ENCFF142GIM)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF142GIM.bed.gz" ${REGION} \
  > "${OUTDIR}/BipolarNeuron_DNase_narrowpeak_rep2.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/BipolarNeuron_DNase_narrowpeak_rep2.bed") intervals"

echo "  [30] Bipolar neuron DNase narrowpeak rep3 (ENCFF635QFN)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/narrowpeak/hg38/ENCFF635QFN.bed.gz" ${REGION} \
  > "${OUTDIR}/BipolarNeuron_DNase_narrowpeak_rep3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/BipolarNeuron_DNase_narrowpeak_rep3.bed") intervals"

echo "  [31] Bipolar neuron DNase footprints rep1 (ENCFF545WGR)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF545WGR.bed.gz" ${REGION} \
  > "${OUTDIR}/BipolarNeuron_DNase_footprints_rep1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/BipolarNeuron_DNase_footprints_rep1.bed") intervals"

echo "  [32] Bipolar neuron DNase footprints rep2 (ENCFF624WEZ)"
${TABIX} "${FILER}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF624WEZ.bed.gz" ${REGION} \
  > "${OUTDIR}/BipolarNeuron_DNase_footprints_rep2.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/BipolarNeuron_DNase_footprints_rep2.bed") intervals"

echo ""

# =============================================================================
# GENERIC NEURON — Mint-ChIP-seq
# =============================================================================

echo "=== NEURON (generic) ==="

echo "  [33] Neuron H3K4me1 (ENCFF705VHY) — enhancers"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF705VHY.bed.gz" ${REGION} \
  > "${OUTDIR}/Neuron_H3K4me1.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Neuron_H3K4me1.bed") intervals"

echo "  [34] Neuron H3K9me3 (ENCFF549SPA) — heterochromatin"
${TABIX} "${FILER}/ENCODE/data/Mint-ChIP-seq/narrowpeak/hg38/ENCFF549SPA.bed.gz" ${REGION} \
  > "${OUTDIR}/Neuron_H3K9me3.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${OUTDIR}/Neuron_H3K9me3.bed") intervals"

echo ""
echo "============================================"
echo "Done. All files saved to: ${OUTDIR}"
echo "============================================"
