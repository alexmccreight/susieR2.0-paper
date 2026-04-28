#!/usr/bin/env bash
# =============================================================================
# Download FILER annotation tracks for ENSG00000163431 locus (chr1:201750000-201970000)
# Tissue: Dorsolateral prefrontal cortex (DLPFC) — Middle frontal area BA46 / E073
# Genome build: hg38
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TABIX="/Users/alexmccreight/.pixi/bin/tabix"
REGION="chr1:201750000-201970000"

FILER_BASE="https://tf.lisanwanglab.org/GADB/Annotationtracks"
HISTONE="${FILER_BASE}/ENCODE/data/Histone_ChIP-seq/narrowpeak/hg38/1"

echo "Downloading FILER tracks for ${REGION} (DLPFC)..."
echo ""

# --- Active marks ---

echo "  [1/8] H3K27ac — active enhancers + promoters (ENCFF183BLP)"
${TABIX} "${HISTONE}/ENCFF183BLP.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/H3K27ac_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/H3K27ac_DLPFC.bed") intervals"

echo "  [2/8] H3K4me3 — active promoters (ENCFF916WSX)"
${TABIX} "${HISTONE}/ENCFF916WSX.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/H3K4me3_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/H3K4me3_DLPFC.bed") intervals"

echo "  [3/8] H3K4me1 — enhancers (ENCFF538IHT)"
${TABIX} "${HISTONE}/ENCFF538IHT.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/H3K4me1_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/H3K4me1_DLPFC.bed") intervals"

echo "  [4/8] H3K9ac — active promoters (ENCFF370ANN)"
${TABIX} "${HISTONE}/ENCFF370ANN.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/H3K9ac_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/H3K9ac_DLPFC.bed") intervals"

echo "  [5/8] H3K36me3 — transcribed gene bodies (ENCFF117QYN)"
${TABIX} "${HISTONE}/ENCFF117QYN.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/H3K36me3_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/H3K36me3_DLPFC.bed") intervals"

# --- Open chromatin ---

echo "  [6/8] DNase-seq — open chromatin footprints (ENCFF988POL)"
${TABIX} "${FILER_BASE}/ENCODE/data/DNase-seq/bed5-FDR_footprints/hg38/ENCFF988POL.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/DNase_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/DNase_DLPFC.bed") intervals"

# --- Enhancer-specific tracks ---

echo "  [7/8] ROADMAP Enhancers — E073 ChromHMM enhancer states"
${TABIX} "${FILER_BASE}/ROADMAP_Enhancers/ChIP-seq/bed4/hg38/E073_15_coreMarks_hg38lift_mnemonics_enhancer.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/ROADMAP_enhancers_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/ROADMAP_enhancers_DLPFC.bed") intervals"

echo "  [8/8] EpiMap Enhancers — BSS01271 ChromHMM enhancer calls"
${TABIX} "${FILER_BASE}/EpiMap_enhancers/ChromHMM/bed9/hg38/formatted_output_BSS01271_Roadmap_2015_Enh_hg38.bed.gz" ${REGION} \
  > "${SCRIPT_DIR}/EpiMap_enhancers_DLPFC.bed" 2>/dev/null || true
echo "    -> $(wc -l < "${SCRIPT_DIR}/EpiMap_enhancers_DLPFC.bed") intervals"

echo ""
echo "Done. Files saved to: ${SCRIPT_DIR}"
