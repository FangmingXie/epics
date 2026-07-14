#!/usr/bin/env bash
# Stage 3a — MACS3 peak calling on the V1 L2/3 ATAC fragments (hg38, UCSC chroms).
# Mirrors v1astro 10.call_peaks_macs3.sh, adapted to: subset fragments to L2/3 barcodes,
# human genome (-g hs), and the UCSC/primary-chrom fragment file from Stage 2b.
# Run with plain `bash` — macs3/bedtools are addressed by explicit env paths below.
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
WORK=${DATA}/work
DB=${DATA}/db

FRAG_UCSC_BGZ=${DATA}/atac_fragments.ucsc.tsv.bgz     # Stage 2b output (UCSC, primary chroms)
L23_BARCODES=${DATA}/v1_l23_barcodes.txt              # Stage 1 output (5,558 barcodes)
HG38_BLACKLIST=${DB}/hg38-blacklist.v2.bed            # Stage 2 output

SAMPLE=l23
L23_BED=${WORK}/${SAMPLE}_fragments.bed               # L2/3 fragment intervals (MACS3 input)
NARROWPEAK=${WORK}/${SAMPLE}_peaks.narrowPeak         # MACS3 output
SUMMITS=${WORK}/${SAMPLE}_summits.bed                 # MACS3 output (--call-summits)
EXTENDED=${WORK}/${SAMPLE}_summits_501.bed            # summits ±250 → 501 bp
PEAKS_FILTERED=${WORK}/${SAMPLE}_summits_501_filtered.bed   # blacklist-filtered consensus peaks

# macs3 + bedtools live in the `macs3` env (not epics).
MACS3=/home/qlyu/miniforge3/envs/macs3/bin/macs3
BEDTOOLS=/home/qlyu/miniforge3/envs/macs3/bin/bedtools

# --- fail-fast checks ---
for f in "${FRAG_UCSC_BGZ}" "${L23_BARCODES}" "${HG38_BLACKLIST}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
for t in "${MACS3}" "${BEDTOOLS}"; do
    [[ -x "${t}" ]] || { echo "ERROR: ${t} not found/executable." >&2; exit 1; }
done
mkdir -p "${WORK}"

# =============================================================================
# STEP 1 — L2/3 fragments → BED (keep fragments whose barcode is in L23_BARCODES)
# =============================================================================
echo "[$(date)] Step 1: extracting L2/3 fragment intervals -> ${L23_BED}"
zcat "${FRAG_UCSC_BGZ}" \
  | awk -F'\t' -v OFS='\t' -v bcfile="${L23_BARCODES}" '
        BEGIN{ while ((getline line < bcfile) > 0) bc[line]=1 }
        /^#/ { next }
        ($4 in bc) { print $1, $2, $3 }' \
  > "${L23_BED}"
echo "  L2/3 fragments: $(wc -l < "${L23_BED}")"

# =============================================================================
# STEP 2 — MACS3 callpeak (human, -g hs)
# =============================================================================
echo "[$(date)] Step 2: MACS3 callpeak"
"${MACS3}" callpeak \
    -f BED -t "${L23_BED}" -g hs -n "${SAMPLE}" --outdir "${WORK}" \
    --nomodel --shift -75 --extsize 150 -q 0.01 --call-summits \
    2> "${WORK}/macs3.log"
echo "  peaks called: $(wc -l < "${NARROWPEAK}")"

# =============================================================================
# STEP 3 — Extend summits ±250 → 501 bp (clamp start ≥ 0), filter blacklist
# =============================================================================
echo "[$(date)] Step 3: extend summits and filter blacklist"
awk 'BEGIN{OFS="\t"} /^chr/ { s=$2-250; if (s<0) s=0; print $1, s, $3+250, $4 }' \
    "${SUMMITS}" > "${EXTENDED}"
"${BEDTOOLS}" intersect -v -a "${EXTENDED}" -b "${HG38_BLACKLIST}" > "${PEAKS_FILTERED}"
echo "  extended regions:  $(wc -l < "${EXTENDED}")"
echo "  after blacklist:   $(wc -l < "${PEAKS_FILTERED}")  -> ${PEAKS_FILTERED}"

echo "[$(date)] Stage 3a complete."