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

FRAG_UCSC_BGZ=${DATA}/atac_fragments.ucsc.tsv.bgz     # Stage 2b output (UCSC, primary chroms)
L23_BARCODES=${DATA}/v1_l23_barcodes.txt              # Stage 1 output (5,558 barcodes)

SAMPLE=l23
L23_BED=${WORK}/${SAMPLE}_fragments.bed               # L2/3 fragment intervals (MACS3 input)
NARROWPEAK=${WORK}/${SAMPLE}_peaks.narrowPeak         # MACS3 output (→ Stage 3a consensus, scripts/12b)

# macs3 lives in the `macs3` env (not epics).
MACS3=/home/qlyu/miniforge3/envs/macs3/bin/macs3

# --- fail-fast checks ---
for f in "${FRAG_UCSC_BGZ}" "${L23_BARCODES}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
[[ -x "${MACS3}" ]] || { echo "ERROR: ${MACS3} not found/executable." >&2; exit 1; }
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

echo "[$(date)] MACS3 done -> ${NARROWPEAK}"
echo "  Next: scripts/12b_consensus_peaks.py (epics env) builds the non-overlapping consensus."