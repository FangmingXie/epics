#!/usr/bin/env bash
# Stage 3a (v2) — MACS3 peak calling on ALL IT-neuron ATAC fragments (hg38, UCSC chroms).
# New copy of 12_call_peaks.sh adapted to the all-IT scope (plan/04): pooled single MACS3 call
# over the 49,283 IT barcodes (PFC+V1, L2/3–L6), output to the isolated work_it/ dir.
# Run with plain `bash` — macs3/bedtools are addressed by explicit env paths below.
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
WORK=${DATA}/work_it                                  # NEW isolated work dir

FRAG_UCSC_BGZ=${DATA}/atac_fragments.ucsc.tsv.gz      # Stage 2b output (shared, UCSC, primary chroms)
IT_BARCODES=${DATA}/it_barcodes.txt                   # Stage 1 v2 output (49,283 barcodes)

SAMPLE=it
IT_BED=${WORK}/${SAMPLE}_fragments.bed                # IT fragment intervals (MACS3 input)
NARROWPEAK=${WORK}/${SAMPLE}_peaks.narrowPeak         # MACS3 output (→ Stage 3a consensus, 12b.v2)

# macs3 lives in the `macs3` env (not epics).
MACS3=/home/qlyu/miniforge3/envs/macs3/bin/macs3

# --- fail-fast checks ---
for f in "${FRAG_UCSC_BGZ}" "${IT_BARCODES}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
[[ -x "${MACS3}" ]] || { echo "ERROR: ${MACS3} not found/executable." >&2; exit 1; }
mkdir -p "${WORK}"

# =============================================================================
# STEP 1 — IT fragments → BED (keep fragments whose barcode is in IT_BARCODES)
# =============================================================================
echo "[$(date)] Step 1: extracting IT fragment intervals -> ${IT_BED}"
zcat "${FRAG_UCSC_BGZ}" \
  | awk -F'\t' -v OFS='\t' -v bcfile="${IT_BARCODES}" '
        BEGIN{ while ((getline line < bcfile) > 0) bc[line]=1 }
        /^#/ { next }
        ($4 in bc) { print $1, $2, $3 }' \
  > "${IT_BED}"
echo "  IT fragments: $(wc -l < "${IT_BED}")"

# =============================================================================
# STEP 2 — MACS3 callpeak (human, -g hs)
# =============================================================================
echo "[$(date)] Step 2: MACS3 callpeak"
"${MACS3}" callpeak \
    -f BED -t "${IT_BED}" -g hs -n "${SAMPLE}" --outdir "${WORK}" \
    --nomodel --shift -75 --extsize 150 -q 0.01 --call-summits \
    2> "${WORK}/macs3.log"
echo "  peaks called: $(wc -l < "${NARROWPEAK}")"

echo "[$(date)] MACS3 done -> ${NARROWPEAK}"
echo "  Next: scripts/12b.v2.consensus_peaks_it.py (epics env) builds the non-overlapping consensus."
