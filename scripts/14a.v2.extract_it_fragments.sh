#!/usr/bin/env bash
# Stage 3c-prep (v2) — extract an IT-only fragment file (all 5 columns) from the full UCSC fragments.
# New copy of 14a_extract_l23_fragments.sh for the all-IT scope (plan/04). create_cistopic_object_
# from_fragments reads the ENTIRE fragment file (all ~243k nuclei, 42 GB) and filters by valid_bc
# internally — slow to iterate. Pre-filter once to the 49,283 IT barcodes → a smaller file (~14 GB,
# ~9× the L2/3 file). We also hard-validate columns (exactly 5 cols, integer count) so no malformed
# row can reach pycisTopic's Score.astype(int32).
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
FRAG_UCSC_BGZ=${DATA}/atac_fragments.ucsc.tsv.gz     # Stage 2b (shared: full, UCSC, primary chroms)
IT_BARCODES=${DATA}/it_barcodes.txt                  # Stage 1 v2 (49,283 barcodes)
IT_FRAG_GZ=${DATA}/it_fragments.ucsc.tsv.gz          # OUTPUT: IT-only fragments (bgzipped)
IT_FRAG_TBI=${IT_FRAG_GZ}.tbi

BGZIP=/home/qlyu/miniforge3/envs/chrombpnet/bin/bgzip
TABIX=/home/qlyu/miniforge3/envs/chrombpnet/bin/tabix

for t in "${BGZIP}" "${TABIX}"; do [[ -x "$t" ]] || { echo "ERROR: $t missing" >&2; exit 1; }; done
for f in "${FRAG_UCSC_BGZ}" "${IT_BARCODES}"; do [[ -s "$f" ]] || { echo "ERROR: $f missing/empty" >&2; exit 1; }; done

echo "[$(date)] extracting IT fragments (barcode∈list, exactly 5 cols, integer count) -> ${IT_FRAG_GZ}"
zcat "${FRAG_UCSC_BGZ}" \
  | awk -F'\t' -v OFS='\t' -v bcfile="${IT_BARCODES}" '
        BEGIN{ while ((getline b < bcfile) > 0) keep[b]=1 }
        /^#/ { next }
        (NF==5 && $5 ~ /^[0-9]+$/ && ($4 in keep)) { print }
        (NF!=5 || $5 !~ /^[0-9]+$/) { dropped++ }
        END{ if (dropped) print "  dropped " dropped " malformed rows" > "/dev/stderr" }' \
  | "${BGZIP}" -c > "${IT_FRAG_GZ}"

echo "[$(date)] indexing -> ${IT_FRAG_TBI}"
"${TABIX}" -p bed "${IT_FRAG_GZ}"

echo "[$(date)] done."
ls -lh "${IT_FRAG_GZ}" "${IT_FRAG_TBI}"
echo "  IT fragment rows: $(zcat "${IT_FRAG_GZ}" | wc -l)"
echo "  chroms: $("${TABIX}" -l "${IT_FRAG_GZ}" | tr '\n' ' ')"
