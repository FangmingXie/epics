#!/usr/bin/env bash
# Stage 3c prep — extract an L2/3-only fragment file (all 5 columns) from the full UCSC fragments.
# Why: create_cistopic_object_from_fragments reads the ENTIRE fragment file (all ~243k nuclei, 40 GB)
# and filters by valid_bc internally — slow to iterate. Pre-filtering to the 5,558 L2/3 barcodes yields
# a small file (fast cisTopic reads) and drops every non-L2/3 row. We also hard-validate columns
# (exactly 5 cols, integer count) so no malformed row can reach pycisTopic's Score.astype(int32).
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
FRAG_UCSC_BGZ=${DATA}/atac_fragments.ucsc.tsv.gz     # Stage 2b (full, UCSC, primary chroms)
L23_BARCODES=${DATA}/v1_l23_barcodes.txt             # Stage 1 (5,558 barcodes)
L23_FRAG_GZ=${DATA}/l23_fragments.ucsc.tsv.gz        # OUTPUT: L2/3-only fragments (bgzipped)
L23_FRAG_TBI=${L23_FRAG_GZ}.tbi

BGZIP=/home/qlyu/miniforge3/envs/chrombpnet/bin/bgzip
TABIX=/home/qlyu/miniforge3/envs/chrombpnet/bin/tabix

for t in "${BGZIP}" "${TABIX}"; do [[ -x "$t" ]] || { echo "ERROR: $t missing" >&2; exit 1; }; done
for f in "${FRAG_UCSC_BGZ}" "${L23_BARCODES}"; do [[ -s "$f" ]] || { echo "ERROR: $f missing/empty" >&2; exit 1; }; done

echo "[$(date)] extracting L2/3 fragments (barcode∈list, exactly 5 cols, integer count) -> ${L23_FRAG_GZ}"
zcat "${FRAG_UCSC_BGZ}" \
  | awk -F'\t' -v OFS='\t' -v bcfile="${L23_BARCODES}" '
        BEGIN{ while ((getline b < bcfile) > 0) keep[b]=1 }
        /^#/ { next }
        (NF==5 && $5 ~ /^[0-9]+$/ && ($4 in keep)) { print }
        (NF!=5 || $5 !~ /^[0-9]+$/) { dropped++ }
        END{ if (dropped) print "  dropped " dropped " malformed rows" > "/dev/stderr" }' \
  | "${BGZIP}" -c > "${L23_FRAG_GZ}"

echo "[$(date)] indexing -> ${L23_FRAG_TBI}"
"${TABIX}" -p bed "${L23_FRAG_GZ}"

echo "[$(date)] done."
ls -lh "${L23_FRAG_GZ}" "${L23_FRAG_TBI}"
echo "  L2/3 fragment rows: $(zcat "${L23_FRAG_GZ}" | wc -l)"
echo "  chroms: $("${TABIX}" -l "${L23_FRAG_GZ}" | tr '\n' ' ')"