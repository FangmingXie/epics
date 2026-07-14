#!/usr/bin/env bash
# Stage 2b — Normalize the Wang25 ATAC fragments to UCSC chrom naming + primary chromosomes.
# CELLxGENE fragments are Ensembl-style (1, 2, MT, GL...) but the fasta/blacklist/tutorials are UCSC
# (chr1, ...). A mismatch silently yields empty region<->fragment overlaps -> empty GRN. This step maps
# 1..22/X/Y -> chr*, MT -> chrM, drops all scaffolds/patches, then re-indexes with tabix.
# Runs ONCE after the Stage 0 download completes. All of Stage 3 reads the output, not the raw file.
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
FRAG_BGZ=${DATA}/atac_fragments.tsv.bgz              # input: raw Ensembl-style fragments (42 GB)
FRAG_UCSC_BGZ=${DATA}/atac_fragments.ucsc.tsv.bgz    # output: UCSC + primary chroms, bgzipped
FRAG_UCSC_TBI=${FRAG_UCSC_BGZ}.tbi

if [[ ! -s "${FRAG_BGZ}" ]]; then
    echo "ERROR: ${FRAG_BGZ} missing/empty — Stage 0 download not complete." >&2
    exit 1
fi

echo "[$(date)] converting fragments to UCSC + primary chroms -> ${FRAG_UCSC_BGZ}"
# rows stay chrom-contiguous and position-sorted (rename only, no reorder) -> tabix-valid.
zcat "${FRAG_BGZ}" \
  | awk -F'\t' 'BEGIN{OFS="\t";
                 split("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y",a," ");
                 for(i in a) keep[a[i]]=1}
                /^#/{print; next}
                { c=$1; if(c=="MT") c="chrM"; else if(c in keep) c="chr" c; else next; $1=c; print }' \
  | bgzip -c > "${FRAG_UCSC_BGZ}"

echo "[$(date)] indexing -> ${FRAG_UCSC_TBI}"
tabix -p bed "${FRAG_UCSC_BGZ}"

echo "[$(date)] done."
ls -lh "${FRAG_UCSC_BGZ}" "${FRAG_UCSC_TBI}"
echo "--- head (expect chr* in col1, barcodes unchanged) ---"
zcat "${FRAG_UCSC_BGZ}" | grep -v '^#' | head -3
echo "--- chroms present ---"
tabix -l "${FRAG_UCSC_BGZ}" | tr '\n' ' '; echo
