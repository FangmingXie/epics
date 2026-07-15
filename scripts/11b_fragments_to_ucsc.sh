#!/usr/bin/env bash
# Stage 2b — Restrict the Wang25 ATAC fragments to primary UCSC chromosomes.
#
# Runtime finding (tabix -l on the raw index): the fragment file is MIXED —
#   primary chroms are ALREADY UCSC (chr1..chr22, chrX, chrY); only the scaffolds are
#   Ensembl-named (GL000009.2, KI270711.1, ...); there is no chrM.
# The on-disk fasta/blacklist/gene-annotation are UCSC with UCSC-named scaffolds
# (chrUn_GL...), so the Ensembl-named fragment scaffolds match nothing downstream and
# are not in our chromsizes. We therefore DROP scaffolds and keep the primary chroms
# as-is (no renaming), then re-index. Runs ONCE after Stage 0; all of Stage 3 reads the output.
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
FRAG_BGZ=${DATA}/atac_fragments.tsv.bgz              # input: raw fragments (40 GB; primary=UCSC, scaffolds=Ensembl)
FRAG_UCSC_BGZ=${DATA}/atac_fragments.ucsc.tsv.gz     # output: primary UCSC chroms only, bgzipped (.gz ext so pycisTopic gzip-opens it)
FRAG_UCSC_TBI=${FRAG_UCSC_BGZ}.tbi

# htslib (bgzip/tabix) is not in the `epics` env; use the chrombpnet env's binaries (htslib 1.22.1).
BGZIP=/home/qlyu/miniforge3/envs/chrombpnet/bin/bgzip
TABIX=/home/qlyu/miniforge3/envs/chrombpnet/bin/tabix

for tool in "${BGZIP}" "${TABIX}"; do
    [[ -x "${tool}" ]] || { echo "ERROR: ${tool} not found/executable." >&2; exit 1; }
done
if [[ ! -s "${FRAG_BGZ}" ]]; then
    echo "ERROR: ${FRAG_BGZ} missing/empty — Stage 0 download not complete." >&2
    exit 1
fi

echo "[$(date)] filtering fragments to primary UCSC chroms -> ${FRAG_UCSC_BGZ}"
# Keep only chr1..chr22, chrX, chrY, chrM (chrM absent here but harmless). Rows stay
# chrom-contiguous + position-sorted (filter only, no reorder) -> tabix-valid.
zcat "${FRAG_BGZ}" \
  | awk -F'\t' 'BEGIN{OFS="\t";
                 split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",a," ");
                 for(i in a) keep[a[i]]=1}
                /^#/{print; next}
                $1 in keep' \
  | "${BGZIP}" -c > "${FRAG_UCSC_BGZ}"

echo "[$(date)] indexing -> ${FRAG_UCSC_TBI}"
"${TABIX}" -p bed "${FRAG_UCSC_BGZ}"

echo "[$(date)] done."
ls -lh "${FRAG_UCSC_BGZ}" "${FRAG_UCSC_TBI}"
echo "--- head (expect chr* in col1, barcodes unchanged) ---"
zcat "${FRAG_UCSC_BGZ}" | grep -v '^#' | head -3
echo "--- chroms present (expect chr1..chr22, chrX, chrY only) ---"
"${TABIX}" -l "${FRAG_UCSC_BGZ}" | tr '\n' ' '; echo
