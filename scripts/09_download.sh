#!/usr/bin/env bash
# Stage 0 — Acquire Wang25 paired 10x snMultiome (CELLxGENE deposit).
# RNA h5ad (2.8 GB) + ATAC fragments (42 GB) + tabix index. ~45 GB total.
# Both files come from the same deposit so RNA<->ATAC barcodes match.
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
GEX_H5AD=${DATA}/GEX_with_raw_counts.h5ad
FRAG_BGZ=${DATA}/atac_fragments.tsv.bgz
FRAG_TBI=${DATA}/atac_fragments.tsv.bgz.tbi

# --- source URLs (CELLxGENE datasets) ---
GEX_URL=https://datasets.cellxgene.cziscience.com/a4310202-4dc8-4e1b-a96d-d9675f5b14d1.h5ad
FRAG_URL=https://datasets.cellxgene.cziscience.com/c9ba15d8-45a1-4f68-ac91-a59ad2d44a2d-fragment.tsv.bgz
FRAG_TBI_URL=https://datasets.cellxgene.cziscience.com/c9ba15d8-45a1-4f68-ac91-a59ad2d44a2d-fragment.tsv.bgz.tbi

mkdir -p "${DATA}"

# -c = resume if partially downloaded; fail-fast on any error (set -e + wget non-zero exit)
echo "[$(date)] downloading RNA h5ad -> ${GEX_H5AD}"
wget -c -O "${GEX_H5AD}" "${GEX_URL}"

echo "[$(date)] downloading ATAC fragments -> ${FRAG_BGZ}"
wget -c -O "${FRAG_BGZ}" "${FRAG_URL}"

echo "[$(date)] downloading tabix index -> ${FRAG_TBI}"
wget -c -O "${FRAG_TBI}" "${FRAG_TBI_URL}"

echo "[$(date)] done. sizes:"
ls -lh "${GEX_H5AD}" "${FRAG_BGZ}" "${FRAG_TBI}"
