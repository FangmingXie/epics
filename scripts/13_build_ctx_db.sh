#!/usr/bin/env bash
# Stage 3b — Build a custom cisTarget motif database on the V1 L2/3 consensus peaks.
# Mirrors v1astro 12.prep_cistopic_db.sh, adapted to hg38 (UCSC) + the aertslab v10nr motif set.
# MUST run inside the `epics` env (provides bedtools + cbust + python/pyarrow on PATH):
#   conda run --no-capture-output -n epics bash scripts/13_build_ctx_db.sh
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
WORK=${DATA}/work
DB=${DATA}/db

DB_BUILDER=/data/qlyu/code/trial/create_cisTarget_databases
HG38_FASTA=/home/qlyu/mydata/data/v1_multiome_ai/chrombpnet_tutorial/data/downloads/hg38.fa
HG38_CHROMSIZES=${DB}/hg38.chrom.sizes                # Stage 2 output (UCSC)
PEAKS_FILTERED=${WORK}/l23_summits_501_filtered.bed   # Stage 3a output
MOTIF_DIR=/data/qlyu/data/common/aertslab_motif_collection/v10nr_clust_public/singletons

OUT_MOTIF_LIST=${DB}/motifs.txt                       # list of motif .cb filenames
OUT_FASTA=${DB}/l23_regions_1kb_bg.fa                 # padded-background region FASTA
CTX_DB_PREFIX=${DB}/wang25_l23                        # → *.regions_vs_motifs.{rankings,scores}.feather

BG_PADDING=1000
N_CPU=20

# --- fail-fast checks (PATH tools come from the epics env) ---
command -v bedtools >/dev/null || { echo "ERROR: bedtools not on PATH — run under: conda run -n epics bash $0" >&2; exit 1; }
command -v cbust    >/dev/null || { echo "ERROR: cbust (Cluster-Buster) not on PATH — run under the epics env." >&2; exit 1; }
for f in "${HG38_FASTA}" "${HG38_CHROMSIZES}" "${PEAKS_FILTERED}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
[[ -d "${MOTIF_DIR}" ]] || { echo "ERROR: motif dir missing: ${MOTIF_DIR}" >&2; exit 1; }
mkdir -p "${DB}"

# --- motif list (one .cb filename per line) ---
ls "${MOTIF_DIR}" > "${OUT_MOTIF_LIST}"
echo "[$(date)] motifs listed: $(wc -l < "${OUT_MOTIF_LIST}") -> ${OUT_MOTIF_LIST}"

# --- padded-background FASTA for the consensus regions ---
echo "[$(date)] building padded-bg FASTA (${BG_PADDING} bp) -> ${OUT_FASTA}"
"${DB_BUILDER}/create_fasta_with_padded_bg_from_bed.sh" \
    "${HG38_FASTA}" "${HG38_CHROMSIZES}" "${PEAKS_FILTERED}" \
    "${OUT_FASTA}" "${BG_PADDING}" yes

# --- score regions × motifs → rankings/scores feathers ---
echo "[$(date)] scoring regions × motifs (t=${N_CPU}) -> ${CTX_DB_PREFIX}.*"
"${DB_BUILDER}/create_cistarget_motif_databases.py" \
    -f "${OUT_FASTA}" \
    -M "${MOTIF_DIR}" -m "${OUT_MOTIF_LIST}" \
    -o "${CTX_DB_PREFIX}" --bgpadding "${BG_PADDING}" -t "${N_CPU}"

echo "[$(date)] Stage 3b complete:"
ls -lh "${CTX_DB_PREFIX}".regions_vs_motifs.* 2>/dev/null